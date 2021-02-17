
#include <CNS.H>
#include <CNS_F.H>

#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EB_LeastSquares_3D_K.H>

using namespace amrex;

Real
CNS::advance (Real time, Real dt, int iteration, int ncycle)
{
    BL_PROFILE("CNS::advance()");
        
    for (int i = 0; i < num_state_data_types; ++i) {
        state[i].allocOldData();
        state[i].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab dSdt(grids,dmap,NUM_STATE,0,MFInfo(),Factory());
    MultiFab Sborder(grids,dmap,NUM_STATE,NUM_GROW,MFInfo(),Factory());

 
    MultiFab& C_new = get_new_data(Cost_Type);
    C_new.setVal(0.0);

    EBFluxRegister* fr_as_crse = nullptr;
    if (do_reflux && level < parent->finestLevel()) {
        CNS& fine_level = getLevel(level+1);
        fr_as_crse = &fine_level.flux_reg;
    }

    EBFluxRegister* fr_as_fine = nullptr;
    if (do_reflux && level > 0) {
        fr_as_fine = &flux_reg;
    }

    if (fr_as_crse) {
        fr_as_crse->reset();
    }

    // RK2 stage 1
    FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, NUM_STATE);
    compute_dSdt(Sborder, dSdt, 0.5*dt, fr_as_crse, fr_as_fine);
    // U^* = U^n + dt*dUdt^n
    MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, dSdt, 0, 0, NUM_STATE, 0);
    computeTemp(S_new,0);
    
    // RK2 stage 2
    // After fillpatch Sborder = U^n+dt*dUdt^n
    FillPatch(*this, Sborder, NUM_GROW, time+dt, State_Type, 0, NUM_STATE);
    compute_dSdt(Sborder, dSdt, 0.5*dt, fr_as_crse, fr_as_fine);
    // S_new = 0.5*(Sborder+S_old) = U^n + 0.5*dt*dUdt^n
    MultiFab::LinComb(S_new, 0.5, Sborder, 0, 0.5, S_old, 0, 0, NUM_STATE, 0);
    // S_new += 0.5*dt*dSdt
    MultiFab::Saxpy(S_new, 0.5*dt, dSdt, 0, 0, NUM_STATE, 0);
    // We now have S_new = U^{n+1} = (U^n+0.5*dt*dUdt^n) + 0.5*dt*dUdt^*
    computeTemp(S_new,0);
    
    return dt;
}

void
CNS::compute_dSdt (const MultiFab& S,  MultiFab& dSdt, Real dt,
                   EBFluxRegister* fr_as_crse, EBFluxRegister* fr_as_fine)
{
    BL_PROFILE("CNS::compute_dSdt()");

    const Real* dx = geom.CellSize();
    const int ncomp = dSdt.nComp();

    int as_crse = (fr_as_crse != nullptr);
    int as_fine = (fr_as_fine != nullptr);


    MultiFab& cost = get_new_data(Cost_Type);


   //////////////////////////////////////////////////////////////////////////
   // 	Define variables required by least squares
   //////////////////////////////////////////////////////////////////////////

    MultiFab Prim(grids,dmap,1,NUM_GROW,MFInfo(),Factory());
    MultiFab Prim_eb(grids,dmap,1,NUM_GROW,MFInfo(),Factory());
    MultiFab grad_eb(grids,dmap,1,NUM_GROW,MFInfo(),Factory());
    MultiFab grad_x(grids,dmap,1,NUM_GROW,MFInfo(),Factory());
    MultiFab grad_y(grids,dmap,1,NUM_GROW,MFInfo(),Factory());

    auto const& fact = dynamic_cast<EBFArrayBoxFactory const&>(S.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();

    const Box& domain_box = geom.Domain();

    const int domlo_x = domain_box.smallEnd(0)-5;
    const int domhi_x = domain_box.bigEnd(0)+5;
    const int domlo_y = domain_box.smallEnd(1)-5;
    const int domhi_y = domain_box.bigEnd(1)+5;
    const int domlo_z = domain_box.smallEnd(2)-5;
    const int domhi_z = domain_box.bigEnd(2)+5;

    const bool on_x_face = !(geom.isPeriodic(0));
    const bool on_y_face = !(geom.isPeriodic(1));
    const bool on_z_face = !(geom.isPeriodic(2));


#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {
        std::array<FArrayBox,AMREX_SPACEDIM> flux;
        FArrayBox dm_as_fine(Box::TheUnitBox(),ncomp);
        FArrayBox fab_drho_as_crse(Box::TheUnitBox(),ncomp);
        IArrayBox fab_rrflag_as_crse(Box::TheUnitBox());

        for (MFIter mfi(S, MFItInfo().EnableTiling(hydro_tile_size).SetDynamic(true));
                        mfi.isValid(); ++mfi)
    //for (MFIter mfi (S); mfi.isValid(); ++mfi)
        {
            auto wt = amrex::second();

            const Box& bx = mfi.tilebox();

            const auto& flag = flags[mfi];

//////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// Compute the required gradients using least squares here //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

	// Compute temperature and velocities and store in a MultiFab named Prim



            //if (1){
            if (flag.getType(bx) == FabType::covered) {
                dSdt[mfi].setVal<RunOn::Host>(0.0, bx, 0, ncomp);
            } else {

                // flux is used to store centroid flux needed for reflux
                for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
                    flux[idim].resize(amrex::surroundingNodes(bx,idim),ncomp);
                }

                if (flag.getType(amrex::grow(bx,1)) == FabType::regular)
                {
                    cns_compute_dudt(BL_TO_FORTRAN_BOX(bx),
                                     BL_TO_FORTRAN_ANYD(dSdt[mfi]),
                                     BL_TO_FORTRAN_ANYD(S[mfi]),
                                     BL_TO_FORTRAN_ANYD(flux[0]),
                                     BL_TO_FORTRAN_ANYD(flux[1]),
                                     BL_TO_FORTRAN_ANYD(flux[2]),
                                     dx, &dt,&level);

                    if (fr_as_crse) {
                        fr_as_crse->CrseAdd(mfi,{&flux[0],&flux[1],&flux[2]},dx,dt,RunOn::Cpu);
                    }

                    if (fr_as_fine) {
                        fr_as_fine->FineAdd(mfi,{&flux[0],&flux[1],&flux[2]},dx,dt,RunOn::Cpu);
                    }
                }
                else
                {
                    FArrayBox* p_drho_as_crse = (fr_as_crse) ?
                        fr_as_crse->getCrseData(mfi) : &fab_drho_as_crse;
                    const IArrayBox* p_rrflag_as_crse = (fr_as_crse) ?
                        fr_as_crse->getCrseFlag(mfi) : &fab_rrflag_as_crse;

                    if (fr_as_fine) {
                        dm_as_fine.resize(amrex::grow(bx,1),ncomp);
                    }

			compute_primitive_variables(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_ANYD(S[mfi]),
                                       BL_TO_FORTRAN_ANYD(Prim[mfi]));

	// Set the eb to a specified Dirichlet value in a MultiFab named Prim_eb

			set_eb_dirichlet_for_prim(BL_TO_FORTRAN_BOX(bx),
				  BL_TO_FORTRAN_ANYD(Prim_eb[mfi]));

	// Initialize eb gradient value to zero

	//set_grad_eb_val(BL_TO_FORTRAN_BOX(bx),
	//		BL_TO_FORTRAN_ANYD(grad_eb[mfi]));
	
	/// All required variables 
			Array4<const Real> const& phi_arr     = Prim.array(mfi);
			Array4<const Real> const& phi_eb_arr  = Prim_eb.array(mfi);
			Array4<      Real> const& grad_eb_arr = grad_eb.array(mfi);
			Array4<      Real> const& grad_x_arr  = grad_x.array(mfi);
			Array4<      Real> const& grad_y_arr  = grad_y.array(mfi);

	//print_grad_eb(BL_TO_FORTRAN_BOX(bx),
	//		BL_TO_FORTRAN_ANYD(flag),
	//		BL_TO_FORTRAN_ANYD(Prim[mfi]));

			Array4<Real const> const& fcx   = (fact.getFaceCent())[0]->const_array(mfi);
			Array4<Real const> const& fcy   = (fact.getFaceCent())[1]->const_array(mfi);
        		Array4<Real const> const& fcz   = (fact.getFaceCent())[2]->const_array(mfi);
        		Array4<Real const> const& vfrac = (fact.getVolFrac()).array(mfi);
        		Array4<Real const> const& ccent = (fact.getCentroid()).array(mfi);
        		Array4<Real const> const& bcent = (fact.getBndryCent()).array(mfi);
        		Array4<Real const> const& norm  = (fact.getBndryNormal()).array(mfi);
        		Array4<Real const> const& apx   = (fact.getAreaFrac())[0]->const_array(mfi);
        		Array4<Real const> const& apy   = (fact.getAreaFrac())[1]->const_array(mfi);
        		Array4<Real const> const& apz   = (fact.getAreaFrac())[2]->const_array(mfi);

			const FabArray<EBCellFlagFab>* iflags = &(fact.getMultiEBCellFlagFab());
			Array4<EBCellFlag const> const& iflag = iflags->const_array(mfi);

			amrex::ParallelFor(amrex::grow(bx,3), 1,
        			[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        		{
				Real nx = norm(i,j,k,0);
				Real ny = norm(i,j,k,1);
				Real nz = norm(i,j,k,2);

				Real yloc_on_xface = fcx(i,j,k,0);
				Real xloc_on_yface = fcy(i,j,k,0);
				Real zloc_on_xface = fcx(i,j,k,1);
				Real zloc_on_yface = fcy(i,j,k,1);
				Real xloc_on_zface = fcz(i,j,k,0);
				Real yloc_on_zface = fcz(i,j,k,1);

				bool is_eb_inhomog  = true;
				bool is_eb_dirichlet  = true;

				bool needs_bdry_stencil = (on_x_face && (i <= domlo_x || i >= domhi_x)) ||
                                      (on_y_face && (j <= domlo_y || j >= domhi_y));
				needs_bdry_stencil = needs_bdry_stencil ||
					(on_z_face && (k <= domlo_z || k >= domhi_z));

				
				if( iflag(i,j,k).isRegular() || iflag(i,j,k).isSingleValued()){

                 				grad_x_arr(i,j,k,n) = (apx(i,j,k) == 0.0) ? 0.0 :
                   				grad_x_of_phi_on_centroids(i, j, k, n, phi_arr, phi_eb_arr,
                                              			iflag, ccent, bcent,
                                              			yloc_on_xface, zloc_on_xface, is_eb_dirichlet, is_eb_inhomog);
	
						 grad_y_arr(i,j,k,n) = (apy(i,j,k) == 0.0) ? 0.0:
			                         grad_y_of_phi_on_centroids(i, j, k, n, phi_arr, phi_eb_arr,
                               			               iflag, ccent, bcent,
                                              			xloc_on_yface, zloc_on_yface, is_eb_dirichlet, is_eb_inhomog);				
			
				}
				if( iflag(i,j,k).isCovered()){
              				grad_x_arr(i,j,k,n)  = 0.0;
              				grad_y_arr(i,j,k,n)  = 0.0;
				}

	      			if (iflag(i,j,k).isSingleValued()){
				// Compute dphidn on the eb here 
              				grad_eb_arr(i,j,k,n) = grad_eb_of_phi_on_centroids_extdir(i, j, k, n, phi_arr, phi_eb_arr,
                        					iflag, ccent, bcent, vfrac, nx, ny, nz, is_eb_inhomog,
                        					on_x_face, domlo_x, domhi_x,
                        					on_y_face, domlo_y, domhi_y,
                        					on_z_face, domlo_z, domhi_z);

				}
			});

                   	 cns_eb_compute_dudt(BL_TO_FORTRAN_BOX(bx),
                                        BL_TO_FORTRAN_ANYD(dSdt[mfi]),
                                        BL_TO_FORTRAN_ANYD(S[mfi]),
					BL_TO_FORTRAN_ANYD(grad_eb[mfi]),
					BL_TO_FORTRAN_ANYD(grad_x[mfi]),
					BL_TO_FORTRAN_ANYD(grad_y[mfi]),
                                        BL_TO_FORTRAN_ANYD(flux[0]),
                                        BL_TO_FORTRAN_ANYD(flux[1]),
                                        BL_TO_FORTRAN_ANYD(flux[2]),
                                        BL_TO_FORTRAN_ANYD(flag),
                                        BL_TO_FORTRAN_ANYD((*volfrac)[mfi]),
                                        BL_TO_FORTRAN_ANYD((*bndrycent)[mfi]),
                                        BL_TO_FORTRAN_ANYD((*areafrac[0])[mfi]),
                                        BL_TO_FORTRAN_ANYD((*areafrac[1])[mfi]),
                                        BL_TO_FORTRAN_ANYD((*areafrac[2])[mfi]),
                                        BL_TO_FORTRAN_ANYD((*facecent[0])[mfi]),
                                        BL_TO_FORTRAN_ANYD((*facecent[1])[mfi]),
                                        BL_TO_FORTRAN_ANYD((*facecent[2])[mfi]),
                                        &as_crse,
                                        BL_TO_FORTRAN_ANYD(*p_drho_as_crse),
                                        BL_TO_FORTRAN_ANYD(*p_rrflag_as_crse),
                                        &as_fine,
                                        BL_TO_FORTRAN_ANYD(dm_as_fine),
                                        BL_TO_FORTRAN_ANYD(level_mask[mfi]),
                                        dx, &dt,&level);

                    if (fr_as_crse) {
                        fr_as_crse->CrseAdd(mfi, {&flux[0],&flux[1],&flux[2]}, dx,dt,
                                            (*volfrac)[mfi],
                                            {&((*areafrac[0])[mfi]),
                                             &((*areafrac[1])[mfi]),
                                             &((*areafrac[2])[mfi])},
                                            RunOn::Cpu);
                    }

                    if (fr_as_fine) {
                        fr_as_fine->FineAdd(mfi, {&flux[0],&flux[1],&flux[2]}, dx,dt,
                                            (*volfrac)[mfi],
                                            {&((*areafrac[0])[mfi]),
                                             &((*areafrac[1])[mfi]),
                                             &((*areafrac[2])[mfi])},
                                            dm_as_fine,
                                            RunOn::Cpu);
                    }
                }
            }

            wt = (amrex::second() - wt) / bx.d_numPts();
            cost[mfi].plus<RunOn::Host>(wt, bx);
        }
    }
}
