#include <Advection_GNN.H>

using namespace amrex;

void CreateAllMasks(const AmrData& amrData, const int ng, Vector<iMultiFab>& allmasks)
{
	const int nLev = amrData.FinestLevel() + 1;
	
	 allmasks.resize(nLev);
    for (int lev = 0; lev < nLev; ++lev) {
        const BoxArray& ba = amrData.boxArray(lev);
        const DistributionMapping dmap(ba);
        allmasks[lev].define(ba,dmap,1,ng);
        allmasks[lev].setVal(1);
    }
    for (int lev = 0; lev < nLev; ++lev) {
        const Box& domain     = amrData.ProbDomain()[lev];
        allmasks[lev].BuildMask(domain,Periodicity::NonPeriodic(),0,1,2,0);
    }
}

void WriteAllMasksIntoVTK(const AmrData& amrData, const int lev, Vector<iMultiFab>& allmasks)
{
	const Vector<Real>& plo = amrData.ProbLo();
    const Vector<Real>& dx0  = amrData.DxLevel()[0];
    const int nLev = amrData.FinestLevel() + 1;

	if(lev > nLev-1){
		std::cout << "There are only " << nLev << " levels and hence the maximum level specified for WriteAllMasksIntoVTK " << 
				     "can only be " << nLev-1 << ". Exiting....""\n";
		exit(0);
	}

	FILE *allmasks_cf_vtk, *mask_physbnd_vtk;
    allmasks_cf_vtk = fopen("allmasks_cf.vtk","w");
    fprintf(allmasks_cf_vtk, "%s\n","# vtk DataFile Version 3.0");
    fprintf(allmasks_cf_vtk, "%s\n","All mask coarse-fine data");
    fprintf(allmasks_cf_vtk, "%s\n","ASCII");
    fprintf(allmasks_cf_vtk, "%s\n","DATASET POLYDATA");
    fprintf(allmasks_cf_vtk, "%s %ld %s\n", "POINTS", 0, "float");

    mask_physbnd_vtk = fopen("mask_physbnd.vtk","w");
    fprintf(mask_physbnd_vtk, "%s\n","# vtk DataFile Version 3.0");
    fprintf(mask_physbnd_vtk, "%s\n","Mask phys bnd data");
    fprintf(mask_physbnd_vtk, "%s\n","ASCII");
    fprintf(mask_physbnd_vtk, "%s\n","DATASET POLYDATA");
    fprintf(mask_physbnd_vtk, "%s %ld %s\n", "POINTS", 0, "float");

        //iMultiFab& allmasks_mf = allmasks[lev];
        iMultiFab& allmasks_mf = allmasks[lev];
        for (MFIter mfi(allmasks_mf, true); mfi.isValid(); ++mfi) {
			const int ng = allmasks_mf.nGrow();
            Array4<int> const& allmasks_array = allmasks_mf.array(mfi);
            Box bx = mfi.growntilebox(ng);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
                if(allmasks_array(i,j,k,0) == 1 and k==0){
                    std::vector<Real> coords = get_coords(lev, i, j, k, dx0, plo);
                    fprintf(allmasks_cf_vtk, "%0.15g %0.15g %0.15g\n", coords[0], coords[1], coords[2]);
                }
				 if(allmasks_array(i,j,k,0) == 2 and k==0){
                    std::vector<Real> coords = get_coords(lev, i, j, k, dx0, plo);
                    fprintf(mask_physbnd_vtk, "%0.15g %0.15g %0.15g\n", coords[0], coords[1], coords[2]);
                }

            });
        }
        fclose(allmasks_cf_vtk);
}



