#include <Advection_GNN.H>

using namespace amrex;


void CreateFineMask(const AmrData& amrData, Vector<iMultiFab>& finemask)
{
	const int nLev = amrData.FinestLevel() + 1;
	const int finest_lev = amrData.FinestLevel();

	Vector<MultiFab> phi;
    phi.resize(nLev);
    for (int lev = 0; lev < nLev; ++lev) {
        const BoxArray ba = amrData.boxArray(lev);
        const DistributionMapping dmap(ba);
        phi[lev].define(ba,dmap,1,0);
    }

    finemask.resize(nLev);
	for (int lev = 0; lev < nLev; ++lev) {
        const BoxArray ba = amrData.boxArray(lev);
        const DistributionMapping dmap(ba);
        finemask[lev].define(ba,dmap,1,3);
		finemask[lev].setVal(-100);
    }


    for (int lev = 0; lev < nLev-1; ++lev) {
        const IntVect ratio{2};

        finemask[lev] = makeFineMask(phi[lev], phi[lev+1], IntVect(3),
                                          ratio,Periodicity::NonPeriodic(),
                                          0, 1);
    }
	
	finemask[finest_lev].setVal(-100);

	/*Vector<Vector<BoxArray> >  grids;
    Vector<Vector<DistributionMapping>> dmapvec;
    grids.resize(nLev);
    dmapvec.resize(nLev);

    for (int lev = 0; lev < nLev; ++lev) {
        grids[lev].push_back(amrData.boxArray(lev));
        const DistributionMapping dmap(amrData.boxArray(lev));
        dmapvec[lev].push_back(dmap);
    }

	for (int lev = 0; lev < nLev-1; ++lev) {
        const IntVect ratio{2};
        const BoxArray ba = amrData.boxArray(lev);
        const DistributionMapping dmap(ba);

        finemask[lev] = makeFineMask(grids[lev][0], dmapvec[lev][0],
                          grids[lev+1][0],
                          ratio, 1, 0);
    }*/
}

void WriteFineMaskIntoVTK(const AmrData& amrData, const int lev, Vector<iMultiFab>& finemask)
{

    const Vector<Real>& plo = amrData.ProbLo();
    const Vector<Real>& dx0  = amrData.DxLevel()[0];
    const int finest_lev = amrData.FinestLevel();

    FILE* finemask_vtk;
    finemask_vtk = fopen("finemask.vtk","w");
    fprintf(finemask_vtk, "%s\n","# vtk DataFile Version 3.0");
    fprintf(finemask_vtk, "%s\n","Fine mask data");
    fprintf(finemask_vtk, "%s\n","ASCII");
    fprintf(finemask_vtk, "%s\n","DATASET POLYDATA");
    fprintf(finemask_vtk, "%s %ld %s\n", "POINTS", 0, "float");

    if(lev > finest_lev){
       std::cout << "There are only " << finest_lev+1 << " levels and hence the maximum level specified for WriteFineMaskIntoVTK " << 
					"can only be " << finest_lev << " but the level specified is " << lev << ". Exiting....." << "\n";
       exit(0);
     }

    iMultiFab& finemask_mf = finemask[lev];
    for (MFIter mfi(finemask_mf); mfi.isValid(); ++mfi) {
        Array4<int> const& finemask_array = finemask_mf.array(mfi);
        Box bx = mfi.growntilebox(1);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
            if(finemask_array(i,j,k,0) == 1 and k==1){
                std::vector<Real> coords = get_coords(lev, i, j, k, dx0, plo);
                fprintf(finemask_vtk, "%0.15g %0.15g %0.15g\n", coords[0], coords[1], coords[2]);
            }
        });
    }
    fclose(finemask_vtk);
}
