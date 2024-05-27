#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

void CreateFineMask(const AmrData& amrData, Vector<iMultiFab>& finemask)
{
	const int nLev = amrData.FinestLevel() + 1;

	Vector<MultiFab> phi;
    phi.resize(nLev);
    for (int lev = 0; lev < nLev; ++lev) {
        const BoxArray ba = amrData.boxArray(lev);
        const DistributionMapping dmap(ba);
        phi[lev].define(ba,dmap,1,0);
    }

    finemask.resize(nLev-1);
    for (int lev = 0; lev < nLev-1; ++lev) {
        const IntVect ratio{2};

        finemask[lev] = makeFineMask(phi[lev], phi[lev+1], IntVect(0),
                                          ratio,Periodicity::NonPeriodic(),
                                          1, 0);
    }

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
