#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>

#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFabUtil.H>

#include <Utils_Adv.H>

using namespace amrex;

void WriteFineMaskIntoVTK(const AmrData& amrData, const int lev, Vector<iMultiFab>& finemask)
{

	const Vector<Real>& plo = amrData.ProbLo();
    const Vector<Real>& dx0  = amrData.DxLevel()[0];
	const int nLev = amrData.FinestLevel() + 1;

	FILE* finemask_vtk;
    finemask_vtk = fopen("finemask.vtk","w");
    fprintf(finemask_vtk, "%s\n","# vtk DataFile Version 3.0");
    fprintf(finemask_vtk, "%s\n","Fine mask data");
    fprintf(finemask_vtk, "%s\n","ASCII");
    fprintf(finemask_vtk, "%s\n","DATASET POLYDATA");
    fprintf(finemask_vtk, "%s %ld %s\n", "POINTS", 0, "float");

    if(lev > nLev-2){
       std::cout << "The number of levels in the grid is " << nLev << " and hence there is fine mask only till level " << nLev - 2 << "\n"; 
	   std::cout << "The level specified is " << lev << "Exiting....." << "\n";
       exit(0);
     }

    iMultiFab& finemask_mf = finemask[lev];
	for (MFIter mfi(finemask_mf); mfi.isValid(); ++mfi) {
     	Array4<int> const& finemask_array = finemask_mf.array(mfi);
     	Box bx = mfi.validbox();
     	ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
      		if(finemask_array(i,j,k,0) == 1 and k==1){
      			std::vector<Real> coords = get_coords(lev, i, j, k, dx0, plo);
      			fprintf(finemask_vtk, "%0.15g %0.15g %0.15g\n", coords[0], coords[1], coords[2]);
     		}
      	});
	}
    fclose(finemask_vtk);
}
