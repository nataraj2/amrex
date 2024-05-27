#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>

using namespace amrex;

void WriteBoxesIntoVTK(const AmrData& amrData)
{

	const int nLev = amrData.FinestLevel() + 1;
	const Vector<Real>& plo = amrData.ProbLo();
	std::vector<Real> xcoord, ycoord, zcoord;
	xcoord.resize(8);
    ycoord.resize(8);
    zcoord.resize(8);

    for (int lev=0; lev<nLev; ++lev) {
        FILE* ba_in_lev_vtk;
        std::string filename = "ba_" + std::to_string(lev) + ".vtk";
        ba_in_lev_vtk = fopen(filename.c_str(),"w");

        //const Box  domain       = amrData.ProbDomain()[lev];
        const BoxArray ba       = amrData.boxArray(lev);
        const Vector<Real>& dx  = amrData.DxLevel()[lev];

        const DistributionMapping dmap(ba);

        fprintf(ba_in_lev_vtk,"%s\n","# vtk DataFile Version 3.0");
        fprintf(ba_in_lev_vtk,"%s\n", "Box layout");
        fprintf(ba_in_lev_vtk,"%s\n","ASCII");
        fprintf(ba_in_lev_vtk,"%s\n","DATASET POLYDATA");
        fprintf(ba_in_lev_vtk,"%s %ld %s\n","POINTS", ba.size()*8, "float");

        std::cout << "Level is " << lev << "\n";
        for (int i = 0; i < ba.size(); ++i) {
            amrex::Box bx = ba[i];
            const int* lo = bx.loVect();
            const int* hi = bx.hiVect();

            xcoord[0] = lo[0]*dx[0];     ycoord[0] = lo[1]*dx[1];     zcoord[0] = lo[2]*dx[2];
            xcoord[1] = (hi[0]+1)*dx[0]; ycoord[1] = lo[1]*dx[1];     zcoord[1] = lo[2]*dx[2];
            xcoord[2] = (hi[0]+1)*dx[0]; ycoord[2] = (hi[1]+1)*dx[1]; zcoord[2] = lo[2]*dx[2];
            xcoord[3] = lo[0]*dx[0];     ycoord[3] = (hi[1]+1)*dx[1]; zcoord[3] = lo[2]*dx[2];
            xcoord[4] = lo[0]*dx[0];     ycoord[4] = lo[1]*dx[1];     zcoord[4] = (hi[2]+1)*dx[2];
            xcoord[5] = (hi[0]+1)*dx[0]; ycoord[5] = lo[1]*dx[1];     zcoord[5] = (hi[2]+1)*dx[2];
            xcoord[6] = (hi[0]+1)*dx[0]; ycoord[6] = (hi[1]+1)*dx[1]; zcoord[6] = (hi[2]+1)*dx[2];
            xcoord[7] = lo[0]*dx[0];     ycoord[7] = (hi[1]+1)*dx[1]; zcoord[7] = (hi[2]+1)*dx[2];

            for(int ipt=0;ipt<8;ipt++){
                xcoord[ipt] = plo[0] + xcoord[ipt];
                ycoord[ipt] = plo[1] + ycoord[ipt];
                zcoord[ipt] = plo[2] + zcoord[ipt];
                fprintf(ba_in_lev_vtk,"%0.15g %0.15g %0.15g\n", xcoord[ipt], ycoord[ipt], zcoord[ipt]);
            }

        }

        fprintf(ba_in_lev_vtk,"%s %ld %ld\n", "LINES", ba.size()*12, ba.size()*12*3);
        for(int ibx=0;ibx<ba.size();ibx++){
            fprintf(ba_in_lev_vtk, "%d %ld %ld\n", 2, static_cast<long int> (ibx*8),     static_cast<long int> (ibx*8 + 1));
            fprintf(ba_in_lev_vtk, "%d %ld %ld\n", 2, static_cast<long int> (ibx*8 + 1), static_cast<long int> (ibx*8 + 2));
            fprintf(ba_in_lev_vtk, "%d %ld %ld\n", 2, static_cast<long int> (ibx*8 + 2), static_cast<long int> (ibx*8 + 3));
            fprintf(ba_in_lev_vtk, "%d %ld %ld\n", 2, static_cast<long int> (ibx*8 + 3), static_cast<long int> (ibx*8     ));
            fprintf(ba_in_lev_vtk, "%d %ld %ld\n", 2, static_cast<long int> (ibx*8 + 4), static_cast<long int> (ibx*8 + 5));
            fprintf(ba_in_lev_vtk, "%d %ld %ld\n", 2, static_cast<long int> (ibx*8 + 5), static_cast<long int> (ibx*8 + 6));
            fprintf(ba_in_lev_vtk, "%d %ld %ld\n", 2, static_cast<long int> (ibx*8 + 6), static_cast<long int> (ibx*8 + 7));
            fprintf(ba_in_lev_vtk, "%d %ld %ld\n", 2, static_cast<long int> (ibx*8 + 7), static_cast<long int> (ibx*8 + 4));
            fprintf(ba_in_lev_vtk, "%d %ld %ld\n", 2, static_cast<long int> (ibx*8    ), static_cast<long int> (ibx*8 + 4));
            fprintf(ba_in_lev_vtk, "%d %ld %ld\n", 2, static_cast<long int> (ibx*8 + 1), static_cast<long int> (ibx*8 + 5));
            fprintf(ba_in_lev_vtk, "%d %ld %ld\n", 2, static_cast<long int> (ibx*8 + 2), static_cast<long int> (ibx*8 + 6));
            fprintf(ba_in_lev_vtk, "%d %ld %ld\n", 2, static_cast<long int> (ibx*8 + 3), static_cast<long int> (ibx*8 + 7));
        }

        fclose(ba_in_lev_vtk);
	}
}
