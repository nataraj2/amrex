/*
  A very simple example of reading a plotfile and calling a function to perform a pointwise transformation
  based on a set of components that are specified by name on the command line.  The transformation is done
  in the accompanying fortran routine.  No grow cells are used, so the transformation cannot involve a
  stencil operation.

  The output is a new plotfile with a single component set to the output of the transform routine.  This
  new plotfile has metadata (number of levels, boxarray, grid spacing, etc) that is identical to the original
  plotfile.
 */
#include <string>
#include <iostream>

#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <AMReX_WritePlotFile.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFabUtil.H>

#include <AMReX_BLFort.H>

#include <Advection_GNN.H>

extern "C" {
  void transform (const int* lo, const int* hi,
                  const amrex_real* sIn, const int* sInlo, const int* sInhi, const int* ncIn,
                  amrex_real* sOut, const int* sOutlo, const int* sOuthi, const int* ncOut);
}

using namespace amrex;

static
void
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=<plotfilename> varNames=v1 v2 ... \n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);

    ParmParse pp;

    const std::string farg = amrex::get_command_argument(1);
    if (farg == "-h" || farg == "--help")
      print_usage(argc,argv);

    std::string infile; pp.get("infile",infile);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int nv = pp.countval("varNames");
    Vector<std::string> varNames(nv); pp.getarr("varNames",varNames,0,nv);

    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    int nCompIn = varNames.size();
    Vector<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i) {
      destFillComps[i] = i;
    }


    for (int i=0; i<nCompIn; ++i) {
      int ivar = -1;
      for (int j=0; j<plotVarNames.size(); ++j) {
        if (plotVarNames[j] == varNames[i]) {ivar = j;}
      }
      if (ParallelDescriptor::IOProcessor() && ivar<0) {
        Abort("Cannot find variable="+varNames[i]+" in pltfile");
      }
    }

    const int nCompOut = 1;
    const int nGrow = 0;
    const int nLev = amrData.FinestLevel() + 1;
	const int finest_lev = amrData.FinestLevel();


	std::vector<Real> xcoord, ycoord, zcoord;
	xcoord.resize(8);
	ycoord.resize(8);
	zcoord.resize(8);

	const Vector<Real>& plo = amrData.ProbLo();
	const Vector<Real>& dx0  = amrData.DxLevel()[0];
	//const Vector<Real>& phi = amrData.ProbHi();
	
	for (int lev=0; lev<nLev; ++lev) {
      	const BoxArray ba 	    = amrData.boxArray(lev);
		for (int ibx = 0; ibx < ba.size(); ++ibx) {
           	amrex::Box bx = ba[ibx];
           	const int* lo = bx.loVect();
            const int* hi = bx.hiVect();	
			ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
				if(i==lo[0] and k==0){
					std::cout << "lev, i, j, k " << lev << " " << i << " " << j << " " << k << "\n";
				}
			});
		}
	}

	// Write the boxes as VTK for visualization
	WriteBoxesIntoVTK(amrData);

	// Create a finemask on all coarse levels ie. all levels 
	//except the finest level
	Vector<iMultiFab> finemask;
	if(nLev > 1){	
		CreateFineMask(amrData, finemask);
		WriteFineMaskIntoVTK(amrData, 0, finemask);
	}

	// Build mask to identify cells at coarse-fine interface
	Vector<iMultiFab> allmasks;
	CreateAllMasks(amrData, allmasks);
	WriteAllMasksIntoVTK(amrData, 1, allmasks);

    Vector<MultiFab> stateout;
	stateout.resize(nLev);
	for (int lev=0; lev <= finest_lev; ++lev) {
      	const BoxArray ba 	    = amrData.boxArray(lev);	
      	const DistributionMapping dmap(ba);

        stateout[lev].define(ba,dmap,1,0);
		MultiFab& stateout_mf = stateout[lev];
	
       // Load input data from pltfile
      	amrData.FillVar(stateout_mf,lev,varNames,destFillComps);

		for (MFIter mfi(stateout_mf); mfi.isValid(); ++mfi) {
        	const Array4<const Real>& stateout_array = stateout_mf.array(mfi);
        	const Box& box = mfi.validbox();

			ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
			{
				//std::cout << stateout_array(i,j,k) << "\n";
			});
		}
    }


	// Trying to write out a multi level VTK file 
	// using the unstructured vtk format. Not successful

	/*FILE* solution;
    std::string file_solution = "solution.vtk";
    solution = fopen(file_solution.c_str(),"w");

    fprintf(solution,"%s\n","# vtk DataFile Version 3.0");
    fprintf(solution,"%s\n", "Multi level solution");
    fprintf(solution,"%s\n","ASCII");
    fprintf(solution,"%s\n","DATASET UNSTRUCTURED_GRID");	

	WriteSolution(amrData, allmasks, stateout, solution);
	exit(0);*/


	FILE* connect;
   	std::string file_connect = "connect.vtk";
    connect = fopen(file_connect.c_str(),"w");

    fprintf(connect,"%s\n","# vtk DataFile Version 3.0");
    fprintf(connect,"%s\n", "Connectivity");
    fprintf(connect,"%s\n","ASCII");
    fprintf(connect,"%s\n","DATASET POLYDATA");


	//WriteGraphForAllLevelsExceptFinest(amrData, finemask, allmasks, stateout, connect);
	if(finest_lev > 0){
		WriteGraphForFinestLev(amrData, allmasks, stateout, connect);
	}
	
	std::cout  << "Reaching here .... " << "\n";
		//exit(0);
    // Write result to new plotfile in local folder
    /*std::string outfile=getFileRoot(infile) + "_tr";
    Vector<std::string> outNames;
    outNames.push_back("transform");
    WritePlotFile(stateout,amrData,outfile,false,outNames);

    for (int lev=0; lev<nLev; ++lev)
    {
        delete stateout[lev];
    }*/
  }
  amrex::Finalize();
  return 0;
}
