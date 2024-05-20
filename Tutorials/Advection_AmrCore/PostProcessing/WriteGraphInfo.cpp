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

#include <AMReX_BLFort.H>

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

Vector<Real> get_coord(const int lev, const int i, const int j, const int k, const Vector<Real>& dx0, const Vector<Real>& plo)
{
	Vector<Real> temp;
	temp.push_back(plo[0] + (i+0.5)*dx0[0]/std::pow(2,lev));
	temp.push_back(plo[1] + (j+0.5)*dx0[1]/std::pow(2,lev));
	temp.push_back(plo[2] + (k+0.5)*dx0[2]/std::pow(2,lev));
	return temp;
}

void WriteBoxesIntoVTK(AmrData& amrData);

//Real getcoord(int i, int j, int k)

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

    Vector<MultiFab*> stateOut(nLev);

	std::vector<Real> xcoord, ycoord, zcoord;
	xcoord.resize(8);
	ycoord.resize(8);
	zcoord.resize(8);

	const Vector<Real>& plo = amrData.ProbLo();
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

	// CHeck cell
	int ilev = 2;
	int i_ind = 200;
	int j_ind = 211;
	int k_ind = 0;

	WriteBoxesIntoVTK(amrData);
	
	const Vector<Real>& dx0  = amrData.DxLevel()[0];

	for (int lev=0; lev<nLev; ++lev) {
	  	//const Box  domain       = amrData.ProbDomain()[lev];
      	const BoxArray ba 	    = amrData.boxArray(lev);
		const Vector<Real>& dx  = amrData.DxLevel()[lev];
      	
      	const DistributionMapping dmap(ba);
      	MultiFab stateIn(ba,dmap,nCompIn,nGrow);
      	stateOut[lev] = new MultiFab(ba,dmap,nCompOut,0);


		FILE* connect;
        std::string file_connect = "connect.vtk";
        connect = fopen(file_connect.c_str(),"w");

        fprintf(connect,"%s\n","# vtk DataFile Version 3.0");
        fprintf(connect,"%s\n", "Connectivity");
        fprintf(connect,"%s\n","ASCII");
        fprintf(connect,"%s\n","DATASET POLYDATA");

	  // Create graph
		for (int ibx = 0; ibx < ba.size(); ++ibx) {
            amrex::Box bx = ba[ibx];
            const int* lo = bx.loVect();
            const int* hi = bx.hiVect();	
			ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
				std::vector<std::vector<double>> vec_coords;
				std::vector<Real> coords1 = get_coord(lev, i, j, k, dx0, plo);
				vec_coords.push_back(coords1);
				if(i > lo[0] and i < hi[0])
				{
					std::vector<Real> coords = get_coord(lev, i+1, j, k, dx0, plo);
					vec_coords.push_back(coords);
					coords = get_coord(lev, i-1, j, k, dx, plo);
					vec_coords.push_back(coords);
				}
				if(j > lo[1] and j < hi[1])
				{
					std::vector<Real> coords = get_coord(lev, i, j+1, k, dx0, plo);
					vec_coords.push_back(coords);
					coords = get_coord(lev, i, j-1, k, dx, plo);
					vec_coords.push_back(coords);
				}
				if(k > lo[2] and k < hi[2])
				{
					std::vector<Real> coords = get_coord(lev, i, j, k+1, dx0, plo);
					vec_coords.push_back(coords);
					coords = get_coord(lev, i, j, k-1, dx, plo);
					vec_coords.push_back(coords);
				}

				if(lev > 0 and (i==lo[0] or i==hi[0] or j==lo[1] or j==hi[1] or k==lo[2] or k==hi[2])){
					for (int ibx1 = 0; ibx1 < ba.size(); ++ibx1) {
						 amrex::Box bx = ba[ibx];
            			 const int* lo = bx.loVect();
            			 const int* hi = bx.hiVect();	
						if(hi[0] == lo[0]-1){
							
						}
					}
				}	
				
				if(lev == ilev and i == i_ind and j == j_ind and k == k_ind){
					fprintf(connect,"%s %ld %s\n","POINTS",vec_coords.size(), "float");
					for(long unsigned int ipt=0;ipt<vec_coords.size();ipt++)
					{
						fprintf(connect, "%0.15g %0.15g %0.15g\n", vec_coords[ipt][0], vec_coords[ipt][1], vec_coords[ipt][2]);		
					}
					fprintf(connect,"%s %ld %ld\n", "LINES",vec_coords.size()-1, (vec_coords.size()-1)*3);
					for(long unsigned int ipt=1;ipt<vec_coords.size();ipt++){	
						fprintf(connect,"%d %d %ld\n", 2, 0, ipt);
					}
					exit(0);
				}
			});	
		}
	
		fclose(connect);


      // Load input data from pltfile
      amrData.FillVar(stateIn,lev,varNames,destFillComps);

		/*for (MFIter mfi(stateIn); mfi.isValid(); ++mfi) {
        	const Array4<const Real>& sIn = stateIn.array(mfi);
        	const Box& box = mfi.validbox();

			ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
			{
				std::cout << sIn(i,j,k) << "\n";
			});
		}*/

		//exit(0);


      // Compute transformation
      for (MFIter mfi(stateIn); mfi.isValid(); ++mfi) {
        const FArrayBox& sIn = stateIn[mfi];
        FArrayBox& sOut = (*stateOut[lev])[mfi];
        const Box& box = mfi.validbox();

        transform(BL_TO_FORTRAN_BOX(box),
                  BL_TO_FORTRAN_ANYD(sIn),&nCompIn,
                  BL_TO_FORTRAN_ANYD(sOut),&nCompOut);

      }
    }

    // Write result to new plotfile in local folder
    std::string outfile=getFileRoot(infile) + "_tr";
    Vector<std::string> outNames;
    outNames.push_back("transform");
    WritePlotFile(stateOut,amrData,outfile,false,outNames);

    for (int lev=0; lev<nLev; ++lev)
    {
        delete stateOut[lev];
    }
  }
  amrex::Finalize();
  return 0;
}
