/*
  A very simple example of reading a plotfile and calling a function to perform a pointwise transformation
  based on a set of components that are specified by name on the command line.  The transformation is done
  in the accompanying fortran routine.  No grow cells are used, so the transformation cannot involve a
  stencil operation.

  The output is a new plotfile with a single component set to the output of the transform routine.  This
  new plotfile has metadata (number of levels, boxarray, grid spacing, etc) that is identical to the original
  plotfile.
 */

#include <Utils.H>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);
    if (argc < 3)
      print_usage(argc,argv);

    ParmParse pp;

    const std::string farg = amrex::get_command_argument(1);
    if (farg == "-h" || farg == "--help")
      print_usage(argc,argv);

    std::string infile; 
	pp.get("infile",infile);
	
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
	std::cout << nCompIn << "\n";
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

    const int nCompOut = nCompIn;
    const int nGrow = 0;
    const int nLev = amrData.FinestLevel() + 1;
	const int finest_lev = amrData.FinestLevel();
	const Vector<Real>& plo = amrData.ProbLo();
	const Vector<Real>& dx0  = amrData.DxLevel()[0];

	FILE* outfile;
	outfile = fopen("solution.csv","w");
	fprintf(outfile,"%s, %s, %s, %s, %s\n", "x", "y", "z", "x_velocity", "y_velocity");
	
    Vector<MultiFab> stateout;
	stateout.resize(nLev);
	for (int lev=0; lev <= finest_lev; ++lev) {
      	const BoxArray ba 	    = amrData.boxArray(lev);	
      	const DistributionMapping dmap(ba);

        stateout[lev].define(ba,dmap,nCompOut,0);
		MultiFab& stateout_mf = stateout[lev];
	
       // Load input data from pltfile
      	amrData.FillVar(stateout_mf,lev,varNames,destFillComps);

		for (MFIter mfi(stateout_mf); mfi.isValid(); ++mfi) {
        	const Array4<const Real>& stateout_array = stateout_mf.array(mfi);
        	const Box& box = mfi.validbox();

			ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
			{
				Vector<Real> centroid = get_coords(lev, i, j, k, dx0, plo);
				fprintf(outfile,"%0.15g, %0.15g, %0.15g, %0.15g, %0.15g\n", centroid[0], centroid[1], 1e-12, stateout_array(i,j,k,0), stateout_array(i,j,k,1) );
			});
		}
    }
  fclose(outfile);
  amrex::Finalize();
  return 0;
}
