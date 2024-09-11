#include <string>
#include <iostream>

#include <AMReX.H>
#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <AMReX_WritePlotFile.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFabUtil.H>

#include <AMReX_BLFort.H>

inline static
void
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=<plotfilename> varNames=v1 v2 ng=<number of layers of cell neighbors for each point for graph> ... \n";
  exit(1);
}

inline std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = amrex::Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

inline amrex::Vector<amrex::Real> get_coords (const int lev, const int i, const int j, const int k, 
									 const amrex::Vector<amrex::Real>& dx0, 
									 const amrex::Vector<amrex::Real>& plo)
{
	amrex::Vector<amrex::Real> temp;
    temp.push_back(plo[0] + (i+0.5)*dx0[0]/std::pow(2,lev));
    temp.push_back(plo[1] + (j+0.5)*dx0[1]/std::pow(2,lev));
    //temp.push_back(plo[2] + (k+0.5)*dx0[2]/std::pow(2,lev));
    return temp;
}
