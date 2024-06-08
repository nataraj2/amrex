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

class Quadruplet {
public:
    int a, b, c, d;

    Quadruplet(int w, int x, int y, int z) : a(w), b(x), c(y), d(z) {}

    // Equality operator to compare two Quadruplets
    bool operator==(const Quadruplet& other) const {
        return std::tie(a, b, c, d) == std::tie(other.a, other.b, other.c, other.d);
    }

    // Less than operator to use Quadruplet in sorted containers
    bool operator<(const Quadruplet& other) const {
        return std::tie(a, b, c, d) < std::tie(other.a, other.b, other.c, other.d);
    }

	int get_lev() {
		return a;
	}
	
	int get_i() {
		return b;
	}

	int get_j() {
		return c;
	}

	int get_k() {
		return d;
	}

};

// Define the QuadrupletStore class
class QuadrupletStore {
private:
    std::vector<Quadruplet> quadruplets;

public:
    // Function to add a triplet if it does not already exist
    bool addQuadruplet(int w, int x, int y, int z) {
        Quadruplet newQuadruplet(w, x, y, z);
        if (std::find(quadruplets.begin(), quadruplets.end(), newQuadruplet) == quadruplets.end()) {
            quadruplets.push_back(newQuadruplet);
			return true;	
        }
		else {
			return false;
		}
    }

	int get_size() {
		return quadruplets.size();
	}

	Quadruplet get_quadruplet(int i) {
		return this->quadruplets[i];
	} 
};

inline static
void
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=<plotfilename> varNames=v1 v2 ... \n";
  exit(1);
}

inline std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = amrex::Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

inline amrex::Vector<amrex::Real> get_node_coords (const int lev, const int i, const int j, const int k, 
									 const amrex::Vector<amrex::Real>& dx0, 
									 const amrex::Vector<amrex::Real>& plo)
{
	amrex::Vector<amrex::Real> temp;
    temp.push_back(plo[0] + i*dx0[0]/std::pow(2,lev));
    temp.push_back(plo[1] + j*dx0[1]/std::pow(2,lev));
    temp.push_back(plo[2] + k*dx0[2]/std::pow(2,lev));
    return temp;
}

inline amrex::Vector<amrex::Real> get_coords (const int lev, const int i, const int j, const int k, 
									 const amrex::Vector<amrex::Real>& dx0, 
									 const amrex::Vector<amrex::Real>& plo)
{
	amrex::Vector<amrex::Real> temp;
    temp.push_back(plo[0] + (i+0.5)*dx0[0]/std::pow(2,lev));
    temp.push_back(plo[1] + (j+0.5)*dx0[1]/std::pow(2,lev));
    temp.push_back(plo[2] + (k+0.5)*dx0[2]/std::pow(2,lev));
    return temp;
}



void WriteBoxesIntoVTK(const amrex::AmrData& amrData);

void CreateFineMask(const amrex::AmrData& amrData, 
					amrex::Vector<amrex::iMultiFab>& finemask);

void WriteFineMaskIntoVTK(const amrex::AmrData& amrData, 
						  const int lev, 
						  amrex::Vector<amrex::iMultiFab>& finemask);

void CreateAllMasks(const amrex::AmrData& amrData, 
					amrex::Vector<amrex::iMultiFab>& allmasks);

void WriteAllMasksIntoVTK(const amrex::AmrData& amrData, 
						  const int lev, 
						  amrex::Vector<amrex::iMultiFab>& allmasks);

void FindNeighbors(const int lev,
				   const int finest_lev,
                   const int i,
                   const int j,
                   const int k,
                   amrex::Array4<const int> const& allmasks_array,
                   amrex::Array4<const int> const& finemask_array,
                   const amrex::Vector<amrex::Real>& dx0,
                   const amrex::Vector<amrex::Real>& plo,
                   amrex::Vector<amrex::Vector<amrex::Real>>& vec_coords,
				   QuadrupletStore& quad);

void WriteGraphForAllLevels(const amrex::AmrData& amrData, 
										const amrex::Vector<amrex::iMultiFab>& finemask, 
										const amrex::Vector<amrex::iMultiFab>& allmasks,
										const amrex::Vector<amrex::MultiFab>& stateout,
										const std::string& file_point_str,
                                        const std::string& file_neighbors_str,
                                        const std::string& file_connect_str);

void WriteSolution(const int lev,
				   const amrex::AmrData& amrData,
                   const amrex::Vector<amrex::iMultiFab>& allmasks,
                   const amrex::Vector<amrex::MultiFab>& stateout,
                   const std::string& sol_file_string);
