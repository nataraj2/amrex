//#include <AMReX_ParallelDescriptor.H>

inline amrex::Vector<amrex::Real> get_coords(const int lev, const int i, const int j, const int k, 
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
void CreateFineMask(const amrex::AmrData& amrData, amrex::Vector<amrex::iMultiFab>& finemask);
void WriteFineMaskIntoVTK(const amrex::AmrData& amrData, const int lev, amrex::Vector<amrex::iMultiFab>& finemask);
void CreateAllMasks(const amrex::AmrData& amrData, amrex::Vector<amrex::iMultiFab>& allmasks);
void WriteAllMasksIntoVTK(const amrex::AmrData& amrData, const int lev, amrex::Vector<amrex::iMultiFab>& allmasks);
