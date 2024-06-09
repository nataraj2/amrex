#include <Advection_GNN.H>

using namespace amrex;

void FindNeighbors(const int lev,
				   const int finest_lev,
				   const int i,
				   const int j,
				   const int k,
				   Array4<const int> const& allmasks_array,
				   Array4<const int> const& finemask_array,
				   const Vector<Real>& dx0,
				   const Vector<Real>& plo,
				   Vector<Vector<Real>>& vec_coords,
				   QuadrupletStore& quad)
{

	if(lev < finest_lev) {
		if(finemask_array(i,j,k,0) == 0) {
			if(allmasks_array(i,j,k,0) == 0) {
				bool isnew = quad.addQuadruplet(lev,i,j,k);
				if(isnew) {
       				Vector<Real> coords = get_coords(lev, i, j, k, dx0, plo);
        			vec_coords.push_back(coords);
				}
        	}
        	if(allmasks_array(i,j,k,0) == 1) {
				bool isnew = quad.addQuadruplet(lev-1,i/2,j/2,k);
				if(isnew) {
        			Vector<Real> coords = get_coords(lev-1, i/2, j/2, k, dx0, plo);
            		vec_coords.push_back(coords);
				}
        	}
    	} else if(finemask_array(i,j,k,0) == 1) {
			Vector<Real> coords;
			bool isnew = quad.addQuadruplet(lev+1, 2*i, 2*j, k);
			if(isnew) {
    			coords = get_coords(lev+1, 2*i, 2*j, k, dx0, plo);
       			vec_coords.push_back(coords);
			}
			isnew = quad.addQuadruplet(lev+1, 2*i+1, 2*j, k);
			if(isnew) {
    			coords = get_coords(lev+1, 2*i+1, 2*j, k, dx0, plo);
       			vec_coords.push_back(coords);
			}
			isnew = quad.addQuadruplet(lev+1, 2*i, 2*j+1, k);
			if(isnew) {
       			coords = get_coords(lev+1, 2*i, 2*j+1, k, dx0, plo);
       			vec_coords.push_back(coords);
			}
			isnew = quad.addQuadruplet(lev+1, 2*i+1, 2*j+1, k);
			if(isnew) {
       			coords = get_coords(lev+1, 2*i+1, 2*j+1, k, dx0, plo);
       			vec_coords.push_back(coords);
			}

    	}
	} else if(lev == finest_lev) {

		if(allmasks_array(i,j,k,0) == 0) {
			bool isnew = quad.addQuadruplet(lev,i,j,k);
			if(isnew) {
    			Vector<Real> coords = get_coords(lev, i, j, k, dx0, plo);
        		vec_coords.push_back(coords);
			}
    	}
    	if(allmasks_array(i,j,k,0) == 1)
    	{
			bool isnew = quad.addQuadruplet(lev-1,i/2,j/2,k);
			if(isnew){
    			Vector<Real> coords = get_coords(lev-1, i/2, j/2, k, dx0, plo);
    			vec_coords.push_back(coords);
			}
    	}
	}
}
	
void FindNeighborsFinestLev(const int lev,
				   const int i,
				   const int j,
				   const int k,
				   Array4<const int> const& allmasks_array,
				   const Vector<Real>& dx0,
				   const Vector<Real>& plo,
				   Vector<Vector<Real>>& vec_coords,
				   QuadrupletStore& quad)
{
	if(allmasks_array(i,j,k,0) == 0)
    {
		bool isnew = quad.addQuadruplet(lev,i,j,k);
		if(isnew) {
    		Vector<Real> coords = get_coords(lev, i, j, k, dx0, plo);
        	vec_coords.push_back(coords);
		}
    }
    if(allmasks_array(i,j,k,0) == 1)
    {
		bool isnew = quad.addQuadruplet(lev-1,i/2,j/2,k);
		if(isnew){
    		Vector<Real> coords = get_coords(lev-1, i/2, j/2, k, dx0, plo);
    		vec_coords.push_back(coords);
		}
    }
}


			   
