#include <Advection_GNN.H>

using namespace amrex;	
	
void WriteGraphForAllLevels(const AmrData& amrData, 
										const Vector<iMultiFab>& finemask, 
										const Vector<iMultiFab>& allmasks, 
										const Vector<MultiFab>& stateout,
										const std::string& file_point_str,
										const std::string& file_neighbors_str,
										const std::string& file_connect_str)
{	
    FILE* file_point = fopen(file_point_str.c_str(),"w");

    fprintf(file_point,"%s\n","# vtk DataFile Version 3.0");
    fprintf(file_point,"%s\n", "Connectivity");
    fprintf(file_point,"%s\n","ASCII");
    fprintf(file_point,"%s\n","DATASET POLYDATA");

    FILE* file_neighbors = fopen(file_neighbors_str.c_str(),"w");

    fprintf(file_neighbors,"%s\n","# vtk DataFile Version 3.0");
    fprintf(file_neighbors,"%s\n", "Neighbors");
    fprintf(file_neighbors,"%s\n","ASCII");
    fprintf(file_neighbors,"%s\n","DATASET POLYDATA");

	 FILE* file_connect = fopen(file_connect_str.c_str(),"w");

    fprintf(file_connect,"%s\n","# vtk DataFile Version 3.0");
    fprintf(file_connect,"%s\n", "Neighbors");
    fprintf(file_connect,"%s\n","ASCII");
    fprintf(file_connect,"%s\n","DATASET POLYDATA");

	const Vector<Real>& plo = amrData.ProbLo();
    const Vector<Real>& dx0  = amrData.DxLevel()[0];
    int finest_lev = amrData.FinestLevel();

	for (int lev=0; lev <= finest_lev; ++lev) {
	
    	const iMultiFab& allmasks_mf = allmasks[lev];	
		const iMultiFab& finemask_mf = finemask[lev];
	
	    for (MFIter mfi(finemask_mf); mfi.isValid(); ++mfi) {
	    	Array4<const int> const& allmasks_array = allmasks_mf.const_array(mfi);
	    	Array4<const int> const& finemask_array = finemask_mf.const_array(mfi);
	        Box bx = mfi.validbox();
	        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
	       		if(k == 0 and allmasks_array(i,j,k,0) == 0) {
	            	Vector<Vector<Real>> vec_coords;
	            	Vector<Real> coords = get_coords(lev, i, j, k, dx0, plo);
	                vec_coords.push_back(coords);
	
	                QuadrupletStore quad;
	                for(int i_shift=-3; i_shift<=3; i_shift++){
	                	for(int j_shift=-3; j_shift<=3; j_shift++){
	                		if(!(i_shift == 0  and j_shift == 0)){
	                        	FindNeighbors(lev, finest_lev, i+i_shift, j+j_shift, k, 
											  allmasks_array, finemask_array, dx0, plo, 
											  vec_coords, quad);
	                        }
	                    }
	                }
	
	                AMREX_ALWAYS_ASSERT(quad.get_size() == vec_coords.size()-1);
	
	                if(lev == 2 and allmasks_array(i+1,j,k,0) == 1 and allmasks_array(i,j+1,k,0) == 1 ){
	                	fprintf(file_point,"%s %ld %s\n","POINTS", static_cast<long int>(1), "float");
	                	fprintf(file_connect,"%s %ld %s\n","POINTS",vec_coords.size(), "float");

	                    fprintf(file_point, "%0.15g %0.15g %0.15g\n", vec_coords[0][0], vec_coords[0][1], vec_coords[0][2]);
	                    fprintf(file_connect, "%0.15g %0.15g %0.15g\n", vec_coords[0][0], vec_coords[0][1], vec_coords[0][2]);

	                	fprintf(file_neighbors,"%s %ld %s\n","POINTS",vec_coords.size()-1, "float");
	                	for(int ipt=1;ipt<vec_coords.size();ipt++)
	                	{
	                    	fprintf(file_neighbors, "%0.15g %0.15g %0.15g\n", vec_coords[ipt][0], vec_coords[ipt][1], vec_coords[ipt][2]);
	                    	fprintf(file_connect, "%0.15g %0.15g %0.15g\n", vec_coords[ipt][0], vec_coords[ipt][1], vec_coords[ipt][2]);
	                	}
	                    fprintf(file_connect,"%s %ld %ld\n", "LINES",vec_coords.size()-1, (vec_coords.size()-1)*3);
	                    for(long unsigned int ipt=1;ipt<vec_coords.size();ipt++){
	                    	fprintf(file_connect,"%d %d %ld\n", 2, 0, ipt);
	                    }
	                    fclose(file_point);
	                    fclose(file_neighbors);
	                    fclose(file_connect);
	                    exit(0);
	                }
	            }
	        });
	    }
	}
}

void WriteSolution(const int lev, 
				   const AmrData& amrData,
				   const Vector<iMultiFab>& allmasks,  
                   const Vector<MultiFab>& stateout,
                   const std::string& sol_file_string)

{

    FILE* solution = fopen(sol_file_string.c_str(),"w");

    fprintf(solution,"%s\n","# vtk DataFile Version 3.0");
    fprintf(solution,"%s\n", "Multi level solution");
    fprintf(solution,"%s\n","ASCII");
    fprintf(solution,"%s\n","DATASET UNSTRUCTURED_GRID");

	const Vector<Real>& plo = amrData.ProbLo();
    const Vector<Real>& dx0  = amrData.DxLevel()[0];

	int cell_count = 0;

	QuadrupletStore quad;


	 const MultiFab& stateout_mf = stateout[lev];

     	for (MFIter mfi(stateout_mf, true); mfi.isValid(); ++mfi) {
        	Box bx = mfi.tilebox();
        	ParallelFor(bx, [=, &quad, &cell_count] AMREX_GPU_DEVICE (int i, int j, int k){
				if(k == 0){
					 cell_count++;
					 //amrex::Gpu::Atomic::Add(cell_count, 1);
					 quad.addQuadruplet(lev,i,j,k);
				}
			});
		}
					
	 
	fprintf(solution, "%s %ld %s\n", "POINTS", static_cast<long int> (cell_count*4), "float");

	for(int ipt=0;ipt<quad.get_size();ipt++) {

		Quadruplet iquad = quad.get_quadruplet(ipt);

		int lev = iquad.get_lev();
		int i   = iquad.get_i();
		int j   = iquad.get_j();
		int k   = iquad.get_k();

		Vector<Vector<Real>> node_coords;
		node_coords.push_back(get_node_coords(lev, i, j, k, dx0, plo)); 
		node_coords.push_back(get_node_coords(lev, i+1, j, k, dx0, plo)); 
		node_coords.push_back(get_node_coords(lev, i+1, j+1, k, dx0, plo)); 
		node_coords.push_back(get_node_coords(lev, i, j+1, k, dx0, plo));

		for(int inode=0;inode<node_coords.size();inode++){
			fprintf(solution, "%0.15g %0.15g %0.15g\n", node_coords[inode][0], node_coords[inode][1], 1e-12);
		}
	}

	fprintf(solution,"%s %ld %ld\n", "CELLS", static_cast<long int>(cell_count), static_cast<long int>(cell_count*5));

	for(int ipt=0;ipt<quad.get_size();ipt++) {
		fprintf(solution,"%ld %ld %ld %ld %ld\n", static_cast<long int> (4), static_cast<long int> (ipt*4), static_cast<long int> (ipt*4 + 1), 
											      static_cast<long int> (ipt*4 + 2), static_cast<long int> (ipt*4 + 3));
	}

	fprintf(solution,"%s %ld\n", "CELL_TYPES", static_cast<long int> (cell_count));

	for(int ipt=0;ipt<quad.get_size();ipt++) {
		fprintf(solution,"%ld\n", static_cast<long int> (9));
	}
		
	fprintf(solution,"%s %ld\n", "CELL_DATA", static_cast<long int>(cell_count));
	fprintf(solution,"%s %s %s\n", "SCALARS", "phi", "double");
	fprintf(solution,"%s %s\n","LOOKUP_TABLE", "default");

	for (MFIter mfi(stateout_mf, true); mfi.isValid(); ++mfi) {
    	Array4<const Real> const& stateout_array = stateout_mf.array(mfi);
    	Box bx = mfi.tilebox();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
       		if(k == 0){
				fprintf(solution, "%0.15g ", stateout_array(i,j,k,0));
			}
		});
	}
	fclose(solution);
}

