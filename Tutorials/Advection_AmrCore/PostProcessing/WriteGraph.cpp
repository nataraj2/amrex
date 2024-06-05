#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>

#include <Advection_GNN.H>

using namespace amrex;	
	
void WriteGraphForAllLevelsExceptFinest(const AmrData& amrData, 
										const Vector<iMultiFab>& finemask, 
										const Vector<iMultiFab>& allmasks, 
										const Vector<MultiFab>& stateout,
										FILE* connect)
{
	const Vector<Real>& plo = amrData.ProbLo();
    const Vector<Real>& dx0  = amrData.DxLevel()[0];
    int finest_lev = amrData.FinestLevel();

	for (int lev=0; lev < finest_lev; ++lev) {
	
    	const iMultiFab& allmasks_mf = allmasks[lev];
		const iMultiFab& finemask_mf = finemask[lev];
	
	    for (MFIter mfi(allmasks_mf); mfi.isValid(); ++mfi) {
	    	Array4<const int> const& allmasks_array = allmasks_mf.const_array(mfi);
	    	Array4<const int> const& finemask_array = finemask_mf.const_array(mfi);
	        Box bx = mfi.validbox();
	        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
	       		if(k == 0 and allmasks_array(i,j,k,0) == 0 and finemask_array(i,j,k,0) == 0) {
	            	Vector<Vector<Real>> vec_coords;
	            	Vector<Real> coords = get_coords(lev, i, j, k, dx0, plo);
	                vec_coords.push_back(coords);
	
	                QuadrupletStore quad;
	                for(int i_shift=-3; i_shift<=3; i_shift++){
	                	for(int j_shift=-3; j_shift<=3; j_shift++){
	                		if(!(i_shift == 0  and j_shift == 0)){
	                        	FindNeighbors(lev, i+i_shift, j+j_shift, k, allmasks_array, finemask_array, dx0, plo, vec_coords, quad);
	                        }
	                    }
	                }
	
	                AMREX_ALWAYS_ASSERT(quad.get_size() == vec_coords.size()-1);
	
	                if(finemask_array(i,j+1,k,0) == 0 and allmasks_array(i,j+1,k,0) == 2){
	                	fprintf(connect,"%s %ld %s\n","POINTS",vec_coords.size(), "float");
	                	for(long unsigned int ipt=0;ipt<vec_coords.size();ipt++)
	                	{
	                    	fprintf(connect, "%0.15g %0.15g %0.15g\n", vec_coords[ipt][0], vec_coords[ipt][1], vec_coords[ipt][2]);
	                	}
	                    fprintf(connect,"%s %ld %ld\n", "LINES",vec_coords.size()-1, (vec_coords.size()-1)*3);
	                    for(long unsigned int ipt=1;ipt<vec_coords.size();ipt++){
	                    	fprintf(connect,"%d %d %ld\n", 2, 0, ipt);
	                    }
	                    fclose(connect);
	                    exit(0);
	                }
	            }
	        });
	    }
	}
}

void WriteGraphForFinestLev(const AmrData& amrData, 
							const Vector<iMultiFab>& allmasks,
							const Vector<MultiFab>& stateout,
							FILE* connect)
{
	 int lev = amrData.FinestLevel(); 
	 const Vector<Real>& plo = amrData.ProbLo();
     const Vector<Real>& dx0  = amrData.DxLevel()[0];

     const iMultiFab& allmasks_mf = allmasks[lev];
	 const MultiFab& stateout_mf = stateout[lev];

     for (MFIter mfi(allmasks_mf); mfi.isValid(); ++mfi) {
     	Array4<const int> const& allmasks_array  = allmasks_mf.array(mfi);
     	Array4<const Real> const& stateout_array = stateout_mf.array(mfi);
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
                    		FindNeighborsFinestLev(lev, i+i_shift, j+j_shift, k, allmasks_array, dx0, plo, vec_coords, quad);
                    	}
                	}
            	}
			
				AMREX_ALWAYS_ASSERT(quad.get_size() == vec_coords.size()-1);

            	if(allmasks_array(i,j+1,k,0) == 0 and
               		allmasks_array(i,j-1,k,0) == 0 and
               		allmasks_array(i+1,j,k,0) == 0 and
               		allmasks_array(i-1,j,k,0) == 0){
             		fprintf(connect,"%s %ld %s\n","POINTS",vec_coords.size(), "float");
             		for(long unsigned int ipt=0;ipt<vec_coords.size();ipt++)
                	{
                		fprintf(connect, "%0.15g %0.15g %0.15g\n", vec_coords[ipt][0], vec_coords[ipt][1], vec_coords[ipt][2]);
                	}
                	fprintf(connect,"%s %ld %ld\n", "LINES",vec_coords.size()-1, (vec_coords.size()-1)*3);
                	for(long unsigned int ipt=1;ipt<vec_coords.size();ipt++){
                		fprintf(connect,"%d %d %ld\n", 2, 0, ipt);
                	}
                	fclose(connect);
                	exit(0);
             	}
         	}
     	});
    }
}

void WriteSolution(const AmrData& amrData,
				   const Vector<iMultiFab>& allmasks,  
                   const Vector<MultiFab>& stateout,
                   FILE* solution)

{

	int finest_lev = amrData.FinestLevel();

	const Vector<Real>& plo = amrData.ProbLo();
    const Vector<Real>& dx0  = amrData.DxLevel()[0];

	Vector<int> vec_cell_count;
	vec_cell_count.resize(finest_lev+1);

	QuadrupletStore quad;

    for (int lev=0; lev <= finest_lev; ++lev) {

        const iMultiFab& allmasks_mf = allmasks[lev];
	 	const MultiFab& stateout_mf = stateout[lev];

	    int *total_cell_count = &vec_cell_count[lev];
		for (MFIter mfi(allmasks_mf); mfi.isValid(); ++mfi) {
            Array4<const int> const& allmasks_array  = allmasks_mf.array(mfi);
			Box bx = mfi.validbox();
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
				if(k == 0){
					amrex::Gpu::Atomic::Add(total_cell_count, 1);
				}
			});
		}
		//std::cout << "total cell count at level " << lev << " is " << vec_cell_count[lev] << "\n";
		//exit(0);

     	for (MFIter mfi(allmasks_mf); mfi.isValid(); ++mfi) {
        	Array4<const int> const& allmasks_array  = allmasks_mf.array(mfi);
        	Array4<const Real> const& stateout_array = stateout_mf.array(mfi);
        	Box bx = mfi.validbox();
        	ParallelFor(bx, [=, &quad] AMREX_GPU_DEVICE (int i, int j, int k){
				if(k == 0 and allmasks_array(i,j,k,0) == 0){
					 amrex::Gpu::Atomic::Add(total_cell_count, 1);
					 quad.addQuadruplet(lev,i,j,k);
				}
			});
		}
	}
					 
	int cell_count = std::accumulate(vec_cell_count.begin(), vec_cell_count.end(), 0);

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

		for(int ipt=0;ipt<node_coords.size();ipt++){
			fprintf(solution, "%0.15g %0.15g %0.15g\n", node_coords[ipt][0], node_coords[ipt][1], node_coords[ipt][2]);
		}
	}

	fprintf(solution,"%s %ld %ld\n", "CELLS", cell_count, cell_count*5);

	for(int ipt=0;ipt<quad.get_size();ipt++) {
		fprintf(solution,"%ld %ld %ld %ld %ld\n", static_cast<long int> (4), static_cast<long int> (ipt*4), static_cast<long int> (ipt*4 + 1), 
											      static_cast<long int> (ipt*4 + 2), static_cast<long int> (ipt*4 + 3));
	}

	fprintf(solution,"%s %ld\n", "CELL_TYPES", cell_count);

	for(int ipt=0;ipt<quad.get_size();ipt++) {
		fprintf(solution,"%ld\n", static_cast<long int> (9));
	}
		
	/*fprintf(solution,"%s %ld\n", "CELL_DATA", cell_count);
	fprintf(solution,"%s %s %s\n", "SCALARS", "phi", "double");
	fprintf(solution,"%s %s\n","LOOKUP_TABLE", "default");

	for(int lev = 0;lev<=finest_lev;lev++) {
        const MultiFab& stateout_mf = stateout[lev];
		const iMultiFab& allmasks_mf = allmasks[lev];
		for (MFIter mfi(stateout_mf); mfi.isValid(); ++mfi) {
            Array4<const Real> const& stateout_array = stateout_mf.array(mfi);
            Array4<const int> const& allmasks_array  = allmasks_mf.array(mfi);
            Box bx = mfi.validbox();
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
                if(k == 0 and allmasks_array(i,j,k,0) == 0){
					fprintf(solution, "%0.15g ", stateout_array(i,j,k,0));
				}
			});
		}
	}*/
}

