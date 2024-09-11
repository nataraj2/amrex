# Post processor to read in plot files 

1. `git clone --recursive https://github.com/nataraj2/amrex.git`
2. `git checkout ReadPltFiles`
3. `cd Tools/Postprocessing/ReadPltFiles` 
4. In the GNUMakefile, set `NDIM=2` for 2D and `NDIM=3` for 3D    
   `make -j`
5. `./<executable> infile=plt_00000 varNames=x_velocity y_velocity`  
Note that currently the `fprintf` to write the `csv` files are hard coded for 2D 
and 2 variables. 
