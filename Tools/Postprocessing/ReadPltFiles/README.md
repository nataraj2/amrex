# Post processor to read in plot files 
The plot files to be read in should be placed in this directory.
1. `git clone --recursive https://github.com/nataraj2/amrex.git`
2. `cd amrex`
3. `git checkout ReadPltFiles`
4. `cd Tools/Postprocessing/ReadPltFiles` 
5. In the `GNUMakefile`, set `NDIM=2` for 2D and `NDIM=3` for 3D    
   `make -j`
6. `./<executable> infile=<plt filename> varNames=var1 var2 var3`  
So for eg.     
`./main2d.gnu.ex infile=plt_00000 varNames=x_velocity y_velocity`
Note that currently the `fprintf` to write the `csv` files are hard coded for 2D 
and 2 variables. 
