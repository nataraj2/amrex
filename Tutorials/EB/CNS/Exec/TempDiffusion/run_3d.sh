sh run_clean.sh
rm -rf CNS3d.gnu.MPI.ex
make -j 
mpirun -n 4 ./CNS3d.gnu.MPI.ex inputs
