rm -rf WriteGraphInfo3d.gnu.DEBUG.ex
rm -rf *vtk
make -j
./WriteGraphInfo3d.gnu.DEBUG.ex infile=../Exec/plt00050 varNames=phi
