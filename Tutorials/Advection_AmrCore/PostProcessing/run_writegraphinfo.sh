rm -rf main3d.gnu.ex
rm -rf *vtk
make -j
./main3d.gnu.ex infile=../Exec/plt00030 varNames=phi
