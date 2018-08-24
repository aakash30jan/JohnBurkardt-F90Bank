#!/bin/bash
#
mpif90 wave_mpi.f90
#
if [ $? -ne 0 ]; then
  echo "Errors compiling wave_mpi.f90"
  exit
fi
#
#  Rename the executable.
#
mv a.out wave
#
#  Run the program under MPI.
#
mpirun -np 4 ./wave > wave_local_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wave"
  exit
fi
#
#  Clean up.
#
rm wave
#
echo "Program output written to wave_local_output.txt"

