#!/bin/bash
#
#  Compile
#
gfortran -c d_hb.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling d_hb.f90"
  exit
fi
#
gcc -c -I/$HOME/include c_bridge_pdgssv.c
if [ $? -ne 0 ]; then
  echo "Errors compiling c_bridge_pdgssv.c"
  exit
fi
#
#  Link and load
#
gfortran -fopenmp d_hb.o c_bridge_pdgssv.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_openmp -lm -lblas
if [ $? -ne 0 ]; then
  echo "Errors linking and loading d_hb.o + c_bridge_pdgssv.o"
  exit
fi
rm d_hb.o
rm c_bridge_pdgssv.o
mv a.out d_hb
#
#  Run with 1 processor.
#
export OMP_NUM_THREADS=1
./d_hb < g20_rua.txt > d_hb_1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running d_hb"
  exit
fi
echo "Program output written to d_hb_1_output.txt"
#
#  Run with 4 processors.
#
export OMP_NUM_THREADS=4
./d_hb < g20_rua.txt > d_hb_4_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running d_hb"
  exit
fi
echo "Program output written to d_hb_4_output.txt"
#
#  Terminate.
#
rm d_hb
