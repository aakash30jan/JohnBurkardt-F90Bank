#!/bin/bash
#
#  Compile
#
gfortran -c z_hb.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling z_hb.f90"
  exit
fi
#
gcc -c -I/$HOME/include c_bridge_pzgssv.c
if [ $? -ne 0 ]; then
  echo "Errors compiling c_bridge_pzgssv.c"
  exit
fi
#
#  Link and load
#
gfortran -fopenmp z_hb.o c_bridge_pzgssv.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_openmp -lm -lblas
if [ $? -ne 0 ]; then
  echo "Errors linking and loading z_hb.o + c_bridge_pzgssv.o"
  exit
fi
rm z_hb.o
rm c_bridge_pzgssv.o
mv a.out z_hb
#
#  Run with 1 processor.
#
export OMP_NUM_THREADS=1
./z_hb < cg20_cua.txt > z_hb_1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running z_hb"
  exit
fi
echo "Program output written to z_hb_1_output.txt"
#
#  Run with 4 processors.
#
export OMP_NUM_THREADS=4
./z_hb < cg20_cua.txt > z_hb_4_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running z_hb"
  exit
fi
echo "Program output written to z_hb_4_output.txt"
#
#  Terminate.
#
rm z_hb
