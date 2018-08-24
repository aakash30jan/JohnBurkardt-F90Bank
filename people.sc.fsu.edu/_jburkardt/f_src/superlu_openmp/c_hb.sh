#!/bin/bash
#
#  Compile
#
gfortran -c c_hb.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling c_hb.f90"
  exit
fi
#
gcc -c -I/$HOME/include c_bridge_pcgssv.c
if [ $? -ne 0 ]; then
  echo "Errors compiling c_bridge_pcgssv.c"
  exit
fi
#
#  Link and load
#
gfortran -fopenmp c_hb.o c_bridge_pcgssv.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_openmp -lm -lblas
if [ $? -ne 0 ]; then
  echo "Errors linking and loading c_hb.o + c_bridge_pcgssv.o"
  exit
fi
rm c_hb.o
rm c_bridge_pcgssv.o
mv a.out c_hb
#
#  Run with 1 processor.
#
export OMP_NUM_THREADS=1
./c_hb < cg20_cua.txt > c_hb_1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running c_hb"
  exit
fi
echo "Program output written to c_hb_1_output.txt"
#
#  Run with 4 processors.
#
export OMP_NUM_THREADS=4
./c_hb < cg20_cua.txt > c_hb_4_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running c_hb"
  exit
fi
echo "Program output written to c_hb_4_output.txt"
#
#  Terminate.
#
rm c_hb
