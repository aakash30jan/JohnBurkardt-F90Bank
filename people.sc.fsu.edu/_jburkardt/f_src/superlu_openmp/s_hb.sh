#!/bin/bash
#
#  Compile
#
gfortran -c s_hb.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling s_hb.f90"
  exit
fi
#
gcc -c -I/$HOME/include c_bridge_psgssv.c
if [ $? -ne 0 ]; then
  echo "Errors compiling c_bridge_psgssv.c"
  exit
fi
#
#  Link and load
#
gfortran -fopenmp s_hb.o c_bridge_psgssv.o -L$HOME/lib/$ARCH \
  -L/$HOME/libc/$ARCH -lsuperlu_openmp -lm -lblas
if [ $? -ne 0 ]; then
  echo "Errors linking and loading s_hb.o + c_bridge_psgssv.o"
  exit
fi
rm s_hb.o
rm c_bridge_psgssv.o
mv a.out s_hb
#
#  Run with 1 processor.
#
export OMP_NUM_THREADS=1
./s_hb < g20_rua.txt > s_hb_1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running s_hb"
  exit
fi
echo "Program output written to s_hb_1_output.txt"
#
#  Run with 4 processors.
#
export OMP_NUM_THREADS=4
./s_hb < g20_rua.txt > s_hb_4_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running s_hb"
  exit
fi
echo "Program output written to s_hb_4_output.txt"
#
#  Terminate.
#
rm s_hb
