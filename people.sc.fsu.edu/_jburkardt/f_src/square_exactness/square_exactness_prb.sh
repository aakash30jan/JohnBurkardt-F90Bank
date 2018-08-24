#!/bin/bash
#
gfortran -c square_exactness_prb.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling square_exactness_prb.f90"
  exit
fi
#
gfortran square_exactness_prb.o -L$HOME/lib/$ARCH -lsquare_exactness
if [ $? -ne 0 ]; then
  echo "Errors linking and loading square_exactness_prb.o"
  exit
fi
rm square_exactness_prb.o
#
mv a.out square_exactness_prb
./square_exactness_prb > square_exactness_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running square_exactness_prb"
  exit
fi
rm square_exactness_prb
#
echo "Test program output written to square_exactness_prb_output.txt."
