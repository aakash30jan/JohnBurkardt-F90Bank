#!/bin/bash
#
gfortran -c square_hex_grid_prb.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling square_hex_grid_prb.f90"
  exit
fi
#
gfortran square_hex_grid_prb.o -L$HOME/lib/$ARCH -lsquare_hex_grid
if [ $? -ne 0 ]; then
  echo "Errors linking and loading square_hex_grid_prb.o"
  exit
fi
rm square_hex_grid_prb.o
#
mv a.out square_hex_grid_prb
./square_hex_grid_prb > square_hex_grid_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running square_hex_grid_prb"
  exit
fi
rm square_hex_grid_prb
#
echo "Test program output written to square_hex_grid_prb_output.txt."
