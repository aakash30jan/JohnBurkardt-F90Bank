#! /bin/bash
#
mkdir temp
cd temp
rm *
~/bin/f90split ../r8lib.f90
#
for FILE in `ls -1 *.f90`;
do
  gfortran -c $FILE
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
done
rm *.f90
#
ar qc libr8lib.a *.o
rm *.o
#
mv libr8lib.a ~/lib
cd ..
rmdir temp
#
echo "Library installed as ~/lib/libr8lib.a"
