#!/bin/sh
cd globalfiles
rm *
cd ..
rm *.mod Glue.x log Glue.o*
NETCDFROOT=/lustre/software/pnetcdf-1.12.2_intel2016/

mpiifort -O3 -r8 -fpp -convert big_endian -I/lustre/software/pnetcdf-1.12.2_intel2016/include Glue.f90 -o Glue.x -L/lustre/software/pnetcdf-1.12.2_intel2016/lib -lpnetcdf  
mpirun -n 120 ./Glue.x
