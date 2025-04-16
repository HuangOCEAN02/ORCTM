#!/bin/sh
rm *.nc
rm *.o *.mod
rm GI* INI*

ifort -r8 -convert big_endian chpar.f90
# Make Grid 
ifort -r8 -convert big_endian -o mkanta.o mkanta.f90
./mkanta.o -Wl,-t

ifort -r8 -convert big_endian -o mkgrid.o mkgrid.f90
./mkgrid.o -Wl,-T
# Make Initial Fields
ifort -r8 -convert big_endian -o mkini.o mkini.f90
./mkini.o -Wl,-T
# Make Forcing Fields
ifort -r8 -convert big_endian -o mkGI.o mkGI.f90
./mkGI.o -Wl,-T


# output
cdo -f nc copy anta_region anta_region.nc
cdo -f nc copy anta anta.nc
cdo -f nc copy arcgri arcgri.nc

#cp anta arcgri GI* INI* ../

