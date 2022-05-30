#!/bin/sh
rm anta.nc arcgri.nc
rm a.out anta arcgri GI* INI*
ifort -r8 -convert big_endian mkgrid.f90
./a.out -Wl,-T
cdo -f nc copy anta anta.nc
cdo -f nc copy arcgri arcgri.nc
cp anta arcgri GI* INI* ../
