#!/bin/sh
rm anta.nc arcgri.nc
rm mkgrid.o anta arcgri GI* INI*
ifort -r8 -convert big_endian mkgrid.f90 -o mkgrid.o
./mkgrid.o -Wl,-T
cp anta arcgri GI* INI* ../

