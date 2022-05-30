#!/bin/sh
cdo -f nc selindexbox,3,1002,3,3 -setgrid,r1004x5 ../fort.82 eta.nc
cdo -f nc selindexbox,3,1002,3,3 -setgrid,r1004x5 ../fort.72 sal.nc
cdo -f nc selindexbox,3,1003,3,3 -setgrid,r1005x5 ../fort.73 uuu.nc
cdo -f nc selindexbox,3,1002,3,3 -setgrid,r1004x5 ../fort.146 www.nc
cdo -f nc copy ../GRID_INFO.ext GRID_INFO.nc


