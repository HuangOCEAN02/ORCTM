#!/bin/sh
cdo -f nc selindexbox,3,502,3,3 -setgrid,r504x5 ../fort.72 sal.nc
cdo -f nc selindexbox,3,503,3,3 -setgrid,r505x5 ../fort.73 uuu.nc
cdo -f nc selindexbox,3,502,3,3 -setgrid,r504x5 ../fort.146 www.nc
cdo -f nc selindexbox,3,502,3,3 -setgrid,r504x5 ../fort.82 eta.nc
cdo -f nc copy ../GRID_INFO.ext GRID_INFO.nc


