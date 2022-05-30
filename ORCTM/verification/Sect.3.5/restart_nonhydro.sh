#!/bin/sh

## compile source code
cd src_orctm
make clean
rm make.log
make > make.log

## prepare for initial files
cd ../BT_forc/
rm mkgrid.o anta arcgri GI* INI*
ifort -r8 -convert big_endian mkgrid.f90 -o mkgrid.o
./mkgrid.o -Wl,-T

## prepare for Open boundary files
cd ../OBCS/
matlab -nosplash -nodesktop -r "Gen_Boundary.m;exit"

## copy to run folder
cd ..
rm oceout* Z3* topo GRID_INFO.ext
rm anta arcgri GI* INI* orctm_nh.*
rm -r fort* diagnosis* restart_file

cp ./src_orctm/orctm_nh.x ./
cp ./BT_forc/anta    ./
cp ./BT_forc/arcgri  ./
cp ./BT_forc/GI*     ./
cp ./BT_forc/INI*    ./
cp ./obcs/OB*.bin ./


## start simulation
mpirun -n 120 ./orctm_nh.x -cg3d_pc_type mg -cg3d_pc_mg_type full -cg3d_ksp_type fgmres -cg3d_ksp_monitor_short -cg3d_pc_mg_levels 1

# merge the fragmental file
cd ./post_process/Glue
bash ./Make_Glue.sh