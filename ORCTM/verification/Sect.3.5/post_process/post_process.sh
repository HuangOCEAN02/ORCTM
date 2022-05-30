#!/bin/sh
rm *.nc
cdo -f nc copy ../GRID_INFO.ext GRID_INFO.nc

cdo -f nc copy ../diagnosis/diagnosis.71_007 tem_01.nc
cdo -f nc copy ../diagnosis/diagnosis.71_015 tem_02.nc
cdo -f nc copy ../diagnosis/diagnosis.71_022 tem_03.nc
cdo -f nc copy ../diagnosis/diagnosis.71_030 tem_04.nc
cdo -f nc copy ../diagnosis/diagnosis.71_037 tem_05.nc
cdo -f nc copy ../diagnosis/diagnosis.71_045 tem_06.nc
cdo -f nc copy ../diagnosis/diagnosis.71_067 tem_07.nc
cdo -f nc copy ../diagnosis/diagnosis.71_112 tem_08.nc

cdo -f nc copy ../diagnosis/diagnosis.72_007 sal_01.nc
cdo -f nc copy ../diagnosis/diagnosis.72_015 sal_02.nc
cdo -f nc copy ../diagnosis/diagnosis.72_022 sal_03.nc
cdo -f nc copy ../diagnosis/diagnosis.72_030 sal_04.nc
cdo -f nc copy ../diagnosis/diagnosis.72_037 sal_05.nc
cdo -f nc copy ../diagnosis/diagnosis.72_045 sal_06.nc
cdo -f nc copy ../diagnosis/diagnosis.72_067 sal_07.nc
cdo -f nc copy ../diagnosis/diagnosis.72_112 sal_08.nc


cdo -f nc copy ../diagnosis/diagnosis.73_007 uko_01.nc
cdo -f nc copy ../diagnosis/diagnosis.73_015 uko_02.nc
cdo -f nc copy ../diagnosis/diagnosis.73_022 uko_03.nc
cdo -f nc copy ../diagnosis/diagnosis.73_030 uko_04.nc
cdo -f nc copy ../diagnosis/diagnosis.73_037 uko_05.nc
cdo -f nc copy ../diagnosis/diagnosis.73_045 uko_06.nc
cdo -f nc copy ../diagnosis/diagnosis.73_067 uko_07.nc
cdo -f nc copy ../diagnosis/diagnosis.73_112 uko_08.nc

#cdo -f nc copy ../diagnosis/diagnosis.74_007 vke_01.nc
#cdo -f nc copy ../diagnosis/diagnosis.74_015 vke_02.nc
#cdo -f nc copy ../diagnosis/diagnosis.74_022 vke_03.nc
#cdo -f nc copy ../diagnosis/diagnosis.74_030 vke_04.nc
#cdo -f nc copy ../diagnosis/diagnosis.74_037 vke_05.nc
#cdo -f nc copy ../diagnosis/diagnosis.74_045 vke_06.nc
#cdo -f nc copy ../diagnosis/diagnosis.74_067 vke_07.nc
#cdo -f nc copy ../diagnosis/diagnosis.74_112 vke_08.nc


cdo -f nc copy ../diagnosis/diagnosis.63_007 uso_01.nc
cdo -f nc copy ../diagnosis/diagnosis.63_015 uso_02.nc
cdo -f nc copy ../diagnosis/diagnosis.63_022 uso_03.nc
cdo -f nc copy ../diagnosis/diagnosis.63_030 uso_04.nc
cdo -f nc copy ../diagnosis/diagnosis.63_037 uso_05.nc
cdo -f nc copy ../diagnosis/diagnosis.63_045 uso_06.nc
cdo -f nc copy ../diagnosis/diagnosis.63_067 uso_07.nc
cdo -f nc copy ../diagnosis/diagnosis.63_112 uso_08.nc

#cdo -f nc copy ../diagnosis/diagnosis.64_007 vze_01.nc
#cdo -f nc copy ../diagnosis/diagnosis.64_015 vze_02.nc
#cdo -f nc copy ../diagnosis/diagnosis.64_022 vze_03.nc
#cdo -f nc copy ../diagnosis/diagnosis.64_030 vze_04.nc
#cdo -f nc copy ../diagnosis/diagnosis.64_037 vze_05.nc
#cdo -f nc copy ../diagnosis/diagnosis.64_045 vze_06.nc
#cdo -f nc copy ../diagnosis/diagnosis.64_067 vze_07.nc
#cdo -f nc copy ../diagnosis/diagnosis.64_112 vze_08.nc


cdo -f nc copy ../diagnosis/diagnosis.82_007 eta_01.nc
cdo -f nc copy ../diagnosis/diagnosis.82_015 eta_02.nc
cdo -f nc copy ../diagnosis/diagnosis.82_022 eta_03.nc
cdo -f nc copy ../diagnosis/diagnosis.82_030 eta_04.nc
cdo -f nc copy ../diagnosis/diagnosis.82_037 eta_05.nc
cdo -f nc copy ../diagnosis/diagnosis.82_045 eta_06.nc
cdo -f nc copy ../diagnosis/diagnosis.82_067 eta_07.nc
cdo -f nc copy ../diagnosis/diagnosis.82_112 eta_08.nc


cdo -f nc copy ../diagnosis/diagnosis.146_007 woo_01.nc
cdo -f nc copy ../diagnosis/diagnosis.146_015 woo_02.nc
cdo -f nc copy ../diagnosis/diagnosis.146_022 woo_03.nc
cdo -f nc copy ../diagnosis/diagnosis.146_030 woo_04.nc
cdo -f nc copy ../diagnosis/diagnosis.146_037 woo_05.nc
cdo -f nc copy ../diagnosis/diagnosis.146_045 woo_06.nc
cdo -f nc copy ../diagnosis/diagnosis.146_067 woo_07.nc
cdo -f nc copy ../diagnosis/diagnosis.146_112 woo_08.nc

#cdo -f nc copy ../diagnosis/
#cdo -f nc copy ../diagnosis/
#cdo -f nc copy ../diagnosis/
#cdo -f nc copy ../diagnosis/
#cdo -f nc copy ../diagnosis/
#cdo -f nc copy ../diagnosis/

