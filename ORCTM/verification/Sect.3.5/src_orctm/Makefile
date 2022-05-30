PROG =	orctm_nh.x	

VPATH = . 

SRCS =	absturz.f90 adisit1.f90 adisitj.f90 aufr.f90 aufw.f90 aufr_fragment.f90 aufw_fragment.f90 \
        bartim.f90 beleg.f90 beleg_zero.f90 boden.f90 budget.f90 \
        calcgmvis.f90 coriol.f90 cppoutput.f90 findalfa.f90 \
        findbot.f90 growth.f90 itprep.f90 levire.f90 mo_commo1.f90 \
        mo_commo2.f90 mo_commoau1.f90 mo_commoau2.f90 mo_commoau3.f90 \
        mo_commobbl.f90 mo_diagnosis.f90 mo_elicom.f90 mo_kind.f90 mo_levitus.f90 \
        mo_mean.f90 mo_mpi.f90 mo_obcs.f90 mo_octdiff.f90 \
        mo_para2.f90 mo_parallel.f90 mo_param1.f90 mo_param3.f90 \
        mo_petsc.f90 mo_units.f90 monlen_feb.f90 orctm.f90 \
        nlopps.f90 obudget.f90 ocadfs.f90 ocadpo.f90 ocbarp.f90 occlit.f90 \
        ocice.f90 ocjitr.f90 ocmodmom.f90 ocschep.f90 octdiff_base.f90 \
        octdiff_ts.f90 octher.f90 octimf.f90 ocuad.f90 ocvad.f90 \
        ocvisc.f90 ocvtot.f90 ocvtot_nh.f90 ocvtro.f90 ocwad.f90 \
        ocwind.f90 open_stdout.f90 rho.f90 rho1j.f90 rho2.f90 \
        setunits.f90 slopetrans.f90 sqrtrnd.f90 trian.f90 troneu.f90 \
        trotest.f90 vardrag.f90 wrte_diagnosis.f90 wrte_gridinfo.f90 wrte_mean.f90



OBJSM =	absturz.o adisit1.o adisitj.o aufr.o aufw.o aufr_fragment.o aufw_fragment.o \
        bartim.o beleg.o beleg_zero.o boden.o budget.o \
        calcgmvis.o coriol.o cppoutput.o findalfa.o \
        findbot.o growth.o itprep.o levire.o mo_commo1.o \
        mo_commo2.o mo_commoau1.o mo_commoau2.o mo_commoau3.o \
        mo_commobbl.o mo_diagnosis.o mo_elicom.o mo_kind.o mo_levitus.o \
        mo_mean.o mo_mpi.o mo_obcs.o mo_octdiff.o \
        mo_para2.o mo_parallel.o mo_param1.o mo_param3.o \
        mo_petsc.o mo_units.o monlen_feb.o orctm.o \
        nlopps.o obudget.o ocadfs.o ocadpo.o ocbarp.o occlit.o \
        ocice.o ocjitr.o ocmodmom.o ocschep.o octdiff_base.o \
        octdiff_ts.o octher.o octimf.o ocuad.o ocvad.o \
        ocvisc.o ocvtot.o ocvtot_nh.o ocvtro.o ocwad.o \
        ocwind.o open_stdout.o rho.o rho1j.o rho2.o \
        setunits.o slopetrans.o sqrtrnd.o trian.o troneu.o \
        trotest.o vardrag.o wrte_diagnosis.o wrte_gridinfo.o wrte_mean.o


DEF =   -DZZNOMPI -DVERSIONGR03 -DLEVELS40 -DMEAN -DRESYEAR \
        -DZZDEBUG_ONEDAY -DZZQLOBERL -DBULK_KARA -DZZYEAR360 \
        -DEISREST -DREDWMICE -DALBOMIP -DZZISOPYK -DZZGMBOLUS \
        -DBOLK05 -DADPO -DZZSLOPECON_ADPO -DOPEND55 -DNURDIF \
        -DZZNOSKIP -DRYEAR -DZZSW089 -DZZLEADCLOSE -DSOR \
        -DFORCE_DWLW -DZZTIMECHECK -DRESTORE_MON -DFREESLIP \
        -DZZCYCLICNS -DZZCYCLICWE -DNON_HYDROSTATIC -DZZCG2D4NONHY \
        -DZZDIFFDIAG -DOBC_UV_FORC -DZZOBC_ETA_FORC -DOBC_TS_FORC \
        -DOBC_TIDE -DZZOBC_SPONGE_EXP -DBAROTROPIC -DDIAGNOSIS \
        -DFRAGMENT
        


PETSC_ARCH = arch-linux2-c-debug
PETSC_DIR = ../../../petsc-3.10.2
#$(info PETSC_DIR is $(PETSC_DIR))

include ${PETSC_DIR}/lib/petsc/conf/variables
PETSC_INCLUDE=-I${PETSC_DIR}/src

INCLUDES = ${MPI_INCLUDE}
LIBS = ${MPI_LIB}

#CC = mpiicc
CFLAGS = -g


#FOR LINUX INTEL
#F90 = mpiifort
#F90 = mpif90
#F90FLAGS = -O3 -ipo -static -tpp7 -r8 -cpp -Vaxlib -convert big_endian ${INCLUDES} ${LIBS} 
#F90FLAGS = -O3 -ipo -static -r8 -cpp -Vaxlib -convert big_endian ${INCLUDES} ${LIBS} 
#F90FLAGS = -O3 -static -r8 -cpp -Vaxlib -convert big_endian ${INCLUDES} ${LIBS} 
F90FLAGS = -O3 -r8 -fpp -convert big_endian ${INCLUDES} ${LIBS}


LDFLAGS = $(F90FLAGS)

all: $(PROG)
#
P_FC    = $(PETSC_FCOMPILE)
P_CC    = $(PETSC_COMPILE)
$(PROG): $(OBJSM)
	$(FLINKER) $(LDFLAGS) -o $@ $(OBJSM) $(LIBS) $(PETSC_LIB)

clean:
	rm -f $(PROG) $(OBJSM) *.mod *.o *.L

.SUFFIXES: $(SUFFIXES) .F90 .f90 .f .F .c
%.o: %.c
	$(P_CC) $(CFLAGS) -c $(DEF) $<
%.o: %.F90
	$(P_FC) $(F90FLAGS) -c $(DEF) $<
%.o: %.f90
	$(P_FC) $(F90FLAGS) -c $(DEF) $<
%.o: %.F
	$(P_FC) $(F90FLAGS) -c $(DEF) $<
%.o: %.f
	$(P_FC) $(F90FLAGS) -c $(DEF) $<
#
#--------------------------------------------------------------------------------
# Dependencies
#

absturz.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o
adisit1.o : mo_param3.o
adisitj.o :mo_param1.o mo_param3.o
aufr.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o mo_mpi.o
aufw.o : mo_commo1.o mo_commo2.o mo_commoau2.o mo_commobbl.o mo_param1.o mo_units.o mo_parallel.o mo_mpi.o
aufr_fragment.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o mo_mpi.o
aufw_fragment.o : mo_commo1.o mo_commo2.o mo_commoau2.o mo_commobbl.o mo_param1.o mo_units.o mo_parallel.o mo_mpi.o
bartim.o : mo_commo1.o mo_elicom.o mo_param1.o mo_units.o mo_parallel.o mo_mpi.o
beleg.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o mo_mpi.o
beleg_zero.o : mo_commo1.o mo_commo2.o mo_commoau2.o mo_commoau3.o mo_commobbl.o mo_elicom.o mo_mean.o mo_para2.o mo_param1.o
boden.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o mo_mpi.o
budget.o : mo_commoau1.o mo_param1.o
calcgmvis.o : mo_commo1.o mo_commoau1.o mo_param1.o mo_units.o mo_parallel.o
coriol.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o mo_mpi.o
cppoutput.o : mo_commo1.o mo_param1.o mo_units.o
findalfa.o : mo_commo1.o mo_commobbl.o mo_param1.o mo_parallel.o mo_units.o mo_mpi.o
findbot.o : mo_commo1.o mo_param1.o
growth.o : mo_commoau1.o mo_commoau3.o mo_param1.o mo_units.o
itprep.o : mo_commo1.o mo_para2.o mo_param1.o mo_parallel.o
levire.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o mo_levitus.o mo_mpi.o
mo_commo1.o : mo_kind.o mo_param1.o
mo_commo2.o : mo_param1.o
mo_commoau1.o :
mo_commoau2.o : mo_param1.o
mo_commoau3.o : mo_param1.o
mo_commobbl.o : mo_param1.o
mo_diagnosis.o : mo_param1.o mo_mpi.o mo_parallel.o mo_commo1.o mo_kind.o
mo_elicom.o : mo_param1.o
mo_kind.o :
mo_levitus.o : mo_param1.o
mo_mean.o : mo_param1.o
mo_mpi.o : mo_kind.o
mo_obcs.o : mo_commo1.o mo_param1.o mo_mpi.o mo_parallel.o mo_units.o
mo_octdiff.o : mo_param1.o
mo_para2.o : mo_param1.o
mo_parallel.o : mo_param1.o mo_mpi.o mo_kind.o mo_units.o
mo_param1.o :
mo_param3.o :
mo_petsc.o : mo_commo1.o mo_param1.o mo_parallel.o mo_obcs.o
mo_units.o :
monlen_feb.o :
orctm.o : mo_commo1.o mo_commo2.o mo_commoau1.o mo_commoau2.o mo_commoau3.o mo_commobbl.o mo_mean.o mo_param1.o mo_units.o mo_parallel.o mo_mpi.o mo_obcs.o mo_diagnosis.o mo_octdiff.o mo_para2.o mo_elicom.o mo_levitus.o
nlopps.o : mo_param1.o
obudget.o : mo_commoau1.o mo_param1.o
ocadfs.o : mo_commo1.o mo_commo2.o mo_commoau1.o mo_param1.o mo_parallel.o
ocadpo.o : mo_commo1.o mo_commo2.o mo_commoau1.o mo_commobbl.o mo_param1.o mo_units.o mo_parallel.o
ocbarp.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o
occlit.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o
ocice.o : mo_commo1.o mo_commo2.o mo_commoau1.o mo_commoau2.o mo_param1.o mo_units.o mo_parallel.o
ocjitr.o : mo_commo1.o mo_commo2.o mo_commoau1.o mo_param1.o mo_units.o mo_parallel.o
ocmodmom.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o
ocschep.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o
octdiff_base.o : mo_commo1.o mo_commoau1.o mo_param1.o mo_units.o mo_parallel.o mo_octdiff.o
octdiff_ts.o : mo_commo1.o mo_commoau1.o mo_param1.o mo_units.o mo_parallel.o mo_octdiff.o mo_levitus.o
octher.o : mo_commo1.o mo_commoau1.o mo_commoau2.o mo_mpi.o mo_param1.o mo_units.o mo_parallel.o
octimf.o : mo_commo1.o mo_commo2.o mo_param1.o
ocuad.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o
ocvad.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o
ocvisc.o : mo_commo1.o mo_commo2.o mo_param1.o mo_parallel.o
ocvtot.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o mo_mpi.o
ocvtot_nh.o : mo_commo1.o  mo_param1.o mo_parallel.o mo_petsc.o
ocvtro.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o
ocwad.o : mo_commo1.o mo_commo2.o mo_param1.o mo_units.o mo_parallel.o
ocwind.o : mo_commo1.o mo_commo2.o mo_commoau1.o mo_param1.o
open_stdout.o : mo_parallel.o
rho.o :
rho1j.o : mo_param1.o
rho2.o :
setunits.o : mo_units.o
slopetrans.o : mo_commo1.o mo_commobbl.o mo_param1.o mo_parallel.o
sqrtrnd.o :
trian.o : mo_commo1.o mo_commo2.o mo_elicom.o mo_param1.o mo_units.o mo_parallel.o mo_mpi.o
troneu.o : mo_commo1.o mo_para2.o mo_param1.o mo_units.o mo_parallel.o mo_obcs.o
trotest.o : mo_commo1.o mo_para2.o mo_param1.o mo_units.o mo_parallel.o
vardrag.o : mo_param1.o
wrte_diagnosis.o : mo_param1.o mo_commo1.o mo_diagnosis.o
wrte_gridinfo.o : mo_param1.o mo_commo1.o mo_units.o mo_parallel.o mo_kind.o
wrte_mean.o : mo_param1.o mo_commo1.o mo_units.o mo_parallel.o mo_commoau1.o mo_commoau2.o mo_mean.o



