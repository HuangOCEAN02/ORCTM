      SUBROUTINE SETUNITS
      USE MO_UNITS
!     SET THE INPUT UNITS
!HH   29.11.99
!
      IO_STDOUT=6
      IO_IN_OCTL=10
!
!     MODELGRID AND TOPOFILE
      IO_IN_ARCG=21
      IO_IN_GIGI=91
      IO_IN_ANTA=81
!
!     RESTART FILES
      IO_IN_Z370=27
      IO_IN_Z380=28
!
!     FORCING
      IO_IN_SURS=38
      IO_IN_SURT=338
      IO_IN_INIS=33
      IO_IN_INIT=32
      IO_IN_GTEM=37
      IO_IN_GTDE=54
      IO_IN_GSLP=57
      IO_IN_GWIX=40
      IO_IN_GWIY=41
      IO_IN_GU10=55
      IO_IN_GSWR=53
      IO_IN_GLWR=56
      IO_IN_GCLO=51
      IO_IN_GPRE=52
      IO_IN_RVAL=20
!
#ifdef MEAN
      IO_OU_MTHO=71
      IO_OU_MSAO=72
      IO_OU_MUKO=73
      IO_OU_MVKE=74
      IO_OU_MWO=146
      IO_OU_MPO=76
      IO_OU_MRO=880
      IO_OU_MMZO=82
#ifdef BAROTROPIC      
      IO_OU_MUSO=63
      IO_OU_MVSE=64  
#endif           
      IO_OU_EMIP=79
      IO_OU_FLUM=84
      IO_OU_MPEM=85
      IO_OU_SICT=86
      IO_OU_SICO=87
      IO_OU_SICU=88
      IO_OU_SICV=89
      IO_OU_SICTRU=147
      IO_OU_SICTRV=148
      IO_OU_SICS=136
      IO_OU_QSWO=137
      IO_OU_QLWO=138
      IO_OU_QSEO=140
      IO_OU_QLAO=139
      IO_OU_PREC=141
      IO_OU_RIVRUN=305
#ifdef GMBOLUS
      IO_OU_MWGO=246
      IO_OU_MBOLX=159
      IO_OU_MBOLY=160
#endif /*GMBOLUS*/

#ifdef DIFFDIAG
      IO_OU_MAVO=144
      IO_OU_MDVO=145
!
      IO_OU_WTMI=245 ! Wind mixing
#endif /*DIFFDIAG*/


#ifdef NON_HYDROSTATIC
      IO_OU_NHPO=77
      IO_OU_NHDIV=78
#ifdef CG2D4NONHY        
      IO_OU_NHPOBOT=177
      IO_OU_NHDIVBOT=178
#endif /*CG2D4NONHY*/      
#endif /*NON_HYDROSTATIC*/

#ifdef ORLANSKI
     IO_IN_O370 =370
     IO_IN_O380 =380
#endif
#endif /*MEAN*/
      IO_OU_GRID=999
!
      RETURN
      END
