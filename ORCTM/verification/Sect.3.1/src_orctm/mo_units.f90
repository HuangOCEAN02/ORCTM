      MODULE MO_UNITS
! ---------------------------------------------------------------------
!
!*    *COMMON* *UNITS*      - Logical units for I/O of the ocean model.
!
!                           - included in all subroutines.
!                               
!
!     S. LEGUTKE          *DKRZ*           07.01.00
!
!*    VARIABLE     TYPE        PURPOSE.
!     --------     ----        --------
!     *IO_STDERR*  *INTEGER*   Logical unit for ocean std eror
!     *IO_STDOUT*  *INTEGER*   Logical unit for ocean stdout
!     *IO_IN_OCTL* *INTEGER*   Logical unit for ocean namelist (OCECTL)
!     ** *INTEGER*   
!
! ---------------------------------------------------------------------
!
      COMMON /UNITS/                                                    &
     &     IO_STDOUT , IO_STDERR , IO_IN_OCTL                           &
     &   , IO_IN_ARCG, IO_IN_GIGI, IO_IN_ANTA                           &
     &   , IO_IN_Z370, IO_IN_Z380                                       &
     &   , IO_IN_SURS, IO_IN_SURT, IO_IN_INIS, IO_IN_INIT               &
     &   , IO_IN_GTEM, IO_IN_GTDE, IO_IN_GSLP                           &
     &   , IO_IN_GWIX, IO_IN_GWIY, IO_IN_GU10                           &
     &   , IO_IN_GSWR, IO_IN_GLWR, IO_IN_GCLO                           &
     &   , IO_IN_GPRE, IO_IN_RVAL                                       &
#ifdef MEAN
     &   , IO_OU_MTHO,   IO_OU_MSAO,   IO_OU_MUKO,   IO_OU_MVKE         &
     &   , IO_OU_MWO,    IO_OU_MPO                                      &
     &   , IO_OU_MMZO,   IO_OU_EMIP,   IO_OU_FLUM,   IO_OU_MPEM         &
     &   , IO_OU_SICT,   IO_OU_SICO,   IO_OU_SICU,   IO_OU_SICV         &
     &   , IO_OU_SICTRU, IO_OU_SICTRV, IO_OU_SICS                       &
     &   , IO_OU_QSWO,   IO_OU_QLWO,   IO_OU_QSEO,   IO_OU_QLAO         &
     &   , IO_OU_PREC,   IO_OU_RIVRUN                                   &
#ifdef GMBOLUS
     &   , IO_OU_MWGO, IO_OU_MBOLX, IO_OU_MBOLY                         &
#endif /*GMBOLUS*/
#ifdef BAROTROPIC 
     &   , IO_OU_MUSO, IO_OU_MVSE                                       &
#endif /*BAROTROPIC*/
#ifdef DIFFDIAG
     &   , IO_OU_MAVO, IO_OU_MDVO, IO_OU_WTMI                           &
#endif /*DIFFDIAG*/
#ifdef NON_HYDROSTATIC
     &   , IO_OU_NHPO, IO_OU_NHDIV, IO_OU_NHPOBOT, IO_OU_NHDIVBOT       &
#endif
#endif /*MEAN*/
     &   , IO_OU_GRID

      END MODULE MO_UNITS
