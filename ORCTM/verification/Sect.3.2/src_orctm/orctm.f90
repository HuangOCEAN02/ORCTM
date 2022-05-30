      PROGRAM ORCTM   
!
!    ORCTM: AN OCEANIC REGIONAL CIRCULATION AND TIDE MODEL FOR INTERNAL SOLITARY WAVES SIMULTAIONS   
!
!    NONHYDROSTATIC VERSION OF THE MAX-PLANCK-INSTITUTE OCEAN MODEL ACROSS A REGIONAL OCEAN.
!
!    DEVELOPED BY Pengyang Song(SPY), Jiaqi Guo(GJQ), Hao Huang(H.H), and Shi Qiu(QS) 2022
!
!    INCLUDES CONFORMAL MAPPING     
!
!    MODIFIED :
!    ORCTM:           Nov 20 2020  SPY         - Orlanski radiation boundary conditions
!    ORCTM:           Nov 20 2020  SPY&H.H     - Barotropical signals Sponge Relaxation boundary conditions
!    ORCTM:           FRI 20 2021  GJQ         - Flather boundary and 2-D Orlanski boundary conditions for tidal currents
!    ORCTM:           JUL 03 2021  QS          - Boundary Volume Balance condition
!    ORCTM:                                      - above all created mo_obcs.f90 & octher.f90
!    ORCTM:           MAR 21 2021  H.H         - N-S/W-E Periodic boundary conditions 
!                     JUL 05 2021  SPY&H.H     - Nonhydrostatic Module(3D/2D-Possion) 
!                                                - created ocwad.f90, mo_petsc.f90 & ocvtot_nh.f90
!                                                - updated ocvisc.f90, ocmodmom.f90 & IO Moudles.
!    ORCTM:           JUL 30 2021  H.H         - Diagnositic Output Module 
!                                                - created diagnosis.f90 & mo_diagnosis.f90
!    ORCTM:           SEP 17 2021  H.H         - Fragmentization IO 
!                                                - created aufr_fragment.f90 & aufw_fragment.f90
!*************************************************************************
!
!
!                       MP ORCTM
!                      SBR BELEG
!                          BODEN
!                          CORIOL
!                          ITPREP
!                            !
!                 !  --->  THERMODYNAMIC FORCING(Surface and Open Boundary)
!                 !        WIND FORCING
!                 !        OPEN BOUNDARY CONDITIONS
!      TIME       !        DECOMPOSITION INTO BAROTROPIC AND BAROCLINIC FIELD
!      STEPPING   !        BAROTROPIC SYSTEM
!                 !        BAROCLINIC SYSTEM
!                 !        MOMENTUM ADVECTION
!                 !        TRACER ADVECTION
!                 !        TRACER DIFFUSION
!                 !        MOMENTUM DIFFUSION
!                 !        NONHYDROSTAIC DYNAMICS
!                 ! <---     !
!                            !
!                          OUTPUT-ROUTINES

!          ------------------------------
! NUMBER OF VERTICAL LAYERS IS KE.
!
!************************************************************
! PARAMETER :  IE   NUMBER OF GRID POINTS IN X 
!              JE                            Y
!              KE                            Z
!
! SOME IMPORTANT VARIABLES AND FIELDS :
!             DT        TIME STEP
!             TIESTU    DEPTH OF HORIZONTAL VELOCITY POINTS
!             TIESTW    DEPTH OF VERTICAL VELOCITY POINTS
!             ZO        SEA SURFACE ELEVATION
!             AMSUE/O   LAND/SEA-MASK FOR VECTORFIELD (LAND=0/SEA=1)
!             WETO                    FOR SCALARFIELD       " 
!             UKO       ZONAL VELOCITY COMPONENT
!             VKE       MERIDIONAL "       "
!             WO        VERTICAL VELOCITY COMPONENT(HYDROSTATIC VERSION)
!             WNO       VERTICAL VELOCITY COMPONENT(NONHYDROSTATIC VERSION)
!             THO       TEMPERATURE
!             SAO       SALINITY
!             PO        PRESSURE
!             PNH       NONHYDROSTATIC PRESSURE
!             TXO       ZONAL WIND STRESS
!             TYE       MERIDIONAL WIND STRESS
!             AVO       VARIABLE VERTICAL EDDY VISCOSITY
!                       (DEFINED ON SCALAR POINTS!!)
!             DVO       VARIABLE VERTICAL EDDY DIFFUSIVITY
!
!
! E=MERIDIONAL VELOCITY POINTS, O=ZONAL VELOCITY POINTS/SCALAR POINTS
!***********************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_LEVITUS
      USE MO_COMMOAU1
      USE MO_COMMOAU2
      USE MO_COMMOAU3
!
!
#ifdef SLOPECON_ADPO
      USE MO_COMMOBBL
#endif /*SLOPECON_ADPO*/
!
#ifdef MEAN
      USE MO_MEAN
#endif /*MEAN*/
!
      USE MO_ELICOM
      USE MO_PARA2
!
      USE MO_UNITS
      USE MO_OCTDIFF
      USE MO_OBCS

#ifdef NON_HYDROSTATIC
      use mo_petsc
#endif

#ifdef DIAGNOSIS
      use MO_DIAGNOSIS
#endif

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     DECLARATIONS
      INTEGER*8 IDATE
      INTEGER NACTYEAR
      REAL ABACK,CAULAPTS,CAULAPUV,CRELSAL,CRELTEM,CWT,DBACK,GFDL_DIFF
      REAL PIBOG,RRELSAL,RRELTEM
      REAL AWERT, DDPOMAX, DDPOMIN, EWERT
      REAL SUGG, SUZZ, SUZZ1, ZQQ, ZQQ1
      INTEGER I,ICONVA,IMEAN,IWIRB,J,K,M,N,NANF,NDAYS,NDTDAY
      INTEGER ICOU,IJJU,JMANF,JMEND,JMM,JMMM,LDTRUN
      INTEGER IERR, LDTDAY, LEN_FEB, LMON1, LMON2, LMONT, LREAD
      INTEGER LREADMAX, MONMON, NACYEAR, NDTYEAR, IOCAD, IMAL
      INTEGER L,IREADC
      REAL*8 ttts, tttr, ttt
!
! OtB  cannot use variables in more than one module
! ---------------------------------------------------------------------
!
!*    *NAMELIST* *OCECTL*   - Defines namelist parameters.
!                           - included in *ORCTM*.
!
!*     VARIABLE  TYPE        PURPOSE.
!      --------  ----        --------
!      *DT*      *REAL*      Ocean time step in sec
!
! ---------------------------------------------------------------------
!
      NAMELIST /OCECTL/                                                 &
     &         DT                                                       &
     &        ,CAULAPTS, CAULAPUV, CAH00, AUS                           &
     &        ,AV0,DV0,CWT,CSTABEPS,DBACK,ABACK,CRELSAL,CRELTEM         &
     &        ,ICONVA, IWIRB                                            &
     &        ,RRELTEM, RRELSAL                                         &
     &        ,CDVOCON,CAVOCON                                          &
     &        ,ISNFLG                                                   &
     &        ,H0, HMIN, ARMIN, ARMAX, HSNTOICE, SICTHMIN, SICE         &
     &        ,D3                                                       &
     &        ,IAUFR, IAUFW                                             &
     &        ,ISTART,I3DREST                                           &
     &        ,NYEARS, NMONTS, IMEAN, LY_END, LY_START,LM_START
      NAMELIST /OCEDZW/ DZW

      NAMELIST /NPROCS/ nprocx, nprocy ! declared in mo_parallel

! H.H--Open Boundary Condtion in Tidal Forcing Module
      NAMELIST /OCEOBC/ spongew, spongee, spongen, sponges              &
     &        ,Urelaxinner, Urelaxbound, Vrelaxinner, Vrelaxbound       &  
     &        ,TAO,tidalPeriod  
! H.H--Diagnoistic location index in Diagnoisis Module
      NAMELIST /OCEDIAG/ DMEAN, I_POINT_G, J_POINT_G 
!     Initialize MPI

      CALL p_start

!     SPECIFY LOGICAL I/O UNITS

      CALL SETUNITS

!  OPEN STD OUT FILE
!  This is needed for the MPI version since we do not want
!  the output of all processors intermixed (and basically
!  we want only the ouput of processor 0)

      CALL OPEN_STDOUT(io_stdout,'oceout')

!     Open input file and read the number of processors along
!     x- and y-direction

      IF(p_pe==p_io) THEN
        OPEN(IO_IN_OCTL,FILE='OCECTL',STATUS='UNKNOWN',                 &
     &          ACCESS='SEQUENTIAL',FORM='FORMATTED')
!
        READ(IO_IN_OCTL,NPROCS)
      ENDIF

!     Domain decomposition, setting local IE and JE

      CALL p_deco

!     Set some dependent parameters and allocate arrays

      CALL set_param1
      CALL alloc_mem_commo1
      CALL alloc_mem_commo2
      CALL alloc_mem_octdiff
      CALL alloc_mem_commoau2
      CALL alloc_mem_commoau3
#ifdef SLOPECON_ADPO
      CALL alloc_mem_commobbl
#endif /*SLOPECON_ADPO*/
#ifdef MEAN
      CALL alloc_mem_mean
#endif /* MEAN */

      CALL alloc_mem_elicom
      CALL alloc_mem_para2

      CALL alloc_mem_obcs
#ifdef DIAGNOSIS  
      CALL alloc_mem_diagnosis
#endif     
      CALL BELEG_ZERO
!  DEFAULT LENGTH OF INTEGRATION
!
       NYEARS=0
       NMONTS=1
!

      NANF=0
      NNNDT = 0
      NDAYS = 0
      DH = 0.
      AH = 0.
      FOUR=4.0
      TWO=2.
      ONE=1.0
      ZERO=0.0
      HALF=0.5
      FOURTH=0.25
      EIGHTH=0.125
      TENM4=1.E-4
      ALMZER=1.E-19
!
      monlen(:) = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
#ifdef YEAR360
      monlen(:) = 30
#endif /*YEAR360*/

#ifdef DEBUG_ONEDAY
      monlen(:) = 1
#endif /*DEBUG_ONEDAY*/

!
      WINTUR(1)=0.
      WINTUR(2)=1.E-3
      DO K=3,KEP
         WINTUR(K)=0.3*WINTUR(K-1)
      ENDDO
!
!     DZW(K) DECOFTES THE THICKNESS OF LAYER K IN METERS
!     SUM OF DZW(1-KE) IS THE BOTTOM OF THE MODEL, I.E. SOLID GROUND !!!
!     DZ(K) DECOFTES THE DISTANCE BETWEEN VECTOR POINTS BUT WILL BE
!     COMPUTED IN SBR BODEN FOR CONSISTENCY WITH W-POINTS
      DO K=1,KE
         SAF(K)=34.8
         TAF(K)=1.
         DZW(K)=0.
      ENDDO  
!
!----------------------------------------------------------------------
! DEFAULT PARAMETER SETTINGS - CAN BE OVERWRITEN BY THE NAMELIST
! ISTART :: START OPTIONS
! ISTART == 0: COMPLETELY NEW SETUP, 
!              topography read from anta and written to topo (!!! WARNING !!!)
!              start from climatology 
! ISTART == 1: new run, topography read from topo
!              start from horizontally uniform ts-profile (taf,saf)
! ISTART == 2: new run, topography read from topo
!              start from climatology
! ISTART == 3: continuing run (default)
         ISTART =3
         Z3DOC = 1  ! H.H --default for no restart
! I3DREST OPTIONS FOR 3-D RESTORING
! I3DREST == 0: NO RESTORING (DEFAULT)
! I3DREST == 1:  RESTORING to annual climatology
! I3DREST == 2:  RESTORING to monthly climatology
         I3DREST=0
!
!HH    DEFAULT TIME STEPS
         DT=1800.
!
!UWE   CONSTANTS FOR BIHARMONIC DIFFUSION
         CAULAPTS=0.0002
!
!UWE   CONSTANTS FOR HARMONIC DIFFUSION
         CAH00=0.
!
#ifdef ISOPYK
         CAULAPTS=0.
#endif
!
!UWE   CONSTANTS FOR BIHARMONIC FRICTION
         CAULAPUV=0.0045
!
!UWE   CONSTANTS FOR HARMONIC FRICTION
         AUS=3.E-6
!
#ifdef ISOPYK
         CAH00=1000.
#endif
!
!HH    CONSTANT FOR DIFFUSION IN OCTHER  
        DV0=0.5E-2
        AV0=0.5E-2
        CDVOCON=20.
        CAVOCON=0.
!
!HH    CONSTANT FOR WINDMIXING IN OCTHER 
        CWT=0.5E-3
        CSTABEPS=0.05
!
!HH    BACKGROUNDDIFFUSION
        ABACK=5.E-5
        DBACK=5.E-5
!
!HH    DEFAULT MEAN OUTPUT
       IMEAN=2
!
!HH    RELAXATION TIME SALINITY 
        CRELSAL = 3.E-7
!HH    RELAXATION TIME TEMERATURE
        CRELTEM = 0.
!JJ    DEFAULT END OF RUN IN YEAR
       LY_END=5000
!JJ    DEFAULT OFFSET FOR YEAR COUNTER (NEG. VALUE == NO ACTION TAKEN)
       LY_START=-999
       LM_START=-999
!----------------------------------------------------------------------
!     ICONVA : 1   WITH CONVECTIVE ADJUSTMENT
!              0    NO      "          "
      ICONVA=1
!-------------------------------------------------------
!     IWIRB  : 1   COMPUTE VARIABLE EDDY VISCOSITY
!              0   CONSTANT EDDY VISCOSITY
      IWIRB=1
!----------------------------------------------------------------------
! SWITCH FOR RESPECTIVE RESTART FILES
!
      IFLAG=1
!
!      Parameterization of ocean self-attraction and loading in
!      numerical models of the ocean circulation (0.08 0.10)
! cxe  rlsa = 0.085 for foreph effect and 0.00 for no foreph effect
       rlsa = 0.085
       rlsa1 = 1 - rlsa
!----------------------------------------------------------------------
!  READ OCEAN NAMELIST
!
!::
      IF(p_pe==p_io) THEN
!
        READ(IO_IN_OCTL,OCECTL)
        READ(IO_IN_OCTL,OCEDZW)
        READ(IO_IN_OCTL,OCEOBC)
#ifdef DIAGNOSIS 
        READ(IO_IN_OCTL,OCEDIAG)
#endif
        CLOSE(IO_IN_OCTL)
      ENDIF

      CALL p_bcast(DT,p_io)
      CALL p_bcast(CAULAPTS,p_io)
      CALL p_bcast(CAULAPUV,p_io)
      CALL p_bcast(CAH00,p_io)
      CALL p_bcast(AUS,p_io)
      CALL p_bcast(AV0,p_io)
      CALL p_bcast(DV0,p_io)
      CALL p_bcast(CWT,p_io)
      CALL p_bcast(CSTABEPS,p_io)
      CALL p_bcast(DBACK,p_io)
      CALL p_bcast(ABACK,p_io)
      CALL p_bcast(CRELSAL,p_io)
      CALL p_bcast(CRELTEM,p_io)
      CALL p_bcast(ICONVA,p_io)
      CALL p_bcast(IWIRB,p_io)
      CALL p_bcast(RRELTEM,p_io)
      CALL p_bcast(RRELSAL,p_io)
      CALL p_bcast(CDVOCON,p_io)
      CALL p_bcast(CAVOCON,p_io)
      CALL p_bcast(ISNFLG,p_io)
      CALL p_bcast(H0,p_io)
      CALL p_bcast(HMIN,p_io)
      CALL p_bcast(ARMIN,p_io)
      CALL p_bcast(ARMAX,p_io)
      CALL p_bcast(HSNTOICE,p_io)
      CALL p_bcast(SICTHMIN,p_io)
      CALL p_bcast(SICE,p_io)
      CALL p_bcast(D3,p_io)
      CALL p_bcast(IAUFR,p_io)
      CALL p_bcast(IAUFW,p_io)
      CALL p_bcast(NYEARS,p_io)
      CALL p_bcast(NMONTS,p_io)
      CALL p_bcast(IMEAN,p_io)
      CALL p_bcast(LY_END,p_io)
      CALL p_bcast(LY_START,p_io)
      CALL p_bcast(LM_START,p_io)
      CALL p_bcast(ISTART,p_io)
      CALL p_bcast(I3DREST,p_io)
      CALL p_bcast(dzw,p_io)
      CALL p_bcast(spongew,p_io)
      CALL p_bcast(spongee,p_io)
      CALL p_bcast(spongen,p_io)
      CALL p_bcast(sponges,p_io)
      CALL p_bcast(Urelaxinner,p_io)
      CALL p_bcast(Urelaxbound,p_io)
      CALL p_bcast(Vrelaxinner,p_io)
      CALL p_bcast(Vrelaxbound,p_io)
      CALL p_bcast(TAO,p_io)                              
      CALL p_bcast(tidalPeriod,p_io)
#ifdef DIAGNOSIS
      CALL p_bcast(DMEAN,p_io)
!      CALL p_bcast(POINT,p_io)
      CALL p_bcast(I_POINT_G,p_io)
      CALL p_bcast(J_POINT_G,p_io)          
#endif
!
!
      DTI=1./DT
      NDTDAY=NINT(86400./DT)
! Write final namelist parameters
!
      WRITE(IO_STDOUT,*)' SECONDS PER TIMESTEP   (DT): ',DT
      WRITE(IO_STDOUT,*)' TIME STEPS PER DAY (NDTDAY): ',NDTDAY
      WRITE(IO_STDOUT,*)'                  (CAULAPUV): ',CAULAPUV
      WRITE(IO_STDOUT,*)'                  (CAULAPTS): ',CAULAPTS
      WRITE(IO_STDOUT,*)'                       (AUS): ',AUS
      WRITE(IO_STDOUT,*)'                     (CAH00): ',CAH00
      WRITE(IO_STDOUT,*)'                       (AV0): ',AV0
      WRITE(IO_STDOUT,*)'                       (DV0): ',DV0
      WRITE(IO_STDOUT,*)'                       (CWT): ',CWT
      WRITE(IO_STDOUT,*)'                  (CSTABEPS): ',CSTABEPS
      WRITE(IO_STDOUT,*)'                     (DBACK): ',DBACK
      WRITE(IO_STDOUT,*)'                     (ABACK): ',ABACK
      WRITE(IO_STDOUT,*)'                   (CRELSAL): ',CRELSAL
      WRITE(IO_STDOUT,*)'                   (CRELTEM): ',CRELTEM
      WRITE(IO_STDOUT,*)'                   (CDVOCON): ',CDVOCON
      WRITE(IO_STDOUT,*)'                   (CAVOCON): ',CAVOCON
      WRITE(IO_STDOUT,*)'                     (IMEAN): ',IMEAN
      WRITE(IO_STDOUT,*)'                  (LY_START): ',LY_START
      WRITE(IO_STDOUT,*)'                  (LM_START): ',LM_START
      WRITE(IO_STDOUT,*)'                    (LY_END): ',LY_END
      WRITE(IO_STDOUT,*)'                    (ISTART): ',ISTART
      WRITE(IO_STDOUT,*)'                    (I3DREST): ',I3DREST
#ifdef OBC_TIDE           
      WRITE(IO_STDOUT,*)'                    (Spongew): ',Spongew
      WRITE(IO_STDOUT,*)'                    (Spongee): ',Spongee
      WRITE(IO_STDOUT,*)'                    (Spongen): ',Spongen
      WRITE(IO_STDOUT,*)'                    (Sponges): ',Sponges
      WRITE(IO_STDOUT,*)'                (Urelaxinner): ',Urelaxinner
      WRITE(IO_STDOUT,*)'                (Urelaxbound): ',Urelaxbound
      WRITE(IO_STDOUT,*)'                (Vrelaxinner): ',Vrelaxinner
      WRITE(IO_STDOUT,*)'                (Vrelaxbound): ',Vrelaxbound
#endif      
#ifdef OBC_SPONGE_EXP      
      WRITE(IO_STDOUT,*)'                        (TAO): ',TAO
#endif      
#ifdef OBC_TIDE                
      WRITE(IO_STDOUT,*)'                (tidalPeriod): ',tidalPeriod
#endif
#ifdef DIAGNOSIS
      WRITE(IO_STDOUT,*)'                      (DMEAN): ',DMEAN
!      WRITE(IO_STDOUT,*)'                      (Point): ',POINT
      WRITE(IO_STDOUT,*)'                    (I_POINT_G): ',I_POINT_G      
      WRITE(IO_STDOUT,*)'                    (J_POINT_G): ',J_POINT_G 
#endif                
!
      IF (ISTART .LT. 0 .OR. ISTART .GT. 3)                             &
          WRITE(IO_STDOUT,*)'ISTART NOT SUPPORTED!!!!!'
      IF (I3DREST .LT. 0 .OR. I3DREST .GT. 2)                           &
          WRITE(IO_STDOUT,*)'I3DREST NOT SUPPORTED!!!!!'     
!
      IF (I3DREST .GT. 0 .OR. ISTART .LT. 3) THEN
       CALL init_levitus(ie,je,ke)
      ENDIF
#ifdef OBC_ETA_FORC
      WRITE(IO_STDOUT,*)'OPEN BOUNDARY FORCED BY WATER LEVEL !!!!!'
#endif     
 
#ifdef OBC_UV_FORC
      WRITE(IO_STDOUT,*)'OPEN BOUNDARY FORCED BY WATER VELOCITY !!!!!'
#endif     
  
#ifdef DIAGNOSIS  
      WRITE(IO_STDOUT,*)'DIAGNOSTIC MODULE SWITCH ON !!!!!'
#endif     
!
!  DEFAULT ...
       IAUFR=1
       IAUFW=1
      IF (ISTART .LT. 3) IAUFR=0

      AULAPTS=CAULAPTS*DT/3600.
      AULAPUV=CAULAPUV*DT/3600.
      AH00=CAH00/4.E5
      WT=CWT/(6.**3)
!
!HH   CHECK LAYER THICKNESS
      DO K=1,KE
        IF(DZW(K).EQ.0.) THEN
          WRITE(IO_STDOUT,*)' LAYER: ',K,' THICKNESS IS ZERO !!!'
          CALL ABSTURZ 
        ELSE
          WRITE(IO_STDOUT,*) K, ' LAYERTHICKNESS    (DZW): ',DZW(K)
        ENDIF
      ENDDO
!
!HH   BACKGROUNDDIFFUSION
      DO K=1,KEP
        ABACKV(K) = ABACK
        DBACKV(K) = DBACK 
      ENDDO
!
      TIESTW(1) = 0.
      DO K=1,KE
      TIESTW(K+1)  = TIESTW(K) + DZW(K)
      ENDDO
      PI=4.*ATAN(1.)
      PIBOG=180./PI
!
#ifdef DBACKPROFIL
      IDBACK=0
      DO K=1,KEP
         ABACKV(K) = ABACK
         DBACKV(K) = DBACK + (1-IDBACK) * 1.E-4 *                       &
     &               (0.8 + 1.05/PI*ATAN(4.5*1.E-3*(TIESTW(K)-2500.)))* &
     &               (0.5+SIGN(0.5,TIESTW(K)-500.)) *                   &
     &               SQRT(ABS((TIESTW(K)-500.)/(3500.-500.)))
      ENDDO
#endif
!
       DO K=1,KEP
         GFDL_DIFF = 1.E-4 *                                            &
     &            (0.8 + 1.05/PI*ATAN(4.5*1.E-3*(TIESTW(K)-2500.)))
#ifdef DBACKGFDL
         DBACKV(K) = GFDL_DIFF
#endif
       ENDDO
!
#ifdef DBACKGFDL2
         DBACKV(K) = 0.5*(DBACK+GFDL_DIFF)
#endif
       DO K=1,KEP
         WRITE(IO_STDOUT,6002)'BACKGROUND DIFFUSIVITY AT '              &
     &                        ,INT(TIESTW(K))                           &
     &         ,'M : HOPE : ',DBACKV(K),' GFDL : ',GFDL_DIFF
         WRITE(IO_STDOUT,6002)'BACKGROUND VISCOSITY AT ',INT(TIESTW(K)) &
     &         ,'M : HOPE : ',ABACKV(K)
       ENDDO
!
6002   FORMAT(1X,A27,I5,A10,E12.3,A8,E12.3)
!
!HH    RELAXATION TIME SALINITY 
      RELSAL = CRELSAL*20./DZW(1)
      IF (RELSAL.GT. ALMZER) then 
         WRITE(IO_STDOUT,26668)1./(RELSAL*24.*3600.)
      ELSE
         WRITE(IO_STDOUT,*) 'SSS relaxation switched off !!'
      ENDIF
26668 FORMAT('  RELAXATION TIME SALINITY COUPLING : ',F10.2,' DAYS')
!HH    RELAXATION TIME TEMPERATURE
      RELTEM = CRELTEM*20./DZW(1)
      IF (reltem.GT. ALMZER) then 
         WRITE(IO_STDOUT,26669)1./(RELTEM*24.*3600.)
      ELSE
         WRITE(IO_STDOUT,*) 'SST relaxation switched off !!'
      ENDIF
26669 FORMAT('  RELAXATION TIME SST COUPLING : ',F10.2,' DAYS')
!
!
!----------------------------------------------------------------------
! OPEN FILES
!
      IF(p_pe==p_io) THEN
        OPEN(IO_IN_ARCG,FILE='arcgri',                                  &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
!
        OPEN(IO_IN_INIT,FILE='./INITEM',STATUS='UNKNOWN',               &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        OPEN(IO_IN_INIS,FILE='./INISAL',STATUS='UNKNOWN',               &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        IF (RELSAL .GT. ALMZER) THEN
        OPEN(IO_IN_SURS,FILE='./SURSAL',STATUS='UNKNOWN',               &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        ENDIF
        IF (RELTEM .GT. ALMZER) THEN
        OPEN(IO_IN_SURT,FILE='./SURTEM',STATUS='UNKNOWN',               &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        ENDIF
!
        OPEN(IO_IN_GWIX,FILE='./GIWIX',STATUS='UNKNOWN',                &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        OPEN(IO_IN_GWIY,FILE='./GIWIY',STATUS='UNKNOWN',                &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        OPEN(IO_IN_GTEM,FILE='./GITEM',STATUS='UNKNOWN',                &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        OPEN(IO_IN_GPRE,FILE='./GIPREC',STATUS='UNKNOWN',               &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        OPEN(IO_IN_GSWR,FILE='./GISWRAD',STATUS='UNKNOWN',              &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        OPEN(IO_IN_GTDE,FILE='./GITDEW',STATUS='UNKNOWN',               &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        OPEN(IO_IN_GU10,FILE='./GIU10',STATUS='UNKNOWN',                &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
        OPEN(IO_IN_RVAL,FILE='./GIRIV',STATUS='UNKNOWN',                &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
#ifdef FORCE_SLP
        OPEN(IO_IN_GSLP,FILE='./GIPRESS',STATUS='UNKNOWN',              &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
#endif
#ifdef FORCE_DWLW
        OPEN(IO_IN_GLWR,FILE='./GIDWLW',STATUS='UNKNOWN',               &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
#else
        OPEN(IO_IN_GCLO,FILE='./GICLOUD',STATUS='UNKNOWN',              &
     &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
#endif /*FORCE_DWLW*/
!
#ifdef FRAGMENT
     call system("mkdir -p fort.71")
     call system("mkdir -p fort.72")
     call system("mkdir -p fort.73")
     call system("mkdir -p fort.74")     
     call system("mkdir -p fort.76")      
     call system("mkdir -p fort.82")     
     call system("mkdir -p fort.146") 
     call system("mkdir -p restart_file")       
#ifdef BAROTROPIC
     call system("mkdir -p fort.63")
     call system("mkdir -p fort.64") 
#endif    
 
#ifdef NON_HYDROSTATIC
     call system("mkdir -p fort.77") 
#endif
  
#endif /*FRAGMENT*/

#ifdef DIAGNOSIS
     call system("mkdir -p diagnosis") 
#endif
      ENDIF ! p_pe == p_io

!-------------------------------------------------------------------
!
!     DU/DT - F*V = -G*( STABN* (DZ/DX)(NEW) + STABO* (DZ/DX)(OLD) )
!
!     DZ/DT = H*( CONN* (DU/DX)(NEW) + CONO* (DU/DX)(OLD) + ... )
!
      STABN=0.6
      STABO=1.-STABN
!
      CONN=0.50
      CONO=1.-CONN
!
      WRITE(IO_STDOUT,*)' STABN = ',STABN,' CONN = ',CONN
!
!-------------------------------------------------------------------
!  OMEGA  :  ANGULAR VELOCITY OF OUR NICE BLUE PLANET
!  RADIUS :  RADIUS OF THE ABOVE MENTIONED PLANET
!  G      :  GRAVITATIONAL ACCELERATION ( CONSTANT 9.806 M/S**2
!                TO ACCOUNT FOR THE EXISTENCE OF MINI BLACK HOLES )
!  ROCP   :  HEAT CAPACITY PER CUBICMETER
!  ROCD   :  WIND STRESS DRAG COEFFICIENT
!------------------------------------------------------------------
      PI=4.*ATAN(1.)
! cxe siday= 23.*3600.+56.*60.+4.091
! cxe OMEGA= 2.*PI/siday
! cxe gamma0=0.69 for foreph effect and 1.00 for no foreph effect
      gamma0=0.69
      OMEGA=7.292E-5
      RADIUS=6371.E+3
      G=9.806
      GHN=G*STABN
      GHO=G*STABO
      ROCP=4.E06
      ROCD=1.2*1.5E-3
!
!-------------------------------------------------------------------
! PARAMETER ICE MODEL 
!
      DO 139 J=1,JE
      DO 139 I=1,IE
         TICEO(I,J)=0.
         ACLO(I,J)=.7
         PAO(I,J)=101300.
         RPRECO(I,J)=0.
         FRSE(I,J)=0.
139   CONTINUE
!
      ISNFLG=1
!
      ALBI=0.75
      ALBM=0.66
      ALBW=0.10
      ALBSN=0.85
      ALBSNM=0.75
#ifdef ALBMSN07
      ALBSN=0.82
      ALBSNM=0.7
      ALBM=0.63
#endif /*ALBMSN07*/
#ifdef ALBOMIP
!SJM TUNE ALBEDOS FOR OMIP WITH BERYLIAND
      ALBI=0.75    !ICE < 0
      ALBM=0.70    !ICE > 0
      ALBW=0.10    !WATER
      ALBSN=0.85   !SNOW < 0
      ALBSNM=0.70  !SNOW > 0
#endif /*ALBOMIP*/
#ifdef ALBNCEP
!HH   TUNE ALBEDOS FOR NCEP WITH KOCH
      ALBI=0.76
      ALBM=0.71
      ALBW=0.1
      ALBSN=0.86
      ALBSNM=0.72
#endif
#ifdef ALBMELTHI
      ALBSNM=0.8
      ALBSN=0.87
      ALBI=0.78
      ALBM=0.7
#endif /*ALBMELTHI*/
      TMELT=273.16
      TFREZ=-1.9
      CC=4.2E6
      CW=0.0045
      CLO=3.02E8
      CLB=2.70E8
      RHOAIR=1.3E+00
      RHOWAT=1.025E+03
      RHOICE=0.91E+03
      RHOSNO=0.33E+03
      RHOICWA=RHOICE/RHOWAT
      RHOSNWA=RHOSNO/RHOWAT
      RHOSNIC=RHOSNO/RHOICE
      CON=2.1656
      CONSN=0.31
      H0=0.5
      ARMIN=0.15
      ARMAX=1.0
      HMIN=0.05
!UWE  MAXIMUM ICE THICKNESS FOR SNOW-->ICE CONVERSION
!      HSNTOICE = 17.
       HSNTOICE = 0.45 * DZW(1)
!UWE  MINIMUM ICE THICKNESS IN NEW ICE GROWTH
      SICTHMIN=0.5
!
      VAPL=2.5E6
      SUBL=2.834E6
      SICE=5.0
      D3=5.5E-08
!
!JJ OLD VALUES UP TO HOPS 62!
      D1=RHOAIR*1004.*1.75E-3
      D2I=RHOAIR*SUBL*1.75E-3
      D2W=RHOAIR*VAPL*1.75E-3
!
#ifdef DRAGGILL
!HH   VALUES FROM GILL ATMOSPHERE OCEAN DYNAMICS
      D1=RHOAIR*1004.*1.1E-3
      D2I=RHOAIR*SUBL*1.5E-3
      D2W=RHOAIR*VAPL*1.5E-3 
#endif
#ifdef DASILVA
      D1=RHOAIR*1004.
      D2I=RHOAIR*SUBL
      D2W=RHOAIR*VAPL
#endif
#ifdef BULK_KARA
      D1=RHOAIR*1004.67
      D2I=RHOAIR*SUBL
      D2W=RHOAIR*VAPL
#endif
!
      TFREEZ=TFREZ
      ENTMEL=320.E6
      SICHEC=2.
      HICCE=2.
      HICCP=20.
      HTH00=0.5
      STEBOL=5.67E-8
!
      KB=KBB
      KM=KB+1
      KBM=KB+KM
!
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG,DEPTO)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, DLXP)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, DLXU,1)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, DLXV,2)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, DLYP)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, DLYU,1)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, DLYV,2)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, FTWOU,1)
       IF(p_pe==p_io) READ(IO_IN_ARCG) IBLA
       CALL read_slice(IO_IN_ARCG, FTWOV,2)
       IF(p_pe==p_io) CLOSE(IO_IN_ARCG)

!SL 
!SL GRID DEFORMATION
!SL 
!SL 
       CALL bounds_exch(DLXP,DLYP)
       CALL bounds_exch(DLXV,DLYV)
       CALL bounds_exch(DLXU,DLYU)
!
      CALL BELEG
!
      CALL BODEN
!
!-------------------------------------------------------------------
!OtB LMONTS used in LEVIRE but not defined yet
      LMONTS=0
      IF (I3DREST .GT. 0) THEN
! READ 3D LEVITUS DATA FOR USE IN 3D RESTORING
!
      IF(p_pe==p_io) THEN
        REWIND(IO_IN_INIT)
        REWIND(IO_IN_INIS)
      ENDIF
      CALL LEVIRE(-1)
      ENDIF ! I3DREST
!
!-------------------------------------------------------------------
! SPECIFY CORIOLIS PARAMETER ETC.
!
      CALL CORIOL
!
!-------------------------------------------------------------------
!
      IF (ISTART .LT. 3) THEN
      DO K=1,KE
       DO J=1,JE
        DO I=1,IE
         THO(I,J,K)=TLEVI(I,J,K)
         SAO(I,J,K)=SLEVI(I,J,K)
        ENDDO
       ENDDO
      ENDDO
      ENDIF !ISTART
!
#ifdef SLOPECON_ADPO
!JJ PREPARE FOR BOTTOM BOUNDARY PARAMETRIZATIONS
      CALL FINDBOT
      CALL FINDALFA
#endif /*SLOPECON_ADPO*/
!
!-------------------------------------------------------------------
!  PROVIDE MATRIX FOR BAROTROPIC MODE
!     UWE  ADD ITERATION OF MATRIX
!
      IF(IELIMI.GE.1) THEN
         CALL TRIAN
      ELSE
         CALL ITPREP
         CALL TROTEST
      ENDIF
!
!
!-------------------------------------------------------------------
!
!:: INCLUDE WRITEOUT OF CPP OPTIONS AND PARAMETER SETTINGS
!OtB #include "CPPOUTPUT
      CALL CPPOUTPUT
!
!-----------------------------------------------------------------------
!     
!     TIMESTEPPING IS BASED ON 3 LOOPS
!          OUTER LOOP : DO 1000  LYEAR = LYEAR1, LYEAR2
!        1.INNER LOOP : DO 1100  LMONT = 1 OR LMONT1, LMONT2 OR 12
!        2.INNER LOOP : DO 1010  LDAY  = 1,MONLEN(LMONTS)
!        3.INNER LOOP : DO 1001  LDTDAY  = 1,NDTDAY
! ADDITIONAL COUNTERS
!        LYEARS = ACTUAL YEAR
!        LMONTS = ACTUAL MONTH
!        LDAYS  = ACTUAL DAY
!     
!     TIME COUNTER ARE :
!               LDT       : TIME STEP COUNTER FOR EXPERIMENT
!               LDTRUN    : TIME STEP COUNTER FOR RUN
!               LDTYEAR   : TIME STEP COUNTER FOR YEAR
!               LDTMONTH  : TIME STEP COUNTER FOR MONTH
!               NDTMONTH  : TIME STEPS IN ACTUAL MONTH
!               NDTDAY  : TIME STEPS PER DAY
      WRITE(IO_STDOUT,*)'START FROM PREVIOUS CALCULATION YES=1/NO=0 : ' &
     &                  ,IAUFR
      WRITE(IO_STDOUT,*)'WRITE RESTART FILE              YES=1/NO=0 : ' &
     &                  ,IAUFW
      WRITE(IO_STDOUT,*)'THIS JOB WILL TRY TO INTEGRATE '               &
     &                  ,NYEARS,' YEARS AND ',NMONTS,' MONTHS'
!
!
!-----------------------------------------------------------------------
!
!     START FROM STATUS LEVITUS OR HORIZONTAL STRATIFICATION
! 
!  CALL LEVIRE WITH ARGUMENT =  0 : HORIZONTAL STRATIFICATION
!  CALL LEVIRE WITH ARGUMENT = 13 : 3D STRATIFICATION
!  CALL LEVIRE WITH ARGUMENT = -2 : SURFACE SALINITY
!
      IF (RELSAL .GT. ALMZER) THEN
      CALL LEVIRE(-2)
      ENDIF
      IF (RELTEM .GT. ALMZER) THEN
      CALL LEVIRE(-3)
      ENDIF
! 
      IF(IAUFR.EQ.0) THEN
        IF (ISTART .EQ. 0) CALL LEVIRE(13)
        IF (ISTART .EQ. 1) CALL LEVIRE(0)
        IF (ISTART .EQ. 2) CALL LEVIRE(13)
!
!  SET LOGICAL UNIT FOR RESTART FILE; FOR FIRST RUN WRITE TO IO_IN_Z370
        IUNIT = IO_IN_Z380
        IFLAG = -1
!     WRITE RESTART FILES WITH INITIAL STRATIFICATION                   
!     MAKE SURE THAT BOTH RESTART FILES EXIST ON EXIT
#ifndef FRAGMENT
        CALL AUFW
        CALL AUFW
#else
        CALL AUFW_FRAGMENT
        CALL AUFW_FRAGMENT        
#endif /*FRAGMENT*/        
      ENDIF
 
!
!-----------------------------------------------------------------------
!
!     START FROM RESTART FILES Z37000 OR Z38000
!
      IF (IAUFR .EQ. 1) THEN
        Z3DOC = 3 ! for read_slice      
#ifndef FRAGMENT
        CALL AUFR  
#else
        CALL AUFR_FRAGMENT
#endif /*FRAGMENT*/      
        DO J=1,JE
          DO I=1,IE
            SICUDO(I,J)=SICUO(I,J)
            SICVDE(I,J)=SICVE(I,J)
          ENDDO
        ENDDO
      ENDIF
      
      Z3DOC = 1 ! read_slice-no restart      
!-----------------------------------------------------------------------
!
      CALL WRTE_GRIDINFO
      CALL INI_OBCS
#ifdef NON_HYDROSTATIC
      call petsc_start
#endif

#ifdef DIAGNOSIS 
      call INI_diagnosis
#endif      
!
!-----------------------------------------------------------------------
! SET LIMITS OF YEAR/MONTH TIME STEPPING LOOP
!
      LYEARS = 0
      LMONTS = 0
      LDAYS=0 
      LDT = 0
      LDTRUN = 0
!
      WRITE(IO_STDOUT,*)'TIME COUNTER INITIALIZED ' 
!
      LYEAR1 = LY_START
      LMONT1 = LM_START
!
      IF (LMONT1 .GT. 12) THEN
         LMONT1=1
         LYEAR1=LYEAR1+1
      ENDIF
!
      LMONT2 = LMONT1 + 12*NYEARS + NMONTS - 1
      LYEAR2 = LYEAR1 + (LMONT2-1)/12
      LMONT2 = MOD(LMONT2-1,12)+1
!
!--------------------------------------------------------------------
!
#ifdef RYEAR
      WRITE(IO_STDOUT,*)                                                &
     &'THIS RUN ASSUMES REALISTIC YEAR WITH LEAP FEBRUARIES'
#else
!
#ifdef YEAR360
      WRITE(IO_STDOUT,*)                                                &
     &'THIS RUN ASSUMES IDEALIZED YEARS (12*30 DAYS EACH)'
#else
      WRITE(IO_STDOUT,*)                                                &
     &'THIS RUN ASSUMES IDEALIZED YEARS WITH 365 DAYS'
#endif /*YEAR360*/
!
#endif /*RYEAR*/
!
      WRITE(IO_STDOUT,*)                                                &
     &    ' INTEGRATION PERIOD IS YEAR ',LYEAR1,' TO ',LYEAR2
!
      DO 1000 LYEAR=LYEAR1,LYEAR2
!
         LYEARS=LYEAR
         LDTYEAR = 0 
!
         WRITE(IO_STDOUT,*) ' TO BE CALCULATED NOW : YEAR = ',LYEARS
         WRITE(IO_STDOUT,*) ' Grids scheme : IE -> ',IE, 'JE -> ', JE
!
#ifdef RYEAR
         NACTYEAR=LYEAR
         CALL MONLEN_FEB(NACTYEAR,LEN_FEB)
         MONLEN(2) = LEN_FEB
         WRITE(IO_STDOUT,*)                                             &
     & 'FEBRUARY OF YEAR ',NACTYEAR,' HAS ',MONLEN(2),' DAYS'
#endif /*RYEAR*/
!
!:: LREADMAX IS THE MAX. NUMBER OF READ OPERATION FOR SURFACE FORCING
!:: LREAD COUNTS READS OF THE FORCING
       LREADMAX=0
       LREAD=0
!
#ifdef NOSKIP
! initial file start from Jan 1st ! peterspy
       DO I=1,12
         LREADMAX=LREADMAX+MONLEN(I)
       ENDDO
! POSITIONING NEAR-SURFACE FORCING DATA AT BEGINNIG OF ACTUAL MONTH
! 
       IF ( LMONT1 .GT. 1 ) THEN
!
       DO LMONT=1,LMONT1-1
          DO LDAY=1,MONLEN(LMONT)
             LREAD=LREAD+1
             IF(p_pe==p_io) READ(IO_IN_GTEM)IBLA
             CALL read_slice(IO_IN_GTEM,TAFO1)
             IF(p_pe==p_io) READ(IO_IN_GWIX)IBLA
             CALL read_slice(IO_IN_GWIX,TXO1,1)
             IF(p_pe==p_io) READ(IO_IN_GWIY)IBLA
             CALL read_slice(IO_IN_GWIY,TYE1,2)
             IF(p_pe==p_io) READ(IO_IN_GPRE)IBLA
             CALL read_slice(IO_IN_GPRE,FPREC1)
             IF(p_pe==p_io) READ(IO_IN_GSWR)IBLA
             CALL read_slice(IO_IN_GSWR,FSWR1)
             IF(p_pe==p_io) READ(IO_IN_GTDE)IBLA
             CALL read_slice(IO_IN_GTDE,FTDEW1)
             IF(p_pe==p_io) READ(IO_IN_GU10)IBLA
             CALL read_slice(IO_IN_GU10,FU101)
             IF(p_pe==p_io) READ(IO_IN_RVAL)IBLA
             CALL read_slice(IO_IN_RVAL,GIRIV)
#ifdef FORCE_SLP
             IF(p_pe==p_io) READ(IO_IN_GSLP)IBLA
             CALL read_slice(IO_IN_GSLP,FSLP1)
#endif /*FORCE_SLP*/
#ifdef FORCE_DWLW
             IF(p_pe==p_io) READ(IO_IN_GLWR)IBLA
             CALL read_slice(IO_IN_GLWR,FCLOU1)
#else
             IF(p_pe==p_io) READ(IO_IN_GCLO)IBLA
             CALL read_slice(IO_IN_GCLO,FCLOU1)
#endif /*FORCE_DWLW*/
             WRITE(IO_STDOUT,*)'SKIP LESEN: ',LMONT,LDAY,IBLA
          ENDDO
       ENDDO
       DO J=1,JE
          DO I=1,IE
             TAFO(I,J)=TAFO1(I,J)
             FPREC(I,J)=FPREC1(I,J)
             FSWR(I,J)=FSWR1(I,J)
             FTDEW(I,J)=FTDEW1(I,J)
             FU10(I,J)=FU101(I,J)
             FCLOU(I,J)=FCLOU1(I,J)
#ifdef FORCE_SLP
             PAO(I,J)=FSLP1(I,J)
#endif /*FORCE_SLP*/
          ENDDO
       ENDDO
       DO J=1,JE
          DO I=I_start,IE
             TXO(I,J)=TXO1(I,J)
          ENDDO
       ENDDO
       DO J=J_start,JE
          DO I=1,IE
             TYE(I,J)=TYE1(I,J)
          ENDDO
       ENDDO

       WRITE(IO_STDOUT,*) 'POINTER IN FORCING DATA'                     &
     &    ,' IS SET TO  START OF MONTH ',LMONT1
 
       ENDIF
#endif /*NOSKIP*/
!
!--------------------------------------------------------------------
!  SET NUMBER OF TIME STEPS PER YEAR
!
       ndtyear=0
       DO i=1,12
         ndtyear=ndtyear+monlen(i)*ndtday
       ENDDO
!
!  SET MONTH LOOP LIMITS FOR EACH YEAR
!
       LMON1=1
       LMON2=12
       IF(LYEAR1.EQ.LYEAR2) THEN
          LMON1=LMONT1
          LMON2=LMONT2
       ELSE
          IF(LYEARS.EQ.LYEAR2) THEN
             LMON1=1
             LMON2=LMONT2
          ENDIF          
          IF(LYEARS.EQ.LYEAR1) THEN
             LMON1=LMONT1
             LMON2=12
          ENDIF
       ENDIF
       WRITE(IO_STDOUT,*) 'MONTHLY LOOP FROM MONTH ',LMON1,             &
     &                        ' TO',LMON2
!
!:: MONTH LOOP
       DO 1100 LMONT=LMON1,LMON2
!
          LMONTS=LMONT
          LDTMONTH = 0
          NDTMONTH = MONLEN(LMONTS) * NDTDAY
!
          WRITE(IO_STDOUT,18823) LYEARS,LMONTS
18823     FORMAT(' TO BE CALCULATED NOW : YEAR ',I4,' MONTH ',I4)
!
!
#ifdef RESTORE_MON
!:: REWIND SURSAL TO BEGIN OF YEAR
       IF(p_pe==p_io) THEN
         IF (RELSAL .GT. ALMZER) REWIND(IO_IN_SURS)
         IF (RELTEM .GT. ALMZER) REWIND(IO_IN_SURT)
       ENDIF
!:: READ APPROPRIATE SURTEM and SURSAL FOR MONTHLY RESTORING
       WRITE(IO_STDOUT,*)                                               &
     &        'READING MONTHLY SURFACE TEM and SAL IN MONTH ',LMONT
       DO IREADC=1,LMONT
         IF (RELSAL .GT. ALMZER) CALL LEVIRE(-2)
         IF (RELTEM .GT. ALMZER) CALL LEVIRE(-3)
       ENDDO
#endif /*RESTORE_MON*/
!::
       IF (I3DREST .GT. 1) THEN
!:: REWIND INISAL/INITEM TO BEGIN OF YEAR
         IF(p_pe==p_io) THEN
           REWIND(IO_IN_INIT)
           REWIND(IO_IN_INIS)
         ENDIF
!:: READ APPROPRIATE INISAL/INITEM FOR MONTHLY RESTORING
       WRITE(IO_STDOUT,*)                                               &
     &        'READING MONTHLY LEVITUS FIELDS IN MONTH ',LMONT
       DO IREADC=1,LMONT
         CALL LEVIRE(-1)
       ENDDO
       ENDIF ! I3DREST
!
      DO 1010 LDAY=1,MONLEN(LMONTS)
      LDAYS=LDAY
      DO 1001 LDTDAY = 1,NDTDAY
      LDTDAYC = LDTDAY
#if defined (TIMECHECK)
      ttts = MPI_Wtime()
#endif
!
!*       6.  TIMESTEP THE OCEAN.
!            -------- --- ------
      LDT = LDT + 1
      LDTRUN = LDTRUN + 1
      LDTYEAR = LDTYEAR + 1
      LDTMONTH =LDTMONTH + 1

      if (p_pe==p_io) then
        write(0,'(a,/,2(a,i10,/),a,i2,a,i2.2,a,i4.4)') &
             '--------------------------------------------', &
             'current timestep (run): ', ldtrun, &
             'current timestep (day): ', ldtday, &
             'current date:    ', lday, '.', lmonts, '.', lyear 
      endif

!JJ READING OF FORCING FIELDS EVERY DAY
! FOR SURFACE FORCING
      IF (MOD(LDTDAY-1,NDTDAY) .EQ. 0) THEN
          WRITE(IO_STDOUT,*)'READ SURFACE FORCING AT TIMESTEP ',LDTDAY  &
     &              ,LDAYS,LMONTS,LYEARS
! COUNT READ OPERATIONS
          LREAD=LREAD+1
          IF (LREAD.GT.LREADMAX) THEN
           WRITE(IO_STDOUT,*)'POSITION THE FORCING TO 1 JAN ...'
           IF(p_pe==p_io) THEN
             REWIND (IO_IN_GTEM)
             REWIND (IO_IN_GWIX)
             REWIND (IO_IN_GWIY)
             REWIND (IO_IN_GPRE)
             REWIND (IO_IN_GSWR)
             REWIND (IO_IN_GTDE)
             REWIND (IO_IN_GU10)
             REWIND (IO_IN_RVAL)
#ifdef FORCE_SLP
             REWIND (IO_IN_GSLP)
#endif /*FORCE_SLP*/
#ifdef FORCE_DWLW
             REWIND (IO_IN_GLWR)
#else
             REWIND (IO_IN_GCLO)
#endif /*FORCE_DWLW*/
           ENDIF
           LREAD=1
          ENDIF
          WRITE(IO_STDOUT,*)'STEP LESEN: ',LMONTS,LDAYS,IBLA

          IF(p_pe==p_io) READ(IO_IN_GTEM)IBLA
          CALL read_slice(IO_IN_GTEM,TAFO1)
          IF(p_pe==p_io) READ(IO_IN_GWIX)IBLA
          CALL read_slice(IO_IN_GWIX,TXO1,1)
          IF(p_pe==p_io) READ(IO_IN_GWIY)IBLA
          CALL read_slice(IO_IN_GWIY,TYE1,2)
          IF(p_pe==p_io) READ(IO_IN_GPRE)IBLA
          CALL read_slice(IO_IN_GPRE,FPREC1)
          IF(p_pe==p_io) READ(IO_IN_GSWR)IBLA
          CALL read_slice(IO_IN_GSWR,FSWR1)
          IF(p_pe==p_io) READ(IO_IN_GTDE)IBLA
          CALL read_slice(IO_IN_GTDE,FTDEW1)
          IF(p_pe==p_io) READ(IO_IN_GU10)IBLA
          CALL read_slice(IO_IN_GU10,FU101)
          IF(p_pe==p_io) READ(IO_IN_RVAL)IBLA
          CALL read_slice(IO_IN_RVAL,GIRIV)
#ifdef FORCE_SLP
          IF(p_pe==p_io) READ(IO_IN_GSLP)IBLA
          CALL read_slice(IO_IN_GSLP,FSLP1)
#endif /*FORCE_SLP*/
#ifdef FORCE_DWLW
          IF(p_pe==p_io) READ(IO_IN_GLWR)IBLA
          CALL read_slice(IO_IN_GLWR,FCLOU1)
          WRITE(IO_STDOUT,*) 'ATTN: OPTION DW forcing DWLW'
#else
          IF(p_pe==p_io) READ(IO_IN_GCLO)IBLA
          CALL read_slice(IO_IN_GCLO,FCLOU1)
          WRITE(IO_STDOUT,*) 'ATTN: OPTION DW forcing CLOUD'

#endif /*FORCE_DWLW*/
!
!      UPDATE FORCING ONCE PER DAY
!
#ifdef SW089
          WRITE(IO_STDOUT,*) 'ACHTUNG!!! NCEP SW-RAD TIMES 0.89'
#endif
!
          DO J=1,JE
          DO I=1,IE
            TAFO(I,J)=TAFO1(I,J)
            FPREC(I,J)=FPREC1(I,J)
#ifdef SW089
            FSWR(I,J)=FSWR1(I,J)*0.89
#else
            FSWR(I,J)=FSWR1(I,J)
#endif
            FTDEW(I,J)=FTDEW1(I,J)
            FU10(I,J)=FU101(I,J)
            FCLOU(I,J)=FCLOU1(I,J)
#ifdef FORCE_SLP
            PAO(I,J)=FSLP1(I,J)
#endif /*FORCE_SLP*/
          ENDDO
          ENDDO

          DO J=1,JE
          DO I=I_start,IE
             TXO(I,J)=TXO1(I,J)
          ENDDO
          ENDDO
  
          DO J=J_start,JE
          DO I=1,IE
             TYE(I,J)=TYE1(I,J)
          ENDDO
          ENDDO
!
!      For periodic boundaries
!
          CALL bounds_exch(TXO)
          CALL bounds_exch(TYE)
          CALL bounds_exch(TAFO,FCLOU,FSWR,FU10,FPREC,FTDEW)
#ifdef FORCE_SLP
          CALL bounds_exch(PAO)
#endif /*FORCE_SLP*/
!
          DO 335 J=1,JE
          DO 335 I=1,IE
            ACLO(I,J)=FCLOU(I,J)
            RPRECO(I,J)=FPREC(I,J)
            TAIRO(I,J)=TAFO(I,J)
            TDO(I,J)=FTDEW(I,J)
335       CONTINUE
      ENDIF  ! READING FLUXES AT BEGINNING OF DAY
!
!  UPDATE FORCING ONCE PER DAY END
!
!--------------------------------------------------------------------
! OCTHER : THERMODYNAMIC FORCING
!
#ifdef TIMECHECK
      tttr = MPI_Wtime()
#endif
      CALL OCTHER(LDTRUN)
#ifdef TIMECHECK
      ttt = MPI_Wtime()
      write(0,123)'Time OCTHER',ttt-tttr
#endif
123   format(a,f10.3)
#ifdef GMVISETAL
! UPDATE GM COEFFICIENTS ONCE/DAY !
      IF (LDTDAY .EQ. 1) CALL CALCGMVIS
#endif /*GMVISETAL*/
!
!--------------------------------------------------------------------
! OCWIND : 
!
      CALL OCWIND
!
!--------------------------------------------------------------------
!  OBCS : BOUNDARY BT TIDAL FORCINGS ---- ETA
!         Forced ETA 
#ifdef OBC_ETA_FORC
      CALL FORC_OBCS_ETA(LDTRUN) 
         
#ifdef OBC_SPONGE_EXP
      CALL SPONGE_OBCS_ETA_EXP(ZO) 
#else
      CALL SPONGE_OBCS_ETA(ZO) 
#endif      
         
#endif
!
!  New OBCS : BOUNDARY FORCINGS ---- UV
!         Forced Total Flow
#ifdef OBC_UV_FORC
      CALL FORC_OBCS_UV(LDTRUN)
!      CALL PRED_OBCS_UV(UOO,VOE)

#ifdef OBC_SPONGE_EXP
      CALL SPONGE_OBCS_UV_EXP(UOO,VOE)
#else
      CALL SPONGE_OBCS_UV(UOO,VOE)
#endif

#endif
!--------------------------------------------------------------------
!  OCTIMF : UPDATE VELOCITY
!
      CALL OCTIMF
!
!--------------------------------------------------------------------
!  OCMODMOM : DECOMPOSITION INTO BAROTROPIC AND BAROCLINIC FIELD
!
#ifdef TIMECHECK
      tttr = MPI_Wtime()
#endif
      CALL OCMODMOM
#ifdef TIMECHECK
      ttt = MPI_Wtime()
      write(0,123)'Time OCMODMOM',ttt-tttr
#endif
!
!--------------------------------------------------------------------
!    
      CALL OCBARP
!
!--------------------------------------------------------------------
!  Original OBCS : BOUNDARY BT TIDAL FORCINGS ---- UV
!         Forced BT Flow & Radiative BC Flow
! #ifdef OBC_UV_FORC
!       CALL PRED_OBCS_UV(UKO,VKE)
!       CALL FORC_OBCS_UV(LDTRUN)
! #endif
!
!--------------------------------------------------------------------
!     
#ifdef TIMECHECK
      tttr = MPI_Wtime()
#endif
      CALL OCCLIT
#ifdef TIMECHECK
      ttt = MPI_Wtime()
      write(0,123)'Time OCCLIT',ttt-tttr
#endif
!
!--------------------------------------------------------------------
!
#ifdef TIMECHECK
      tttr = MPI_Wtime()
#endif
      IF(IELIMI.GE.1) THEN
        CALL BARTIM
      ELSE
        CALL TRONEU
      ENDIF
#ifdef TIMECHECK
      ttt = MPI_Wtime()
      write(0,123)'Time Solver',ttt-tttr
#endif
60023 FORMAT(4E20.12)
!
!--------------------------------------------------------------------
!    
      CALL OCVTRO
!
!--------------------------------------------------------------------
!
#ifdef TIMECHECK
      tttr = MPI_Wtime()
#endif
      CALL OCVTOT

#ifdef TIMECHECK
      ttt = MPI_Wtime()
      write(0,123)'Time OCVTOT',ttt-tttr
#endif
!
!--------------------------------------------------------------------
!
!UWE   RESET STABIO TO DRHODZ
       DO K=1,KE
        DO J=1,JE
         DO I=1,IE
          STABIO(I,J,K)=(1000./DZ(K))*STABIO(I,J,K)
         ENDDO
        ENDDO
       ENDDO
!
! UWE :NEW MOMENTUM ADVECTION
!
#ifdef TIMECHECK
      tttr = MPI_Wtime()
#endif
      CALL OCUAD(UOO)
      CALL OCVAD(VOE)
#ifdef NON_HYDROSTATIC
      CALL OCWAD(WNO)
#endif
#ifdef TIMECHECK
      ttt = MPI_Wtime()
      write(0,123)'Time OCUAD/OCVAD',ttt-tttr
#endif
!
#ifdef SLOPECON_ADPO
      CALL SLOPETRANS
#endif /*SLOPECON_ADPO*/
!--------------------------------------------------------------------
!  TRACER ADVECTION ROUTINES : SEVERAL OPTIONS 
!             1: 
!             2: UPWIND
!             3: 
!             4: 
!     
       IOCAD=4
!
#ifdef ADPO
       IOCAD=3
#endif /*ADPO*/
!
#ifdef TIMECHECK
      tttr = MPI_Wtime()
#endif
#ifdef ADPO
       IF(IOCAD.EQ.3)THEN
        CALL OCADPO(THO)
!:: SALT ADVECTION
        CALL OCADPO(SAO)
       ENDIF 
#endif /*ADPO*/
!
       IF(IOCAD.EQ.4) CALL OCADFS
#ifdef TIMECHECK
      ttt = MPI_Wtime()
      write(0,123)'Time OCADPO/OCADFS',ttt-tttr
#endif
!
!--------------------------------------------------------------------
! OCJITR : ISOPYCNAL DIFFUSION (W/O GM PARAMETERISATION)
!
#ifdef GMBOLUS
#ifdef TIMECHECK
       tttr = MPI_Wtime()
#endif
       CALL OCJITR
#ifdef TIMECHECK
       ttt = MPI_Wtime()
       write(0,123)'Time OCJITR',ttt-tttr
#endif
#endif /*GMBOLUS*/
!
!--------------------------------------------------------------------
!  Tracer diffusion
!
#ifdef TIMECHECK
       tttr = MPI_Wtime()
#endif
       DO IMAL=1,1
          CALL OCTDIFF_BASE
          CALL OCTDIFF_TS
       ENDDO
#ifdef TIMECHECK
       ttt = MPI_Wtime()
       write(0,123)'Time OCTDIFF',ttt-tttr
#endif
! cxe calculate the surface integrals
       CALL OCTIMF
!
!--------------------------------------------------------------------
!  MOMENTUM DIFFUSION
!
#ifdef TIMECHECK
       tttr = MPI_Wtime()
#endif
       IF (AUS .GT. 1.E-12) CALL OCSCHEP
#ifdef TIMECHECK
       ttt = MPI_Wtime()
       write(0,123)'Time OCSCHEP',ttt-tttr
#endif
!
!--------------------------------------------------------------------
!  BIHARMONIC MOMENTUM DIFFUSION, BOTTOM FRICTION, 
!  SHEAR DEPENDENT DIFFUSION
!
#ifdef TIMECHECK
       tttr = MPI_Wtime()
#endif
       CALL OCVISC
#ifdef TIMECHECK
       ttt = MPI_Wtime()
       write(0,123)'Time OCVISC',ttt-tttr
#endif

!--------------------------------------------------------------------
!  NONHYDROSTATIC PRESSURE CORRECTION
!  VERTICAL VELOCITY(WNO) UPDATE
!
#ifdef NON_HYDROSTATIC
      CALL OCVTOT_NH
#endif

!
       CALL OCTIMF
!

!--------------------------------------------------------------------
! OUTPUT SETTINGS
!
#ifdef MEAN
       IF ( MOD( LDTRUN,NINT(0.1/DT) ) .EQ.  1 ) THEN ! Global output interval: 0.1 s  
        CALL WRTE_MEAN(LDTDAY,LDAYS,LMONTS,LYEARS,IMEAN,NANF)
       ENDIF
#ifdef DIAGNOSIS
       IF ( MOD( LDTRUN,NINT(0.05/DT) ) .EQ.  1 ) THEN ! Diagnostic output interval: 0.05 s  
        CALL WRTE_DIAGNOSIS(LDTDAY,LDAYS,LMONTS,LYEARS,DMEAN,NANF)
       ENDIF
#endif /*DIAGNOSIS*/     
       if (LDTRUN .GE. NINT((150)/DT) )  then ! ISWs simulation period: 150 s 
#ifndef FRAGMENT
          CALL AUFW
#else
          CALL AUFW_FRAGMENT
#endif /*FRAGMENT*/
          call sleep(5)
#ifdef NON_HYDROSTATIC
          call petsc_end
#endif
          CALL ABSTURZ
      endif
! added by peterspy, hourly instant output, if climate, then comment them.
#endif /*MEAN*/

!
#if defined (TIMECHECK)
      ttt = MPI_Wtime()
      if (p_pe==p_io) then
        write(0,'(a,f7.3,a)')'Time for timestep: ', ttt-ttts, ' s'
        write(0,123)'Time for 1 timestep',ttt-ttts
      endif
#endif

! END OF ONE STEP
1001  CONTINUE
!
! END OF ONE DAY
!----------------------------------------------------------------------
!
1010  CONTINUE
!
#ifndef RESYEAR
!  WRITE RESTART FILE AT END OF MONTH
      IF(IAUFW.EQ.1)THEN
#ifndef FRAGMENT
        CALL AUFW
#else
        CALL AUFW_FRAGMENT
#endif /*FRAGMENT*/
      ELSE
        WRITE(IO_STDOUT,*)'STOPPED WITHOUT WRITING RESTART FILE,        &
     &                     IAUFW= ',IAUFW
      ENDIF
#endif
!
! END OF ONE MONTH
!----------------------------------------------------------------------
!
1100  CONTINUE
!
#ifdef RESYEAR
!  WRITE RESTART FILE AT END OF YEAR
      IF(IAUFW.EQ.1)THEN
#ifndef FRAGMENT
        CALL AUFW
#else
        CALL AUFW_FRAGMENT
#endif /*FRAGMENT*/
      ELSE
        WRITE(IO_STDOUT,*)'STOPPED WITHOUT WRITING RESTART FILE         &
     &     AT END OF YEAR,',LYEARS,'IAUFW= ',IAUFW
      ENDIF
#endif /*RESYEAR*/
!
! END OF ONE YEAR
!----------------------------------------------------------------------
!
1000  CONTINUE
!
! END OF TIME STEPPING
!----------------------------------------------------------------------
!
#ifdef NON_HYDROSTATIC
      call petsc_end
#endif
!
!     Branch target for finishing cleanly
!
      CALL print_stats
!
      WRITE(nerr,*) 'NORMAL END OF ORCTM'
!
      CALL p_stop
!
      END


