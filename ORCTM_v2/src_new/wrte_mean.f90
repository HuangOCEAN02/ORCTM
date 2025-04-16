      SUBROUTINE WRTE_MEAN(kdtday,kdays,kmonts,kyears,kmean,nanf)
!C****************************************************************
!C
!C**** *WRTE_MEAN* - save mean diagnostic output.
!
!C     CHH,    *MPI-Met, HH*   14.01.99
!C
!C     Modified
!C     --------
!C     S.Legutke,        *MPI-MaD, HH*    01.10.01
!C     - separate routine extracted from OLLIE (MAIN)
!C
!C     Purpose
!C     -------
!C     Accumulate fields, average, and save.
!C
!C     Method
!C     -------
!CHH   NO OUTPUT     : KMEAN=0
!!CHH   MONTLY AVERAGE: KMEAN=2
!CHH   YEARLY AVERAGE: KMEAN=3 
!C     
!C
!C**   Interface.
!C     ----------
!C!
!C     *CALL*       *WRTE_MEAN(kdtday,kdays,kmonts,kmean)*
!C
!C     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!C     *COMMON*     *COMMO1.h*     - ocean/sediment tracer arrays.
!C     *COMMON*     *UNITS.h*      - std I/O logical units.
!C
!C**   Interface to calling routine (parameter list):
!C     ----------------------------------------------
!C
!C     *INTEGER* *KYEARS*   - actual year.
!C     *INTEGER* *KMONT*   - actual month.
!C     *INTEGER* *KDAYS*    - actual day.
!C
!C
!C     Externals
!C     ---------
!C     none.
!C
!C**************************************************************************
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMOAU1
      USE MO_COMMOAU2
      USE MO_UNITS
      USE MO_MEAN
      USE MO_GLUE
      USE MO_KIND
      NDTDAY=NINT(86400./DT)
!
! Total heat flux and total fresh water flux
      DO I=1,IE
        DO J=1,JE
          FLUM(I,J)=QSWO(I,J)+QLWO(I,J)+QLAO(I,J)+QSEO(I,J)
          PEM(I,J)=PRECH(I,J)+EMINPO(I,J)
        ENDDO
      ENDDO
      DO I=1,IE1
        DO J=1,JE1
! Sea ice transport (x-direction)
          SICTRU(I,J)=( (ABS(SICUO(I,J))+SICUO(I,J))                    &
     &                 *(SICTHO(I,J)+SICSNO(I,J)*RHOSNIC)               &
     &               +( SICUO(I,J)-ABS(SICUO(I,J)))                     &
     &                 *(SICTHO(I+1,J)+SICSNO(I+1,J)*RHOSNIC) )*0.5
! Sea ice transport (y-direction)
          SICTRV(I,J)=( (ABS(SICVE(I,J))+SICVE(I,J))                    &
     &                 *(SICTHO(I,J+1)+SICSNO(I,J+1)*RHOSNIC)           &
     &               +( SICVE(I,J)-ABS(SICVE(I,J)))                     &
     &                 *(SICTHO(I,J)+SICSNO(I,J)*RHOSNIC) )*0.5         
        ENDDO
      ENDDO
       CALL bounds_exch(FLUM,PEM,SICTRU,SICTRV)
!
! Set 'actual' time step and output time step for cases KMEAN=1/2/3:
!
!spy 5-day average
      IF (KMEAN.EQ.5) THEN
        NSTEP=0
        DO I=1,(KMONTS-1)
          NSTEP=NSTEP+MONLEN(I)*NDTDAY
        ENDDO
        NSTEP=NSTEP+(KDAYS-1)*NDTDAY+KDTDAY
        NEND=NDTDAY*5
        NANF=MOD((NSTEP-1),NEND)+1

!HH step by step
      ELSEIF (KMEAN.EQ.4) THEN

        NANF=1
        NEND=1
      
!HH Daily average
      ELSEIF (KMEAN.EQ.1) THEN

        NANF=KDTDAY
        NEND=NDTDAY
      
!HH Monthly average
      ELSEIF (KMEAN.EQ.2) THEN

        IF (KDAYS.EQ.1) THEN
          NANF=KDTDAY
        ELSEIF (KDAYS.GT.1) THEN
          NANF=KDTDAY+((KDAYS-1)*NDTDAY)
        ENDIF
        NEND=NDTDAY*MONLEN(KMONTS)
 
!HH Yearly average
      ELSEIF (KMEAN.EQ.3) THEN

        IF ((KDAYS.EQ.1).AND.(KMONTS.EQ.1)) THEN
          NANF=KDTDAY
        ELSE
          NANF=NANF+1
        ENDIF
        NEND=0
        DO I=1,12
          NEND=NEND+MONLEN(I)
        ENDDO
        NEND=NEND*NDTDAY

!HH No diagnostic output
      ELSEIF (KMEAN.EQ.0) THEN
        CONTINUE
      ELSE
         STOP 'Stop in WRTE_MEAN due to wrong parameter.'
      ENDIF
 
      IF (KMEAN.NE.0) THEN
      WRITE(IO_STDOUT,*) 'nanf,nend',nanf,nend  

! spy: ivec: hrizontal grid -1 -> w point;  0 -> p point; 1 -> u point; 2 -> v point
!HH 3D temporal means
     !   CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MTHO            &
     ! &              ,THO,SUM_THO,2,0)
     !  CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MSAO            &
     ! &              ,SAO,SUM_SAO,5,0)
       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MUKO            &
      &              ,UKO,SUM_UKO,3,1)
      CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MVKE            &
      &              ,VKE,SUM_VKE,4,2)
#ifdef NON_HYDROSTATIC
      IF (FIRST_STEP) THEN
       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MWO             &
     &              ,WNO,SUM_WO,7,-1)
       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_NHPO            &
     &              ,PNH,SUM_PNH,16,0)
       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_NHDIV           &
     &              ,DIVG,SUM_DIVG,26,-1)
      ELSE
       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MWO             &
     &              ,WO,SUM_WO,7,-1)
      ENDIF
#else
!       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MWO             &
!     &              ,WO,SUM_WO,7,-1)
#endif
       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MPO             &
     &              ,PO,SUM_PO,6,0)
#ifdef DENSITYIO
! Variables of those are in mo_glue.f90
       DO K=1,KE
             PPO(:,:) = PO(:,:,K)/100._sp
             TTO = THO(:,:,K)
             SSO = SAO(:,:,K)
            CALL INSITU_THO(TTO,SSO,PPO)
            CALL INSITU_ROO(TTO,SSO,PPO,RRO(:,:))
            ROO(:,:,K) = RRO
       ENDDO
       
       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MRO             &
     &              ,ROO,SUM_ROO,8,0)
#endif
#ifdef GMBOLUS
!       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MWGO            &
!     &              ,WGO,SUM_WGO,207,-1)
#endif
#ifdef DIFFDIAG
!       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_WTMI            &
!     &              ,WTMIX,SUM_WTMIX,612,0)
       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MAVO            &
     &              ,AVO,SUM_AVO,110,-1)
       CALL MMEAN3D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MDVO            &
     &              ,DVO,SUM_DVO,111,-1)
#endif
!HH     2D ZEITMITTEL
!        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_EMIP           &
!     &              ,EMINPO,SUM_EMINPO,67,0)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MMZO           &
     &              ,ZO,SUM_ZO,1,0)
#ifdef BAROTROPIC
! H.H   
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MUSO           &
     &              ,USO,SUM_USO,203,1)
        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MVSE           &
     &              ,VSE,SUM_VSE,204,2) 
#endif       
!        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_FLUM           &
!     &              ,FLUM,SUM_FLUM,70,0)
!        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_MPEM           &
!     &              ,PEM,SUM_PEM,79,0)
!        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICT           &
!     &              ,SICTHO,SUM_SICTHO,13,0)
!        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICO           &
!     &              ,SICOMO,SUM_SICOMO,15,0)
!        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICU           &
!     &              ,SICUO,SUM_SICUO,35,0)
!        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICV           &
!     &              ,SICVE,SUM_SICVE,36,0)
!        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICS           &
!    &              ,SICSNO,SUM_SICSNO,141,0)
!        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICTRU         &
!     &              ,SICTRU,SUM_SICTRU,142,0)
!        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_SICTRV         &
!     &              ,SICTRV,SUM_SICTRV,143,0)
#ifdef NON_HYDROSTATIC
#ifdef CG2D4NONHY 
!        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_NHPOBOT        &
!     &              ,PNH_bot,SUM_PNHBOT,17,0)
!        CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_NHDIVBOT       &
!     &              ,DIVG_bot_surf,SUM_DIVGBOT,18,0)          
#endif /*CG2D4NONHY*/      
#endif /*NON_HYDROSTATIC*/    
! spy: note that obcs for sea ice has not been developed yet,
! spy: ivec for SICUO SICVE SICTRU SICTRV are set to zero now.
#ifdef GMBOLUS
     !    CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_MBOLX          &
     ! &              ,BOLX,SUM_BOLX,159,0)
     !   CALL MMEAN2D(LDAYS,LMONTS,LYEARS,NANF,NEND,IO_OU_MBOLY          &
     ! &              ,BOLY,SUM_BOLY,160,0)
! spy: no need to use vector dimensions for GMBOLUS (BOLX/BOLY)
#endif /*GMBOLUS*/
!       CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_QSWO            &
!     &              ,QSWO,SUM_QSWO,176,0)
!       CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_QLWO            &
!     &              ,QLWO,SUM_QLWO,177,0)
!       CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_QLAO            &
!     &              ,QLAO,SUM_QLAO,147,0)
!       CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_QSEO            &
!     &              ,QSEO,SUM_QSEO,146,0)
!       CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_PREC            &
!     &              ,PRECH,SUM_PRECH,65,0)
!       CALL MMEAN2D(KDAYS,KMONTS,KYEARS,NANF,NEND,IO_OU_RIVRUN          &
!     &              ,RIVRUN,SUM_RIVRUN,305,0)
      ENDIF

      RETURN
      END
