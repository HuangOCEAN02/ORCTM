      SUBROUTINE OCADFS
#ifdef ADFS
!**********************************************************************
!
!
!      OOO    CCCCC  AAAAAA   DDDDD   FFFFFF   SSS
!     O   O  C       A    A   D    D  F       S   
!     O   O  C       AAAAAA   D    D  FFF      SSS
!     O   O  C       A    A   D    D  F           S
!      OOO    CCCCC  A    A   DDDDD   F        SSS
!
!
!=======================================================================
!
!     SBR OCADFS
!
!     
!     PURPOSE :
!
!     A) 3-D ADVECTION OF TEMPERATURE AND SALINITY
!        PREDICTOR-CORRECTOR SCHEME
!        IF COMPILE OPTION QUICK THEN
!          USING QUICK-SCHEME AS PROPOSED BY FARROW AND STEVENS 1995
!          JPO 25, 1731...
!          3RD ORDER ACCURATE IN SPACE, 2ND IN TIME
!        ELSEIF COMPILE OPTION QUICK2
!          USE TEMPERATURE OF THE INFLOWING WATER IN I-1 INSTEAD OF 
!          TEMPERATURE IN I-2
!          SEEMS TO BE MORE EFFICIENT IN REDUCING OVERSHOOTINGS.
!          BETTER CONDITIONED AT BOUNDARIES
!        ELSE
!          ORDINARY CENTERED DIFFERENCES
! 
!
!      BY UWE MIKOLAJEWICZ    1/99
!      OPTION QUICK2 ADDED    3/99
!      DIFFUSION REMOVED      5/99
!      Changed for MPI-Parallelization R. Johanni 11/03
!
      USE MO_PARAM1
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_COMMOAU1

      USE MO_PARALLEL

      DIMENSION CURVUPT(IE,JE),CURVLOT(IE,JE),CURVEAT(IE,JE)            &
     &         ,CURVWET(IE,JE),CURVNOT(IE,JE),CURVSOT(IE,JE)
      DIMENSION CURVUPS(IE,JE),CURVLOS(IE,JE),CURVEAS(IE,JE)            &
     &         ,CURVWES(IE,JE),CURVNOS(IE,JE),CURVSOS(IE,JE)
      DIMENSION ZSURF(KE)
      DIMENSION DLXP2(IE,JE),DLYP2(IE,JE)
#ifdef QUICK2
      DIMENSION TINF(IE,JE,KE),SINF(IE,JE,KE)
#endif
#ifdef QUICK
      REAL DSX(IE,JE,KE),DSY(IE,JE,KE),DTX(IE,JE,KE),DTY(IE,JE,KE)
#endif
      DO K=1,KE
       ZSURF(K)=0.
      ENDDO
      ZSURF(1)=1.
!
!C      WRITE(IO_STDOUT,*)'FS ADVEKTION!!!!'
!---------------------------------------------------------------------
!
!     A.1)
!
!     1. HALF STEP IN TIME
!
!
!---------------------------------------------------------------------
!
!     A.1.1)
!
       FOURTH=0.25
!----------------------------------------------------------------------
!
!     A.1.2)
!
!     ODD FIELDS
!
!
      DO 175 K=1,KE
!
!
      KUP = MAX(1,K-1)
      KLO = MIN(KE,K+1)
!
!
      DO 175 J=2,JE1
      DO 175 I=2,IE1
!
      TM = THO(I,J,K)
      SM = SAO(I,J,K)
!
      WUN = WO(I,J,K+1)*DLXP(I,J)*DLYP(I,J)                             &
     &     /(DZW(KLO)+DZW(K))
      WOB = WO(I,J,K)*DLXP(I,J)*DLYP(I,J)                               &
     &     /(DZW(KUP)+DZW(K))
      UWE = UKO(I-1,J,K)*DDUO(I-1,J,K)*DLYU(I-1,J)                      &
     &     /(DLXP(I,J)+DLXP(I-1,J))
      UOS = UKO(I,J,K)*DDUO(I,J,K)*DLYU(I,J)                            &
     &     /(DLXP(I,J)+DLXP(I+1,J))
      VSU = VKE(I,J,K)*DDUE(I,J,K)*DLXV(I,J)                            &
     &     /(DLYP(I,J)+DLYP(I,J+1))
      VNO = VKE(I,J-1,K)*DLXV(I,J-1)*DDUE(I,J-1,K)                      &
     &     /(DLYP(I,J)+DLYP(I,J-1))
!
      T1O(I,J,K) = THO (I,J,K) + 0.5*DT*WETO(I,J,K)*(                   &
     &         UWE * ( DLXP(I,J)*THO(I-1,J,K) +DLXP(I-1,J)*TM )         &
     &     -   UOS * ( DLXP(I+1,J)*TM + DLXP(I,J)*THO(I+1,J,K) )        &
     &     +   VSU * ( DLYP(I,J)*THO(I,J+1,K) + DLYP(I,J+1)*TM )        &
     &     -   VNO * ( DLYP(I,J-1)*TM + DLYP(I,J)*THO(I,J-1,K) )        &
     &     +   WUN * ( DZW(K)*THO(I,J,KLO) + DZW(KLO)*TM )              &
     &     -   WOB * ( DZW(KUP)*TM + DZW(K)*THO(I,J,KUP) )   )          &
!CUWE  INCLUDE ZETA IN UPPERMOST LEVEL THICKNESS 12/99
     &  /((DDPO(I,J,K)+ZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)            &
     &       -RHOSNWA*SICSNO(I,J))+ALMZER)                              &
     &      *DLXP(I,J)*DLYP(I,J))
!
!
      S1O(I,J,K) = SAO (I,J,K) + 0.5*DT*WETO(I,J,K)*(                   &
     &         UWE * ( DLXP(I,J)*SAO(I-1,J,K) +DLXP(I-1,J)*SM )         &
     &     -   UOS * ( DLXP(I+1,J)*SM + DLXP(I,J)*SAO(I+1,J,K) )        &
     &     +   VSU * ( DLYP(I,J)*SAO(I,J+1,K) + DLYP(I,J+1)*SM )        &
     &     -   VNO * ( DLYP(I,J-1)*SM + DLYP(I,J)*SAO(I,J-1,K) )        &
     &     +   WUN * ( DZW(K)*SAO(I,J,KLO) + DZW(KLO)*SM )              &
     &     -   WOB * ( DZW(KUP)*SM + DZW(K)*SAO(I,J,KUP) )   )          &
     &  /((DDPO(I,J,K)+ZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)            &
     &                          -RHOSNWA*SICSNO(I,J))+ALMZER)           &
     &      *DLXP(I,J)*DLYP(I,J))
!
175    CONTINUE
!
      CALL bounds_exch(T1O)
      CALL bounds_exch(S1O)
!
       FOURTH=0.25
!
!
!----------------------------------------------------------------------
!
!     A.2)
!
!     2. HALF STEP IN TIME
!
!
!----------------------------------------------------------------------
!
      TURDIF = 1.E-3
!
#ifdef QUICK2
!
!UWE
!     COMPUTE T AND S OF INFLOWING WATER
!     
      DO K=1,KE
       KLO=MIN(K+1,KE)
       KUP=MAX(K-1,1)
       DO J=2,JE-1
        DO I=2,IE-1
         WUN=HALF*(WO(I,J,K+1)+ABS(WO(I,J,K+1)) )*DLXP(I,J)*DLYP(I,J)
         WOB=HALF*(ABS(WO(I,J,K))-WO(I,J,K))*DLXP(I,J)*DLYP(I,J)
         UWE=HALF*(UKO(I-1,J,K)+ABS(UKO(I-1,J,K)))                      &
     &           *DLYU(I-1,J)*DDUO(I-1,J,K)
         UOS=HALF*(ABS(UKO(I,J,K))-UKO(I,J,K))*DLYU(I,J)*DDUO(I,J,K)
         VSU=HALF*(VKE(I,J,K)+ABS(VKE(I,J,K)))*DLYV(I,J)*DDUE(I,J,K)
         VNO=HALF*(ABS(VKE(I,J-1,K))-VKE(I,J-1,K))                      &
     &           *DLYV(I,J-1)*DDUE(I,J-1,K)
         GEW=ALMZER+WUN+WOB+UWE+UOS+VSU+VNO
         TINF(I,J,K)=(ALMZER*T1O(I,J,K)                                 &
     &        +WUN*T1O(I,J,KLO)+WOB*T1O(I,J,KUP)                        &
     &        +UWE*T1O(I-1,J,K)+UOS*T1O(I+1,J,K)                        &
     &        +VSU*T1O(I,J+1,K)+VNO*T1O(I,J-1,K))/GEW
         SINF(I,J,K)=(ALMZER*S1O(I,J,K)                                 &
     &        +WUN*S1O(I,J,KLO)+WOB*S1O(I,J,KUP)                        &
     &        +UWE*S1O(I-1,J,K)+UOS*S1O(I+1,J,K)                        &
     &        +VSU*S1O(I,J+1,K)+VNO*S1O(I,J-1,K))/GEW
        ENDDO
       ENDDO
       DO I=1,IE
        TINF(I,1,K)=0.
        SINF(I,1,K)=0.
        TINF(I,JE,K)=0.
        SINF(I,JE,K)=0.
       ENDDO
      ENDDO
      CALL bounds_exch(TINF)
      CALL bounds_exch(SINF)
!
#endif
!
!
!     Calculate some auxiliary arrays

      DLXP2(:,:) = 0.
      DLYP2(:,:) = 0.
      DO J=1,JE-1
      DO I=2,IE-1
        DLXP2(I,J) = DLXP(I,J)+DLXP(I+1,J)
        DLYP2(I,J) = DLYP(I,J)+DLYP(I,J+1)
      ENDDO
      ENDDO

      CALL bounds_exch(DLXP2,DLYP2)

#ifdef QUICK
      DO K=1,KE
        DO J=1,JE-1
        DO I=2,IE-1
          DSX(I,J,K) = (S1O(I+1,J,K)-S1O(I,J,K))/(DLXP(I,J)+DLXP(I+1,J))
          DTX(I,J,K) = (T1O(I+1,J,K)-T1O(I,J,K))/(DLXP(I,J)+DLXP(I+1,J))
          DSY(I,J,K) = (S1O(I,J,K)-S1O(I,J+1,K))/(DLYP(I,J)+DLYP(I,J+1))
          DTY(I,J,K) = (T1O(I,J,K)-T1O(I,J+1,K))/(DLYP(I,J)+DLYP(I,J+1))
        ENDDO
        ENDDO
        DO I=2,IE-1
          DSX(I,JE,K) = 0.
          DTX(I,JE,K) = 0.
          DSY(I,JE,K) = 0.
          DTY(I,JE,K) = 0.
        ENDDO
      ENDDO

      CALL bounds_exch(DSX)
      CALL bounds_exch(DTX)
      CALL bounds_exch(DSY)
      CALL bounds_exch(DTY)
#endif

!
      DO 185 K=1,KE
      DO J=1,JE
       DO I=1,IE
        CURVNOT(I,J)=0.
        CURVSOT(I,J)=0.
        CURVWET(I,J)=0.
        CURVEAT(I,J)=0.
        CURVNOS(I,J)=0.
        CURVSOS(I,J)=0.
        CURVWES(I,J)=0.
        CURVEAS(I,J)=0.
       ENDDO
      ENDDO
!
      KUP = MAX(1,K-1)
      KLO = MIN(KE,K+1)
!
      KUPUP=MAX(1,K-2)
      KLOLO=MIN(KE,K+2)
      IF(K.EQ.1)THEN
      DO J=1,JE
       DO I=1,IE
         CURVUPT(I,J)=0.
         CURVUPS(I,J)=0.
         CURVLOT(I,J)=0.
         CURVLOS(I,J)=0.
       ENDDO
      ENDDO
      ELSE
      DO J=1,JE
        DO I=1,IE
         CURVUPT(I,J)=CURVLOT(I,J)
         CURVUPS(I,J)=CURVLOS(I,J)
        ENDDO
      ENDDO
      ENDIF
#ifdef QUICK2
!
      DO J=1,JE
       DO I=2,IE-1
        CURVLOT(I,J)=0.
        CURVLOS(I,J)=0.
        IF(WO(I,J,K+1).GT.0.)THEN
!
!       UPWARD FLOW
!HH UM SWICH OFF ?
!CWET        IF(K.LE.KE-1.AND.WETO(I,J,KLO).GT.0.5)THEN
!C        IF(K.LE.KE-2.AND.WETO(I,J,KLOLO).GT.0.5)THEN
        CURVLOT(I,J)=((T1O(I,J,K)-T1O(I,J,KLO))/(DZW(K)+DZW(KLO))       &
     &          -(T1O(I,J,KLO)-TINF(I,J,KLO))/(DZW(KLO)+DZW(KLOLO)))    &
     &  *DZW(K)*DZW(KLO)/(DZW(KLOLO)+2.*DZW(KLO)+DZW(K))
        CURVLOS(I,J)=((S1O(I,J,K)-S1O(I,J,KLO))/(DZW(K)+DZW(KLO))       &
     &          -(S1O(I,J,KLO)-SINF(I,J,KLO))/(DZW(KLO)+DZW(KLOLO)))    &
     &  *DZW(K)*DZW(KLO)/(DZW(KLOLO)+2.*DZW(KLO)+DZW(K))
!CWET        ENDIF
        ELSE
!
!       DOWNWARD FLOW
!
!C        IF(K.GT.1.AND.K.LT.KE)THEN
!C        IF(K.LT.KE)THEN
!C        CURVLOT(I,J)=((T1O(I,J,KUP)-T1O(I,J,K))/(DZW(K)+DZW(KUP))
        CURVLOT(I,J)=((TINF(I,J,K)-T1O(I,J,K))/(DZW(K)+DZW(KUP))        &
     &         -(T1O(I,J,K)-T1O(I,J,KLO))/(DZW(KLO)+DZW(K)))            &
     &  *WETO(I,J,KLO)*DZW(KLO)*DZW(K)/(DZW(KLO)+2.*DZW(K)+DZW(KUP))
        CURVLOS(I,J)=((SINF(I,J,K)-S1O(I,J,K))/(DZW(K)+DZW(KUP))        &
     &         -(S1O(I,J,K)-S1O(I,J,KLO))/(DZW(KLO)+DZW(K)))            &
     &  *WETO(I,J,KLO)*DZW(KLO)*DZW(K)/(DZW(KLO)+2.*DZW(K)+DZW(KUP))
!C        ENDIF
        ENDIF
       ENDDO
      ENDDO
#endif
#ifdef QUICK
!
      DO J=1,JE
       DO I=2,IE-1
        CURVLOT(I,J)=0.
        CURVLOS(I,J)=0.
        IF(WO(I,J,KLO).GE.0.)THEN
!
!       UPWARD FLOW
!
        IF(K.LT.KE-1.AND.WETO(I,J,KLOLO).GT.0.5)THEN
        CURVLOT(I,J)=((T1O(I,J,K)-T1O(I,J,KLO))/(DZW(K)+DZW(KLO))       &
     &          -(T1O(I,J,KLO)-T1O(I,J,KLOLO))/(DZW(KLO)+DZW(KLOLO)))   &
     &  *DZW(K)*DZW(KLO)/(DZW(KLOLO)+2.*DZW(KLO)+DZW(K))
        CURVLOS(I,J)=((S1O(I,J,K)-S1O(I,J,KLO))/(DZW(K)+DZW(KLO))       &
     &          -(S1O(I,J,KLO)-S1O(I,J,KLOLO))/(DZW(KLO)+DZW(KLOLO)))   &
     &  *DZW(K)*DZW(KLO)/(DZW(KLOLO)+2.*DZW(KLO)+DZW(K))
        ENDIF
        ELSE
!
!       DOWNWARD FLOW
!
!CC        IF(K.GT.1.AND.K.LT.KE)THEN
        CURVLOT(I,J)=((T1O(I,J,KUP)-T1O(I,J,K))/(DZW(K)+DZW(KUP))       &
     &         -(T1O(I,J,K)-T1O(I,J,KLO))/(DZW(KLO)+DZW(K)))            &
     &  *WETO(I,J,KLO)*DZW(KLO)*DZW(K)/(DZW(KLO)+2.*DZW(K)+DZW(KUP))
        CURVLOS(I,J)=((S1O(I,J,KUP)-S1O(I,J,K))/(DZW(K)+DZW(KUP))       &
     &         -(S1O(I,J,K)-S1O(I,J,KLO))/(DZW(KLO)+DZW(K)))            &
     &  *WETO(I,J,KLO)*DZW(KLO)*DZW(K)/(DZW(KLO)+2.*DZW(K)+DZW(KUP))
!CC        ENDIF
        ENDIF
       ENDDO
      ENDDO
#endif
!
      CALL bounds_exch(CURVLOT,CURVLOS)
!
#ifdef QUICK2
      DO J=2,JE1
       DO I=2,IE1
!
!       SOUTHERLY POINT
!
        IF(VKE(I,J,K).GT.0.)THEN
!
!        INFLOW FROM SOUTH
!
         CURVSOT(I,J)=                                                  &
     &     ((T1O(I,J,K)-T1O(I,J+1,K))/DLYP2(I,J)                        &
     &      -(T1O(I,J+1,K)-TINF(I,J+1,K))/DLYP2(I,J+1))                 &
     &     *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J)+DLYP2(I,J+1))
!
         CURVSOS(I,J)=                                                  &
     &     ((S1O(I,J,K)-S1O(I,J+1,K))/DLYP2(I,J)                        &
     &      -(S1O(I,J+1,K)-SINF(I,J+1,K))/DLYP2(I,J+1))                 &
     &     *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J)+DLYP2(I,J+1))
!
        ELSE
!
!        OUTFLOW TO SOUTH
!
         CURVSOT(I,J)=                                                  &
     &     ((TINF(I,J,K)-T1O(I,J,K))/DLYP2(I,J-1)                       &
     &      -(T1O(I,J,K)-T1O(I,J+1,K))/DLYP2(I,J))                      &
     &     *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J-1)+DLYP2(I,J))
!
         CURVSOS(I,J)=                                                  &
     &     ((SINF(I,J,K)-S1O(I,J,K))/DLYP2(I,J-1)                       &
     &      -(S1O(I,J,K)-S1O(I,J+1,K))/DLYP2(I,J))                      &
     &     *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J-1)+DLYP2(I,J))
        ENDIF
       ENDDO
      ENDDO

      CALL bounds_exch(CURVSOT,CURVSOS)

      DO J=2,JE
       DO I=2,IE-1
        CURVNOT(I,J)=CURVSOT(I,J-1)
        CURVNOS(I,J)=CURVSOS(I,J-1)
       ENDDO
      ENDDO

      CALL bounds_exch(CURVNOT,CURVNOS)

      DO J=2,JE1
       DO I=2,IE1
        IF(UKO(I,J,K).GT.0.)THEN
!
!        INFLOW
!
         CURVEAT(I,J)=                                                  &
     &     ((T1O(I+1,J,K)-T1O(I,J,K))/DLXP2(I,J)                        &
     &      -(T1O(I,J,K)-TINF(I,J,K))/DLXP2(I-1,J))                     &
     &     *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I-1,J)+DLXP2(I,J))
!
         CURVEAS(I,J)=                                                  &
     &     ((S1O(I+1,J,K)-S1O(I,J,K))/DLXP2(I,J)                        &
     &      -(S1O(I,J,K)-SINF(I,J,K))/DLXP2(I-1,J))                     &
     &     *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I-1,J)+DLXP2(I,J))
        ELSE
!
!        OUTFLOW
!
         CURVEAT(I,J)=                                                  &
     &     ((TINF(I+1,J,K)-T1O(I+1,J,K))/DLXP2(I+1,J)                   &
     &      -(T1O(I+1,J,K)-T1O(I,J,K))/DLXP2(I,J))                      &
     &     *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I,J)+DLXP2(I+1,J))
!
         CURVEAS(I,J)=                                                  &
     &     ((SINF(I+1,J,K)-S1O(I+1,J,K))/DLXP2(I+1,J)                   &
     &      -(S1O(I+1,J,K)-S1O(I,J,K))/DLXP2(I,J))                      &
     &     *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I,J)+DLXP2(I+1,J))
        ENDIF
!
       ENDDO
      ENDDO

      CALL bounds_exch(CURVEAT,CURVEAS)

      DO J=1,JE
       DO I=2,IE-1
        CURVWET(I,J)=CURVEAT(I-1,J)
        CURVWES(I,J)=CURVEAS(I-1,J)
       ENDDO
      ENDDO

      CALL bounds_exch(CURVWET,CURVWES)
#endif
#ifdef QUICK
      DO J=2,JE1
       DO I=2,IE1
!
!       SOUTHERLY POINT
!
        IF(AMSUE(I,J,K).GT.0.5)THEN
         IF(VKE(I,J,K).GT.0.)THEN
!
!         INFLOW FROM SOUTH
!
          CURVSOT(I,J)=AMSUE(I,J+1,K)*(DTY(I,J,K)-DTY(I,J+1,K))         &
     &                *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J)+DLYP2(I,J+1))
!
          CURVSOS(I,J)=AMSUE(I,J+1,K)*(DSY(I,J,K)-DSY(I,J+1,K))         &
     &                *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J)+DLYP2(I,J+1))
!
         ELSE
!
!         OUTFLOW TO SOUTH
!
          CURVSOT(I,J)=AMSUE(I,J-1,K)*(DTY(I,J-1,K)-DTY(I,J,K))         &
     &                *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J-1)+DLYP2(I,J))
!
          CURVSOS(I,J)=AMSUE(I,J-1,K)*(DSY(I,J-1,K)-DSY(I,J,K))         &
     &                *DLYP(I,J)*DLYP(I,J+1)/(DLYP2(I,J-1)+DLYP2(I,J))
         ENDIF
        ENDIF
       ENDDO
      ENDDO

      CALL bounds_exch(CURVSOT,CURVSOS)

      DO J=2,JE
       DO I=2,IE-1
        CURVNOT(I,J)=CURVSOT(I,J-1)
        CURVNOS(I,J)=CURVSOS(I,J-1)
       ENDDO
      ENDDO

      CALL bounds_exch(CURVNOT,CURVNOS)

      DO J=2,JE1
       DO I=2,IE1
        IF(AMSUO(I,J,K).GT.0.5)THEN
         IF(UKO(I,J,K).GT.0.)THEN
!
!         INFLOW
!
          CURVEAT(I,J)=AMSUO(I-1,J,K)*(DTX(I,J,K)-DTX(I-1,J,K))         &
     &                *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I-1,J)+DLXP2(I,J))
!
          CURVEAS(I,J)=AMSUO(I-1,J,K)*(DSX(I,J,K)-DSX(I-1,J,K))         &
     &                *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I-1,J)+DLXP2(I,J))
         ELSE
!
!         OUTFLOW
!
          CURVEAT(I,J)=AMSUO(I+1,J,K)*(DTX(I+1,J,K)-DTX(I,J,K))         &
     &                *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I,J)+DLXP2(I+1,J))
!
          CURVEAS(I,J)=AMSUO(I+1,J,K)*(DSX(I+1,J,K)-DSX(I,J,K))         &
     &                *DLXP(I,J)*DLXP(I+1,J)/(DLXP2(I,J)+DLXP2(I+1,J))
         ENDIF
        ENDIF
       ENDDO
      ENDDO

      CALL bounds_exch(CURVEAT,CURVEAS)

      DO J=1,JE
       DO I=2,IE-1
        CURVWET(I,J)=CURVEAT(I-1,J)
        CURVWES(I,J)=CURVEAS(I-1,J)
       ENDDO
      ENDDO

      CALL bounds_exch(CURVWET,CURVWES)
#endif
!
!
      DO 185 J=2,JE1
      DO 185 I=2,IE1
!
      TM = T1O(I,J,K)
      SM = S1O(I,J,K)
!
      TA=THO(I,J,K)
      SA=SAO(I,J,K)
      WUN = WO(I,J,K+1)*DLXP(I,J)*DLYP(I,J)                             &
     &     /(DZW(KLO)+DZW(K))
      WOB = WO(I,J,K)*DLXP(I,J)*DLYP(I,J)                               &
     &     /(DZW(KUP)+DZW(K))
      UWE = UKO(I-1,J,K)*DDUO(I-1,J,K)*DLYU(I-1,J)                      &
     &     /(DLXP(I,J)+DLXP(I-1,J))
      UOS = UKO(I,J,K)*DDUO(I,J,K)*DLYU(I,J)                            &
     &     /(DLXP(I,J)+DLXP(I+1,J))
      VSU = VKE(I,J,K)*DDUE(I,J,K)*DLXV(I,J)                            &
     &     /(DLYP(I,J)+DLYP(I,J+1))
      VNO = VKE(I,J-1,K)*DLXV(I,J-1)*DDUE(I,J-1,K)                      &
     &     /(DLYP(I,J)+DLYP(I,J-1))
!
      UK1O(I,J,K) = THO (I,J,K) + DT*WETO(I,J,K)*(                      &
     &         UWE * ( DLXP(I,J)*T1O(I-1,J,K) +DLXP(I-1,J)*TM           &
     &         -(DLXP(I,J)+DLXP(I+1,J))*CURVWET(I,J))                   &
     &     -   UOS * ( DLXP(I+1,J)*TM + DLXP(I,J)*T1O(I+1,J,K)          &
     &         -(DLXP(I,J)+DLXP(I+1,J))*CURVEAT(I,J))                   &
     &     +   VSU * ( DLYP(I,J)*T1O(I,J+1,K) + DLYP(I,J+1)*TM          &
     &         -(DLYP(I,J)+DLYP(I,J-1))*CURVSOT(I,J))                   &
     &     -   VNO * ( DLYP(I,J-1)*TM + DLYP(I,J)*T1O(I,J-1,K)          &
     &         -(DLYP(I,J)+DLYP(I,J-1))*CURVNOT(I,J))                   &
     &     +WUN*(DZW(K)*T1O(I,J,KLO)+DZW(KLO)*TM                        &
     &         -(DZW(KLO)+DZW(K))*CURVLOT(I,J))                         &
     &     -WOB*(DZW(KUP)*TM + DZW(K)*T1O(I,J,KUP)                      &
     &         -(DZW(KUP)+DZW(K))*(1.-ZSURF(K))*CURVUPT(I,J)))          &
     &  /((DDPO(I,J,K)+ZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)            &
     &                  -RHOSNWA*SICSNO(I,J))+ALMZER)                   &
     &      *DLXP(I,J)*DLYP(I,J))
!
      VK1E(I,J,K) = SAO (I,J,K) + DT*WETO(I,J,K)*(                      &
     &         UWE * ( DLXP(I,J)*S1O(I-1,J,K) +DLXP(I-1,J)*SM           &
     &         -(DLXP(I,J)+DLXP(I-1,J))*CURVWES(I,J))                   &
     &     -   UOS * ( DLXP(I+1,J)*SM + DLXP(I,J)*S1O(I+1,J,K)          &
     &         -(DLXP(I,J)+DLXP(I+1,J))*CURVEAS(I,J))                   &
     &     +   VSU * ( DLYP(I,J)*S1O(I,J+1,K) + DLYP(I,J+1)*SM          &
     &         -(DLYP(I,J)+DLYP(I,J+1))*CURVSOS(I,J))                   &
     &     -   VNO * ( DLYP(I,J-1)*SM + DLYP(I,J)*S1O(I,J-1,K)          &
     &         -(DLYP(I,J)+DLYP(I,J-1))*CURVNOS(I,J))                   &
     &     +WUN*(DZW(K)*S1O(I,J,KLO)+DZW(KLO)*SM                        &
     &         -(DZW(KLO)+DZW(K))*CURVLOS(I,J))                         &
     &     -WOB*(DZW(KUP)*SM + DZW(K)*S1O(I,J,KUP)                      &
     &         -(DZW(KUP)+DZW(K))*CURVUPS(I,J)))                        &
     & /((DDPO(I,J,K)+ZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)             &
     &       -RHOSNWA*SICSNO(I,J))+ALMZER)                              &
     &       *DLXP(I,J)*DLYP(I,J))
!
185    CONTINUE
!
      CALL bounds_exch(UK1O)
      CALL bounds_exch(VK1E)
!
       DO K=1,KE
        DO J=2,JE1
         DO I=2,IE1
          THO(I,J,K)=UK1O(I,J,K)
          SAO(I,J,K)=VK1E(I,J,K)
         ENDDO
        ENDDO
       ENDDO
!
#endif /*ADFS*/
      RETURN
      END
