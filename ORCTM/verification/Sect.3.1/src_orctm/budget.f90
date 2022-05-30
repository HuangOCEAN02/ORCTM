      SUBROUTINE BUDGET(ALAT,FHS,FHB,T,LRHS,A,                        &
     &         TAIR,TD,ACL,PA,UG,SLN,OM,HICE,ALB,A2,QSW,QLW,QSE,QLA)
!
!=======================================================================
!     A. STOESSEL               MPI, HAMBURG                         1989
!
!
!  MODIFIED
!  --------
!  S. LEGUTKE           *DKRZ*               16.04.95 
!     - SEPARATE CALCULATION OF SURFACE/BOTTOM MELT/GROWTH.
!
!  UWE    15.1.2000
!  - REMOVE ERROR IN SW CALCULATION (DIAGNOSTIC)
!  - INCLUDE UPDATE OF TEMPERATURE IN SW-FLUXES IN ITERATION
!  - MAKE ALBEDO CONTINUOUSLY DEPENDENT OF TEMEPRATURE
!  HELMUTH 19.04.2000
!  - LONGWAVE-FLUX FROM BERLIAND(1952)
!
!  PURPOSE:
!  --------
!      CALCULATION OF GROWTH RATES FOR THE ICE COVERED PART OF A GRID
!      CELL WITH STANDARD BULK FORMULAS.
!
!  METHOD:
!  -------
!      ICE OR SNOW SURFACE TEMPERATURES, RESPECTIVELY, ARE CALCULATED BY
!      ITERATION (REGULA FALSI)
!
!  INTERFACE:
!  ----------
!     INPUT:
!     *A*        ICE COMPACTNESS             [FRAC.]
!     *TAIR*     AIR TEMPERATURES            [DEG C]
!     *TD*       DEW POINT TEMPERATURES      [K]
!     *ACL*      FRACTIONAL CLOUD COVER      [FRAC.]
!     ifdef FORCE_DWLW ::: ACL CONTAINS DOWNWARD LONGWAVE [W/M**2] !!!
!     *PA*       ATMOSPHERIC SURFACE PRESURE [PA]
!     *UG*       2 M WIND SPEED              [M/SEC]
!     *SLN*      INCOMING SURFACE SOLAR RAD. [W/M**2]
!     *OM*       LAND/SEA MASK               [0/1]
!     *HICE*     ICE THICKNESS               [M]
!     *ALB*      SURFACE ALBEDO
!     *A2*       SNOW FLAG
!
!     CHANGED:
!     *T*        SURFACE SNOW/ICE TEMPERATURE[DEG C]
!
!     OUTPUT:
!     *FHS*      SURFACE MELT THICK ICE               [M/SEC]
!     *FHB*      BOTTOM MELT/GROWTH OF THICK ICE      [M/SEC]
!     *QSW*      ABSORBED SOLAR RADIATION             [W/M**2]
!     *QLW*      OUTGOING LONGWAVE HEAT FLUX          [W/M**2]
!     *QSE*      SENSIBLE HEAT FLUX                   [W/M**2]
!     *QLA*      LATENT HEAT FLUX                     [W/M**2]      
!
!  EXTERNALS:
!  ----------
!     *VAPOR*    CALCULATES VAPOR PRESSURE
!     
!=======================================================================
!
      USE MO_PARAM1
      USE MO_COMMOAU1
      DIMENSION T(IE,JE),FHS(IE,JE),FHB(IE,JE),SLN(IE,JE)
      DIMENSION TAIR(IE,JE),TD(IE,JE),ACL(IE,JE),PA(IE,JE),UG(IE,JE)
      DIMENSION OM(IE,JE)
      DIMENSION A(IE,JE,2),HICE(IE,JE),ALB(IE,JE),A2(IE,JE)
      DIMENSION TA1(IE*JE),TD1(IE*JE),ACL1(IE*JE),PA1(IE*JE),UG1(IE*JE),&
     &SLN1(IE*JE),FLSE1(IE*JE),FLLA1(IE*JE),FHS1(IE*JE),FHB1(IE*JE),    &
     &T1(IE*JE),ESTA(IE*JE),ESTI(IE*JE),FAKTS(IE*JE)
      DIMENSION HICE1(IE*JE),ALB1(IE*JE),A21(IE*JE),FLAI(IE*JE)
      DIMENSION QSW(IE,JE),QLW(IE,JE),QSE(IE,JE),QLA(IE,JE)
      DIMENSION Q1(IE*JE),Q2(IE*JE),Q3(IE*JE)
      DIMENSION XLAT(IE*JE),ALAT(IE,JE)
!
      DIMENSION STP(IE*JE),STPP(IE*JE),FP(IE*JE),FPP(IE*JE),DIFF(IE*JE) &
     &      ,DRAGL(IE*JE),DRAGS(IE*JE),SPHUMIDA(IE*JE),SPHUMIDO(IE*JE)
!
#ifdef BULK_KARA
      DIMENSION RHOAIR_K(IE*JE)
#endif

      L = IE
      M = JE
!
!-----------------------------------------------------------------------
!  SET MAXIMUM NUMBER OF ITERATION STEPS
!  SET ICE TEMPERATURE AT SEA ICE INTERFACE
!-----------------------------------------------------------------------
      RGAS=287.1
      CPA=1004.67
      ALMZER=1.E-19
! Frank Roeske's budget closing factor:
      FR_FAC=1.1925
#ifdef FORCE_DWLW
      D3=0.97*5.67E-8
#endif
!
!
      IMAX=10
      TB=TFREZ+TMELT
!UWE   WIDTH OF ALBEDOTRANSITION MELT TO FREEZ
!
!     LARGE VALE: SHARP TRANSITION
!     LOW VALUE:  SMOOTHE TRANSITION
      ALBTRANS=0.5
!-----------------------------------------------------------------------
!  STORE EXTERNAL VARIABLES INTO ONE-DIMENSIONAL ARRAY
!-----------------------------------------------------------------------
!
#ifdef QLOBERL
#ifdef FORCE_DWLW
      WRITE(IO_STDOUT,*) 'WARNING: OPTIONS QLOBERL &
	  &  AND FORCE_DWLW MAKE NO SENSE'
#endif
#endif
      K=0
      DO 82 J=1,M
      DO 82 I=1,L
!
        IF(OM(I,J) .LT. 0.5) GOTO 82
        IF(A(I,J,LRHS) .LT. 1.E-6) GOTO 82
!
         K=K+1
         HICE1(K) = HICE(I,J)
         ALB1 (K) = ALB (I,J)
         A21  (K) = A2  (I,J)
         T1   (K) = T   (I,J)+TMELT
         TA1  (K) = TAIR(I,J)+TMELT
         TD1  (K) = TD  (I,J)
         ACL1 (K) = ACL (I,J)
         PA1  (K) = PA  (I,J)
         UG1  (K) = MAX(UG(I,J),1.E-6)
         SLN1 (K) = SLN (I,J)
!
         XLAT(K)=MIN(ABS(ALAT(I,J)),60.)
!
 82   CONTINUE
!
      IF(K.EQ.0) GOTO 87
!
!-----------------------------------------------------------------------
!  COMPUTE WATER VAPOR PRESSURE AT 2M AND AT ICE SURFACE
!-----------------------------------------------------------------------
!HH      CALL VAPOR(TD1,ESTA,1,K)
!HH      CALL VAPOR(T1 ,ESTI,2,K)
!HH   2M
      DO N=1,K
      ESTA(N)=                                                          &
     &         611.21*EXP((18.729-(TD1(N)-TMELT)/227.3)*                &
     &        (TD1(N)-TMELT)/(TD1(N)-TMELT+257.87))
!HH   ICE SURFACE
      ESTI(N)=                                                          &
     &         611.15*EXP((23.036-(MIN(T1(N),TMELT)-TMELT)/333.7)*      &
     &     (MIN(T1(N),TMELT)-TMELT)/(MAX(T1(N),200.)-TMELT+279.82))
      ENDDO
!
!
!
      DO N=1,K
       DRAGL(N)=1.0
       DRAGS(N)=1.0
       SPHUMIDA(N)=0.
       SPHUMIDO(N)=0.
      ENDDO
#ifdef DASILVA
      DO N=1,K
       SPHUMIDA(N)=0.622*ESTA(N)/(1.E5-0.378*ESTA(N))
       SPHUMIDO(N)=0.622*ESTI(N)/(1.E5-0.378*ESTI(N))
      ENDDO

      CALL VARDRAG(DRAGL,DRAGS,K,TD1,T1,SPHUMIDA,SPHUMIDO,UG1)

!
#endif /*DASILVA*/
!
#ifdef BULK_KARA
      DO N=1,K
       SPHUMIDA(N)=0.62197*ESTA(N)/(PA1(N)-0.378*ESTA(N))
       SPHUMIDO(N)=0.62197*ESTI(N)/(PA1(N)-0.378*ESTI(N))
       RHOAIR_K(N)=PA1(N)/(RGAS*ta1(n)*(1.0+0.61*SPHUMIDA(N)))
       UGG=MAX(2.5,MIN(32.5,UG1(N)) )
       DRAGL(N)=MAX(0.5e-3,1.e-3*( (0.8195+0.0506*ugg-0.0009*ugg*ugg)  &
     &    + (-0.0154+0.5698/ugg-0.6743/(ugg*ugg))*(t1(n)-ta1(n))) )
       DRAGL(N)=MIN(DRAGL(N),3.0e-3)
       DRAGS(N)=0.96*DRAGL(N)
      ENDDO
#endif
!
!     DO N=1,K
!      WRITE(IO_STDOUT,*) DRAGL(N),T1(N),TA1(N)
!     ENDDO
!-----------------------------------------------------------------------
!  COMPUTE CLOUDINESS FACTOR FOR DOWNGOING LONGWAVE RADIATION
!               (KOCH 1988; BEITR.PHYS.ATMOSPH.,61(4),344-354);
!                ABSORBED SURFACE RADIATION;
!          OR    BERLIAND, 1952
!-----------------------------------------------------------------------
!
      DO 31 N=1,K
!
!UWE   MAKE ALBEDO TEMPERATURE DEPENDENT
!
!         FLAG=(0.5+SIGN(0.5,T1(N)-TMELT))
       FLAG=1./(1.+ALBTRANS*(T1(N)-TMELT)**2)
       ALB1(N)=FLAG*(A21(N)*ALBSNM+(1.-A21(N))*ALBM)                    &
     &        +(1.-FLAG)*(A21(N)*ALBSN+(1.-A21(N))*ALBI)
       FLAI (N) = (1.-ALB1(N))*SLN1(N)
!
#ifndef FORCE_DWLW
#ifndef QLOBERL
       FAKTS(N) = 1. + 0.3*ACL1(N)**2
#else
       FAKTS(N)=(1.-(0.5+0.4/90.*XLAT(N))*(ACL1(N)**2))
#endif /*QLOBERL*/
#endif /*FORCE_DWLW*/
!
!
   31 CONTINUE
!
!-----------------------------------------------------------------------
!  MAKE FIRST GUESS FOR SURFACE TEMPERATURE AND HEAT BALANCE
!-----------------------------------------------------------------------
!
      DO 33 N=1,K
!
         STP(N)=T1(N)
       FLAG=1./(1.+ALBTRANS*(STP(N)-TMELT)**2)
       ALB1(N)=FLAG*(A21(N)*ALBSNM+(1.-A21(N))*ALBM)                    &
     &        +(1.-FLAG)*(A21(N)*ALBSN+(1.-A21(N))*ALBI)
       FLAI (N) = (1.-ALB1(N))*SLN1(N)
!
#ifndef QLOBERL
#ifdef FORCE_DWLW
         Q33=ACL1(N)
#else
         FEU   = 0.601 + 5.95*1.E-7*ESTA(N)*EXP(1500./TA1(N))
         Q33 =FAKTS(N)*FEU*D3*TA1(N)**4
#endif /* FORCE_DWLW*/
!
         FP(N) = D3*STP(N)**4-FLAI(N)-Q33                               &
#ifdef BULK_KARA
     &  -RHOAIR_K(N)*CPA*UG1(N)*(TA1(N)-STP(N))*DRAGS(N)*FR_FAC         &
     &  -RHOAIR_K(N)*SUBL*UG1(N)*(SPHUMIDA(N)-SPHUMIDO(N))              &
     &   *FR_FAC*DRAGL(N)                                               &
#else
     &          -D1*UG1(N)*(TA1(N)-STP(N))*DRAGS(N)                     &
     &          -D2I*UG1(N)*(ESTA(N)-ESTI(N))*0.623/PA1(N)              &
     &                *DRAGL(N)                                         &
#endif /*BULK_KARA*/
     &          +(STP(N)-TB)/HICE1(N)*CON
!
#else /*QLOBERL*/
         FEU   = 0.39-0.05*SQRT(ESTA(N)/100.)
         FP(N) = -FLAI(N)+(FAKTS(N)*FEU*D3*TA1(N)**4                    &
     &           +4.*D3*(TA1(N)**3)*(STP(N)-TA1(N)))                    &
#ifdef BULK_KARA
     &  -RHOAIR_K(N)*CPA*UG1(N)*(TA1(N)-STP(N))*DRAGS(N)*FR_FAC         &
     &  -RHOAIR_K(N)*SUBL*UG1(N)*(SPHUMIDA(N)-SPHUMIDO(N))              &
     &   *FR_FAC*DRAGL(N)                                               &
#else
     &          -D1*UG1(N)*(TA1(N)-STP(N))*DRAGS(N)                     &
     &          -D2I*UG1(N)*(ESTA(N)-ESTI(N))*0.623/PA1(N)              &
     &                *DRAGL(N)                                         &
#endif /*BULK_KARA*/
     &          +(STP(N)-TB)/HICE1(N)*CON
#endif /*QLOBERL*/
         T1(N) = T1(N)+1.
!
 33   CONTINUE
!
!-----------------------------------------------------------------------
!  CALCULATE THE SURFACE TEMPERATURE (START OF ITERATION PROCEDURE)
!-----------------------------------------------------------------------
!
!
!
      DO 3 ITER=1,IMAX
!
      DO N=1,K
        ESTI(N)=                                                        &
     &  611.15*EXP((23.036-(MIN(T1(N),TMELT)-TMELT)/333.7)*             &
     &  (MIN(T1(N),TMELT)-TMELT)/(MAX(T1(N),200.)-TMELT+279.82))
#ifdef BULK_KARA
       SPHUMIDO(N)=0.62197*ESTI(N)/(1.E5-0.378*ESTI(N))
#endif
      ENDDO
!
             DO 34 N=1,K
!
            STPP(N)=STP(N)
            FPP (N)=FP(N)
            STP (N)=T1(N)
       FLAG=1./(1.+ALBTRANS*(STP(N)-TMELT)**2)
       FLAI (N) = (1.-ALB1(N))*SLN1(N)
#ifndef QLOBERL
#ifdef FORCE_DWLW
        Q33=ACL1(N)
#else
        FEU    = 0.601+5.95*1.E-7*ESTA(N)*EXP(1500./TA1(N))
        Q33=FAKTS(N)*FEU*D3*TA1(N)**4
#endif /*FORCE_DWLW*/
            FP(N)  = D3*STP(N)**4-FLAI(N)-Q33                           &
#ifdef BULK_KARA
     &  -RHOAIR_K(N)*CPA*UG1(N)*(TA1(N) -STP(N) )*DRAGS(N)*FR_FAC       &
     &  -RHOAIR_K(N)*SUBL*UG1(N)*(SPHUMIDA(N)-SPHUMIDO(N))              &
     &   *FR_FAC*DRAGL(N)                                               &
#else 
     &              -D1 *UG1(N)*(TA1(N) -STP(N) )*DRAGS(N)              &
     &              -D2I*UG1(N)*(ESTA(N)-ESTI(N))*0.623/PA1(N)          &
     &                *DRAGL(N)                                         &
#endif /*BULK_KARA*/
     &              +(STP(N)-TB)/HICE1(N)*CON
#else /*QLOBERL*/
            FEU    = 0.39-0.05*SQRT(ESTA(N)/100.)
            FP(N)  = -FLAI(N)+(FAKTS(N)*FEU*D3*TA1(N)**4                &
     &               +4.*D3*(TA1(N)**3)*(STP(N)-TA1(N)))                &
#ifdef BULK_KARA
     & -RHOAIR_K(N)*CPA*UG1(N)*(TA1(N)-STP(N))*DRAGS(N)*FR_FAC          &
     & -RHOAIR_K(N)*SUBL*UG1(N)*(SPHUMIDA(N)-SPHUMIDO(N))               &
     &  *FR_FAC*DRAGL(N)                                                &
#else 
     &               -D1*UG1(N)*(TA1(N)-STP(N))*DRAGS(N)                &
     &               -D2I*UG1(N)*(ESTA(N)-ESTI(N))* 0.623/PA1(N)        &
     &                *DRAGL(N)                                         &
#endif /*BULK_KARA*/
     &               +(STP(N)-TB)/HICE1(N)*CON
#endif /*QLOBERL*/
!
            FDIFF  = FP(N)-FPP(N)
            T1(N)  = MAX(STP(N)-(STP(N)-STPP(N))*FP(N)/                 &
     &               MAX(ABS(FDIFF),1.E-10)*SIGN(1.,FDIFF),100.)
            DIFF(N)= T1(N)-STP(N)
!
 34      CONTINUE
    3 CONTINUE
!
!-----------------------------------------------------------------------
!  CALCULATE GROWTH RATES WITH UPDATED HEAT BALANCE EQUATION:
!            FHS IS MELTING RATE AT THE SNOW/ICE SURFACE
!            FHB IS MELTING/GROWTH RATE AT THE BOTTOM OF THE ICE.
!-----------------------------------------------------------------------
!
      DO 83 N=1,K
!
         T1(N)=MIN(T1(N),TMELT)
#ifndef QLOBERL 
         FEU   = 0.601+5.95*1.E-7*ESTA(N)*EXP(1500./TA1(N))
#else
         FEU   = 0.39-0.05*SQRT(ESTA(N)/100.)
#endif
         A1       = (0.5+SIGN(0.5,T1(N)-TMELT))
!
#ifdef BULK_KARA
       FLSE1(N) = RHOAIR_K(N)*CPA*UG1(N)*(TA1(N)-T1(N))*DRAGS(N)*FR_FAC
       FLLA1(N) = RHOAIR_K(N)*SUBL*UG1(N)*(SPHUMIDA(N)-SPHUMIDO(N))     &
     &            *FR_FAC*DRAGL(N)
#else
         FLSE1(N) = D1*UG1(N)*(TA1(N)-T1(N))*DRAGS(N)
         FLLA1(N) = D2I*UG1(N)*(ESTA(N)-ESTI(N))*0.623/PA1(N)*DRAGL(N)
#endif /*BULK_KARA*/
!
#ifndef QLOBERL 
         Q1   (N) = D3*T1(N)**4
#else
         Q1   (N) = 0.
#endif /*QLOBERL*/
!
         Q2   (N) =  (1.-ALB1(N))*SLN1(N)
!
#ifdef FORCE_DWLW
         Q3(N)=ACL1(N)
#else /*FORCE_DWLW*/
#ifndef QLOBERL 
         Q3   (N) = FAKTS(N)*FEU*D3*TA1(N)**4
#else
         Q3   (N) = -(FAKTS(N)*FEU*D3*TA1(N)**4                         &
     &              +4.*D3*(TA1(N)**3)*(T1(N)-TA1(N)))
#endif /*QLOBERL*/
#endif /*FORCE_DWLW*/
!
         FHS1(N)  = (Q1(N)-Q2(N)-Q3(N)-FLSE1(N)-FLLA1(N)                &
     &             -(TB-T1(N))/HICE1(N)*CON)/CLB
!
         FHB1(N)  =( (TB-T1(N))/HICE1(N)*CON )/CLB
!
 83   CONTINUE
!-----------------------------------------------------------------------
!  UNSCRAMBLE FOR TWO-DIMENSIONAL FIELD
!-----------------------------------------------------------------------
!
      K=0
      DO 84 J=1,M
      DO 84 I=1,L
!
         FHS(I,J)=0.0
         FHB(I,J)=0.0
         QSW(I,J)=0.0
         QLW(I,J)=0.0
         QSE(I,J)=0.0
         QLA(I,J)=0.0
!
         IF(OM(I,J).LT. 0.5) GOTO 84
         IF(A (I,J,LRHS).LT.1.E-6) GOTO 84
         K = K+1
!
         FHS(I,J) = FHS1 (K)
         FHB(I,J) = FHB1 (K)
         T  (I,J) = T1   (K)-TMELT
         QSW(I,J) = Q2   (K)
         QLW(I,J) = Q3   (K)-Q1(K)
         QSE(I,J) = FLSE1(K)
         QLA(I,J) = FLLA1(K)
!
 84   CONTINUE
!
 87   RETURN
      END
