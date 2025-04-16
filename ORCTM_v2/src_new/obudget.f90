      SUBROUTINE OBUDGET(ALAT,QT,TAIR,TD,ACL,PA,UG,SLN,FO,OM,           &
     &                   QSW,QLW,QSE,QLA)
!
!=======================================================================
!  PROGRAMMED BY:
!  --------------
!     A.STOESSEL                 MPI, HAMBURG                      1989
!
!  MODIFIED
!  --------
!  S.LEGUTKE                     *DKRZ*          13.8.95
!    - NET ATMOSPHERIC HEAT FLUX USED FOR SURFACE TEMP. > TFREZ
!
!  HELMUTH 19.04.2000
!  - LONGWAVE-FLUX FROM BERLIAND(1952)
!
!
!  PURPOSE:
!  --------
!     CALCULATES GROWTH RATES OF NEW ICE IN THE ICE FREE PART OF A GRID
!     CELL WITH BULK FORMULAS OR NET ATM. HEAT FLUX
!
!  METHOD:
!  -------
!     HEAT BUDGET EQUATION FOR OPEN WATER
!
!  INTERFACE:
!  ----------
!     INPUT:
!     *QT*       SEA SURFACE=OML TEMPERATURE [DEG C]
!     *TAIR*     2M AIR TEMPERATURES         [DEG C]
!UWE     *SAIR*     SURFACE TEMPERATURES        [K]
!     *TD*       DEW POINT TEMPERATURES      [K]
!     *ACL*      FRACTIONAL CLOUD COVER      [FRAC.]
!     FOR CPP OPTION FORCE_DWLW ACL CONTAINS DOWNWARD LONGWAVE RADIATION!!!!!
!     *PA*       ATMOSPHERIC SURFACE PRESURE [PA]
!     *UG*       2 M WIND SPEED              [M/SEC]
!     *SLN*      INCOMING SURFACE SOLAR RAD. [W/M**2]
!     *OM*       LAND/SEA MASK               [0/1]
!     *A*        ICE COMPACTNESS             [FRAC.]
!UWE     *AOHFLX*   ATMOSPHERIC HEAT FLUX       [W/M**2]
!UWE     *AODFLX*   HEAT FLUX CORRECTION        [W/M**2]
!
!     OUTPUT:
!     *FO*       GROWTH RATE OF THIN ICE     [M/SEC]
!     *QSW*      ABSORBED SOLAR RADIATION    [W/M**2]
!     *QLW*      OUTGOING LONGWAVE HEAT FLUX [W/M**2]
!     *QSE*      SENSIBLE HEAT FLUX          [W/M**2]
!     *QLA*      LATENT HEAT FLUX            [W/M**2]
!
!
!  EXTERNALS:
!  ----------
!     VAPOR: CALCULATES VAPOR PRESSURE
!
!=======================================================================
      USE MO_PARAM1
      USE MO_COMMOAU1
      DIMENSION QT(IE,JE),FO(IE,JE),SLN(IE,JE)
!UWE      DIMENSION AOHFLX(IE,JE)
      DIMENSION TAIR(IE,JE),TD(IE,JE),ACL(IE,JE),PA(IE,JE),UG(IE,JE)
      DIMENSION OM(IE,JE)
!UWE     1         ,SAIR(IE,JE)
      DIMENSION TA1(IE*JE),TD1(IE*JE),ACL1(IE*JE),PA1(IE*JE),UG1(IE*JE) &
     &         ,SLN1(IE*JE),FLSE1(IE*JE),FLLA1(IE*JE),FH1(IE*JE)        &
     &         ,QT1(IE*JE),ESTA(IE*JE),ESTW(IE*JE),FAKTS(IE*JE)
!UWE     3         ,SA1(IE*JE),AODFLX1(IE*JE),AOHFLX1(IE*JE)
      DIMENSION QSW(IE,JE),QLW(IE,JE),QSE(IE,JE),QLA(IE,JE)             &
     &            ,DRAGS(IE*JE),DRAGL(IE*JE)                            &
     &            ,SPHUMIDO(IE*JE),SPHUMIDA(IE*JE)
      DIMENSION Q1(IE*JE),Q2(IE*JE),Q3(IE*JE)
      DIMENSION XLAT(IE*JE),ALAT(IE,JE)
!
!-----------------------------------------------------------------------
!  SELECT THE GRID CELLS AND STORE THEM IN SELECTIVE 1-DIMENSIONAL ARRAY
!-----------------------------------------------------------------------
!
#ifdef FORCE_DWLW
#ifdef QLOBERL
      WRITE(IO_STDOUT,*) 'WARNING: CPP OPTIONS &
	  & QLOBERL AND FORCE_DWLW TOGETHER'
#endif /*QLOBERL*/
#endif /*FORCE_DWLW*/
!
      L = IE
      M = JE
      LMDP = IE*JE

      RGAS=287.1
      CPA=1004.67
!:: Frank Roeske's close budget factor:
      FR_FAC=1.1925
      ALMZER=1.E-19
#ifdef FORCE_DWLW
      D3=0.97*5.67E-8
#endif
!OtB initialize fields to 0.
      DO I=1,LMDP
        Q1(I) = .0
        Q2(I) = .0
        Q3(I) = .0
        SLN1(I) = .0
      ENDDO
!
      K=0
      DO 93 J=1,M
      DO 93 I=1,L
!
         IF(OM(I,J).LT.0.5) GOTO 93
!
         K = K+1
!
         QT1    (K) = QT    (I,J)+TMELT
         TA1    (K) = TAIR  (I,J)+TMELT
         TD1    (K) = TD    (I,J)
!UWE         SA1    (K) = SAIR  (I,J)
         ACL1   (K) = ACL   (I,J)
         PA1    (K) = PA    (I,J)
         UG1    (K) = MAX(UG(I,J),1.E-6)
         SLN1   (K) = SLN   (I,J)
!UWE         AOHFLX1(K) = AOHFLX(I,J)
!UWE         AODFLX1(K) = AODFLX(I,J)
         XLAT(K)=MIN(ABS(ALAT(I,J)),60.)
 93   CONTINUE
!
!-----------------------------------------------------------------------
!  COMPUTE WATER VAPOR PRESSURE AT 2M AND THE SATURATION VAPOR PRESSURE 
!   AT THE SEA SURFACE
!-----------------------------------------------------------------------
!hh      CALL VAPOR(TD1,ESTA,1,K)
!hh      CALL VAPOR(QT1,ESTW,3,K)
!HH   2M
      DO N=1,K
      ESTA(N)=                                                          &
     &    611.21*EXP((18.729-(TD1(N)-TMELT)/227.3)*                     &
     &    (TD1(N)-TMELT)/(TD1(N)-TMELT+257.87))
!HH   SEA SURFACE
      ESTW(N)=                                                          &
     &  0.9815*611.21*EXP((18.729-(QT1(N)-TMELT)/227.3)*                &
     &    (QT1(N)-TMELT)/(QT1(N)-TMELT+257.87))
      ENDDO
!
!
!-----------------------------------------------------------------------
!  COMPUTE CLOUDINESS FACTOR FOR DOWNGOING LONGWAVE RADIATION
!               (KOCH 1988; BEITR.PHYS.ATMOSPH.,61(4),344-354);
!-----------------------------------------------------------------------
      TFREZK=TMELT+TFREZ

      DO n=1,k
       dragl(n)=1.0
       drags(n)=1.0
       sphumida(n)=0.0
       sphumido(n)=0.0
      ENDDO
#ifdef DASILVA
      DO N=1,K
       SPHUMIDA(N)=0.622*ESTA(N)/(1.E5-0.378*ESTA(N))
       SPHUMIDO(N)=0.622*ESTW(N)/(1.E5-0.378*ESTW(N))
      ENDDO
!
      CALL VARDRAG(DRAGL,DRAGS,K,QT1,TA1,SPHUMIDA,SPHUMIDO,UG1)

#endif
#ifdef BULK_KARA
      DO N=1,K
       SPHUMIDA(N)=0.62197*ESTA(N)/(PA1(N)-0.378*ESTA(N))
       SPHUMIDO(N)=0.62197*ESTW(N)/(PA1(N)-0.378*ESTW(N))
!
       UGG=MAX(2.5,MIN(32.5,UG1(N)) )
! BRACKET DRAG COEFF. BETWEEN 0.5e-3 AND 2e-3!
!      DRAGL(N)=1.e-3*( (0.8195+0.0506*ugg-0.0009*ugg*ugg)              &
!    &    + (-0.0154+0.5698/ugg-0.6743/(ugg*ugg))*(qt1(n)-ta1(n)))
       DRAGL(N)=MAX(0.5e-3,1.e-3*( (0.8195+0.0506*ugg-0.0009*ugg*ugg)  &
     &    + (-0.0154+0.5698/ugg-0.6743/(ugg*ugg))*(qt1(n)-ta1(n))) )
       DRAGL(N)=MIN(DRAGL(N),3.00e-3)
       DRAGS(N)=DRAGL(N)*0.96
      ENDDO

#endif /*BULK_KARA*/
!
      DO 25 N=1,K
!
#ifndef FORCE_DWLW
#ifndef QLOBERL
         FAKTS(N)=1.+0.3*ACL1(N)**2
#else
         FAKTS(N)=(1.-(0.5+0.4/90.*XLAT(N))*(ACL1(N)**2))
#endif /*QLOBERL*/
#endif /*FORCE_DWLW*/
!
!-----------------------------------------------------------------------
!  CALCULATE HEAT FLUXES AND GROWTH RATES:
!            USE HEAT BALANCE EQUATION WHERE THE ATMOSPHERE HAD ICE
!            USE NET ATMOSPHERIC HEAT FLUX + FLUX CORRECTION ELSEWHERE
!-----------------------------------------------------------------------
!
!-ET      TFREZK=TMELT+TFREZ
!-ET      DO 25 N=1,K
#ifndef FORCE_DWLW
#ifndef QLOBERL
       FEU      = 0.605+5.95*1.E-7*ESTA(N)*EXP(1500./TA1(N))
#else
       FEU      = 0.39-0.05*SQRT(ESTA(N)/100.)
#endif /*QLOBERL*/
#endif /*FORCE_DWLW*/
!
#ifdef BULK_KARA
!
      RHOAIR_K=PA1(N)/(RGAS*TA1(N)*(1.0+0.61*SPHUMIDA(N)) )
      FLSE1(N) = RHOAIR_K*CPA*UG1(N)*(TA1 (N)-QT1(N) )*DRAGS(N)*FR_FAC
!
      FLLA1(N) = (2.501-0.00237*(qt1(n)-TMELT))*1.e+6                   &
     & *RHOAIR_K*UG1(N)*(SPHUMIDA(N)-SPHUMIDO(N))*DRAGL(N)*FR_FAC
#else
       FLSE1(N) = D1 *UG1(N)*(TA1 (N)-QT1(N) )                          &
     &            *DRAGS(N)
!
       FLLA1(N) = D2W*UG1(N)*(ESTA(N)-ESTW(N))*.623/PA1(N)              &
     &            *DRAGL(N)
#endif /*BULK_KARA*/
!
#ifndef QLOBERL
       Q1   (N) = D3*QT1(N)**4
#else
       Q1   (N) = 0.
#endif
       Q2   (N) = (1.-ALBW)*SLN1(N)
!
#ifdef FORCE_DWLW
       Q3(N)=ACL1(N)
#else
#ifndef QLOBERL
       Q3   (N) = FAKTS(N)*FEU*D3*TA1(N)**4
#else
       Q3   (N) = -(FAKTS(N)*FEU*D3*TA1(N)**4                           &
     &              +4.*D3*(TA1(N)**3)*(QT1(N)-TA1(N)))
#endif /*QLOBERL*/
#endif /*FORCE_DWLW*/
       FH1  (N) = (Q1(N)-Q2(N)-Q3(N)-FLSE1(N)-FLLA1(N))/CLB
!
   25 CONTINUE
!-----------------------------------------------------------------------
!  UNSCRAMBLE INTO TWO-DIMENSIONAL FIELD
!-----------------------------------------------------------------------
!
      K=0
      DO 81 J=1,M
      DO 81 I=1,L
!
         QSW(I,J) = 0.0
         QLW(I,J) = 0.0
         QSE(I,J) = 0.0
         QLA(I,J) = 0.0
         FO (I,J) = 0.0
!
       IF(OM(I,J).LT.0.5) GOTO 81
!
       K=K+1
!
       FO (I,J) = FH1(K)
       QSW(I,J) = Q2(K)
       QLW(I,J) = Q3(K)-Q1(K)
       QSE(I,J) = FLSE1(K)
       QLA(I,J) = FLLA1(K)
!
   81 CONTINUE
!
      RETURN
      END
