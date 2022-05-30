      SUBROUTINE OCTHER(i_step)
      USE MO_PARAM1
      USE mo_mpi
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMOAU1
      USE MO_COMMOAU2
      USE MO_UNITS
      USE MO_OBCS
#ifdef MEAN
#ifdef DIFFDIAG
      USE MO_MEAN
#endif /*DIFFDIAG*/
#endif /*MEAN*/      
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!     SBR OCTHER COMPUTES
!
!          A) BOUNDARY FORCING ON SALT, TEMPERATURE AND ZETA
!          B) BAROCLINIC PRESSURE IN EACH LAYER
!          C) CONVECTIVE ADJUSTMENT
!          D) RICHARDSON-NUMBER DEPENDING COEFFICIENTS FOR
!             VERTICAL DIFFUSION OF MOMENTUM (AV) AND
!             TEMPERATURE AND SALINITY (DV)
!
!          LAST MODIFIED:
!          UM 13.12.99 MASS CONSERVATION SURFACE LAYER
!
!
!-----------------------------------------------------------------------
!
      INTEGER I,ILEN,J,K,KO,KU,N,I_STEP
      REAL RHUPPO(IE,JE),ZSURF(KE),CONVEFO(IE,JE)
      REAL CRA,CRD,CTMPOTB,DISTI,DRDZ0,DRDZO,DUDO,DUDZ,GOR,HELPSWI,HHO
      REAL REISCON,RELAX,RELNE,RINUMO,SQ,STABEPS,SWITC,THERMIN,THICK
      REAL TOPNAH,TQ,WPENDEP,WTDECAY,ZZZDZ,DTTS,DRIV

!     A)
!
!     BOUNDARY FORCING ON TEMPERATURE, SALT AND ZETA
!
!     RELAXATION TIME ON SALINITY   : 1./ RELSAL
!
      ZSURF(1)=1.
      DO K=2,KE
       ZSURF(K)=0.
      ENDDO
!
      STEBOL=3.E-8
      THERMIN=4.E6*DZW(1)
!
!      CALL OCICE
!
      DO J=1,JE
      DO I=1,IE
        EMINPO(I,J)=0.
        REISCON=0.
        IF(SICOMO(I,J).LE.0.01) THEN
          REISCON=1.
        ENDIF
!
#ifdef EISREST
        REISCON=1.-SICOMO(I,J)
#endif
!
        OLDSO=SAO(I,J,1)
!
!:: PERFORM SSS RESTORING ONLY in NON COUPLED CASE
        SAO(I,J,1)=SAO(I,J,1)+DT*RELSAL*REISCON*                        &
     &            (RELSAO(I,J)-SAO(I,J,1))
        THO(I,J,1)=THO(I,J,1)+DT*RELTEM*                                &
     &            (RELTHO(I,J)-THO(I,J,1))
!
       EMINPO(I,J)=(DDPO(I,J,1)+ZO(I,J)-SICTHO(I,J)*RHOICWA             &
     &                                 -SICSNO(I,J)*RHOSNWA)            &
     &            *(MAX(OLDSO,1.E-3)/MAX(SAO(I,J,1),1.E-3)-1.)
!
      ENDDO   
      ENDDO
!
      DTI=1./DT
      DO J=1,JE
      DO I=1,IE
      ZO(I,J)=(ZO(I,J)+EMINPO(I,J))*WETO(I,J,1)
      EMINPO(I,J)=EMINPO(I,J)*DTI
      ENDDO    
      ENDDO
!
!UWE
!     RIVER INPUTS

      DO j=1,jE
      DO i=1,ie   
        ZZZDZ=MAX(ALMZER,DDPO(I,J,1)+ZO(I,J)-SICTHO(I,J)*RHOICWA        &
     &            -SICSNO(I,J)*RHOSNWA)
!
!        DRIV=GIRIV(I,J)*DT/(DLXP(I,J)*DLYP(I,J))
        DRIV=GIRIV(I,J)*DT
!  changed by peterspy, wrong unit for GIRIV

!       USE ACTUAL LAYERTHICKNESS FOR MASS/SALT CONSERVATION
        SAO(I,J,1)=SAO(I,J,1)*ZZZDZ/(ZZZDZ+DRIV)
        ZO(I,J)=ZO(I,J)+DRIV
        PRECO(I,J)=PRECO(I,J)+DRIV/DT
        PRECH(I,J)=PRECH(I,J)+DRIV/DT
        RIVRUN(I,J)=DRIV/DT
      ENDDO
      ENDDO
!
!=====================================================================
!
!  OBCS : BOUNDARY TS FORCINGS
!         SPONGE NUDGED CLIM
#ifdef OBC_TS_FORC
      CALL FORC_OBCS_TS(I_STEP)
      
!      CALL PRED_OBCS_TS(THO,SAO)
      
#ifdef OBC_SPONGE_EXP
      CALL SPONGE_OBCS_TS_EXP(THO,SAO) 
#else
      CALL SPONGE_OBCS_TS(THO,SAO) 
#endif
#endif

!     B)
!
!     BAROCLINIC PRESSURE AND STABILITY
!
!
!---------------------------------------------------------------------
!
!     B.1) UPPER LAYER
!

!$OMP PARALLEL PRIVATE(i,j,k,ko,ku, &
!$OMP   thick, disti, tq, sq, switc, helpswi, tupper, tlower, &
!$OMP   supper, slower, ddhelp, ctmpotb, dtts, &
!$OMP   relne, relax, dudz, drdz0, gor, cra, crd, wpendep, wtdecay, &
!$OMP   stabeps, topnah, drdzo, dudo, hho, rinumo)

!$OMP DO
      DO J=1,JE
!
      DO 42 I=1,IE
      SHELP(I,J)=SAO(I,J,1)
      THELP(I,J)=THO(I,J,1)
42    CONTINUE
!
      CALL ADISITJ(THELP,SHELP,PREFF(1),J)
      CALL RHO1J(THELP,SHELP,PREFF(1),RHELP,J)
!
      DO 43 I=1,IE
      STABIO(I,J,1) =0.
      PO(I,J,1)     = G*TIESTU(1)*0.00097*RHELP(I,J)
      S1O(I,J,1)=RHELP(I,J)
43    CONTINUE
!
!---------------------------------------------------------------------
!
!     B.2) ALL OTHER LAYERS ==> K = 2 , KE
!
!
!
!-ET      CONVEREI=0.
!
      DO 457 K=2,KE
!
      DISTI=1./DZ(K)
!
      DO 46 I=1,IE
      SHELP(I,J)=SAO(I,J,K)
      THELP(I,J)=THO(I,J,K)
46    CONTINUE
!
      CALL ADISITJ(THELP,SHELP,PREFF(K),J)
      CALL RHO1J(THELP,SHELP,PREFF(K),RHELP,J)
!
      DO 47 I=1,IE
      SHELP(I,J)=SAO(I,J,K-1)
      THELP(I,J)=THO(I,J,K-1)
47    CONTINUE
!
      CALL ADISITJ(THELP,SHELP,PREFF(K),J)
      CALL RHO1J(THELP,SHELP,PREFF(K),RHUPPO,J)
!
      DO 48 I=1,IE
      S1O(I,J,K)=RHELP(I,J)
      STABIO(I,J,K) = DISTI * ( RHELP(I,J) - RHUPPO(I,J) )
      PO(I,J,K) = PO(I,J,K-1) + G*DZ(K)*0.00049*(RHELP(I,J)             &
     & +RHUPPO(I,J))
!
#ifndef PLUME
#ifndef UMKLAP
!UWE     SALT CONSERVATION 12/99
      TQ=((DDPO(I,J,K-1)+ZSURF(K-1)*(ZO(I,J)-SICTHO(I,J)*RHOICWA        &
     &                                      -SICSNO(I,J)*RHOSNWA))      &
     &      *THO(I,J,K-1)+DDPO(I,J,K)*THO(I,J,K))                       &
     &  /(DDPO(I,J,K)+DDPO(I,J,K-1)                                     &
     &   +ZSURF(K-1)*(ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)  &
     &                    +(1.-WETO(I,J,K)))
      SQ=((DDPO(I,J,K-1)+ZSURF(K-1)*(ZO(I,J)-SICTHO(I,J)*RHOICWA        &
     &                                      -SICSNO(I,J)*RHOSNWA))      &
     &    *SAO(I,J,K-1)+DDPO(I,J,K)*SAO(I,J,K))                         &
     &  /(DDPO(I,J,K)+DDPO(I,J,K-1)                                     &
     &   +ZSURF(K-1)*(ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)  &
     &                    +(1.-WETO(I,J,K)))
#endif /*UMKLAP*/
!
      SWITC=(HALF-SIGN(HALF,STABIO(I,J,K)))*WETO(I,J,K)
      SWITC=MAX(0.,-STABIO(I,J,K)/(1.E-11+ABS(STABIO(I,J,K))))          &
     &     *WETO(I,J,K)
!-ET      CONVEREI=CONVEREI+SWITC
!
#ifdef NURDIF
      HELPSWI=SWITC
      SWITC=0.
#endif /*NURDIF*/
!
#ifndef UMKLAP
      THO(I,J,K-1) = TQ * SWITC + (ONE-SWITC) * THO(I,J,K-1)
      THO(I,J,K)   = TQ * SWITC + (ONE-SWITC) * THO(I,J,K)
      SAO(I,J,K-1) = SQ * SWITC + (ONE-SWITC) * SAO(I,J,K-1)
      SAO(I,J,K)   = SQ * SWITC + (ONE-SWITC) * SAO(I,J,K)
#else /*UMKLAP*/
         IF (SWITC.GE.0.5) THEN
           TUPPER=THO(I,J,K-1)
           TLOWER=THO(I,J,K)
           SUPPER=SAO(I,J,K-1)
           SLOWER=SAO(I,J,K)
           DDHELP=DDPO(I,J,K-1)+ZSURF(K-1)*(ZO(I,J)-SICTHO(I,J)*RHOICWA &
     &             -SICSNO(I,J)*RHOSNWA)
!
           IF(DDPO(I,J,K).GT.DDHELP) THEN
            THO(I,J,K-1)=TLOWER
            SAO(I,J,K-1)=SLOWER
            THO(I,J,K)=TLOWER+(TUPPER-TLOWER)                           &
     &                 *(DDHELP/DDPO(I,J,K))
            SAO(I,J,K)=SLOWER+(SUPPER-SLOWER)                           &
     &                 *(DDHELP/DDPO(I,J,K))
           ELSE
            THO(I,J,K)=TUPPER
            SAO(I,J,K)=SUPPER
            THO(I,J,K-1)=TUPPER+(TLOWER-TUPPER)                         &
     &                 *(DDPO(I,J,K)/DDHELP)
            SAO(I,J,K-1)=SUPPER+(SLOWER-SUPPER)                         &
     &                 *(DDPO(I,J,K)/DDHELP)
           ENDIF
!
         ENDIF
#endif /*UMKLAP*/
!
      STABIO(I,J,K)=(ONE-SWITC)*STABIO(I,J,K)
      STABIO(I,J,K)=MAX(STABIO(I,J,K),0.)
!
#ifdef NURDIF
      SWITC=HELPSWI
#endif /*NURDIF*/
!
      IF(KCONDEP(I,J).EQ.K-1)KCONDEP(I,J)=KCONDEP(I,J)+NINT(SWITC)
#endif /*PLUME*/
48    CONTINUE
457   CONTINUE

      ENDDO ! J-Loop
!
!$OMP SINGLE
#ifdef PLUME
!SJM PLUME CONVECTION
      DTTS=DT
      CALL NLOPPS(THO,SAO,DZW,DDPO,PREFF,KCONDEP,DTTS)
#endif /*PLUME*/
!
!$OMP END SINGLE
!
!
!
!======================================================================
!
!     D)
!
!                   CALCULATION OF RICHARDSON NUMBER DEPENDENT
!                      VERTICAL EDDY VISCOSITY   AV     AND
!                      VERTICAL EDDY DIFFUSIVITY DV
!
!
!         RINUM : ( (G/RHO)*D(RHO)/DZ ) / ( (D(VELOCITY)/DZ)**2 )
!                  RICHARDSON NUMBER (EVEN,ODD ==> RINUME,RINUMO)
!         GOR   : G/RHO
!
!         AV0   : NUMERICAL VALUE OF VERTICAL EDDY VISCOSITY IN CASE
!                 OF NEUTRAL STABILITY, I.E. FREE TURBULENCE
!
!---------------------------------------------------------------------
!
!
      RELNE=0.4
      RELAX=1.-RELNE
!
      DUDZ=1.E4
      DRDZ0=1.E-3
!
      GOR=G/1025.
!
!--------------------------------------------------------------------
!
!     C.1)
!
!     VERTICAL EDDY VISCOSITY  (FOR MOMENTUM EQUATION)
!
!--------------------------------------------------------------------
!
!     D.1)    MIXED-LAYER TURBULENCE
!
!  AMPLITUDE OF TURBULENCE DECAYS BY FACTOR WTDECAY EVERY MODEL LEVEL.
!  TURBULENCE STOPS ONCE DENSITY DIFFERENCE REACHES THE EQUIVALENT OF
!   WTDT
!  TEMPERATURE DIFFERENCE.
!  TURBULENCE UNDER ICE IS / IS NOT  ENHANCED.
!                             ==
      CRA=5.
      CRD=5.
!
      WPENDEP=40.
      RELNE=0.4
      RELAX=1.-RELNE
      WTDECAY=EXP(-DZW(1)/WPENDEP)
!
!$OMP DO
      DO 4100 J=1,JE
        DO 4100 I=1,IE
#ifdef REDWMICE
          T1O(I,J,1)=WT*(1.-SICOMO(I,J))**2*FU10(I,J)**3
          S1O(I,J,1)=WT*(1.-SICOMO(I,J))**2*FU10(I,J)**3
#else
          T1O(I,J,1)=WT*(1.-SICOMO(I,J))*FU10(I,J)**3
          S1O(I,J,1)=WT*(1.-SICOMO(I,J))*FU10(I,J)**3
#endif /*REDWMICE*/
#ifdef T1O2
       T1O(I,J,1)=T1O(I,J,1)*2.
#endif
4100  CONTINUE
!
!$OMP DO
      DO J=(J_start+1),JE1
      DO I=(I_start+1),IE1

      DO K=2,KE
      KU=MIN(K+1,KE)
      KO=MAX(K-1,1)
      STABEPS=CSTABEPS/DZW(KO)
      WTDECAY=EXP(-DZW(KO)/WPENDEP)
      TOPNAH=0.0

!
      T1O(I,J,1)=T1O(I,J,1)*WTDECAY*STABEPS                             &
     &  /(STABEPS+0.5*(STABIO(I,J,K)+(1.-ZSURF(KO))*STABIO(I,J,KO)))
!
      S1O(I,J,1)=S1O(I,J,1)*WTDECAY*STABEPS                             &
     & /(STABEPS+0.5*(STABIO(I,J,K)+(1.-ZSURF(KO))*STABIO(I,J,KO)))
!
#ifdef MEAN
#ifdef DIFFDIAG
      WTMIX(I,J,K)=S1O(I,J,1)
#endif /*DIFFDIAG*/
#endif /*MEAN*/

      DRDZO =(STABIO(I,J,K)+STABIO(I+1,J,K))*G*0.001
!
      DRDZO=MAX(0.,DRDZO)
!
      DUDO=ALMZER + WETO(I,J,K) * DI(K)**2                              &
     & * (   ( UKO(I-1,J,K) - UKO(I-1,J,KO) )**2                        &
     &     + ( UKO(I,J,K)   - UKO(I,J,KO)   )**2                        &
     &     + ( VKE(I,J-1,K) - VKE(I,J-1,KO) )**2                        &
     &     + ( VKE(I,J,K)   - VKE(I,J,KO)   )**2   ) *0.5
!
      HHO=WETO(I,J,K)*(AMSUO(I,J,K)+AMSUO(I-1,J,K)+AMSUE(I,J-1,K)       &
     & +AMSUE(I,J,K))*0.25
!
      RINUMO=HHO*MAX(GOR*STABIO(I,J,K)/DUDO,0.)
!
      SWITC=(HALF-SIGN(HALF,STABIO(I,J,K)))*WETO(I,J,K)
!
      DVO(I,J,K)=(RELAX*MIN(DVO(I,J,K),DV0+S1o(i,j,1))+RELNE*(S1O(I,J,1)&
     &    +DV0/((1.+ CRD*RINUMO)**3)+DBACKV(K)+TOPNAH))*WETO(I,J,K)
!
      AVO(I,J,K)=(RELAX*MIN(AVO(I,J,K),AV0+ABACKV(K))+RELNE*(T1O(I,J,1) &
     &    +AV0/((1.+ CRA*RINUMO)**2)+ABACKV(K)+TOPNAH))*WETO(I,J,K)
!
      AVO(I,J,K)=WETO(I,J,K)*MAX(CAVOCON*(1.E-11-STABIO(I,J,K))         &
     &          /(1.E-11+ABS(STABIO(I,J,K))),AVO(I,J,K))
!
#ifndef UMKLAP
#ifndef NURMISCH
      DVO(I,J,K)=WETO(I,J,K)*MAX(CDVOCON*(1.E-11-STABIO(I,J,K))         &
     &          /(1.E-11+ABS(STABIO(I,J,K))),DVO(I,J,K))
#endif /*NURMISCH*/
#endif /*UMKLAP*/
!
      STABIO(I,J,K)=MAX(STABIO(I,J,K),0.)

      ENDDO

      ENDDO
      ENDDO ! J-Loop
!
!$OMP SINGLE
      CALL bounds_exch(S1O(:,:,1),T1O(:,:,1))
      CALL bounds_exch(DVO)
      CALL bounds_exch(AVO)
      CALL bounds_exch(STABIO)
!$OMP END SINGLE
!
!=====================================================================
!
!$OMP DO
       DO J=1,JE

        DO I=1,IE
         AVO(I,J,KEP)=0.
         DVO(I,J,KEP)=0.
         SHELP(I,J)=0.0
         THELP(I,J)=0.0
        ENDDO
!
!UWE  INCLUDE DOWNWARD PROPAGATION OF THO AND SAO IN LAND
!
       DO K=2,KE
        DO I=1,IE
         IF(WETO(I,J,K).LT.0.5)THEN
           THO(I,J,K)=THO(I,J,K-1)
           SAO(I,J,K)=SAO(I,J,K-1)
         ENDIF
        ENDDO
       ENDDO
!
!UWE  COMPUTE PRESSURE AND DENSITY NEW
!
       DO K=1,KE
       DO I=1,IE
        THELP(I,J)=THO(I,J,K)
        SHELP(I,J)=SAO(I,J,K)
       ENDDO
       CALL ADISITJ(THELP,SHELP,PREFF(K),J)
       CALL RHO1J(THELP,SHELP,PREFF(K),RHELP,J)
!
!UWE   PO NOCH SAUBER BESTIMMEN, AUS RHOO NOCH MACHEN!!!
!
       DO I=1,IE
        RHOO(I,J,K)=RHELP(I,J)
       ENDDO
       ENDDO
      ENDDO ! J-Loop
!$OMP END PARALLEL
!
      RETURN
      END
