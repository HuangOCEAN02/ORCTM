      SUBROUTINE OCICE  
!     HIBLERS ICE MODEL AS DESCRIBED IN
!     W.D.hIBLER III, A DYNAMIC THERMODYNAMIC SEA ICE MODEL. J. PHYS. OCEANOGR.
!     9, 815 - 846, 1979
!
!     MODIFICATIONS:
!     UWE 2.2.00
!       INCLUDE ADVECTION OF ICE VELOCITIES
!       CORRECT TAU-W TERM USING MIXTURE OF OLD AND NEW VELOCITIES
!        NEW FILEDS INCLUDED:
!         SPEEDU    DIFFERENCE IN SPEED BETWEEN WATER AND ICE (OLD)
!         SPEEDV
!     UWE 8.3.00
!        INCLUDE TAUWAT, PROPER DETERMINATION OF WATER ICE STRESS
!        TO BE USED IN OCWIND
!     UWE 28.6.00
!        MAKE SW-PENETRATION APPLICABLE UNDER SEA ICE TOO
!     UWE 10/00
!        INCLUDE IMPLICITE TREATMENT OF MASS ADVECTION IN ITERATION 
!          OF ICE VELOCITIES
!
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_COMMOAU1
      USE MO_COMMOAU2
      USE MO_UNITS
!
      IMPLICIT NONE
      REAL SWSUM, SWSUMI
!
      REAL EPS11(IE,JE),EPS22(IE,JE),EPS12(IE,JE)
      REAL SIOUO(IE,JE),SIOVE(IE,JE),UH(IE,JE),VH(IE,JE)
      REAL EFFICO(IE,JE),EFFICE(IE,JE),UCOR(IE,JE),VCOR(IE,JE),         &
     &       SPEEDU(IE,JE),SPEEDV(IE,JE),ZZL(IE,JE),ZZR(IE,JE),         &
     &       ZZO(IE,JE),ZZU(IE,JE),CWEFFO(IE,JE),CWEFFE(IE,JE)
!
      REAL WH(IE,JE),ZDIFF(IE,JE),ZSIC(IE,JE)
      REAL SIOTHO(IE,JE),SIOOMO(IE,JE),SIOSNO(IE,JE)

      REAL SWRAB(KE),HEATABS(IE,JE),HEATABB(IE,JE)
     
      INTEGER I,ITER,J,K
      INTEGER ICHECK1, ICHECK2, JCHECK1, JCHECK2
      REAL ALPALT,ALPNEU,DXDX,DXDY,DYDY,E11,E12,E22,HIBCC
      REAL OPENDEP,PST,REVAL,RHOICSN,RHOWAIC,RS,RST,RV
      REAL SICM,SICOM,SIMAE,SIMAO,SPEEDMI
      REAL UEIN,UEOU,UEPSI,UWIN,UWOU,UX,UY,VNIN,VNOU,VSIN,VSOU,VX,VY
      REAL ZMAE,ZMAO,ZMIE,ZMIO,ZSCHALT,ZSCHALTD
      REAL AUFTRIEB, BEMMA, COMMI, SCHNEE, UEISMAX, VEISMAX

! RJ: The following is not necesssary, but it helps debugging the parallel version
      HEATABS(:,:)=0
      HEATABB(:,:)=0
      WH(:,:)=0
      SIOUO(:,:)=0
      SIOVE(:,:)=0
!
      SPEEDMI=0.01

      RHOWAIC=RHOWAT/RHOICE
      RHOICSN=RHOICE/RHOSNO
!     TAF ATMOSPHERIC TEMPERATURE
!     HTH00: ARBITRARY CONSTANT FOR DISCRIMINATION BETWEEN THIN AND THICK ICE
!          SET TO 1.5 M   NOT IN USE
!     SICTH  ICE THICKNESS
!     HICCE,HICCP THE COEFFICIENTS E AND P FROM (7), (8)
!     HICCP IS DIVIDED BY MEAN DENSITY OF 1000 KG / M**3
!     SICOM  ICE COMPACTNESS
!     SICDI  ICE FLOW DIVERGENCE
!     SICSH  ICE FLOW SHEAR
!     HIBZET,HIBDEL,HIBET ARE THE FIELDS OF (7) - (9)
!     ENTMEL:MELTING ENTHALPY OF ICE = 320 MILLIONS WS/M**3
!      STEBOL: STEFAN BOLTZMANN CONSTANT =5.67 10**(-8)W/M**2/K**4
!     SICHEC: HEAT CONDUCTIVITY OF ICE = 2 W/M/K   CUWE NOT IN USE
!
!      PART 4: INCREASE OF EXISTING ICE
! IN THEORY A DISTINCTION SHOULD BE MADE BEWTEEN THE ICE-COVERED PART OF A 
! GRID-CELL FOR WHICH THIS DO-LOOP APPLIES AND THE ICE FREE PART, WHICH
! SHOULD BE TREATED AS DO LOOPS 20 AND 25. SUCH A DISTINCTION ENHANCES THE
! OVERALL ICE-GROWTH WHICH RESULTS IN IRREALISTICALLY HIGH VALUES IN SOME
! GRID POINTS AND SUBSEQUENT NUMERICAL INSTABILITIES. INCLUSION OF DIFFUSION
! ON ICETHICKNESS AND COMPACTNESS ALSO GAVE RISE TO INSTABILITIES. 
!
!      CUTOFF VALUE FOR PROGNOSTIC CALCULATION OF ICE VELOCITIES
!  
       REVAL=0.01

       ZMAO=0.
       ZMAE=0.
       ZMIO=0.
       ZMIE=0.
       SIMAE=0.
       SIMAO=0.
       UEPSI=1.E-8
!
!     CHECK SEA ICE PARAMETERS
!
!$OMP PARALLEL PRIVATE(i,j,ux,uy,vx,vy,e12,e11,e22,hibcc,alpneu,alpalt, &
!$OMP   pst,rst,sicm,sicom,dxdx,dydy,dxdy,rs,rv)
!$OMP DO
       DO 31 J=1,JE
       DO 31 I=1,IE
       SICTHO(I,J)=MAX(0.,SICTHO(I,J)*WETO(I,J,1))
       IF(SICTHO(I,J).LE.0.) SICOMO(I,J)=0.
       SICOMO(I,J)=MAX(0.,SICOMO(I,J)*WETO(I,J,1))
       SICSNO(I,J)=MAX(0.,SICSNO(I,J)*WETO(I,J,1))
       SICOMO(I,J)=MIN(1.,SICOMO(I,J))
       SPEEDU(I,J)=0.
       SPEEDV(I,J)=0.
       ZZL(I,J)=0.
       ZZR(I,J)=0.
       ZZO(I,J)=0.
       ZZU(I,J)=0.
       ZDIFF(I,J)=0.
       ZSIC(I,J)=0.
       EFFICO(I,J)=0.
       EFFICE(I,J)=0.
       CWEFFO(I,J)=0.
       CWEFFE(I,J)=0.
       UCOR(I,J)=0.
       VCOR(I,J)=0.
31     CONTINUE
!
!     STORE OLD VALUES FOR SR GROWTH
!
!$OMP DO
      DO J=1,JE
       DO I=1,IE
       SIOTHO(I,J)=SICTHO(I,J)
       SIOOMO(I,J)=SICOMO(I,J)
       SIOSNO(I,J)=SICSNO(I,J)
       ENDDO
      ENDDO
!UWE
!$OMP DO
        DO J=1,JE1
        DO I=1,IE1
        SICOMU(I,J)=1.-0.5*(SICOMO(I,J)+SICOMO(I+1,J))
        SICOMV(I,J)=1.-0.5*(SICOMO(I,J)+SICOMO(I,J+1))
        SICOMP(I,J)=SICOMO(I,J)
        ENDDO
        ENDDO
!$OMP SINGLE
      CALL bounds_exch(SICOMU,SICOMV,SICOMP)
!$OMP END SINGLE
!
!      DYNAMICS
!
!$OMP DO
       DO 49 J=1,JE
       DO 49 I=1,IE
       EPS11(I,J)=0.
       EPS12(I,J)=0.
       EPS22(I,J)=0.
       UH(I,J)=0.
       VH(I,J)=0.
 49    CONTINUE
!
!$OMP DO
       DO 50 J=2,JE1
       DO 50 I=2,IE1
       UX=(SICUO(I,J)-SICUO(I-1,J))/DLXP(I,J)
       VY=(SICVE(I,J-1)-SICVE(I,J))/DLYP(I,J)
       VX=(SICVE(I+1,J)-SICVE(I,J))/DLXP(I,J)
       UY=(SICUO(I,J)-SICUO(I,J+1))/DLYP(I,J)
       EPS11(I,J)=UX
       EPS22(I,J)=VY
       EPS12(I,J)=0.5*(VX+UY)
       IF(AMSUO(I,J,1).GT.0.5)                                          &
     &  SPEEDU(I,J)=MAX(SQRT((UEPSI+SICUO(I,J)-UKO(I,J,1))**2           &
     &   +(0.25*(SICVE(I,J)+SICVE(I+1,J)+SICVE(I,J-1)+SICVE(I+1,J-1)    &
     &     -VKE(I,J,1)-VKE(I+1,J,1)-VKE(I,J-1,1)-VKE(I+1,J-1,1)))**2)   &
     &   ,SPEEDMI)
       IF(AMSUE(I,J,1).GT.0.5)                                          &
     &  SPEEDV(I,J)=MAX(SQRT((UEPSI+SICVE(I,J)-VKE(I,J,1))**2           &
     &   +(0.25*(SICUO(I,J)+SICUO(I-1,J)+SICUO(I,J+1)+SICUO(I-1,J+1)    &
     &     -UKO(I,J,1)-UKO(I-1,J,1)-UKO(I,J+1,1)-UKO(I-1,J+1,1)))**2)   &
     &   ,SPEEDMI)
50     CONTINUE
!$OMP SINGLE
       CALL bounds_exch(EPS11,EPS22,EPS12,SPEEDU,SPEEDV)
!$OMP END SINGLE
!
!$OMP DO
       DO 51 J=2,JE1
       DO 51 I=2,IE1
       E12=0.25*(EPS12(I-1,J-1)+EPS12(I,J-1)+EPS12(I-1,J)+EPS12(I,J))
!      ARGU=((EPS11(I,J)-EPS22(I,J))**2+4.*E12**2)/HICCE**2             &
!    & +(EPS11(I,J)+EPS22(I,J))**2
       HIBDELO(I,J)=SQRT(((EPS11(I,J)-EPS22(I,J))**2+4.*E12**2) /       &
     &               HICCE**2+(EPS11(I,J)+EPS22(I,J))**2)
51     CONTINUE
!
!$OMP DO
       DO 52 J=2,JE1
       DO 52 I=2,IE1
       E11=0.25*(EPS11(I,J)+EPS11(I+1,J)+EPS11(I,J+1)+EPS11(I+1,J+1))
       E22=0.25*(EPS22(I,J)+EPS22(I+1,J)+EPS22(I,J+1)+EPS22(I+1,J+1))
       HIBDELE(I,J)=SQRT(  ((E11-E22)**2+4.*EPS12(I,J)**2)/HICCE**2     &
     & +(E11+E22)**2)
52     CONTINUE
!$OMP SINGLE
       CALL bounds_exch(HIBDELO,HIBDELE)
!$OMP END SINGLE
!C       DO I=1,IE
!C        WRITE(6,*)I,' HIBDEL ', HIBDELO(I,JE/2),HIBDELE(I,JE/2)
!C     X     ,SICTHO(I,JE/2)
!C       ENDDO
!C       CALL ABSTURZ
!
!UWE
!
!      INTRODUCE HIBCC, CORRESPONDS TO C IN HIBLER, JPO 9,825
!
       HIBCC=20.
!
       ALPNEU=1.
!
       ALPALT=1.-ALPNEU
!$OMP DO
       DO 53 J=2,JE1
       DO 53 I=2,IE1
!
!UWE    HIBLER (17), 
       PST=HICCP*SICTHO(I,J)*EXP(-HIBCC*(1.-SICOMO(I,J)))*WETO(I,J,1)
       SICM=(SICTHO(I,J)+SICTHO(I+1,J)+SICTHO(I,J+1)+SICTHO(I+1,J+1))   &
     &/(WETO(I,J,1)+WETO(I+1,J,1)+WETO(I,J+1,1)+WETO(I+1,J+1,1)+1.E-20)
       SICOM=(SICOMO(I,J)+SICOMO(I+1,J)+SICOMO(I,J+1)+SICOMO(I+1,J+1))  &
     &/(WETO(I,J,1)+WETO(I+1,J,1)+WETO(I,J+1,1)+WETO(I+1,J+1,1)+1.E-20)
       SICOM=MIN(SICOM,1.)
       SICOM=MAX(SICOM,0.)
       SICM=MAX(SICM,0.)
       RST=HICCP*SICM*EXP(-HIBCC*(1.-SICOM))
       UH(I,J)=SICUO(I,J)
       VH(I,J)=SICVE(I,J)
       SIOUO(I,J)=SICUO(I,J)
       SIOVE(I,J)=SICVE(I,J)
       IF(HIBDELO(I,J).GT.ALMZER)THEN
         HIBZETO(I,J)=ALPNEU*PST*0.5/MAX(HIBDELO(I,J),ALMZER)           &
     &      +ALPALT*HIBZETO(I,J)
         HIBETO(I,J)=ALPNEU*HIBZETO(I,J)/HICCE**2+ALPALT*HIBETO(I,J)
       ELSE
         HIBZETO(I,J)=0.
         HIBETO(I,J)=0.
       ENDIF
       IF(HIBDELE(I,J).GT.ALMZER) THEN
         HIBZETE(I,J)=ALPNEU*RST*0.5/MAX(HIBDELE(I,J),ALMZER)           &
     &      +ALPALT*HIBZETE(I,J)
         HIBETE(I,J)=ALPNEU*HIBZETE(I,J)/HICCE**2+ALPALT*HIBETE(I,J)
       ELSE
         HIBZETE(I,J)=0.
         HIBETE(I,J)=0.
       ENDIF
!C       SPUR=EPS11(I,J)+EPS22(I,J)
!
       EFFICO(I,J)=0.5*(SICTHO(I,J)+SICTHO(I+1,J)                       &
     &     +RHOICSN*(SICSNO(I,J)+SICSNO(I+1,J)))*AMSUO(I,J,1)
       EFFICE(I,J)=0.5*(SICTHO(I,J)+SICTHO(I,J+1)                       &
     &     +RHOICSN*(SICSNO(I,J)+SICSNO(I,J+1)))*AMSUE(I,J,1)
53     CONTINUE
!$OMP SINGLE
       CALL bounds_exch(UH,VH,SIOUO,SIOVE,EFFICO,EFFICE)
       CALL bounds_exch(HIBZETO,HIBETO,HIBZETE,HIBETE)
!$OMP END SINGLE
!
!      ENHANCE FRICTION COEFFICIENT FOR VERY THIN ICE
!      ==> SMOOTHE TRANSITION TOWARDS WATER VELOCITIES FOR EFFIC ==> 0.
!
!$OMP DO
       DO J=1,JE
        DO I=1,IE
!C         CWEFFO(I,J)=CW*MAX(1.,1./(25.*(EFFICO(I,J)+1.E-6)**2))
!C         CWEFFE(I,J)=CW*MAX(1.,1./(25.*(EFFICE(I,J)+1.E-6)**2))
         CWEFFO(I,J)=CW*MAX(1.,1./(10.*EFFICO(I,J)+1.E-6))**1
         CWEFFE(I,J)=CW*MAX(1.,1./(10.*EFFICE(I,J)+1.E-6))**1
         SICUDO(I,J)=(0.8*SICUDO(I,J)+0.2*SICUO(I,J))*AMSUO(I,J,1)
         SICVDE(I,J)=(0.8*SICVDE(I,J)+0.2*SICVE(I,J))*AMSUE(I,J,1)
         ZSIC(I,J)=SICTHO(I,J)+RHOSNIC*SICSNO(I,J)
        ENDDO
       ENDDO
!
!      OLD ICE VELOCITIES, U ON V-POINT AND V ON U-POINT
!
!      ZSCHALTD: UPWIND TYPE DIFFUSION
!
       ZSCHALTD=1.
!$OMP DO
       DO 54 J=2,JE1
       DO 54 I=2,IE1
       UCOR(I,J)=0.25*(UH(I,J)+UH(I-1,J)+UH(I,J+1)+UH(I-1,J+1))
       VCOR(I,J)=0.25*(VH(I,J)+VH(I+1,J)+VH(I,J-1)+VH(I+1,J-1))
!
!     SEA LEVEL CHANGE DUE TO DIFFUSION OF SEA ICE
!
       ZDIFF(I,J)=ZSCHALTD*0.5*DT*RHOICWA*                              &
     &    (ABS(SICUDO(I-1,J))*DLYU(I-1,J)*(ZSIC(I-1,J)-ZSIC(I,J))       &
     &    +ABS(SICUDO(I,J))*DLYU(I,J)*(ZSIC(I+1,J)-ZSIC(I,J))           &
     &    +ABS(SICVDE(I,J))*DLXV(I,J)*(ZSIC(I,J+1)-ZSIC(I,J))           &
     &    +ABS(SICVDE(I,J-1))*DLXV(I,J-1)*(ZSIC(I,J-1)-ZSIC(I,J)))      &
     &        /(DLXP(I,J)*DLYP(I,J))
!
!     THICKNESS PROTECTION, AVOID NEGATIVE LAYER THICKNESS
!       INDUCE DIVERGENT ICE TRANSPORTS
!
!C     X  +  MAX(0.,ZSIC(I,J)-0.6667*(DZW(1)+ZO(I,J)))
!C     X     /MAX(0.01,DZW(1)+ZO(I,J)-ZSIC(I,J))
!
54     CONTINUE
!$OMP SINGLE
       CALL bounds_exch(UCOR,VCOR,ZDIFF)
       ZSCHALT=1.
!$OMP END SINGLE
!
!      MAIN ITERATION
!
       DO 100 ITER=1,40
!
!$OMP DO
       DO 601  J=1,JE
       DO 601 I=1,IE
       UH(I,J)=0.
       VH(I,J)=0.
601    CONTINUE
!
!$OMP DO
       DO J=2,JE1
        DO I=2,IE1
       ZZL(I,J)=DT*RHOICWA*(                                            &
     &     DLYU(I-1,J)*SICUO(I-1,J)*EFFICO(I-1,J)                       &
     &    +DLXV(I,J)*SICVE(I,J)*EFFICE(I,J)                             &
     &    -DLXV(I,J-1)*SICVE(I,J-1)*EFFICE(I,J-1))                      &
     &   /(DLXP(I,J)*DLYP(I,J))
       ZZR(I,J)=DT*RHOICWA*(                                            &
     &    -DLYU(I+1,J)*SICUO(I+1,J)*EFFICO(I+1,J)                       &
     &    +DLXV(I+1,J)*SICVE(I+1,J)*EFFICE(I+1,J)                       &
     &    -DLXV(I+1,J-1)*SICVE(I+1,J-1)*EFFICE(I+1,J-1))                &
     &   /(DLXP(I+1,J)*DLYP(I+1,J))
       SPEEDU(I,J)=0.5*(MAX(SQRT((UEPSI+SICUO(I,J)-UKO(I,J,1))**2       &
     &   +(0.25*(SICVE(I,J)+SICVE(I+1,J)+SICVE(I,J-1)+SICVE(I+1,J-1)    &
     &     -VKE(I,J,1)-VKE(I+1,J,1)-VKE(I,J-1,1)-VKE(I+1,J-1,1)))**2)   &
     &   ,SPEEDMI)+SPEEDU(I,J))
       SPEEDV(I,J)=0.5*(MAX(SQRT((UEPSI+SICVE(I,J)-VKE(I,J,1))**2       &
     &   +(0.25*(SICUO(I,J)+SICUO(I-1,J)+SICUO(I,J+1)+SICUO(I-1,J+1)    &
     &     -UKO(I,J,1)-UKO(I-1,J,1)-UKO(I,J+1,1)-UKO(I-1,J+1,1)))**2)   &
     &   ,SPEEDMI)+SPEEDV(I,J))
        ENDDO
       ENDDO
!$OMP SINGLE
       CALL bounds_exch(ZZL,ZZR,SPEEDU,SPEEDV)
!$OMP END SINGLE
!
!JJ    DO 61 J0=0,15
!JJ    DO 61 J=2+J0,JE1,16
!JJ RESTORE OLD LOOPS      
!$OMP DO
       DO 61 J=2,JE1
       DO 61 I=2,IE1
!
       IF (EFFICO(I,J).GT.REVAL)THEN
       DXDX=DLXP(I,J)**2
       DYDY=DLYP(I,J)**2
       DXDY=DLXP(I,J)*DLYP(I,J)
       RS=EFFICO(I,J)*(SIOUO(I,J)+0.5*DT*FTWOU(I,J)*(VCOR(I,J)          &
     &+0.25*(SICVE(I,J)+SICVE(I+1,J)+SICVE(I,J-1)+SICVE(I+1,J-1))))     &
     & +CWEFFO(I,J)*DT*UKO(I,J,1)*SPEEDU(I,J)                           &
     &+ DT*HICCP*(SICTHO(I,J)*EXP(-HIBCC*(1.-SICOMO(I,J)))              &
     &   -SICTHO(I+1,J)*EXP(-HIBCC*(1.-SICOMO(I+1,J))))/DLXU(I,J)       &
     &+ DT*EFFICO(I,J)*G*(ZO(I,J)-ZO(I+1,J)+ZDIFF(I,J)-ZDIFF(I+1,J)     &
     &    +ZSCHALT*(ZZL(I,J)-ZZR(I,J)))/DLXU(I,J)                       &
     & +DT*RHOWAIC*TXO(I,J)*(SICOMO(I,J)+SICOMO(I+1,J))*0.5

!SVX 31.08.99
       RV=                                                              &
     &+(HIBZETO(I,J)+HIBETO(I,J))    *SICUO(I-1,J)/DXDX                 &
     &+(HIBZETO(I+1,J)+HIBETO(I+1,J))*SICUO(I+1,J)/DXDX                 &
     &+(HIBZETO(I+1,J)-HIBETO(I+1,J))*(SICVE(I+1,J-1)-SICVE(I+1,J))     &
     & /DXDY                                                            &
     &-(HIBZETO(I,J)-HIBETO(I,J))*(SICVE(I,J-1)-SICVE(I,J))/DXDY        &
     &+(HIBETE(I,J-1)*SICUO(I,J-1)+HIBETE(I,J)*SICUO(I,J+1))/DYDY       &
     &+HIBETE(I,J-1)*(SICVE(I+1,J-1)-SICVE(I,J-1))/DXDY                 &
     &-HIBETE(I,J)*(SICVE(I+1,J)-SICVE(I,J))/DXDY
!
       UH(I,J)=(RS+RV*DT)/(EFFICO(I,J)                                  &
     &      +CWEFFO(I,J)*DT*SPEEDU(I,J)                                 &
     &   +DT*EFFICO(I,J)*G*DT*RHOICWA*EFFICO(I,J)*ZSCHALT*              &
     &  (1./(DLXP(I,J)*DLYP(I,J))+1./(DLXP(I+1,J)*DLYP(I+1,J)))         &
     &*DLYU(I,J)/DLXU(I,J)                                              &
     & +DT*(HIBZETO(I,J)+HIBETO(I,J)+HIBZETO(I+1,J)+HIBETO(I+1,J))/DXDX &
     & +DT*(HIBETE(I,J-1)+HIBETE(I,J))/DYDY)
      ENDIF
61    CONTINUE
!$OMP SINGLE
      CALL bounds_exch(UH)
!$OMP END SINGLE
!OtB JE-2 changed to JE1
!     DO J=2,JE-2
!$OMP DO
      DO J=2,JE1
       DO I=2,IE1
        ZZO(I,J)=DT*RHOICWA*(                                           &
     &      DLYU(I-1,J)*SICUO(I-1,J)*EFFICO(I-1,J)                      &
     &     -DLYU(I,J)*SICUO(I,J)*EFFICO(I,J)                            &
     &     -DLXV(I,J-1)*SICVE(I,J-1)*EFFICE(I,J-1))                     &
     &     /(DLXP(I,J)*DLYP(I,J))
        ZZU(I,J)=DT*RHOICWA*(                                           &
     &      DLYU(I-1,J+1)*SICUO(I-1,J+1)*EFFICO(I-1,J+1)                &
     &     -DLYU(I,J+1)*SICUO(I,J+1)*EFFICO(I,J+1)                      &
     &     +DLXV(I,J+1)*SICVE(I,J+1)*EFFICE(I,J+1))                     &
     &     /(DLXP(I,J+1)*DLYP(I,J+1))
       ENDDO
      ENDDO
!$OMP SINGLE
      CALL bounds_exch(ZZO,ZZU)
!$OMP END SINGLE
!JJ   DO 62 J0=0,15
!JJ   DO 62 J=2+J0,JE1,16
!$OMP DO
      DO 62 J=2,JE1      
      DO 62 I=2,IE1
      IF(EFFICE(I,J).GT.REVAL) THEN
       DXDX=DLXP(I,J)**2
       DYDY=DLYP(I,J)**2
       DXDY=DLXP(I,J)*DLYP(I,J)
      RS=EFFICE(I,J)*(SIOVE(I,J)-0.5*DT*FTWOV(I,J)*(UCOR(I,J)           &
     &+0.25*(SICUO(I,J)+SICUO(I-1,J)+SICUO(I,J+1)+SICUO(I-1,J+1))))     &
     & +CWEFFE(I,J)*DT*VKE(I,J,1)*SPEEDV(I,J)                           &
     & +DT*HICCP*(SICTHO(I,J+1)*EXP(-HIBCC*(1.-SICOMO(I,J+1)))          &
     &  -SICTHO(I,J)*EXP(-HIBCC*(1.-SICOMO(I,J))))/DLYV(I,J)            &
     & +EFFICE(I,J)*G*DT*(ZO(I,J+1)-ZO(I,J)+ZDIFF(I,J+1)-ZDIFF(I,J)     &
     &     +ZSCHALT*(ZZU(I,J)-ZZO(I,J)))/DLYV(I,J)                      &
     & +RHOWAIC*TYE(I,J)*DT*(SICOMO(I,J)+SICOMO(I,J+1))*0.5 

       RV=(HIBZETO(I,J)+HIBETO(I,J))*SICVE(I,J-1)/DYDY                  &
     & +   (HIBZETO(I,J+1)+HIBETO(I,J+1))*SICVE(I,J+1)/DYDY             &
     & +(HIBETE(I,J)*SICVE(I+1,J)+HIBETE(I-1,J)*SICVE(I-1,J))/DXDX      &
     & +(HIBZETO(I,J)-HIBETO(I,J))*(SICUO(I,J)-SICUO(I-1,J))/DXDY       &
     & -(HIBZETO(I,J+1)-HIBETO(I,J+1))*(SICUO(I,J+1)-SICUO(I-1,J+1))    &
     & /DXDY                                                            &
     & +HIBETE(I,J)*(SICUO(I,J)-SICUO(I,J+1))/DXDY                      &
     & -HIBETE(I-1,J)*(SICUO(I-1,J)-SICUO(I-1,J+1))/DXDY    
!
      VH(I,J)=(RS+RV*DT)/(EFFICE(I,J)                                   &
     &        +CWEFFE(I,J)*DT*SPEEDV(I,J)                               &
     &   +DT*EFFICE(I,J)*G*DT*ZSCHALT*EFFICE(I,J)*RHOICWA*              &
     & (1./(DLXP(I,J)*DLYP(I,J))+1./(DLXP(I,J+1)*DLYP(I,J+1)))          &
     &*DLXV(I,J)/DLYV(I,J)                                              &
     & +DT*(HIBZETO(I,J)+HIBETO(I,J)+HIBZETO(I,J+1)+HIBETO(I,J+1))/DYDY &
     & +DT*(HIBETE(I,J)+HIBETE(I-1,J))/DXDX)
!
      ENDIF
62    CONTINUE
!$OMP SINGLE
      CALL bounds_exch(VH)
!$OMP END SINGLE
!$OMP DO
      DO 63 J=2,JE1
      DO 63 I=2,IE1
!
      IF(EFFICO(I,J).GT.REVAL)THEN
      SICUO(I,J)=0.5*(UH(I,J)+SICUO(I,J))*AMSUO(I,J,1)
      ELSE
      SICUO(I,J)=UKO(I,J,1)
      ENDIF
      IF(EFFICE(I,J).GT.REVAL)THEN
      SICVE(I,J)=0.5*(VH(I,J)+SICVE(I,J))*AMSUE(I,J,1)
      ELSE
      SICVE(I,J)=VKE(I,J,1)
      ENDIF
63    CONTINUE
!$OMP SINGLE
      CALL bounds_exch(SICUO,SICVE)
!      WRITE(IO_STDOUT,*) ' ITERATION   U ',ITER
!      WRITE(IO_STDOUT,681)((AMSUO(I,J,1)*SICUO(I,J),J=1,25),I=1,IE)
!      WRITE(IO_STDOUT,*) ' ITERATION   V ',ITER
!      WRITE(IO_STDOUT,681)((AMSUE(I,J,1)*SICVE(I,J),J=1,25),I=1,IE)
!$OMP END SINGLE

100   CONTINUE

!$OMP END PARALLEL
!      WRITE(IO_STDOUT,*) ' FELD '
!      WRITE(IO_STDOUT,699)((VH(I,J),I=46,65),J=1,JE)
699   FORMAT(20F6.2)
!
!     TEST IF ICE ADVECTION VIOLATES CFL CRITERION
!
      DTI=1./DT
      DO J=1,JE
       DO I=1,IE
        IF(ABS(SICUO(I,J)).GT.DLXU(I,J)*DTI)THEN
          WRITE(6,*)'ALARM!  EISADVEKT U: ',I+p_ioff,J+p_joff,SICUO(I,J)&
     &       ,SICUO(I,J)*DT/DLXU(I,J),EFFICO(I,J),EFFICO(I+1,J)         &
     &     ,SPEEDU(I,J),UKO(I,J,1)
        ENDIF
        IF(ABS(SICVE(I,J)).GT.DLYV(I,J)*DTI)THEN
          WRITE(6,*)'ALARM!  EISADVEKT V: ',I+p_ioff,J+p_joff,SICVE(I,J)&
     &       ,SICVE(I,J)*DT/DLYV(I,J),EFFICO(I,J+1),EFFICO(I,J)         &
     &     ,SPEEDV(I,J),VKE(I,J,1)
        ENDIF
       ENDDO
      ENDDO
!
!
!     ADVECTION OF SEA ICE AND SNOW
!
!CCCC      IF (ZSCHALTD.EQ.0)THEN
!
!     UPWIND
!
!$OMP PARALLEL PRIVATE(UWIN,UWOU,UEIN,UEOU,VSIN,VSOU,VNIN,VNOU)
!$OMP DO
      DO 65 J=2,JE1
      DO 65 I=2,IE1
      UWIN=0.5*(SICUO(I-1,J)+ABS(SICUO(I-1,J)))
      UWOU=0.5*(ABS(SICUO(I-1,J))-SICUO(I-1,J))
      UEIN=0.5*(ABS(SICUO(I,J))-SICUO(I,J))
      UEOU=0.5*(ABS(SICUO(I,J))+SICUO(I,J))
      VSIN=0.5*(ABS(SICVE(I,J))+SICVE(I,J))
      VSOU=0.5*(ABS(SICVE(I,J))-SICVE(I,J))
      VNIN=0.5*(ABS(SICVE(I,J-1))-SICVE(I,J-1))
      VNOU=0.5*(ABS(SICVE(I,J-1))+SICVE(I,J-1))
      UH(I,J)=SICTHO(I,J)*(1.-(UWOU+UEOU)*DT/DLXP(I,J)                  &
     & -DT*(VSOU*DLXV(I,J)+VNOU*DLXV(I,J-1))/(DLYP(I,J)*DLXP(I,J)))     &
     & + DT*(UWIN*SICTHO(I-1,J)+UEIN*SICTHO(I+1,J))/DLXP(I,J)           &
     & +DT*( VSIN*DLXV(I,J)*SICTHO(I,J+1)                               &
     &  +    VNIN*DLXV(I,J-1)*SICTHO(I,J-1))/(DLYP(I,J)*DLXP(I,J))
      VH(I,J)=SICOMO(I,J)*(1.-(UWOU+UEOU)*DT/DLXP(I,J)                  &
     & -DT*(VSOU*DLXV(I,J)+VNOU*DLXV(I,J-1))/(DLYP(I,J)*DLXP(I,J)))     &
     & + DT*(UWIN*SICOMO(I-1,J)+UEIN*SICOMO(I+1,J))/DLXP(I,J)           &
     & +DT*( VSIN*DLXV(I,J)*SICOMO(I,J+1)                               &
     &  +    VNIN*DLXV(I,J-1)*SICOMO(I,J-1))/(DLYP(I,J)*DLXP(I,J))
!AS***ADD CONTINUITY OF SNOW:
      WH(I,J)=SICSNO(I,J)*(1.-(UWOU+UEOU)*DT/DLXP(I,J)                  &
     & -DT*(VSOU*DLXV(I,J)+VNOU*DLXV(I,J-1))/(DLYP(I,J)*DLXP(I,J)))     &
     & + DT*(UWIN*SICSNO(I-1,J)+UEIN*SICSNO(I+1,J))/DLXP(I,J)           &
     & +DT*( VSIN*DLXV(I,J)*SICSNO(I,J+1)                               &
     &  +    VNIN*DLXV(I,J-1)*SICSNO(I,J-1))/(DLYP(I,J)*DLXP(I,J))
65    CONTINUE
! No bounds exchange necessary on UH, VH, WH since they are used only
! in the inner domain below
!
!$OMP DO
      DO 66 J=2,JE1
      DO 66 I=2,IE1
!UWE  INCLUDE EFFECT ON SEALEVEL ZO
      ZO(I,J)=ZO(I,J)+(UH(I,J)-SICTHO(I,J))*RHOICWA*ZSCHALT             &
     &               +(WH(I,J)-SICSNO(I,J))*RHOSNWA*ZSCHALT
      SICTHO(I,J)=UH(I,J)
      SICOMO(I,J)=VH(I,J)
!AS***ADD SNOW:
      SICSNO(I,J)=WH(I,J)
66    CONTINUE
!$OMP END PARALLEL
      CALL bounds_exch(ZO,SICTHO,SICOMO,SICSNO)
!AS***CALL SEA-ICE THERMODYNAMICS ROUTINES ACC. TO SIBL
!AS***   (DKRZ TECHN. REP. NO.3)
!
!
!     LENGTH SCALE FOR PENETRATION OF SW-RADIATION [M]
!       	  
#ifdef OPEND55
      OPENDEP=5.5
#else
      OPENDEP=11.
#endif
!
      SWSUM=EXP(-TIESTW(2)/OPENDEP)
      SWSUMI=1./SWSUM
!
      CALL GROWTH(ALAT,SIOTHO,SIOOMO,SIOSNO,SICTHO,SICOMO,SICSNO,       &
     &            TICEO,RPRECO,TAIRO,TDO,ACLO,PAO,FU10,FSWR,            &
     &            SAO,THO,DZW(1),ZO,DT,FWO,SHO,WETO,                    &
     &            QSWO,QLWO,QSEO,QLAO,PRECO,PRECH,HEATABS,SWSUM)
!
      DO K=1,KE
       SWRAB(K)=SWSUMI*(EXP(-TIESTW(K)/OPENDEP)                         &
     &    -EXP(-TIESTW(K+1)/OPENDEP))
      ENDDO
!
!     OUTER LOOP OVER J
!
      DO J=1,JE
!
       DO I=1,IE
        HEATABS(I,J)=WETO(I,J,1)*HEATABS(I,J)*DT/CC
        HEATABB(I,J)=0.
       ENDDO
!

       DO K=KE,2,-1
        DO I=1,IE
         HEATABB(I,J)=HEATABB(I,J)+SWRAB(K)*HEATABS(I,J)
         THO(I,J,K)=THO(I,J,K)+HEATABB(I,J)*DPIO(I,J,K)
!
         HEATABB(I,J)=HEATABB(I,J)*(1.-WETO(I,J,K))
        ENDDO
       ENDDO

!
!    CLOSE LOOP OVER J
!
      ENDDO   
!
       DO J=2,JE1
        DO I=2,IE1
         TAUWATU(I,J)=CW*RHOICWA*SPEEDU(I,J)*(SICUO(I,J)-UKO(I,J,1))
         TAUWATV(I,J)=CW*RHOICWA*SPEEDV(I,J)*(SICVE(I,J)-VKE(I,J,1))
        ENDDO
       ENDDO

       DO I=1,IE
        IF(have_g_js) THEN
          TAUWATU(I,1)=0.
          TAUWATV(I,1)=0.
        ENDIF
        IF(have_g_je) THEN
          TAUWATU(I,JE)=0.
          TAUWATV(I,JE)=0.
        ENDIF
       ENDDO

       CALL bounds_exch(TAUWATU,TAUWATV)
       RETURN
       END
