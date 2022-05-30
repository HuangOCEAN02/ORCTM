       SUBROUTINE OCADPO(TRF)
#ifdef ADPO
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_COMMOAU1
#ifdef SLOPECON_ADPO
      USE MO_COMMOBBL
#endif /* SLOPECON_ADPO */
      USE MO_UNITS
!
!     ROUTINE OCADPO
!
!     COMPUTES ADVECTION OF ONE TRACERFIELD
!
!     BY ERNST MAIER-REIMER 1/1999
!     MODIFIED 1/2000 UWE MIKOLAJEWICZ
!      MAKE MASS CONSERVING WITH FREE SURFACE LAYER
!     MODIFIED 3/2000 UWE MIKOLAJEWICZ
!      INCLUDE BOTTOM BOUNDARY LAYER TRANSPORTS IN
!      ADVECTION SCHEME (OPTION SLOPECON_ADPO)
!      NEEDS TO BE RUN TOGETHER WITH SR SLOPETRANS
!
!     INPUT/OUTPUT  
!     TRF(IE,JE,KE)     TRACER FIELD 
!
!      USES VELOCITIES (UKO,VKE,WO)
!
      DIMENSION TRP(IE,JE,KEP),TRM(IE,JE,0:KE)
      DIMENSION TRF(IE,JE,KE)
!OtB      DIMENSION TRF(IE,JE,KE),TRVOL(IE,JE,KE)
!OtB      DIMENSION TR1(IE,JE,KE)
      DIMENSION WTP(IE,JE,KEP),WTM(IE,JE,0:KE)
#ifdef SLOPECON_ADPO
      DIMENSION WOBACK(IE,JE,KEP),BBU(IE,JE),BBV(IE,JE)
#endif
!
!UWE  ONLY PURPOSE OF EQUIVALENCE: SAVE MEMORY, NO DATA TRANSFER
!
!OtB      EQUIVALENCE (TR1(1,1,1),T1O(1,1,1))
!OtB      EQUIVALENCE (TRVOL(1,1,1),S1O(1,1,1))
!
!$OMP PARALLEL PRIVATE(i,j,k,klo,kup,suminf,abl,zwabl,up,um,scha,sch,schi,rn,bbu,bbv)
#ifdef SLOPECON_ADPO
!$OMP DO
      DO K=1,KEP
       DO J=1,JE
        DO I=1,IE
         WOBACK(I,J,K)=WO(I,J,K)
        ENDDO
       ENDDO
      ENDDO
!
!     COMPUTE NEW VERTICAL VELOCITIES INCLUDING BBL TRANSPORTS
!
!$OMP DO
      DO J=2,JE1
       DO K=KE,1,-1
        DO I=2,IE1
         SUMINF=0.
!
         IF(KDWUBBL(I-1,J).EQ.K)SUMINF=SUMINF+UBBL(I-1,J)
         IF(KUPUBBL(I-1,J).EQ.K)SUMINF=SUMINF-UBBL(I-1,J)
         IF(KDWUBBL(I,J).EQ.K)SUMINF=SUMINF-UBBL(I,J)
         IF(KUPUBBL(I,J).EQ.K)SUMINF=SUMINF+UBBL(I,J)
!
         IF(KDWVBBL(I,J).EQ.K)SUMINF=SUMINF+VBBL(I,J)
         IF(KUPVBBL(I,J).EQ.K)SUMINF=SUMINF-VBBL(I,J)
         IF(KDWVBBL(I,J-1).EQ.K)SUMINF=SUMINF-VBBL(I,J-1)
         IF(KUPVBBL(I,J-1).EQ.K)SUMINF=SUMINF+VBBL(I,J-1)
!
         WO(I,J,K)=WO(I,J,K+1)+SUMINF/(DLXP(I,J)*DLYP(I,J))
!
        ENDDO
       ENDDO
      ENDDO
!
!$OMP DO
      DO K=1,KE
       DO J=2,JE1
        DO I=2,IE1
         WO(I,J,K)=WOBACK(I,J,K)+WO(I,J,K)
        ENDDO
       ENDDO
      ENDDO
!$OMP SINGLE
      CALL bounds_exch(WO)
!$OMP END SINGLE
!
!
#endif
!
!$OMP DO
      DO 1 K=1,KE
      DO 1 J=1,JE
      DO 1 I=1,IE
      S1O(I,J,K)=WETO(I,J,K)*DLXP(I,J)*DLYP(I,J)*DDPO(I,J,K)
!CC     X  +(1.-WETO(I,J,K))*1.
      T1O(I,J,K)=TRF(I,J,K)
      TRP(I,J,K)=0.
      TRM(I,J,K)=0.
      WTP(I,J,K)=0.
      WTM(I,J,K)=0.
1     CONTINUE
!
!$OMP DO
       DO 102 J=1,JE
       DO 102 I=1,IE
       TRM(I,J,0)=0.
       WTM(I,J,0)=0.
       TRP(I,J,KEP)=0.
       WTP(I,J,KEP)=0.
       S1O(I,J,1)=S1O(I,J,1)+DLXP(I,J)*DLYP(I,J)*WETO(I,J,1)*           &
     & (ZO(I,J)-WO(I,J,1)*DT-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)
102   CONTINUE
!
!      VERTICAL TRANSPORTS
!
!$OMP DO
       DO 23 K=1,KE
       KLO=MIN(K+1,KE)
       KUP=MAX(K-1,1)
       DO 23 J=2,JE1
       DO 23 I=2,IE1
       ABL=ABS(TRF(I,J,KUP)-TRF(I,J,KLO))
       ZWABL=ABS(TRF(I,J,KUP)+TRF(I,J,KLO)-2.*TRF(I,J,K))
       UP=0.5*DT*DLYP(I,J)*DLXP(I,J)*(WO(I,J,K)+ABS(WO(I,J,K)))
       UM=0.5*DT*DLYP(I,J)*DLXP(I,J)*(ABS(WO(I,J,K+1))-WO(I,J,K+1))
       SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))
!CUWE       SCH=MIN(1.,SCHA*S1O(I,J,K)/(UP+UM+1.E-20))
#ifdef SMOADV
       SCH=MIN(WETO(I,J,KLO),SCHA)
#else
       SCH=MIN(WETO(I,J,KLO),SCHA*S1O(I,J,K)/(UP+UM+1.E-20))
#endif
       IF(K.EQ.1.OR.K.EQ.KE)SCH=0.
       SCHI=1.-SCH
       TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J,KUP)))
       TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J,KLO)))
       WTP(I,J,K)=UP
       WTM(I,J,K)=UM
23      CONTINUE
!
!$OMP DO
       DO 27 J=1,JE
       DO 27 I=1,IE
       TRP(I,J,KEP)=0.
       WTP(I,J,KEP)=0.
!
       TRM(I,J,0)=0.
       WTM(I,J,0)=0.
!
27     CONTINUE
!
! RJ: no boundary exchange necessary here, since TRP, WTP, TRM, WTM
!     are used only in the inner domain in the following loop
!$OMP DO
       DO 26 K=1,KE
       DO 24 J=2,JE1
       DO 24 I=2,IE1
       IF(WETO(I,J,K).GT.0.5) THEN
       RN=S1O(I,J,K)+WTP(I,J,K+1)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J,K-1)
       T1O(I,J,K)=(TRF(I,J,K)*S1O(I,J,K)+TRP(I,J,K+1)-TRP(I,J,K)        &
     & -TRM(I,J,K)+TRM(I,J,K-1))/RN
       S1O(I,J,K)=RN
        ENDIF
24      CONTINUE
       DO 25 J=2,JE1
       DO 25 I=2,IE1
       TRF(I,J,K)=T1O(I,J,K)
25      CONTINUE 
26     CONTINUE
!
!
!      TRANSPORTS IN X-DIRECTION
!
!$OMP SINGLE
      CALL bounds_exch(TRF)
!$OMP END SINGLE
!$OMP DO
       DO 113 K=1,KE
       DO 113 J=1,JE
       DO 113 I=1,IE
       TRP(I,J,K)=0.
       TRM(I,J,K)=0.
       WTP(I,J,K)=0.
       WTM(I,J,K)=0.
113    CONTINUE
!
!$OMP DO
       DO 3 K=1,KE
#ifdef SLOPECON_ADPO
        DO J=1,JE
         DO I=1,IE
          BBU(I,J)=0.
          IF(KDWUBBL(I,J).EQ.K)BBU(I,J)=UBBL(I,J)
          IF(KUPUBBL(I,J).EQ.K)BBU(I,J)=-UBBL(I,J)
         ENDDO
        ENDDO
#endif
!
       DO 3 J=2,JE1
       DO 3 I=2,IE1
       ABL=ABS(TRF(I+1,J,K)-TRF(I-1,J,K))
       ZWABL=ABS(TRF(I+1,J,K)+TRF(I-1,J,K)-2.*TRF(I,J,K))
#ifdef SLOPECON_ADPO
       UP=0.5*DT*DDUO(I,J,K)*DLYU(I,J)*(UKO(I,J,K)+ABS(UKO(I,J,K)))     &
     &    +0.5*DT*(BBU(I,J)+ABS(BBU(I,J)))
#else
       UP=0.5*DT*DDUO(I,J,K)*DLYU(I,J)*(UKO(I,J,K)+ABS(UKO(I,J,K)))
#endif
#ifdef SLOPECON_ADPO
       UM=0.5*DT*DDUO(I-1,J,K)*(ABS(UKO(I-1,J,K))-UKO(I-1,J,K))         &
     & *DLYU(I-1,J)                                                     &
     &     +0.5*DT*(ABS(BBU(I-1,J))-BBU(I-1,J))
#else
       UM=0.5*DT*DDUO(I-1,J,K)*(ABS(UKO(I-1,J,K))-UKO(I-1,J,K))         &
     & *DLYU(I-1,J)
#endif
       SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))*WETO(I,J,K)                &
     &       *WETO(I+1,J,K)*WETO(I-1,J,K)
#ifdef SMOADH
       SCH=MIN(1.,SCHA)
#else
       SCH=MIN(1.,SCHA*S1O(I,J,K)/(UP+UM+1.E-20))
#endif
       SCHI=1.-SCH
       TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I+1,J,K)))
       TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I-1,J,K)))
       WTP(I,J,K)=UP
       WTM(I,J,K)=UM
3      CONTINUE
!
!$OMP SINGLE
      CALL bounds_exch(TRP)
      CALL bounds_exch(TRM)
      CALL bounds_exch(WTP)
      CALL bounds_exch(WTM)
!$OMP END SINGLE
!
!$OMP DO
       DO 6 K=1,KE
       DO 4 J=2,JE1
       DO 4 I=2,IE1
       IF(WETO(I,J,K).GT.0.5) THEN
!
       RN=S1O(I,J,K)+WTP(I-1,J,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I+1,J,K)
       T1O(I,J,K)=(TRF(I,J,K)*S1O(I,J,K)+TRP(I-1,J,K)-TRP(I,J,K)        &
     & -TRM(I,J,K)+TRM(I+1,J,K))/RN
       S1O(I,J,K)=RN
       ENDIF
4      CONTINUE
       DO 5 J=2,JE1
       DO 5 I=2,IE1
       TRF(I,J,K)=T1O(I,J,K)
5      CONTINUE 
6     CONTINUE
!
!$OMP SINGLE
      CALL bounds_exch(TRF)
!$OMP END SINGLE
!
!      TRANSPORTS IN Y-DIRECTION
!
!$OMP DO
       DO 114 K=1,KE
       DO 114 J=1,JE
       DO 114 I=1,IE
       TRP(I,J,K)=0.
       TRM(I,J,K)=0.
       WTP(I,J,K)=0.
       WTM(I,J,K)=0.
114    CONTINUE
!
!$OMP DO
       DO 13 K=1,KE
#ifdef SLOPECON_ADPO
!      INITIALIZE BBL TRANSPORTS FOR LEVEL
!
       DO J=1,JE
        DO I=1,IE
         BBV(I,J)=0.
         IF(K.EQ.KDWVBBL(I,J))BBV(I,J)=VBBL(I,J)
         IF(K.EQ.KUPVBBL(I,J))BBV(I,J)=-VBBL(I,J)
        ENDDO
       ENDDO
#endif
!
       DO 13 J=2,JE1
       DO 13 I=2,IE1
       ABL=ABS(TRF(I,J-1,K)-TRF(I,J+1,K))
       ZWABL=ABS(TRF(I,J-1,K)+TRF(I,J+1,K)-2.*TRF(I,J,K))
#ifdef SLOPECON_ADPO
       UP=0.5*DT*DDUE(I,J-1,K)*DLXV(I,J-1)                              &
     & *(VKE(I,J-1,K)+ABS(VKE(I,J-1,K)))                                &
     &           +0.5*DT*(BBV(I,J-1)+ABS(BBV(I,J-1)))
#else
       UP=0.5*DT*DDUE(I,J-1,K)*DLXV(I,J-1)                              &
     & *(VKE(I,J-1,K)+ABS(VKE(I,J-1,K)))
#endif
#ifdef SLOPECON_ADPO
       UM=0.5*DT*DDUE(I,J,K)*DLXV(I,J)*(ABS(VKE(I,J,K))-VKE(I,J,K))     &
     &           +0.5*DT*(ABS(BBV(I,J))-BBV(I,J))
#else
       UM=0.5*DT*DDUE(I,J,K)*DLXV(I,J)*(ABS(VKE(I,J,K))-VKE(I,J,K))
#endif
       SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))*WETO(I,J,K)                &
     &     *WETO(I,J-1,K)*WETO(I,J+1,K)
#ifdef SMOADH
       SCH=MIN(1.,SCHA)
#else
       SCH=MIN(1.,SCHA*S1O(I,J,K)/(UP+UM+1.E-20))
#endif
       SCHI=1.-SCH
       TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J-1,K)))
       TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J+1,K)))
       WTP(I,J,K)=UP
       WTM(I,J,K)=UM
13      CONTINUE
!
!$OMP SINGLE
      CALL bounds_exch(TRP)
      CALL bounds_exch(TRM)
      CALL bounds_exch(WTP)
      CALL bounds_exch(WTM)
!$OMP END SINGLE
!
!$OMP DO
       DO 16 K=1,KE
       DO 14 J=2,JE1
       DO 14 I=2,IE1
       IF(WETO(I,J,K).GT.0.5) THEN
!
       RN=S1O(I,J,K)+WTP(I,J+1,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J-1,K)
       T1O(I,J,K)=(TRF(I,J,K)*S1O(I,J,K)+TRP(I,J+1,K)-TRP(I,J,K)        &
     & -TRM(I,J,K)+TRM(I,J-1,K))/RN
       S1O(I,J,K)=RN
       ENDIF
14      CONTINUE
       DO 15 J=2,JE1
       DO 15 I=2,IE1
       TRF(I,J,K)=T1O(I,J,K)
15      CONTINUE 
16     CONTINUE
!
!$OMP SINGLE
      CALL bounds_exch(TRF)
!$OMP END SINGLE
!
!     RESET VERTICAL VELOCITIES TO INITIAL VALUES WITHOUT BBL TRANSPORTS
!
#ifdef SLOPECON_ADPO
!$OMP DO
      DO K=1,KE
       DO J=1,JE
        DO I=1,IE
         WO(I,J,K)=WOBACK(I,J,K)
        ENDDO
       ENDDO
      ENDDO
#endif /*SLOPECON_ADPO*/
!
#endif /*ADPO*/
!$OMP END PARALLEL

      RETURN
      END
