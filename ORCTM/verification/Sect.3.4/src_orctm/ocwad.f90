      SUBROUTINE OCWAD(TRF)
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
!
!
      DIMENSION TRP(IE,JE,KEP),TRM(IE,JE,0:KE)
      DIMENSION TRF(IE,JE,KEP)
      DIMENSION WTP(IE,JE,KEP),WTM(IE,JE,0:KE)
!
!
!$OMP DO
      DO 1 K=1,KE
      DO 1 J=1,JE
      DO 1 I=1,IE
      T1O(I,J,K)=TRF(I,J,K)
      TRP(I,J,K)=0.
      TRM(I,J,K)=0.
      WTP(I,J,K)=0.
      WTM(I,J,K)=0.
1     CONTINUE
      DO 111 K=1,KE
      DO 111 J=1,JE
      DO 111 I=1,IE
      IF (K .EQ. 1) THEN
      S1O(I,J,K)=WETO(I,J,K)*DLXP(I,J)*DLYP(I,J)*(                      &
     &           DDPO(I,J,K)+ZO(I,J))/2.0
      ELSEIF (K .EQ. 2) THEN
      S1O(I,J,K)=WETO(I,J,K)*DLXP(I,J)*DLYP(I,J)*(DDPO(I,J,K-1)         &
     &          +DDPO(I,J,K)+ZO(I,J))/2.0
      ELSE
      S1O(I,J,K)=WETO(I,J,K)*DLXP(I,J)*DLYP(I,J)*(DDPO(I,J,K-1)         &
     &          +DDPO(I,J,K))/2.0
      ENDIF
111   CONTINUE
!
!$OMP DO
       DO 102 J=1,JE
       DO 102 I=1,IE
       TRM(I,J,0)=0.
       WTM(I,J,0)=0.
       TRP(I,J,KEP)=0.
       WTP(I,J,KEP)=0.
102   CONTINUE
!
!      VERTICAL TRANSPORTS
!
!$OMP DO
       DO 23 K=2,KE
       KLO=K+1
       KUP=K-1
       DO 23 J=2,JE1
       DO 23 I=2,IE1
       ABL=ABS(TRF(I,J,KUP)-TRF(I,J,KLO))
       ZWABL=ABS(TRF(I,J,KUP)+TRF(I,J,KLO)-2.*TRF(I,J,K))
       WUP=0.5*(WO(I,J,K)+WO(I,J,KUP))
       WUM=0.5*(WO(I,J,K)+WO(I,J,KLO))
       UP=0.5*DT*DLYP(I,J)*DLXP(I,J)*(WUP+ABS(WUP))
       UM=0.5*DT*DLYP(I,J)*DLXP(I,J)*(ABS(WUM)-WUM)
       SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))
       SCH=MIN(WETO(I,J,KLO),SCHA*S1O(I,J,K)/(UP+UM+1.E-20))
       SCHI=1.-SCH
       TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J,KUP)))
       TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J,KLO)))
       WTP(I,J,K)=UP
       WTM(I,J,K)=UM
23     CONTINUE
!
!$OMP DO
       DO 27 J=1,JE
       DO 27 I=1,IE
       WUP=0.5*(WO(I,J,KEP)+WO(I,J,KE))
       UP=0.5*DT*DLYP(I,J)*DLXP(I,J)*(WUP+ABS(WUP))
       TRP(I,J,KEP)=UP*(0.5*(TRF(I,J,KE)+TRF(I,J,KEP)))
       WTP(I,J,KEP)=UP
       
       WUM=0.5*(WO(I,J,1)+WO(I,J,2))       
       UM=0.5*DT*DLYP(I,J)*DLXP(I,J)*(ABS(WUM)-WUM)
       TRM(I,J,1)=UM*(0.5*(TRF(I,J,1)+TRF(I,J,2)))
       WTM(I,J,1)=UM
27     CONTINUE
!
! RJ: no boundary exchange necessary here, since TRP, WTP, TRM, WTM
!     are used only in the inner domain in the following loop
!$OMP DO
       DO 26 K=2,KE
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
!$OMP SINGLE
      CALL bounds_exch(TRF)
!$OMP END SINGLE
!
!      TRANSPORTS IN X-DIRECTION
!
!$OMP DO
       DO 3 K=2,KE
       DO 3 J=2,JE1
       DO 3 I=2,IE1
!        IF (K .EQ. 1) THEN
!          DDWU=(DDUO(I,J,K)+ZO(I,J))/2.0
!        ELSEIF (K .EQ. 2) THEN
!          DDWU=(DDUO(I,J,K-1)+DDUO(I,J,K)+ZO(I,J))/2.0
!        ELSE
!          DDWU=(DDUO(I,J,K-1)+DDUO(I,J,K))/2.0
!        ENDIF
       ABL=ABS(TRF(I+1,J,K)-TRF(I-1,J,K))
       ZWABL=ABS(TRF(I+1,J,K)+TRF(I-1,J,K)-2.*TRF(I,J,K))
!       UOS=0.5*(UKO(I,J,K)+UKO(I,J,K-1))*DLYU(I,J)
!       UWE=0.5*(UKO(I-1,J,K)+UKO(I-1,J,K-1))*DLYU(I-1,J)
!       UP=0.5*DT*DDWU*(UOS+ABS(UOS))
!       UM=0.5*DT*DDWU*(ABS(UWE)-UWE)
       UOS=0.5*(UKO(I,J,K)*DDUO(I,J,K)+UKO(I,J,K-1)*DDUO(I,J,K-1))
       UWE=0.5*(UKO(I-1,J,K)*DDUO(I-1,J,K)+UKO(I-1,J,K-1)               &
     & *DDUO(I-1,J,K-1))
       UP=0.5*DT*DLYU(I,J)*(UOS+ABS(UOS))
       UM=0.5*DT*DLYU(I-1,J)*(ABS(UWE)-UWE)
       SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))*WETO(I,J,K)                &
     &       *WETO(I+1,J,K)*WETO(I-1,J,K)
       SCH=MIN(1.,SCHA*S1O(I,J,K)/(UP+UM+1.E-20))
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
       DO 6 K=2,KE
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
       DO 13 K=2,KE
       DO 13 J=2,JE1
       DO 13 I=2,IE1
!        IF (K .EQ. 1) THEN
!          DDWV=(DDUE(I,J,K)+ZO(I,J))/2.0
!        ELSEIF (K .EQ. 2) THEN
!          DDWV=(DDUE(I,J,K-1)+DDUE(I,J,K)+ZO(I,J))/2.0
!        ELSE
!          DDWV=(DDUE(I,J,K-1)+DDUE(I,J,K))/2.0
!        ENDIF
       ABL=ABS(TRF(I,J-1,K)-TRF(I,J+1,K))
       ZWABL=ABS(TRF(I,J-1,K)+TRF(I,J+1,K)-2.*TRF(I,J,K))
!       VSU=0.5*(VKE(I,J,K)+VKE(I,J,K-1))*DLXV(I,J)
!       VNO=0.5*(VKE(I,J-1,K)+VKE(I,J-1,K-1))*DLXV(I,J-1)
!       UP=0.5*DT*DDWV*(VNO+ABS(VNO))
!       UM=0.5*DT*DDWV*(ABS(VSU)-VSU)
       VSU=0.5*(VKE(I,J,K)*DDUE(I,J,K)+VKE(I,J,K-1)*DDUE(I,J,K-1))
       VNO=0.5*(VKE(I,J-1,K)*DDUE(I,J-1,K)+VKE(I,J-1,K-1)               &
     & *DDUE(I,J-1,K-1))
       UP=0.5*DT*DLXV(I,J-1)*(VNO+ABS(VNO))
       UM=0.5*DT*DLXV(I,J)*(ABS(VSU)-VSU)
       SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))*WETO(I,J,K)                &
     &     *WETO(I,J-1,K)*WETO(I,J+1,K)
       SCH=MIN(1.,SCHA*S1O(I,J,K)/(UP+UM+1.E-20))
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
       DO 16 K=2,KE
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

!$OMP SINGLE
      CALL bounds_exch(TRF)
!$OMP END SINGLE
!
!
!$OMP END PARALLEL

      RETURN
      END
