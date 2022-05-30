      SUBROUTINE OCVAD(TRF0)
!
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
!
      DIMENSION TRF0(IE,JE-J_start+1,KE)

      DIMENSION TRP(1:IE,J_start:JE,1:KEP),TRM(1:IE,J_start:JE,0:KE)
      DIMENSION TRF(1:IE,J_start:JE,1:KE),TRVOL(1:IE,J_start:JE,1:KE)
      DIMENSION TR1(1:IE,J_start:JE,1:KE)
      DIMENSION WTP(1:IE,J_start:JE,1:KEP),WTM(1:IE,J_start:JE,0:KE)

      TRF(:,:,:)=TRF0(:,:,:)

!$OMP PARALLEL PRIVATE(i,j,k,klo,kup,abl,zwabl,wup,wlo,up,um,scha,sch,schi,rn,uos,uwe,vsu,vno)
!$OMP DO
      DO 1 K=1,KE
      DO 1 J=J_start,JE
      DO 1 I=1,IE
      TRVOL(I,J,K)=AMSUE(I,J,K)*DLXV(I,J)*DLYV(I,J)*DDUE(I,J,K)
      TR1(I,J,K)=TRF(I,J,K)
1     CONTINUE
604   FORMAT(13F9.5)
!
!$OMP DO
       TRP(:,:,:)=0.
       TRM(:,:,:)=0.
       WTP(:,:,:)=0.
       WTM(:,:,:)=0.

!$OMP DO
       DO 23 K=1,KE
       KLO=MIN(K+1,KE)
       KUP=MAX(K-1,1)
       DO 23 J=(J_start+1),JE1
       DO 23 I=2,IE1
       ABL=ABS(TRF(I,J,KUP)-TRF(I,J,KLO))
       ZWABL=ABS(TRF(I,J,KUP)+TRF(I,J,KLO)-2.*TRF(I,J,K))
       WUP=0.5*(WO(I,J,K)+WO(I,J+1,K))
       WLO=0.5*(WO(I,J,K+1)+WO(I,J+1,K+1))
       UP=0.5*DT*DLYP(I,J)*DLXP(I,J)*(WUP+ABS(WUP))
       UM=0.5*DT*DLYP(I,J)*DLXP(I,J)*(ABS(WLO)-WLO)
!
       SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))
       SCH=MIN(1.,SCHA*TRVOL(I,J,K)/(UP+UM+1.E-20))
       SCHI=1.-SCH
!
       TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J,KUP)))
       TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J,KLO)))
       WTP(I,J,K)=UP
       WTM(I,J,K)=UM
23     CONTINUE

!$OMP DO
       DO 27 J=(J_start+1),JE1
       DO 27 I=2,IE1
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
       DO 24 J=(J_start+1),JE1
       DO 24 I=2,IE1
       IF(AMSUE(I,J,K).GT.0.5) THEN
!
       RN=TRVOL(I,J,K)+WTP(I,J,K+1)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J,K-1)
       TR1(I,J,K)=(TRF(I,J,K)*TRVOL(I,J,K)+TRP(I,J,K+1)-TRP(I,J,K)      &
     & -TRM(I,J,K)+TRM(I,J,K-1))/RN
       TRVOL(I,J,K)=RN
        ENDIF
24      CONTINUE
       DO 25 J=(J_start+1),JE1
       DO 25 I=2,IE1
       TRF(I,J,K)=TR1(I,J,K)
25      CONTINUE 
26     CONTINUE

!$OMP SINGLE
      CALL bounds_exch(TR1)
      CALL bounds_exch(TRVOL)
      CALL bounds_exch(TRF)
!$OMP END SINGLE
!
      TRP=0.0
      TRM=0.0
      WTP=0.0
      WTM=0.0
!
!$OMP DO
       DO 3 K=1,KE
       DO 3 J=(J_start+1),JE1
       DO 3 I=2,IE1
       ABL=ABS(TRF(I+1,J,K)-TRF(I-1,J,K))
       ZWABL=ABS(TRF(I+1,J,K)+TRF(I-1,J,K)-2.*TRF(I,J,K))
       UOS=0.5*(UKO(I,J,K)*DDUO(I,J,K)+UKO(I,J+1,K)*DDUO(I,J+1,K))
       UWE=0.5*(UKO(I-1,J,K)*DDUO(I-1,J,K)+UKO(I-1,J+1,K)*DDUO(I-1,J+1,K))
       UP=0.5*DT*(0.5*(DLYU(I,J)+DLYU(I,J+1)))*(UOS+ABS(UOS))
       UM=0.5*DT*(ABS(UWE)-UWE)*(0.5*(DLYU(I-1,J)+DLYU(I-1,J+1)))
       SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))
!
       SCH=MIN(1.,SCHA*TRVOL(I,J,K)/(UP+UM+1.E-20))
       SCHI=1.-SCH
!
       TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I+1,J,K)))
       TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I-1,J,K)))
       WTP(I,J,K)=UP
       WTM(I,J,K)=UM
3      CONTINUE

      IF (have_g_is) THEN
       I=1
       DO 31 K=1,KE
       DO 31 J=(J_start+1),JE1
       UOS=0.5*(UKO(I,J,K)*DDUO(I,J,K)+UKO(I,J+1,K)*DDUO(I,J+1,K))
       UWE=0.5*(UKO(I-1,J,K)*DDUO(I-1,J,K)+UKO(I-1,J+1,K)*DDUO(I-1,J+1,K))
       UP=0.5*DT*(0.5*(DLYU(I,J)+DLYU(I,J+1)))*(UOS+ABS(UOS))
       UM=0.5*DT*(ABS(UWE)-UWE)*(0.5*(DLYU(I-1,J)+DLYU(I-1,J+1)))
       TRP(I,J,K)=UP*TRF(I,J,K)
       TRM(I,J,K)=UM*TRF(I,J,K)
       WTP(I,J,K)=UP
       WTM(I,J,K)=UM
31     CONTINUE
      ENDIF

      IF (have_g_ie) THEN
       I=IE
       DO 32 K=1,KE
       DO 32 J=(J_start+1),JE1
       UOS=0.5*(UKO(I,J,K)*DDUO(I,J,K)+UKO(I,J+1,K)*DDUO(I,J+1,K))
       UWE=0.5*(UKO(I-1,J,K)*DDUO(I-1,J,K)+UKO(I-1,J+1,K)*DDUO(I-1,J+1,K))
       UP=0.5*DT*(0.5*(DLYU(I,J)+DLYU(I,J+1)))*(UOS+ABS(UOS))
       UM=0.5*DT*(ABS(UWE)-UWE)*(0.5*(DLYU(I-1,J)+DLYU(I-1,J+1)))
       TRP(I,J,K)=UP*TRF(I,J,K)
       TRM(I,J,K)=UM*TRF(I,J,K)
       WTP(I,J,K)=UP
       WTM(I,J,K)=UM
32     CONTINUE
      ENDIF

!$OMP SINGLE
      CALL bounds_exch(TRP)
      CALL bounds_exch(TRM)
      CALL bounds_exch(WTP)
      CALL bounds_exch(WTM)
!$OMP END SINGLE
!$OMP DO
       DO 6 K=1,KE
       DO 4 J=(J_start+1),JE1
       DO 4 I=2,IE1
       IF(AMSUE(I,J,K).GT.0.5) THEN
!
       RN=TRVOL(I,J,K)+WTP(I-1,J,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I+1,J,K)
       TR1(I,J,K)=(TRF(I,J,K)*TRVOL(I,J,K)+TRP(I-1,J,K)-TRP(I,J,K)      &
     & -TRM(I,J,K)+TRM(I+1,J,K))/RN
       TRVOL(I,J,K)=RN
       ENDIF
4      CONTINUE
       DO 5 J=(J_start+1),JE1
       DO 5 I=2,IE1
       TRF(I,J,K)=TR1(I,J,K)
5      CONTINUE 
6     CONTINUE
!$OMP SINGLE
      CALL bounds_exch(TR1)
      CALL bounds_exch(TRVOL)
      CALL bounds_exch(TRF)
!$OMP END SINGLE
!
      TRP=0.0
      TRM=0.0
      WTP=0.0
      WTM=0.0
!
!$OMP DO
       DO 13 K=1,KE
       DO 13 J=(J_start+1),JE1
       DO 13 I=2,IE1
       ABL=ABS(TRF(I,J-1,K)-TRF(I,J+1,K))
       ZWABL=ABS(TRF(I,J-1,K)+TRF(I,J+1,K)-2.*TRF(I,J,K))
       VNO=0.5*(VKE(I,J,K)*DDUE(I,J,K)+VKE(I,J-1,K)*DDUE(I,J-1,K))
       VSU=0.5*(VKE(I,J,K)*DDUE(I,J,K)+VKE(I,J+1,K)*DDUE(I,J+1,K))
       UP=0.5*DT*DLXP(I,J)*(VNO+ABS(VNO))
       UM=0.5*DT*DLXP(I,J+1)*(ABS(VSU)-VSU)
!
       SCHA=MAX(0.,(ABL-ZWABL)/(ABL+1.E-20))
       SCH=MIN(1.,SCHA*TRVOL(I,J,K)/(UP+UM+1.E-20))
!
       SCHI=1.-SCH
       TRP(I,J,K)=UP*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J-1,K)))
       TRM(I,J,K)=UM*(SCHI*TRF(I,J,K)+SCH*0.5*(TRF(I,J,K)+TRF(I,J+1,K)))
       WTP(I,J,K)=UP
       WTM(I,J,K)=UM
13     CONTINUE

      IF (have_g_js) THEN
       J=0
       DO 131 K=1,KE
       DO 131 I=2,IE1
       VNO=VKE(I,J,K)*DDUE(I,J,K)
       VSU=0.5*(VKE(I,J,K)*DDUE(I,J,K)+VKE(I,J+1,K)*DDUE(I,J+1,K))      
       UP=0.5*DT*DLXU(I,J)*(VNO+ABS(VNO))
       UM=0.5*DT*DLXP(I,J+1)*(ABS(VSU)-VSU)
       TRP(I,J,K)=UP*TRF(I,J,K)
       TRM(I,J,K)=UM*TRF(I,J,K)
       WTP(I,J,K)=UP
       WTM(I,J,K)=UM
131    CONTINUE
      ENDIF

      IF (have_g_je) THEN
       J=JE
       DO 132 K=1,KE
       DO 132 I=2,IE1
       VNO=0.5*(VKE(I,J,K)*DDUE(I,J,K)+VKE(I,J-1,K)*DDUE(I,J-1,K))
       VSU=VKE(I,J,K)*DDUE(I,J,K)
       UP=0.5*DT*DLXP(I,J)*(VNO+ABS(VNO))
       UM=0.5*DT*DLXU(I,J)*(ABS(VSU)-VSU)
       TRP(I,J,K)=UP*TRF(I,J,K)
       TRM(I,J,K)=UM*TRF(I,J,K)
       WTP(I,J,K)=UP
       WTM(I,J,K)=UM
132    CONTINUE
      ENDIF

!$OMP SINGLE
      CALL bounds_exch(TRP)
      CALL bounds_exch(TRM)
      CALL bounds_exch(WTP)
      CALL bounds_exch(WTM)
!$OMP END SINGLE
!$OMP DO
       DO 16 K=1,KE
       DO 14 J=(J_start+1),JE1
       DO 14 I=2,IE1
       IF(AMSUE(I,J,K).GT.0.5) THEN
!
       RN=TRVOL(I,J,K)+WTP(I,J+1,K)-WTP(I,J,K)-WTM(I,J,K)+WTM(I,J-1,K)
       TR1(I,J,K)=(TRF(I,J,K)*TRVOL(I,J,K)+TRP(I,J+1,K)-TRP(I,J,K)      &
     & -TRM(I,J,K)+TRM(I,J-1,K))/RN
       TRVOL(I,J,K)=RN
       ENDIF
14      CONTINUE
       DO 15 J=(J_start+1),JE1
       DO 15 I=2,IE1
       TRF(I,J,K)=TR1(I,J,K)
15      CONTINUE 
16     CONTINUE
!$OMP END PARALLEL
      CALL bounds_exch(TR1)
      CALL bounds_exch(TRVOL)
      CALL bounds_exch(TRF)
!
      TRP=0.0
      TRM=0.0
      WTP=0.0
      WTM=0.0
!
      TRF0(:,:,:)=TRF(:,:,:)
!
      RETURN
      END
