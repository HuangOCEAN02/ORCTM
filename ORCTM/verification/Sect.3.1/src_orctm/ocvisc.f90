      SUBROUTINE OCVISC
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
!-=====================================================================
!
!     SBR OCVISC
!
!     COMPUTES MOMENTUM DIFFUSION
!     VERTICALLY (IMPLICITE)
!     RAYLEIGH BOTTOM FRICTION
!     HORIZONTALLY BIHARMONIC
!
!     LAST MODIFIED UWE MIKOLAJEWICZ 12.12.99
!     RAYLEIGH FRICTION NOW QUADRATIC!!!!
!     DX,DY INCLUDED IN BIHARMONIC FRICTION!!
!
!     MODIFIED 4.2.2000 UWE MIKOLAJEWICZ
!      INCLUDE QUADRATIC BOTTOM FRICTION IN VERTICAL DIFFUSION, 
!              REMOVE ERROR
!
!
!     BIHARMONIC DIFFUSION COEFFICIENTS, DX,DY ON PSI-POINTS
!
      implicit none
      integer i,j,k,l
      REAL DLXPSI(IE,J_start:JE),DLYPSI(I_start:IE,JE)
      REAL AULUX(IE,JE),AULUY(IE,JE),AULVX(IE,JE),AULVY(IE,JE),         &
     &   SPEEDU(IE,JE),SPEEDV(IE,JE),AULWX(IE,JE),AULWY(IE,JE)
      REAL DZU(IE,JE,KE),DZV(IE,JE,KE),DZP(IE,JE,KE)
      real TRIDSY(ke,3)
      real avup,avlo,aulap,bofric,rayfric
      real minvalue
!
!
!     VERTICAL DIFFUSION OF MOMENTUM
!
      minvalue=MIN(1.0,DZW(1)/10.0)
!
!     BOFRIC  FRICTION COEFFICENT AT BOTTOM [M**2/S]
      BOFRIC=0.E-3
!
!     COEFFICIENT FOR QUADRATIC BOTTOM FRICTION []
      RAYFRIC=0.E-3
!
!     SPEED IN BOTTOM LAYER AT OLD TIME STEP
!       ON U AND V POINTS
!

!$OMP PARALLEL PRIVATE(i,j,k,l,avup,avlo,TRIDSY,aulap)

!$OMP DO
      DO J=1,JE
       DO I=1,IE
        SPEEDU(I,J)=0.
        SPEEDV(I,J)=0.
       ENDDO
      ENDDO
!
!$OMP DO
      DO J=2,JE1
       DO K=1,KE-1
        DO I=(I_start+1),IE1
         IF(AMSUO(I,J,K)-AMSUO(I,J,K+1).GT.0.5)THEN
          SPEEDU(I,J)=SQRT(UOO(I,J,K)**2+                               &
     &      (0.25*(VOE(I,J,K)+VOE(I+1,J,K)                              &
     &      +VOE(I,J-1,K)+VOE(I+1,J-1,K)))**2)
          SPEEDU(I,J)=MAX(SPEEDU(I,J),1.E-3)
         ENDIF
        ENDDO
       ENDDO
       K=KE
       DO I=(I_start+1),IE1
        IF(AMSUO(I,J,K).GT.0.5)THEN
          SPEEDU(I,J)=SQRT(UOO(I,J,K)**2+                               &
     &      (0.25*(VOE(I,J,K)+VOE(I+1,J,K)                              &
     &      +VOE(I,J-1,K)+VOE(I+1,J-1,K)))**2)
          SPEEDU(I,J)=MAX(SPEEDU(I,J),1.E-3)
         ENDIF
       ENDDO
      ENDDO
  
      DO J=(J_start+1),JE1
       DO K=1,KE-1
        DO I=2,IE1
         IF(AMSUE(I,J,K)-AMSUE(I,J,K+1).GT.0.5)THEN
          SPEEDV(I,J)=SQRT(VOE(I,J,K)**2+                               &
     &      (0.25*(UOO(I,J,K)+UOO(I-1,J,K)                              &
     &      +UOO(I,J+1,K)+UOO(I-1,J+1,K)))**2)
          SPEEDV(I,J)=MAX(SPEEDV(I,J),1.E-3)
         ENDIF
        ENDDO
       ENDDO
       K=KE
       DO I=2,IE1
         IF(AMSUE(I,J,K).GT.0.5)THEN
          SPEEDV(I,J)=SQRT(VOE(I,J,K)**2+                               &
     &      (0.25*(UOO(I,J,K)+UOO(I-1,J,K)                              &
     &      +UOO(I,J+1,K)+UOO(I-1,J+1,K)))**2)
          SPEEDV(I,J)=MAX(SPEEDV(I,J),1.E-3)
         ENDIF
       ENDDO
      ENDDO
!
!$OMP SINGLE
      CALL bounds_exch(SPEEDU)
      CALL bounds_exch(SPEEDV)
!$OMP END SINGLE
!
!     OUTER LOOP
!
!DO V
      DO I=2,IE1
      DO J=(J_start+1),JE1
!
      DO 11 K=1,KE-1
      AVUP=0.5*(AVO(I,J,K)+AVO(I,J+1,K))*AMSUE(I,J,K)
      AVLO=0.5*(AVO(I,J,K+1)+AVO(I,J+1,K+1))*AMSUE(I,J,K+1)             &
     &  +(BOFRIC+RAYFRIC*DDUE(I,J,K)*SPEEDV(I,J))                       &
     &      *(AMSUE(I,J,K)-AMSUE(I,J,K+1))
      TRIDSY(K,1)= - DT * AVUP   * AMSUE(I,J,K) * DI(K)                 &
     &          /(ALMZER+DDUE(I,J,K))
      TRIDSY(K,3)= - DT * AVLO * AMSUE(I,J,K) * DI(K+1)                 &
     &          /(ALMZER+DDUE(I,J,K))
      TRIDSY(K,2)= 1. - TRIDSY(K,1) - TRIDSY(K,3)
      TRIDSY(K,3)=TRIDSY(K,3)*AMSUE(I,J,K+1)
11    CONTINUE
!
      K=KE
      AVUP=0.5*(AVO(I,J,K)+AVO(I,J+1,K))*AMSUE(I,J,K)
      AVLO=(BOFRIC+RAYFRIC*                                             &
     &     DDUE(I,J,K)*SPEEDV(I,J))*AMSUE(I,J,K)
      TRIDSY(K,1)= - DT * AVUP   * AMSUE(I,J,K) * DI(K)                 &
     &        /(ALMZER+DDUE(I,J,K))
      TRIDSY(K,3)= - DT * AVLO * AMSUE(I,J,K) * DI(K+1)                 &
     &        /(ALMZER+DDUE(I,J,K))

      TRIDSY(K,2)= 1. - TRIDSY(K,1) - TRIDSY(K,3)
      TRIDSY(K,3)=0.
!
      DO 12 K=2,KE
      TRIDSY(K-1,1)  = TRIDSY(K,1) / TRIDSY(K-1,2)
      TRIDSY(K,2)= TRIDSY(K,2)-TRIDSY(K-1,3)*TRIDSY(K,1)/TRIDSY(K-1,2)
12    CONTINUE
!
      DO 14 K=1,KE
      VK1E(I,J,K)=VOE(I,J,K)
14    CONTINUE
!
      DO 15 K=2,KE
      VK1E(I,J,K)=VK1E(I,J,K)-TRIDSY(K-1,1)*VK1E(I,J,K-1)
15    CONTINUE
!
      K=KE
      VKE(I,J,K)=AMSUE(I,J,K)*  VK1E(I,J,K) / TRIDSY(K,2)
!
      DO 17 K=1,KE1
      L=KE-K
      VKE(I,J,L)=(VK1E(I,J,L)-TRIDSY(L,3)*VKE(I,J,L+1))*AMSUE(I,J,L)    &
     &                  /TRIDSY(L,2)
17    CONTINUE
!
      ENDDO
      ENDDO
!
!DO U
      DO J=2,JE1
      DO I=(I_start+1),IE1
!
      DO 21 K=1,KE-1
      AVUP=0.5*(AVO(I,J,K)+AVO(I+1,J,K))*AMSUO(I,J,K)
      AVLO=0.5*(AVO(I,J,K+1)+AVO(I+1,J,K+1))*AMSUO(I,J,K+1)             &
     &+(BOFRIC+RAYFRIC*DDUO(I,J,K)*SPEEDU(I,J))                         &
     &    *(AMSUO(I,J,K)-AMSUO(I,J,K+1))
      TRIDSY(K,1)= - DT * AVUP   * AMSUO(I,J,K)* DI(K)                  &
     &        /(almzer+dduo(i,j,k))
      TRIDSY(K,3)= - DT * AVLO * amsuo(I,J,K) * DI(K+1)                 &
     &         /(almzer+dduo(i,j,k))
      TRIDSY(K,2)= 1. - TRIDSY(K,1) - TRIDSY(K,3)
      TRIDSY(K,3)=TRIDSY(K,3)*AMSUO(I,J,K+1)
21    CONTINUE
!
      K=KE
      AVUP=0.5*(AVO(I,J,K)+AVO(I+1,J,K))*AMSUO(I,J,K)
      AVLO=(BOFRIC+RAYFRIC                                              &
     &     *DDUO(I,J,K)*SPEEDU(I,J))*AMSUO(I,J,K)
      TRIDSY(K,1)= - DT * AVUP   * amsuO(I,J,K) * DI(K)                 &
     &         /(almzer+dduo(i,j,k))
      TRIDSY(K,3)= - DT * AVLO * amsuO(I,J,K) * DI(K+1)                 &
     &         /(almzer+dduo(i,j,k))
      TRIDSY(K,2)= 1. - TRIDSY(K,1) - TRIDSY(K,3)
      TRIDSY(K,3)=0.
!
      DO 22 K=2,KE
      TRIDSY(K-1,1) = TRIDSY(K,1) / TRIDSY(K-1,2)
      TRIDSY(K,2)=TRIDSY(K,2)-TRIDSY(K-1,3)*TRIDSY(K,1)/TRIDSY(K-1,2)
22    CONTINUE
!
      DO 24 K=1,KE
      UK1O(I,J,K)=UOO(I,J,K)
24    CONTINUE
!
      DO 25 K=2,KE
      UK1O(I,J,K) = UK1O(I,J,K) - TRIDSY(K-1,1) * UK1O(I,J,K-1)
25    CONTINUE
!
      K=KE
      UKO(I,J,K) =AMSUO(I,J,K)*   UK1O(I,J,K) / TRIDSY(K,2)
!
      DO 27 K=1,KE1
      L=KE-K
      UKO(I,J,L)=(UK1O(I,J,L)-TRIDSY(L,3)*UKO(I,J,L+1))*AMSUO(I,J,L)    &
     &          /TRIDSY(L,2)
27    CONTINUE
!
      ENDDO
      ENDDO
!
      CALL bounds_exch(UKO)
      CALL bounds_exch(VKE)
!
#ifdef NON_HYDROSTATIC
! DO W       k from 2 to ke
      DO J=2,JE1
      DO I=2,IE1
!
      DO 31 K=2,KE-1
      AVUP=0.5*(AVO(I,J,K-1)+AVO(I,J,K))*WETO(I,J,K-1)
      AVLO=0.5*(AVO(I,J,K+1)+AVO(I,J,K))*WETO(I,J,K)
      TRIDSY(K,1)= - 2.0 * DT * AVUP * WETO(I,J,K) * DI(K-1)            &
     &            / (2.0 * almzer + ddpo(i,j,k-1) + ddpo(i,j,k))
      TRIDSY(K,3)= - 2.0 * DT * AVLO * WETO(I,J,K) * DI(K)              &
     &            / (2.0 * almzer + ddpo(i,j,k-1) + ddpo(i,j,k))
      TRIDSY(K,2)= 1. - TRIDSY(K,1) - TRIDSY(K,3)
31    CONTINUE
!
      K=KE
      AVUP=0.5*(AVO(I,J,K-1)+AVO(I,J,K))*WETO(I,J,K-1)
      AVLO=0.0
      TRIDSY(K,1)= - 2.0 * DT * AVUP * WETO(I,J,K) * DI(K-1)            &
     &            / (2.0 * almzer + ddpo(i,j,k-1) + ddpo(i,j,k))
      TRIDSY(K,3)= - 2.0 * DT * AVLO * WETO(I,J,K) * DI(K)              &
     &            / (2.0 * almzer + ddpo(i,j,k-1) + ddpo(i,j,k))
      TRIDSY(K,2)= 1. - TRIDSY(K,1) - TRIDSY(K,3)
      TRIDSY(K,3)=0.
!
      DO 32 K=3,KE
      TRIDSY(K-1,1) = TRIDSY(K,1) / TRIDSY(K-1,2)
      TRIDSY(K,2)=TRIDSY(K,2)-TRIDSY(K-1,3)*TRIDSY(K,1)/TRIDSY(K-1,2)
32    CONTINUE
!
      DO 34 K=2,KE
      WN1O(I,J,K)=WNO(I,J,K)
34    CONTINUE
!
      DO 35 K=3,KE
      WN1O(I,J,K) = WN1O(I,J,K) - TRIDSY(K-1,1) * WN1O(I,J,K-1)
35    CONTINUE
!
      K=KE
      WNO(I,J,K) = WETO(I,J,K) * WN1O(I,J,K) / TRIDSY(K,2)
!
      DO 37 K=KE1,2,-1
      WNO(I,J,K)=(WN1O(I,J,K)-TRIDSY(K,3)*WNO(I,J,K+1))*WETO(I,J,K)     &
     &          /TRIDSY(K,2)
37    CONTINUE
!
      ENDDO
      ENDDO
!
      CALL bounds_exch(WNO)
#endif
!
!      END VERTICAL DIFFUSION
!
!
!
!     HORIZONTAL DIFFUSION OF MOMENTUM
!
!     COEFFICIENTS FOR BIHARMONIC HORIZONTAL MOMENTUM DIFFUSION,
!         DEFINED ON P- AND PSI-POINTS
!
!$OMP DO

      DO J=1,JE
       DO I=1,IE1
        DLYPSI(I,J)=0.5*(DLYV(I,J)+DLYV(I+1,J))
       ENDDO
      ENDDO
      DLYPSI(I_start,1:JE)=DLYV( 1,1:JE)
      DLYPSI(IE,     1:JE)=DLYV(IE,1:JE)

      CALL bounds_exch(DLYPSI) 

      DO I=1,IE
       DO J=1,JE1
        DLXPSI(I,J)=0.5*(DLXU(I,J)+DLXU(I,J+1))
       ENDDO
      ENDDO
      DLXPSI(1:IE,J_start)=DLXU(1:IE, 1)
      DLXPSI(1:IE,     JE)=DLXU(1:IE,JE)

      CALL bounds_exch(DLXPSI)

#ifdef AULREDUV
      DO J=1,JE
       DO I=1,IE
        AULUX(I,J)=AULAPUV*1.E4*MAX(DLXP(I,J),minvalue)**3
        AULUY(I,J)=AULAPUV*1.E4*MAX(minvalue,DLYPSI(I,J))**3
        AULVY(I,J)=AULAPUV*1.E4*MAX(minvalue,DLYP(I,J))**3
        AULVX(I,J)=AULAPUV*1.E4*MAX(DLXPSI(I,J),minvalue)**3
        AULWX(I,J)=AULAPUV*1.E4*MAX(DLXU(I,J),minvalue)**3
        AULWY(I,J)=AULAPUV*1.E4*MAX(DLYV(I,J),minvalue)**3
       ENDDO
      ENDDO
#else
      DO J=1,JE
       DO I=1,IE
        AULUX(I,J)=AULAPUV*MAX(DLXP(I,J),minvalue)**4
        AULUY(I,J)=AULAPUV*MAX(minvalue,DLYPSI(I,J))**4
        AULVY(I,J)=AULAPUV*MAX(minvalue,DLYP(I,J))**4
        AULVX(I,J)=AULAPUV*MAX(DLXPSI(I,J),minvalue)**4
        AULWX(I,J)=AULAPUV*MAX(DLXU(I,J),minvalue)**4
        AULWY(I,J)=AULAPUV*MAX(DLYV(I,J),minvalue)**4
       ENDDO
      ENDDO
#endif /*AULREDUV*/
!
!$OMP SINGLE
      CALL bounds_exch(AULUX,AULUY,AULVX,AULVY,AULWX,AULWY)
!$OMP END SINGLE
!
      DO 99 K=1,KE
      DO 99 J=1,JE
      DO 99 I=1,IE
        DZU(I,J,K)=MAX(DDUO(I,J,K),minvalue)
        DZV(I,J,K)=MAX(DDUE(I,J,K),minvalue)
        DZP(I,J,K)=MAX(DDPO(I,J,K),minvalue)
99    CONTINUE
!
!     MAIN LOOP OVER LEVELS
!
      AULAP=AULAPUV
!
#ifdef FREESLIP
!     BIHARMONIC HORIZONTAL MOMENTUM DIFFUSION
!
!$OMP DO
      DO K=1,KE
! 
      DO J=2,JE1
      DO I=(I_start+1),IE1
        UK1O(I,J,K)=AMSUO(I,J,K)*                                       &
     &((AMSUO(I+1,J,K)*(UKO(I+1,J,K)-UKO(I,J,K))/DLXP(I+1,J)            &
     &-AMSUO(I-1,J,K)*(UKO(I,J,K)-UKO(I-1,J,K))/DLXP(I,J))/DLXU(I,J)    &
     &+(AMSUO(I,J-1,K)*(UKO(I,J-1,K)-UKO(I,J,K))/DLYPSI(I,J-1)          &
     &-AMSUO(I,J+1,K)*(UKO(I,J,K)-UKO(I,J+1,K))/DLYPSI(I,J))/DLYU(I,J))
      ENDDO
      ENDDO
!
      DO J=(J_start+1),JE1
      DO I=2,IE1
        VK1E(I,J,K)=AMSUE(I,J,K)*                                       &
     &((AMSUE(I+1,J,K)*(VKE(I+1,J,K)-VKE(I,J,K))/DLXPSI(I,J)            &
     &-AMSUE(I-1,J,K)*(VKE(I,J,K)-VKE(I-1,J,K))/DLXPSI(I-1,J))/DLXV(I,J)&
     &+(AMSUE(I,J-1,K)*(VKE(I,J-1,K)-VKE(I,J,K))/DLYP(I,J)              &
     &-AMSUE(I,J+1,K)*(VKE(I,J,K)-VKE(I,J+1,K))/DLYP(I,J+1))/DLYV(I,J))
      ENDDO
      ENDDO
!
      if(have_g_is)then
      I=0
      DO J=2,JE1
        UK1O(I,J,K)=AMSUO(I,J,K)*                                       &
     &((AMSUO(I+1,J,K)*(UKO(I+1,J,K)-UKO(I,J,K))/DLXP(I+1,J))/DLXU(I,J) &
     &+(AMSUO(I,J-1,K)*(UKO(I,J-1,K)-UKO(I,J,K))/DLYPSI(I,J-1)          &
     &-AMSUO(I,J+1,K)*(UKO(I,J,K)-UKO(I,J+1,K))/DLYPSI(I,J))/DLYU(I,J))
      ENDDO
      I=1
      DO J=(J_start+1),JE1
        VK1E(I,J,K)=AMSUE(I,J,K)*                                       &
     &((AMSUE(I+1,J,K)*(VKE(I+1,J,K)-VKE(I,J,K))/DLXPSI(I,J))/DLXV(I,J) &
     &+(AMSUE(I,J-1,K)*(VKE(I,J-1,K)-VKE(I,J,K))/DLYP(I,J)              &
     &-AMSUE(I,J+1,K)*(VKE(I,J,K)-VKE(I,J+1,K))/DLYP(I,J+1))/DLYV(I,J))
      ENDDO
      endif
!
      if(have_g_ie)then
      I=IE
      DO J=2,JE1
        UK1O(I,J,K)=AMSUO(I,J,K)*                                       &
     &((-AMSUO(I-1,J,K)*(UKO(I,J,K)-UKO(I-1,J,K))/DLXP(I,J))/DLXU(I,J)  &
     &+(AMSUO(I,J-1,K)*(UKO(I,J-1,K)-UKO(I,J,K))/DLYPSI(I,J-1)          &
     &-AMSUO(I,J+1,K)*(UKO(I,J,K)-UKO(I,J+1,K))/DLYPSI(I,J))/DLYU(I,J))
      ENDDO
      DO J=(J_start+1),JE1
        VK1E(I,J,K)=AMSUE(I,J,K)*                                       &
     &((-AMSUE(I-1,J,K)*(VKE(I,J,K)-VKE(I-1,J,K))/DLXPSI(I-1,J))/DLXV(I,J)&
     &+(AMSUE(I,J-1,K)*(VKE(I,J-1,K)-VKE(I,J,K))/DLYP(I,J)              &
     &-AMSUE(I,J+1,K)*(VKE(I,J,K)-VKE(I,J+1,K))/DLYP(I,J+1))/DLYV(I,J))
      ENDDO
      endif
!
      if(have_g_js)then
      J=1
      DO I=(I_start+1),IE1
        UK1O(I,J,K)=AMSUO(I,J,K)*                                       &
     &((AMSUO(I+1,J,K)*(UKO(I+1,J,K)-UKO(I,J,K))/DLXP(I+1,J)            &
     &-AMSUO(I-1,J,K)*(UKO(I,J,K)-UKO(I-1,J,K))/DLXP(I,J))/DLXU(I,J)    &
     &+(-AMSUO(I,J+1,K)*(UKO(I,J,K)-UKO(I,J+1,K))/DLYPSI(I,J))/DLYU(I,J))
      ENDDO
      J=0
      DO I=2,IE1
        VK1E(I,J,K)=AMSUE(I,J,K)*                                       &
     &((AMSUE(I+1,J,K)*(VKE(I+1,J,K)-VKE(I,J,K))/DLXPSI(I,J)            &
     &-AMSUE(I-1,J,K)*(VKE(I,J,K)-VKE(I-1,J,K))/DLXPSI(I-1,J))/DLXV(I,J)&
     &+(-AMSUE(I,J+1,K)*(VKE(I,J,K)-VKE(I,J+1,K))/DLYP(I,J+1))/DLYV(I,J))
      ENDDO
      endif
!
      if(have_g_je)then
      J=JE
      DO I=(I_start+1),IE1
        UK1O(I,J,K)=AMSUO(I,J,K)*                                       &
     &((AMSUO(I+1,J,K)*(UKO(I+1,J,K)-UKO(I,J,K))/DLXP(I+1,J)            &
     &-AMSUO(I-1,J,K)*(UKO(I,J,K)-UKO(I-1,J,K))/DLXP(I,J))/DLXU(I,J)    &
     &+(AMSUO(I,J-1,K)*(UKO(I,J-1,K)-UKO(I,J,K))/DLYPSI(I,J-1))/DLYU(I,J))
      ENDDO
      DO I=2,IE1
        VK1E(I,J,K)=AMSUE(I,J,K)*                                       &
     &((AMSUE(I+1,J,K)*(VKE(I+1,J,K)-VKE(I,J,K))/DLXPSI(I,J)            &
     &-AMSUE(I-1,J,K)*(VKE(I,J,K)-VKE(I-1,J,K))/DLXPSI(I-1,J))/DLXV(I,J)&
     &+(AMSUE(I,J-1,K)*(VKE(I,J-1,K)-VKE(I,J,K))/DLYP(I,J))/DLYV(I,J))
      ENDDO
      endif
!
      ENDDO
!$OMP SINGLE
      CALL bounds_exch(UK1O)
      CALL bounds_exch(VK1E)
!$OMP END SINGLE
!
!$OMP DO
      DO K=1,KE
!
      DO J=2,JE1
      DO I=(I_start+1),IE1
        UOO(I,J,K)=(AMSUO(I,J,K)/(DLXU(I,J)*DLYU(I,J)))                 &
     &            *(DLXU(I,J)*DLYU(I,J)*UKO(I,J,K)                      &
     &             -(AULUX(I+1,J)*DLYP(I+1,J)*AMSUO(I+1,J,K)            &
     &                  *(UK1O(I+1,J,K)-UK1O(I,J,K))/DLXP(I+1,J)        &
     &              -AULUX(I,J)*DLYP(I,J)*AMSUO(I-1,J,K)                &
     &                  *(UK1O(I,J,K)-UK1O(I-1,J,K))/DLXP(I,J)          &
     &              +AULUY(I,J-1)*DLXPSI(I,J-1)*AMSUO(I,J-1,K)          &
     &                  *(UK1O(I,J-1,K)-UK1O(I,J,K))/DLYPSI(I,J-1)      &
     &              -AULUY(I,J)*DLXPSI(I,J)*AMSUO(I,J+1,K)              &
     &                  *(UK1O(I,J,K)-UK1O(I,J+1,K))/DLYPSI(I,J)))
      ENDDO
      ENDDO
!
      DO J=(J_start+1),JE1
      DO I=2,IE1
        VOE(I,J,K)=(AMSUE(I,J,K)/(DLXV(I,J)*DLYV(I,J)))                 &
     &            *(DLXV(I,J)*DLYV(I,J)*VKE(I,J,K)                      &
     &             -(AULVX(I,J)*DLYPSI(I,J)*AMSUE(I+1,J,K)              &
     &                  *(VK1E(I+1,J,K)-VK1E(I,J,K))/DLXPSI(I,J)        &
     &              -AULVX(I-1,J)*DLYPSI(I-1,J)*AMSUE(I-1,J,K)          &
     &                  *(VK1E(I,J,K)-VK1E(I-1,J,K))/DLXPSI(I-1,J)      &
     &              +AULVY(I,J)*DLXP(I,J)*AMSUE(I,J-1,K)                &
     &                  *(VK1E(I,J-1,K)-VK1E(I,J,K))/DLYP(I,J)          &
     &              -AULVY(I,J+1)*DLXP(I,J+1)*AMSUE(I,J+1,K)            &
     &                  *(VK1E(I,J,K)-VK1E(I,J+1,K))/DLYP(I,J+1)))
      ENDDO
      ENDDO
!
      ENDDO
!
#else
!$OMP DO
      DO K=1,KE
!
      DO J=2,JE1
      DO I=(I_start+1),IE1
        UK1O(I,J,K)=AMSUO(I,J,K)*                                       &
     &             (((UKO(I+1,J,K)-UKO(I,J,K))/DLXP(I+1,J)              &
     &              -(UKO(I,J,K)-UKO(I-1,J,K))/DLXP(I,J))/DLXU(I,J)     &
     &             +((UKO(I,J-1,K)-UKO(I,J,K))/DLYPSI(I,J-1)            &
     &              -(UKO(I,J,K)-UKO(I,J+1,K))/DLYPSI(I,J))/DLYU(I,J))
      ENDDO
      ENDDO
!
      DO J=(J_start+1),JE1
      DO I=2,IE1
        VK1E(I,J,K)=AMSUE(I,J,K)*                                       &
     &             (((VKE(I+1,J,K)-VKE(I,J,K))/DLXPSI(I,J)              &
     &              -(VKE(I,J,K)-VKE(I-1,J,K))/DLXPSI(I-1,J))/DLXV(I,J) &
     &             +((VKE(I,J-1,K)-VKE(I,J,K))/DLYP(I,J)                &
     &              -(VKE(I,J,K)-VKE(I,J+1,K))/DLYP(I,J+1))/DLYV(I,J))
      ENDDO
      ENDDO
!
      if(have_g_is)then
      I=0
      DO J=2,JE1
        UK1O(I,J,K)=AMSUO(I,J,K)*                                       &
     &             (((UKO(I+1,J,K)-UKO(I,J,K))/DLXP(I+1,J))/DLXU(I,J)   &
     &             +((UKO(I,J-1,K)-UKO(I,J,K))/DLYPSI(I,J-1)            &
     &              -(UKO(I,J,K)-UKO(I,J+1,K))/DLYPSI(I,J))/DLYU(I,J))
      ENDDO
      I=1
      DO J=(J_start+1),JE1
        VK1E(I,J,K)=AMSUE(I,J,K)*                                       &
     &             (((VKE(I+1,J,K)-VKE(I,J,K))/DLXPSI(I,J))/DLXV(I,J)   &
     &             +((VKE(I,J-1,K)-VKE(I,J,K))/DLYP(I,J)                &
     &              -(VKE(I,J,K)-VKE(I,J+1,K))/DLYP(I,J+1))/DLYV(I,J))
      ENDDO
      endif
!
      if(have_g_ie)then
      I=IE
      DO J=2,JE1
        UK1O(I,J,K)=AMSUO(I,J,K)*                                       &
     &             ((-(UKO(I,J,K)-UKO(I-1,J,K))/DLXP(I,J))/DLXU(I,J)    &
     &             +((UKO(I,J-1,K)-UKO(I,J,K))/DLYPSI(I,J-1)            &
     &              -(UKO(I,J,K)-UKO(I,J+1,K))/DLYPSI(I,J))/DLYU(I,J))
      ENDDO
      DO J=(J_start+1),JE1
        VK1E(I,J,K)=AMSUE(I,J,K)*                                       &
     &             ((-(VKE(I,J,K)-VKE(I-1,J,K))/DLXPSI(I-1,J))/DLXV(I,J)&
     &             +((VKE(I,J-1,K)-VKE(I,J,K))/DLYP(I,J)                &
     &              -(VKE(I,J,K)-VKE(I,J+1,K))/DLYP(I,J+1))/DLYV(I,J))
      ENDDO
      endif
!
      if(have_g_js)then
      J=1
      DO I=(I_start+1),IE1
        UK1O(I,J,K)=AMSUO(I,J,K)*                                       &
     &             (((UKO(I+1,J,K)-UKO(I,J,K))/DLXP(I+1,J)              &
     &              -(UKO(I,J,K)-UKO(I-1,J,K))/DLXP(I,J))/DLXU(I,J)     &
     &             +(-(UKO(I,J,K)-UKO(I,J+1,K))/DLYPSI(I,J))/DLYU(I,J))
      ENDDO
      J=0
      DO I=2,IE1
        VK1E(I,J,K)=AMSUE(I,J,K)*                                       &
     &             (((VKE(I+1,J,K)-VKE(I,J,K))/DLXPSI(I,J)              &
     &              -(VKE(I,J,K)-VKE(I-1,J,K))/DLXPSI(I-1,J))/DLXV(I,J) &
     &             +(-(VKE(I,J,K)-VKE(I,J+1,K))/DLYP(I,J+1))/DLYV(I,J))
      ENDDO
      endif
!
      if(have_g_je)then
      J=JE
      DO I=(I_start+1),IE1
        UK1O(I,J,K)=AMSUO(I,J,K)*                                       &
     &             (((UKO(I+1,J,K)-UKO(I,J,K))/DLXP(I+1,J)              &
     &              -(UKO(I,J,K)-UKO(I-1,J,K))/DLXP(I,J))/DLXU(I,J)     &
     &             +((UKO(I,J-1,K)-UKO(I,J,K))/DLYPSI(I,J-1))/DLYU(I,J))
      ENDDO
      DO I=2,IE1
        VK1E(I,J,K)=AMSUE(I,J,K)*                                       &
     &             (((VKE(I+1,J,K)-VKE(I,J,K))/DLXPSI(I,J)              &
     &              -(VKE(I,J,K)-VKE(I-1,J,K))/DLXPSI(I-1,J))/DLXV(I,J) &
     &             +((VKE(I,J-1,K)-VKE(I,J,K))/DLYP(I,J))/DLYV(I,J))
      ENDDO
      endif
!
      ENDDO
!$OMP SINGLE
      CALL bounds_exch(UK1O)
      CALL bounds_exch(VK1E)
!$OMP END SINGLE
!
!$OMP DO
      DO K=1,KE
!
      DO J=2,JE1
      DO I=(I_start+1),IE1
        UOO(I,J,K)=(AMSUO(I,J,K)/(DLXU(I,J)*DLYU(I,J)))                 &
     &            *(DLXU(I,J)*DLYU(I,J)*UKO(I,J,K)                      &
     &             -(AULUX(I+1,J)*DLYP(I+1,J)                           &
     &                  *(UK1O(I+1,J,K)-UK1O(I,J,K))/DLXP(I+1,J)        &
     &              -AULUX(I,J)*DLYP(I,J)                               &
     &                  *(UK1O(I,J,K)-UK1O(I-1,J,K))/DLXP(I,J)          &
     &              +AULUY(I,J-1)*DLXPSI(I,J-1)                         &
     &                  *(UK1O(I,J-1,K)-UK1O(I,J,K))/DLYPSI(I,J-1)      &
     &              -AULUY(I,J)*DLXPSI(I,J)                             &
     &                  *(UK1O(I,J,K)-UK1O(I,J+1,K))/DLYPSI(I,J)))
      ENDDO
      ENDDO
!
      DO J=(J_start+1),JE1
      DO I=2,IE1
        VOE(I,J,K)=(AMSUE(I,J,K)/(DLXV(I,J)*DLYV(I,J)))                 &
     &            *(DLXV(I,J)*DLYV(I,J)*VKE(I,J,K)                      &
     &             -(AULVX(I,J)*DLYPSI(I,J)                             &
     &                  *(VK1E(I+1,J,K)-VK1E(I,J,K))/DLXPSI(I,J)        &
     &              -AULVX(I-1,J)*DLYPSI(I-1,J)                         &
     &                  *(VK1E(I,J,K)-VK1E(I-1,J,K))/DLXPSI(I-1,J)      &
     &              +AULVY(I,J)*DLXP(I,J)                               &
     &                  *(VK1E(I,J-1,K)-VK1E(I,J,K))/DLYP(I,J)          &
     &              -AULVY(I,J+1)*DLXP(I,J+1)                           &
     &                  *(VK1E(I,J,K)-VK1E(I,J+1,K))/DLYP(I,J+1)))
      ENDDO
      ENDDO
!
      ENDDO
#endif /*FREESLIP*/
!$OMP SINGLE
      CALL bounds_exch(UOO)
      CALL bounds_exch(VOE)
!$OMP END SINGLE
!
#ifdef NON_HYDROSTATIC
#ifdef FREESLIP
      DO K=1,KE
!
      DO J=2,JE1
      DO I=2,IE1
        WN1O(I,J,K)=WETO(I,J,K)*                                        &
     & ((WETO(I+1,J,K)*(WNO(I+1,J,K)-WNO(I,J,K))/DLXU(I,J)              &
     &  -WETO(I-1,J,K)*(WNO(I,J,K)-WNO(I-1,J,K))/DLXU(I-1,J))/DLXP(I,J) &
     & +(WETO(I,J-1,K)*(WNO(I,J-1,K)-WNO(I,J,K))/DLYV(I,J-1)            &
     &  -WETO(I,J+1,K)*(WNO(I,J,K)-WNO(I,J+1,K))/DLYV(I,J))/DLYP(I,J))
      ENDDO
      ENDDO
!
      if(have_g_is)then
      I=1
      DO J=2,JE1
        WN1O(I,J,K)=WETO(I,J,K)*                                        &
     & ((WETO(I+1,J,K)*(WNO(I+1,J,K)-WNO(I,J,K))/DLXU(I,J))/DLXP(I,J)   &
     & +(WETO(I,J-1,K)*(WNO(I,J-1,K)-WNO(I,J,K))/DLYV(I,J-1)            &
     &  -WETO(I,J+1,K)*(WNO(I,J,K)-WNO(I,J+1,K))/DLYV(I,J))/DLYP(I,J))
      ENDDO
      endif
!
      if(have_g_ie)then
      I=IE
      DO J=2,JE1
        WN1O(I,J,K)=WETO(I,J,K)*                                        &
     & ((-WETO(I-1,J,K)*(WNO(I,J,K)-WNO(I-1,J,K))/DLXU(I-1,J))/DLXP(I,J)&
     & +(WETO(I,J-1,K)*(WNO(I,J-1,K)-WNO(I,J,K))/DLYV(I,J-1)            &
     &  -WETO(I,J+1,K)*(WNO(I,J,K)-WNO(I,J+1,K))/DLYV(I,J))/DLYP(I,J))
      ENDDO
      endif
!
      if(have_g_js)then
      J=1
      DO I=2,IE1
        WN1O(I,J,K)=WETO(I,J,K)*                                        &
     & ((WETO(I+1,J,K)*(WNO(I+1,J,K)-WNO(I,J,K))/DLXU(I,J)              &
     &  -WETO(I-1,J,K)*(WNO(I,J,K)-WNO(I-1,J,K))/DLXU(I-1,J))/DLXP(I,J) &
     & +(-WETO(I,J+1,K)*(WNO(I,J,K)-WNO(I,J+1,K))/DLYV(I,J))/DLYP(I,J))
      ENDDO
      endif
!
      if(have_g_je)then
      J=JE
      DO I=2,IE1
        WN1O(I,J,K)=WETO(I,J,K)*                                        &
     & ((WETO(I+1,J,K)*(WNO(I+1,J,K)-WNO(I,J,K))/DLXU(I,J)              &
     &  -WETO(I-1,J,K)*(WNO(I,J,K)-WNO(I-1,J,K))/DLXU(I-1,J))/DLXP(I,J) &
     & +(WETO(I,J-1,K)*(WNO(I,J-1,K)-WNO(I,J,K))/DLYV(I,J-1))/DLYP(I,J))
      ENDDO
      endif
!
      ENDDO
!
      CALL bounds_exch(WN1O)

      DO K=1,KE
      DO J=2,JE1
      DO I=2,IE1
        WNO(I,J,K)=(WETO(I,J,K)/(DLXP(I,J)*DLYP(I,J)*DZP(I,J,K)))       &
     &            *(DLXP(I,J)*DLYP(I,J)*DZP(I,J,K)*WNO(I,J,K)           &
     &            -(+AULWX(I,J)*DLYU(I,J)*DZU(I,J,K)*WETO(I+1,J,K)      &
     &                  *(WN1O(I+1,J,K)-WN1O(I,J,K))/DLXU(I,J)          &
     &              -AULWX(I-1,J)*DLYU(I-1,J)*DZU(I-1,J,K)*WETO(I-1,J,K)&
     &                  *(WN1O(I,J,K)-WN1O(I-1,J,K))/DLXU(I-1,J)        &
     &              +AULWY(I,J-1)*DLXV(I,J-1)*DZV(I,J-1,K)*WETO(I,J-1,K)&
     &                  *(WN1O(I,J-1,K)-WN1O(I,J,K))/DLYV(I,J-1)        &
     &              -AULWY(I,J)*DLXV(I,J)*DZV(I,J,K)*WETO(I,J+1,K)      &
     &                  *(WN1O(I,J,K)-WN1O(I,J+1,K))/DLYV(I,J)))
      ENDDO
      ENDDO
      ENDDO
#else
      DO K=1,KE
!
      DO J=2,JE1
      DO I=2,IE1
        WN1O(I,J,K)=WETO(I,J,K)*                                        &
     &             (((WNO(I+1,J,K)-WNO(I,J,K))/DLXU(I,J)                &
     &              -(WNO(I,J,K)-WNO(I-1,J,K))/DLXU(I-1,J))/DLXP(I,J)   &
     &             +((WNO(I,J-1,K)-WNO(I,J,K))/DLYV(I,J-1)              &
     &              -(WNO(I,J,K)-WNO(I,J+1,K))/DLYV(I,J))/DLYP(I,J))
      ENDDO
      ENDDO
!
      if(have_g_is)then
      I=1
      DO J=2,JE1
        WN1O(I,J,K)=WETO(I,J,K)*                                        &
     &             (((WNO(I+1,J,K)-WNO(I,J,K))/DLXU(I,J))/DLXP(I,J)     &
     &             +((WNO(I,J-1,K)-WNO(I,J,K))/DLYV(I,J-1)              &
     &              -(WNO(I,J,K)-WNO(I,J+1,K))/DLYV(I,J))/DLYP(I,J))
      ENDDO
      endif
!
      if(have_g_ie)then
      I=IE
      DO J=2,JE1
        WN1O(I,J,K)=WETO(I,J,K)*                                        &
     &             ((-(WNO(I,J,K)-WNO(I-1,J,K))/DLXU(I-1,J))/DLXP(I,J)  &
     &             +((WNO(I,J-1,K)-WNO(I,J,K))/DLYV(I,J-1)              &
     &              -(WNO(I,J,K)-WNO(I,J+1,K))/DLYV(I,J))/DLYP(I,J))
      ENDDO
      endif
!
      if(have_g_js)then
      J=1
      DO I=2,IE1
        WN1O(I,J,K)=WETO(I,J,K)*                                        &
     &             (((WNO(I+1,J,K)-WNO(I,J,K))/DLXU(I,J)                &
     &              -(WNO(I,J,K)-WNO(I-1,J,K))/DLXU(I-1,J))/DLXP(I,J)   &
     &             +(-(WNO(I,J,K)-WNO(I,J+1,K))/DLYV(I,J))/DLYP(I,J))
      ENDDO
      endif
!
      if(have_g_je)then
      J=JE
      DO I=2,IE1
        WN1O(I,J,K)=WETO(I,J,K)*                                        &
     &             (((WNO(I+1,J,K)-WNO(I,J,K))/DLXU(I,J)                &
     &              -(WNO(I,J,K)-WNO(I-1,J,K))/DLXU(I-1,J))/DLXP(I,J)   &
     &             +((WNO(I,J-1,K)-WNO(I,J,K))/DLYV(I,J-1))/DLYP(I,J))
      ENDDO
      endif
!
      ENDDO
!
      CALL bounds_exch(WN1O)

      DO K=1,KE
      DO J=2,JE1
      DO I=2,IE1
        WNO(I,J,K)=(WETO(I,J,K)/(DLXP(I,J)*DLYP(I,J)*DZP(I,J,K)))       &
     &            *(DLXP(I,J)*DLYP(I,J)*DZP(I,J,K)*WNO(I,J,K)           &
     &            -(+AULWX(I,J)*DLYU(I,J)*DZU(I,J,K)                    &
     &                  *(WN1O(I+1,J,K)-WN1O(I,J,K))/DLXU(I,J)          &
     &              -AULWX(I-1,J)*DLYU(I-1,J)*DZU(I-1,J,K)              &
     &                  *(WN1O(I,J,K)-WN1O(I-1,J,K))/DLXU(I-1,J)        &
     &              +AULWY(I,J-1)*DLXV(I,J-1)*DZV(I,J-1,K)              &
     &                  *(WN1O(I,J-1,K)-WN1O(I,J,K))/DLYV(I,J-1)        &
     &              -AULWY(I,J)*DLXV(I,J)*DZV(I,J,K)                    &
     &                  *(WN1O(I,J,K)-WN1O(I,J+1,K))/DLYV(I,J)))
      ENDDO
      ENDDO
      ENDDO
#endif /*FREESLIP*/
      CALL bounds_exch(WNO)
#endif
!
!$OMP DO
      CALL OCTIMF
!
!$OMP END PARALLEL
!
      RETURN
      END
