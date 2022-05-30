      SUBROUTINE OCJITR
#ifdef GMBOLUS
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_COMMOAU1
      USE MO_UNITS
!
!     PARAMETRIZATION OF SUBGRID EDDY EFFECTS ACCORDING GENT AND 
!     JIM MCWILLIAMS, 
!      GENT ET AL. 1995, JPO 25,463.
!UWE
!     OUTLINE BY E. MAIER-REIMER
!     IMPLEMENTED BY U. MIKOLAJEWICZ 4/99
!     MODIFIED 6/2000 INCLUDE SEA LEVEL
!
!     BOLX,BOLY  GM COEFFICIENT K
!     BOUNDARY CONDITION IMPLEMENTED IN DETERMINATION OF RHO-GRADIENTS
!     FINAL ADVECTION WITH AN UPWIND SCHEME
!
      implicit none

      real ZSURGM
      integer i,j,k,l
      DIMENSION ZSURGM(KE)
      real tsup_p(ie,je), ssup_p(ie,je), tlow_p(ie,je), slow_p(ie,je)
      real roxo,roxu,royo,royu,rozxo,rozxu,rozyo,rozyu
      real tm,sm,wun,wob,uwe,uos,vsu,vno,dhi
      real bolk,valinf,stabmin
!
!UWE VALUE FOR INFINITY
      VALINF=1.E30
      STABMIN=1.E-3
      DO K=1,KE
       ZSURGM(K)=0.
      ENDDO
      ZSURGM(1)=1.
!
!
#ifdef GMVISETAL
!JJ USE VISBECK ET AL. (JPO, 1997) FORMULATION TO CALCULATE
!   EDDY TRANSFER COEFFICIENTS
!   BOLX, BOLY ARE CALCULATED IN SBR OCTHER
!JJ
#else
!UWE  MAKE BOLUS COEFFICIENT A LINEAR FUNCTION OF DX
!     BOLX IN M**2/S, BOLK IN M/S
!HH   SET GM DIFFUSION TO HARMOIC DIFFCOEF FROM NAMELIST
      BOLK=AH00
#ifdef BOLK025
      BOLK=BOLK*0.25
#endif /*BOLK025*/
#ifdef BOLK05
      BOLK=BOLK*0.5
#endif /*BOLK05*/
      DO J=1,JE
       DO I=1,IE
        BOLX(I,J)=MIN(BOLK*DLXU(I,J),2000.)
        BOLY(I,J)=MIN(BOLK*DLYV(I,J),2000.)
       ENDDO
      ENDDO
!
#endif /*GMVISETAL*/

!$OMP PARALLEL PRIVATE(i,j,k,l,roxo,roxu,royo,royu,rozxo,rozxu,rozyo,rozyu, &
!$OMP    tsup_p,ssup_p,tlow_p,slow_p,tm,sm,wun,wob,uwe,uos,vsu,vno,dhi)

!$OMP DO
      DO K=1,KE
        UPTRT(K)=0.
        DOTRT(K)=0.
        DO J=1,JE
          DO I=1,IE
            WGO(I,J,K)=0.
          ENDDO
        ENDDO
        DO J=1,JE
          DO I=I_start,IE
            UK1O(I,J,K)=0.
          ENDDO
        ENDDO
        DO J=J_start,JE
          DO I=1,IE
            VK1E(I,J,K)=0.
          ENDDO
        ENDDO
      ENDDO
!
!$OMP DO
      DO K=2,KE-1
!
      DO J=1,JE
      DO I=(I_start+1),IE1
      ROXO=0.5*(RHOO(I+1,J,K)-RHOO(I,J,K)+RHOO(I+1,J,K-1)               &
     & -RHOO(I,J,K-1))/DLXP(I,J)
      ROXU=0.5*(RHOO(I+1,J,K)-RHOO(I,J,K)+(RHOO(I+1,J,K+1)              &
     & -RHOO(I,J,K+1))*AMSUO(I,J,K+1))/DLXP(I,J)
!
      ROZXO=0.5*(STABIO(I+1,J,K)+STABIO(I,J,K))
      ROZXU=0.5*(STABIO(I+1,J,K+1)+STABIO(I,J,K+1))
      ROZXO=MAX(ROZXO,STABMIN)
      ROZXU=MAX(ROZXU,STABMIN)
      IF(AMSUO(I,J,K+1).LT.0.5)ROZXU=VALINF
!
      UK1O(I,J,K)=-BOLX(I,J)*(ROXO/ROZXO-ROXU/ROZXU)*AMSUO(I,J,K)       &
     &  /(DDUO(I,J,K)+(1.-AMSUO(I,J,K)))
      ENDDO
      ENDDO
!
      DO J=(J_start+1),JE1
      DO I=1,IE
      ROYO=0.5*(RHOO(I,J,K)-RHOO(I,J+1,K)+RHOO(I,J,K-1)                 &
     & -RHOO(I,J+1,K-1))/DLYP(I,J)
      ROYU=0.5*(RHOO(I,J,K)-RHOO(I,J+1,K)+(RHOO(I,J,K+1)                &
     & -RHOO(I,J+1,K+1))*AMSUE(I,J,K+1))/DLYP(I,J)
!
      ROZYO=0.5*(STABIO(I,J+1,K)+STABIO(I,J,K))
      ROZYU=0.5*(STABIO(I,J+1,K+1)+STABIO(I,J,K+1))
      ROZYO=MAX(ROZYO,STABMIN)
      ROZYU=MAX(ROZYU,STABMIN)
      IF(AMSUE(I,J,K+1).LT.0.5)ROZYU=VALINF
!
      VK1E(I,J,K)=-BOLY(I,J)*(ROYO/ROZYO-ROYU/ROZYU)*AMSUE(I,J,K)       &
     &  /(DDUE(I,J,K)+(1.-AMSUE(I,J,K)))
      ENDDO
      ENDDO

      ENDDO
!
!UWE INCLUDE SURFACE AND BOTTOM LAYERS
!
!     SURFACE LAYER
!
      K=1
!$OMP DO
      DO J=1,JE
      DO I=(I_start+1),IE1
      ROXO=0.
      ROXU=0.5*(RHOO(I+1,J,K)-RHOO(I,J,K)+(RHOO(I+1,J,K+1)              &
     & -RHOO(I,J,K+1))*AMSUO(I,J,K+1))/DLXP(I,J)
      ROZXO=VALINF
      ROZXU=0.5*(STABIO(I+1,J,K+1)+STABIO(I,J,K+1))
      ROZXU=MAX(ROZXU,STABMIN)
      UK1O(I,J,K)=-BOLX(I,J)*DWI(K)*(ROXO/ROZXO-ROXU/ROZXU)*AMSUO(I,J,K)
      ENDDO
      ENDDO

      DO J=(J_start+1),JE1
      DO I=1,IE
      ROYO=0.
      ROYU=0.5*(RHOO(I,J,K)-RHOO(I,J+1,K)+(RHOO(I,J,K+1)                &
     & -RHOO(I,J+1,K+1))*AMSUE(I,J,K+1))/DLYP(I,J)
      ROZYO=VALINF
      ROZYU=0.5*(STABIO(I,J+1,K+1)+STABIO(I,J,K+1))
      ROZYU=MAX(ROZYU,STABMIN)
      VK1E(I,J,K)=-BOLY(I,J)*DWI(K)*(ROYO/ROZYO-ROYU/ROZYU)*AMSUE(I,J,K)
      ENDDO
      ENDDO
!
!     BOTTOM LAYER
!
      K=KE
!$OMP DO
      DO J=1,JE
      DO I=(I_start+1),IE1
      ROXO=0.5*(RHOO(I+1,J,K)-RHOO(I,J,K)+RHOO(I+1,J,K-1)               &
     & -RHOO(I,J,K-1))/DLXP(I,J)
      ROXU=0.
      ROZXO=0.5*(STABIO(I+1,J,K)+STABIO(I,J,K))
      ROZXU=VALINF
      ROZXO=MAX(ROZXO,STABMIN)
      UK1O(I,J,K)=-BOLX(I,J)*(ROXO/ROZXO-ROXU/ROZXU)*AMSUO(I,J,K)       &
     &  /(DDUO(I,J,K)+(1.-AMSUO(I,J,K)))
      ENDDO
      ENDDO

      DO J=(J_start+1),JE1
      DO I=1,IE
      ROYO=0.5*(RHOO(I,J,K)-RHOO(I,J+1,K)+RHOO(I,J,K-1)                 &
     & -RHOO(I,J+1,K-1))/DLYP(I,J)
      ROYU=0.
      ROZYO=0.5*(STABIO(I,J+1,K)+STABIO(I,J,K))
      ROZYU=VALINF
      ROZYO=MAX(ROZYO,STABMIN)
      VK1E(I,J,K)=-BOLY(I,J)*(ROYO/ROZYO-ROYU/ROZYU)*AMSUE(I,J,K)       &
     &  /(DDUE(I,J,K)+(1.-AMSUE(I,J,K)))
      ENDDO
      ENDDO

      DO K=1,KE
        DO J=1,JE
          UK1O(I_start,J,K)=UK1O(I_start+1,J,K)
          UK1O(IE,J,K)=UK1O(IE-1,J,K)
        ENDDO
        DO I=1,IE
          VK1E(I,J_start,K)=VK1E(I,J_start+1,K)
          VK1E(I,JE,K)=VK1E(I,JE-1,K)
        ENDDO
      ENDDO

!$OMP SINGLE
      CALL bounds_exch(UK1O)
      CALL bounds_exch(VK1E)
!$OMP END SINGLE
!
!     VERTICAL VELOCITY = VERTICAL INTEGRAL OF DIVERGENCE OF
!                             HORIZONTAL VELOCITY FIELD
!                         IN OCJITR
!
!$OMP DO
      DO 7001 J=1,JE
      DO 7001 I=1,IE
      WGO(I,J,KEP) = ZERO
 7001 CONTINUE
!
!$OMP DO
      DO J=(J_start+1),JE
      DO I=(I_start+1),IE
      DO K=KE,1,-1
      WGO(I,J,K) = WGO(I,J,K+1)                                         &
     & + DTI*WETO(I,J,K)*(                                              &
     & DTDXPO(I,J)   * (   UK1O(I-1,J,K) * DDUO(I-1,J,K)*DLYU(I-1,J)    &
     &  - UK1O(I,J,K)   * DDUO(I,J,K)*DLYU(I,J)   )/DLYP(I,J)  )        &
     & + DTI*WETO(I,J,K)* (                                             &
     & + DTDYO(I,J) *                                                   &
     & (   VK1E(I,J,K)     * DDUE(I,J,K)*DLXV(I,J)                      &
     & - VK1E(I,J-1,K)   * DDUE(I,J-1,K)*DLXV(I,J-1)  )/DLXP(I,J)  )
        ENDDO
        ENDDO
        ENDDO
!
!$OMP   DO
        DO K=1,KE
          UPTRT(K)=0.
          DOTRT(K)=0.
        DO J=2,JE1
        DO I=2,IE1
          UPTRT(K)=UPTRT(K)+DLXP(I,J)*DLYP(I,J)*                        &
     &                               (WGO(I,J,K)+ABS(WGO(I,J,K)))
          DOTRT(K)=DOTRT(K)-DLXP(I,J)*DLYP(I,J)*                        &
     &                               (WGO(I,J,K)-ABS(WGO(I,J,K)))
        ENDDO
        ENDDO
        UPTRT(K)=UPTRT(K)*1.E-6
        DOTRT(K)=DOTRT(K)*1.E-6
      ENDDO
!$OMP SINGLE
!
      CALL bounds_exch(WGO)
      CALL global_sum(UPTRT)
      CALL global_sum(DOTRT)
!
6626   FORMAT(10F8.2)
!        WRITE(IO_STDOUT,*   ) 'Vertical uptransport by layer'
!        WRITE(IO_STDOUT,6626)UPTRT
!$OMP END SINGLE
!     
!$OMP DO
      DO 175 K=1,KE
      IF(K.EQ.1) THEN
      DO 1751 J=1,JE
      DO 1751 I=1,IE
       TSUP_P(I,J)=THO(I,J,1)
       SSUP_P(I,J)=SAO(I,J,1)
1751  CONTINUE
      ELSE
      DO 1752 J=1,JE
      DO 1752 I=1,IE
       TSUP_P(I,J)=THO(I,J,K-1)
       SSUP_P(I,J)=SAO(I,J,K-1)
1752  CONTINUE
      ENDIF
      IF(K.EQ.KE) THEN
      DO 1753 J=1,JE
      DO 1753 I=1,IE
       TLOW_P(I,J)=THO(I,J,KE)
       SLOW_P(I,J)=SAO(I,J,KE)
1753  CONTINUE
      ELSE
      DO 1754 J=1,JE
      DO 1754 I=1,IE
       TLOW_P(I,J)=THO(I,J,K+1)
       SLOW_P(I,J)=SAO(I,J,K+1)
1754  CONTINUE
      ENDIF
!
!UWE  MAKE ADVECTION MASS CONSERVING
!
      DO 175 J=2,JE1
      DO 175 I=2,IE1
      TM=THO(I,J,K)
      SM=SAO(I,J,K)
      WUN=HALF*(WGO(I,J,K+1)+ABS(WGO(I,J,K+1)))*DLXP(I,J)*DLYP(I,J)
      WOB=HALF*(ABS(WGO(I,J,K))-WGO(I,J,K))*DLXP(I,J)*DLYP(I,J)
      UWE=HALF*(UK1O(I-1,J,K)+ABS(UK1O(I-1,J,K)))                       &
     &          *DLYU(I-1,J)*DDUO(I-1,J,K)
      UOS=HALF*(ABS(UK1O(I,J,K))-UK1O(I,J,K))                           &
     &          *DLYU(I,J)*DDUO(I,J,K)
      VSU=HALF*(VK1E(I,J,K)+ABS(VK1E(I,J,K)))                           &
     &          *DLXV(I,J)*DDUE(I,J,K)
      VNO=HALF*(ABS(VK1E(I,J-1,K))-VK1E(I,J-1,K))                       &
     &          *DLXV(I,J-1)*DDUE(I,J-1,K)
      T1O(I,J,K)=(TM*DLXP(I,J)*DLYP(I,J)*(DDPO(I,J,K)+ALMZER            &
     &   +ZSURGM(K)*(ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA))  &
     &      +DT*                                                        &
     &       (WOB*(TSUP_P(I,J)-TM)+WUN*(TLOW_P(I,J)-TM)                 &
     &       +UWE*(THO(I-1,J,K)-TM)+UOS*(THO(I+1,J,K)-TM)               &
     &       +VNO*(THO(I,J-1,K)-TM)+VSU*(THO(I,J+1,K)-TM)))             &
     &/(DLXP(I,J)*DLYP(I,J)*(DDPO(I,J,K)+ALMZER                         &
     &   +ZSURGM(K)*(ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)))
!
      S1O(I,J,K)=(SM*DLXP(I,J)*DLYP(I,J)*(DDPO(I,J,K)+ALMZER            &
     &    +ZSURGM(K)*(ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)) &
     &       +DT*                                                       &
     &       (WOB*(SSUP_P(I,J)-SM)+WUN*(SLOW_P(I,J)-SM)                 &
     &       +UWE*(SAO(I-1,J,K)-SM)+UOS*(SAO(I+1,J,K)-SM)               &
     &       +VNO*(SAO(I,J-1,K)-SM)+VSU*(SAO(I,J+1,K)-SM)))             &
     &/(DLXP(I,J)*DLYP(I,J)*(DDPO(I,J,K)+ALMZER                         &
     &   +ZSURGM(K)*(ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)))
!
175    CONTINUE
!
!$OMP DO
      DO 76 K=1,KE
      DO 76 J=2,JE1
      DO 76 I=2,IE1
      THO(I,J,K)=T1O(I,J,K)
      SAO(I,J,K)=S1O(I,J,K)
76    CONTINUE
!$OMP SINGLE
      CALL bounds_exch(THO)
      CALL bounds_exch(SAO)
!$OMP END SINGLE

#endif /*GMBOLUS*/
!$OMP END PARALLEL
      RETURN
      END
