      SUBROUTINE OCMODMOM
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
      IMPLICIT NONE
!
!-=====================================================================
!
!    DECOMPOSITION INTO BAROTROPIC AND BAROCLINIC FIELD
!
      REAL UCOR(I_start:IE,1:JE,1:KE),VCOR(1:IE,J_start:JE,1:KE)
      REAL WCOR(IE,JE,KEP)
      REAL UCOR_V,UCOR_W,VCOR_U,VCOR_W,WCOR_U,WCOR_V
      REAL dlon,dlat,dist,FTWOU1,FTWOV1,FTWOW1
      INTEGER i,j,k,iter
!
!---------------------------------------------------------------------
!
!$OMP PARALLEL PRIVATE(i,j,k,iter)

!$OMP DO
      DO J=1,JE
      DO I=I_start,IE
      U1O(I,J)=ZERO
      USO(I,J)=ZERO
      ENDDO
      ENDDO
  
      DO J=J_start,JE
      DO I=1,IE
      V1E(I,J)=ZERO
      VSE(I,J)=ZERO
      ENDDO
      ENDDO
!
!$OMP DO
      DO 822 K=1,KE
  
      DO J=1,JE
      DO I=1,IE
        STABIO(I,J,K)=0.001*DZ(K)*STABIO(I,J,K)*WETO(I,J,K)
        PO(I,J,K)=PO(I,J,K)+G*ZO(I,J)*rlsa1
      ENDDO
      ENDDO
  
      DO J=1,JE
      DO I=I_start,IE
        UKO(I,J,K)=UKO(I,J,K)*DDUO(I,J,K)
        UK1O(I,J,K)=UKO(I,J,K)
      ENDDO
      ENDDO
  
      DO J=J_start,JE
      DO I=1,IE
        VKE(I,J,K)=VKE(I,J,K)*DDUE(I,J,K)
        VK1E(I,J,K)=VKE(I,J,K)
      ENDDO
      ENDDO
  
822   CONTINUE
!
#ifdef NON_HYDROSTATIC
      DO K=1,KEP
       DO J=1,JE
        DO I=1,IE
         WN1O(I,J,K)=WNO(I,J,K)
        ENDDO
       ENDDO
      ENDDO
#endif
!
      DO 820 ITER=1,12
!$OMP DO
      DO K=1,KE
        DO J=1,JE
          DO I=I_start,IE
            UCOR(I,J,K)=STABN*UK1O(I,J,K)+STABO*UKO(I,J,K)
          ENDDO
        ENDDO
        DO J=J_start,JE
          DO I=1,IE
            VCOR(I,J,K)=STABN*VK1E(I,J,K)+STABO*VKE(I,J,K)
          ENDDO
        ENDDO
      ENDDO
#ifdef NON_HYDROSTATIC
      DO K=1,KEP
        DO J=1,JE
          DO I=1,IE
            WCOR(I,J,K)=STABN*WN1O(I,J,K)+STABO*WNO(I,J,K)
          ENDDO
        ENDDO
      ENDDO
#endif
!
!$OMP DO
      DO K=1,KE
  
      DO J=2,JE1
      DO I=(I_start+1),IE1
        VCOR_U=0.25 * (  VCOR(I,J,K)   + VCOR(I+1,J,K)                  &
     &                 + VCOR(I,J-1,K) + VCOR(I+1,J-1,K) )
        UK1O(I,J,K)=UKO(I,J,K)                                          &
     &      +DTDXUO(I,J)*(PO(I,J,K)-PO(I+1,J,K))*DDUO(I,J,K)            &
     &      +DT*FTWOU(I,J)*VCOR_U*AMSUO(I,J,K)
      ENDDO
      ENDDO

      DO J=(J_start+1),JE1
      DO I=2,IE1
        UCOR_V=0.25 * (  UCOR(I,J,K)  + UCOR(I-1,J,K)                   &
     &                 + UCOR(I,J+1,K)+ UCOR(I-1,J+1,K) )
        VK1E(I,J,K)=VKE(I,J,K)                                          &
     &      +DPYE(I,J)*(PO(I,J+1,K)-PO(I,J,K))*DDUE(I,J,K)              &
     &      -DT*FTWOV(I,J)*UCOR_V*AMSUE(I,J,K)
      ENDDO
      ENDDO

#ifdef NON_HYDROSTATIC 
      DO J=2,JE1
      DO I=(I_start+1),IE1
! cos(theta(i,j))
       dlon=(gila(2*I+2,2*J)-gila(2*I+1,2*J))*cos(gila(2*I+1,2*J))
       dlat=giph(2*I+2,2*J)-giph(2*I+1,2*J)
       dist=sqrt(dlon*dlon+dlat*dlat)+almzer
       FTWOU1=2*omega*cos(giph(2*I+1,2*J))

       WCOR_U=0.25*((WCOR(I,J,K)+WCOR(I,J,K+1))*DLXP(I,J)+              &
     &        (WCOR(I+1,J,K)+WCOR(I+1,J,K+1))*DLXP(I+1,J))/DLXU(I,J)
       UK1O(I,J,K)=UK1O(I,J,K)-DT*FTWOU1*WCOR_U*DDUO(I,J,K)*dlon/dist
      ENDDO
      ENDDO

      DO J=(J_start+1),JE1
      DO I=2,IE1
! -sin(theta(i,j))
       dlon=(gila(2*I,2*J)-gila(2*I,2*J+1))*cos(gila(2*I,2*J+1))
       dlat=giph(2*I,2*J)-giph(2*I,2*J+1)
       dist=sqrt(dlon*dlon+dlat*dlat)+almzer
       FTWOV1=2*omega*cos(giph(2*I,2*J+1))

       WCOR_V=0.25*((WCOR(I,J,K)+WCOR(I,J,K+1))*DLYP(I,J)+              &
     &        (WCOR(I,J+1,K)+WCOR(I,J+1,K+1))*DLYP(I,J+1))/DLYV(I,J)
       VK1E(I,J,K)=VK1E(I,J,K)-DT*FTWOV1*WCOR_V*DDUE(I,J,K)*dlon/dist
      ENDDO
      ENDDO

      IF (K .GT. 1) THEN
        DO J=2,JE1
        DO I=2,IE1
! u*cos(theta(i,j))-v*sin(theta(i,j))
          dlon=(gila(2*I+1,2*J)-gila(2*I,2*J))*cos(gila(2*I,2*J))
          dlat=giph(2*I,2*J-1)-giph(2*I,2*J)
          dist=sqrt(dlon*dlon+dlat*dlat)+almzer
          FTWOW1=2*omega*cos(giph(2*I,2*J))

          UCOR_W=0.25*( (UKO(I,J,K)+UKO(I-1,J,K))/(DDPO(I,J,K)+ALMZER)  &
     &           +(UKO(I,J,K-1)+UKO(I-1,J,K-1))/(DDPO(I,J,K-1)+ALMZER) )
          VCOR_W=0.25*( (VKE(I,J,K)+VKE(I,J-1,K))/(DDPO(I,J,K)+ALMZER)  &
     &           +(VKE(I,J,K-1)+VKE(I,J-1,K-1))/(DDPO(I,J,K-1)+ALMZER) )
          WN1O(I,J,K)=WNO(I,J,K)+DT*WETO(I,J,K)*                        &
     &           ( FTWOW1*dlon/dist*UCOR_W - FTWOW1*dlat/dist*VCOR_W )
        ENDDO
        ENDDO
      ENDIF
#endif

      ENDDO
!
!$OMP SINGLE
      CALL bounds_exch(UK1O)
      CALL bounds_exch(VK1E)      
#ifdef NON_HYDROSTATIC
      CALL bounds_exch(WN1O)
      DO K=1,KEP
       DO J=1,JE
        DO I=1,IE
         WNO(I,J,K)=WN1O(I,J,K)
        ENDDO
       ENDDO
      ENDDO
#endif
!$OMP END SINGLE
820   CONTINUE
!
!     CALCULATION OF BAROTROPIC VELOCITIES ON FIELDS U1 AND V1
!
!$OMP DO
! 

      DO 8883 K=1,KE
!
      DO J=1,JE
      DO I=1,IE
      PO(I,J,K)=PO(I,J,K)-G*ZO(I,J)*rlsa1
      ENDDO
      ENDDO
!
      DO J=1,JE
      DO I=I_start,IE
      U1O(I,J) = U1O(I,J)  +   DDUO(I,J,K) * UOO(I,J,K)
      USO(I,J) = USO(I,J)  +   UK1O(I,J,K)
      ENDDO
      ENDDO
!
      DO J=J_start,JE
      DO I=1,IE
      V1E(I,J) = V1E(I,J)  +   DDUE(I,J,K) * VOE(I,J,K)
      VSE(I,J) = VSE(I,J)  +   VK1E(I,J,K)
      ENDDO
      ENDDO
!
8883  CONTINUE


!$OMP DO
      DO 120 K=1,KE
!
      DO J=1,JE
      DO I=I_start,IE
      UKO(I,J,K)=UOO(I,J,K)-AMSUO(I,J,K)*DEUTIO(I,J)*U1O(I,J)
      UK1O(I,J,K)=UK1O(I,J,K)*(amsuo(i,j,k)/(almzer+dduo(i,j,k)))       &
     & -AMSUO(I,J,K)*DEUTIO(I,J)*USO(I,J)
      ENDDO
      ENDDO
!
      DO J=J_start,JE
      DO I=1,IE
      VKE(I,J,K)=VOE(I,J,K)-AMSUE(I,J,K)*DEUTIE(I,J)*V1E(I,J)
      VK1E(I,J,K)=VK1E(I,J,K)*(amsue(i,j,k)/(almzer+ddue(i,j,k)))       &
     & -AMSUE(I,J,K)*DEUTIE(I,J)*VSE(I,J)
      ENDDO
      ENDDO
!
120   CONTINUE

     
!$OMP END PARALLEL

      RETURN
      END
