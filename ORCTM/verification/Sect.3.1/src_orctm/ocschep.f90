      SUBROUTINE OCSCHEP
!
!     HORIZONTAL DIFFUSION OF MOMENTUM
!     
!     BY ERNST MAIER-REIMER
!     MODIFIED BY UWE MIKOLAJEWICZ 2/99
!           MAKE AUSF 2D, ADAPT TO SPATIALLY VARIABLE GRID
!           AUS TIME CONSTANT
!
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
      IMPLICIT NONE

      INTEGER i,j,k
      REAL P11(IE,JE),P22(IE,JE),P12(IE,JE),AUSF(IE,JE),AUSFPSI(IE,JE)
	  REAL DLXPSI(IE,JE),DLYPSI(IE,JE)
      REAL DLXXYY, DLXXYYPSI, UX, UY, VX, VY

!
!UWE  AUS DIMENSIONLESS MOMENTUM DIFFUSION COEFFICIENT
!
!     AUSF     2D-MOMENTUM DIFFUSION COEFFICENT P-POINTS      [M2/S]
!     AUSFPSI  THE SAME FOR PSI-POINTS
!
!     SPATIAL DEPENDENCE ACCORDING TO BRYAN ET AL. 1975
!
      DO J=1,JE
       DO I=1,IE1
        DLYPSI(I,J)=0.5*(DLYV(I,J)+DLYV(I+1,J))
       ENDDO
       DLYPSI(IE,J)=DLYV(IE,J)
      ENDDO
  
      DO I=1,IE
       DO J=1,JE1
        DLXPSI(I,J)=0.5*(DLXU(I,J)+DLXU(I,J+1))
       ENDDO
       DLXPSI(I,JE)=DLXU(I,JE)
      ENDDO
  
      DO 2712 J=1,JE
      DO 2712 I=1,IE
          P11(I,J)=0.
          P12(I,J)=0.
          P22(I,J)=0.
          AUSFPSI(I,J)=0.
          DLXXYY=MAX(DLXP(I,J),DLYP(I,J))
          AUSF(I,J)=AUS*DLXXYY**2
          DLXXYYPSI=MAX(DLXPSI(I,J),DLYPSI(I,J))
          AUSFPSI(I,J)=AUS*DLXXYYPSI**2
2712  CONTINUE
!
      DO 1 K=1,KE
!
      DO J=1,JE
      DO I=I_start,IE
      UK1O(I,J,K)=UKO(I,J,K)*AMSUO(I,J,K)
      ENDDO
      ENDDO
!
      DO J=J_start,JE
      DO I=1,IE
      VK1E(I,J,K)=VKE(I,J,K)*AMSUE(I,J,K)
      ENDDO
      ENDDO
!
1     CONTINUE
!
      CALL bounds_exch(UK1O)
      CALL bounds_exch(VK1E)
!
      DO 8921 K=1,KE

      DO J=(J_start+1),JE
      DO I=(I_start+1),IE
      UX=(DLYU(I,J)*DDUO(I,J,K)*UKO(I,J,K)                              &
     &     -DLYU(I-1,J)*DDUO(I-1,J,K)*UKO(I-1,J,K))                     &
     &      *DPIO(I,J,K)/(DLXP(I,J)*DLYP(I,J))
      VY=(DDUE(I,J-1,K)*DLXV(I,J-1)*VKE(I,J-1,K)                        &
     &    -DDUE(I,J,K)*DLXV(I,J)*VKE(I,J,K))                            &
     &     *DPIO(I,J,K)/(DLXP(I,J)*DLYP(I,J)) 
      IF(I.EQ.IE)THEN
        VX=(VK1E(I,J,K)-VK1E(I-1,J,K))*2./(DLXU(I-1,J)+DLXU(I-1,J+1))
      ELSE
        VX=(VK1E(I+1,J,K)-VK1E(I,J,K))*2./(DLXU(I,J)+DLXU(I,J+1))
      ENDIF
      IF(J.EQ.JE)THEN
        UY=(UK1O(I,J-1,K)-UK1O(I,J,K))*2./(DLYV(I,J-1)+DLYV(I+1,J-1))
      ELSE
          UY=(UK1O(I,J,K)-UK1O(I,J+1,K))*2./(DLYV(I,J)+DLYV(I+1,J))
      ENDIF
      P11(I,J)=AUSF(I,J)*DDPO(I,J,K)*(UX-VY)
      P22(I,J)=-P11(I,J)
      P12(I,J)=AUSFPSI(I,J)*DDPSIO(I,J,K)*(UY+VX)
      ENDDO
      ENDDO
!
      CALL bounds_exch(P11,P22,P12)
!
      DO J=2,JE1
      DO I=(I_start+1),IE1
      UKO(I,J,K)=UOO(I,J,K)+                                            &
     & (amsuo(i,j,k)/(almzer+dduo(i,j,k)))                              &
     &     *(dtdxuo(i,j)*(p11(i+1,j)-p11(i,j))                          &
     & +(dpyo(i,j-1)*p12(i,j-1)-dpyo(i,j)*p12(i,j)))
      ENDDO
      ENDDO
!
      DO J=(J_start+1),JE1
      DO I=2,IE1
      VKE(I,J,K)=VOE(I,J,K)                                             &
     & + (amsue(i,j,k)/(almzer+ddue(i,j,k)))                            &
     &           *((dtdxue(i,j)*p12(i,j)-dtdxue(i-1,j)*p12(i-1,j))      &
     & + dpye(i,j)*(p22(i,j)-p22(i,j+1)))
      ENDDO
      ENDDO
!
      DO J=2,JE1
      DO I=(I_start+1),IE1
      UOO(I,J,K)=UKO(I,J,K)
      ENDDO
      ENDDO
!
      DO J=(J_start+1),JE1
      DO I=2,IE1
      VOE(I,J,K)=VKE(I,J,K)
      ENDDO
      ENDDO
 
8921  CONTINUE

      CALL bounds_exch(UKO)
      CALL bounds_exch(VKE)
      CALL bounds_exch(UOO)
      CALL bounds_exch(VOE)

      RETURN
      END
