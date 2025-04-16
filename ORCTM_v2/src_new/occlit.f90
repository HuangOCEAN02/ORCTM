      SUBROUTINE OCCLIT
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
      IMPLICIT NONE

      INTEGER i,j,k,iter,itermax
      REAL BREMS, CONTRA, OVE, SPEED, UNDER, VORW

      REAL BRUVA(IE,JE),CONTRIJ(IE,JE),CONTRIJ_G(IE_G,JE_G)
!
      REAL WPO(IE,JE,KE)
!
!  ITERATION OF BAROCLINIC SYSTEM DIRECTLY IN BAROCLINIC VELOCITIES
!
!$OMP PARALLEL PRIVATE(i,j,k,under,iter,itermax,ove,speed,brems,vorw)
!$OMP DO
       DO J=1,JE
       DO I=1,IE
          BRUVA(I,J)=0.
          STABIO(I,J,1)=0.
          SLOW(I,J)=0.
       ENDDO
       ENDDO
   
       DO J=1,JE
       DO I=I_start,IE
          PXOIN(I,J)=0.
       ENDDO
       ENDDO

       DO J=J_start,JE
       DO I=1,IE
          PYEIN(I,J)=0.
       ENDDO
       ENDDO
!
!$OMP DO
       DO 2 K=1,KE
!
       DO J=1,JE
       DO I=1,IE
          BRUVA(I,J)=BRUVA(I,J)+TIESTW(K+1)*STABIO(I,J,K)
       ENDDO
       ENDDO
!
       DO J=1,JE
       DO I=(I_start+1),IE1
          PXOIN(I,J)=PXOIN(I,J)+DDUO(I,J,K)*DTDXUO(I,J)*                &
     &        (PO(I+1,J,K)-PO(I,J,K))
       ENDDO
       ENDDO
!
       DO J=(J_start+1),JE1
       DO I=1,IE
          PYEIN(I,J)=PYEIN(I,J)+DDUE(I,J,K)*DPYE(I,J)*                  &
     &        (PO(I,J,K)-PO(I,J+1,K))
       ENDDO
       ENDDO
!
2      CONTINUE
!
!$OMP DO
       DO J=1,JE
       DO I=(I_start+1),IE1
          PXOIN(I,J)=PXOIN(I,J)*DEUTIO(I,J)
       ENDDO
       ENDDO
   
       DO J=(J_start+1),JE1
       DO I=1,IE
          PYEIN(I,J)=PYEIN(I,J)*DEUTIE(I,J)
       ENDDO
       ENDDO
!
!$OMP SINGLE
       CALL bounds_exch(PXOIN)
       CALL bounds_exch(PYEIN)
!      No exchange necessary for BRUVA, used only in inner domain
!$OMP END SINGLE

!$OMP DO
       DO K=1,KE
  
       DO J=1,JE
       DO I=I_start,IE
          UKO(I,J,K)=UKO(I,J,K)*DDUO(I,J,K)
          UK1O(I,J,K)=UK1O(I,J,K)*DDUO(I,J,K)
       ENDDO
       ENDDO
  
       DO J=J_start,JE
       DO I=1,IE
          VKE(I,J,K)=VKE(I,J,K)*DDUE(I,J,K)
          VK1E(I,J,K)=VK1E(I,J,K)*DDUE(I,J,K)
       ENDDO
       ENDDO

       ENDDO
!
!$OMP DO
       DO 8015 J=(J_start+1),JE
       DO 8015 I=(I_start+1),IE
       DO 8015 K=KE,1,-1
!
       WO(I,J,K) = WO(I,J,K+1) + DTI * WETO(I,J,K) * (                  &
     &    DTDXPO(I,J) * (   UKO(I-1,J,K)*DLYU(I-1,J)                    &
     &                    - UKO(I,J,K)*DLYU(I,J)    )/DLYP(I,J)         &
     &  + DTDYO(I,J)  * (   VKE(I,J,K)*DLXV(I,J)                        &
     &                    - VKE(I,J-1,K)*DLXV(I,J-1))/DLXP(I,J)  )
!
8015   CONTINUE
!
!$OMP SINGLE
       CALL bounds_exch(WO)
!$OMP END SINGLE
!
!$OMP DO
       DO 8016 J=1,JE
       DO 8016 K=1,KE
       DO 8016 I=1,IE
         T1O(I,J,K)=WO(I,J,K)
         SLOW(I,J)=SLOW(I,J)+STABIO(I,J,K)*DT*WO(I,J,K)
         WPO(I,J,K)=SLOW(I,J)
8016   CONTINUE
!
      UNDER=0.
!OtB Needs more iterations than 6
!     ITERMAX=6
!     ITERMAX=16
!Uwe says 8 is cheaper
      ITERMAX=8 
!
      DO 1000 ITER=1,ITERMAX
!
!$OMP DO
       DO J=1,JE
       DO I=I_start,IE
          UCOS(I,J)=0.
          PXOIN(I,J)=0.
       ENDDO
       ENDDO

       DO J=J_start,JE
       DO I=1,IE
          VCOS(I,J)=0.
          PYEIN(I,J)=0.
       ENDDO
       ENDDO   
!
!$OMP DO
       DO 12 K=1,KE
!
       DO J=1,JE
       DO I=(I_start+1),IE1
       PXOIN(I,J)=PXOIN(I,J)+DDUO(I,J,K)*DTDXUO(I,J)*                   &
     & (WPO(I+1,J,K)-WPO(I,J,K))*AMSUO(I,J,K)
       ENDDO
       ENDDO
!
       DO J=(J_start+1),JE1
       DO I=1,IE
       PYEIN(I,J)=PYEIN(I,J)+DDUE(I,J,K)*DPYE(I,J)*                     &
     & (WPO(I,J,K)-WPO(I,J+1,K))*AMSUE(I,J,K)
       ENDDO
       ENDDO
!
12     CONTINUE
!
!$OMP DO
       DO J=1,JE
       DO I=(I_start+1),IE1
          PXOIN(I,J)=PXOIN(I,J)*DEUTIO(I,J)
       ENDDO
       ENDDO

       DO J=(J_start+1),JE1
       DO I=1,IE
          PYEIN(I,J)=PYEIN(I,J)*DEUTIE(I,J)
       ENDDO
       ENDDO
!
!$OMP SINGLE
       CALL bounds_exch(PXOIN)
       CALL bounds_exch(PYEIN)
!$OMP END SINGLE
!
!$OMP DO
       DO 20 K=1,KE
!
       DO J=2,JE1
       DO I=(I_start+1),IE1
          UKO(I,J,K)=UK1O(I,J,K)+STABN*(PXOIN(I,J) -                    &
     &        DTDXUO(I,J)*(WPO(I+1,J,K)-WPO(I,J,K)))*DDUO(I,J,K)
          UCOS(I,J)=UCOS(I,J)+UKO(I,J,K)
          UKO(I,J,K)=UKO(I,J,K)*AMSUO(I,J,K)
       ENDDO
       ENDDO
!
       DO J=(J_start+1),JE1
       DO I=2,IE1
          VKE(I,J,K)=VK1E(I,J,K)+STABN*(PYEIN(I,J) -                    &
     &        DPYE(I,J)*(WPO(I,J,K)-WPO(I,J+1,K))) *DDUE(I,J,K)
          VCOS(I,J)=VCOS(I,J)+VKE(I,J,K)
          VKE(I,J,K)=VKE(I,J,K)*AMSUE(I,J,K)
       ENDDO
       ENDDO
!
20     CONTINUE
!
!$OMP DO
       DO 2107 K=1,KE
!
       DO J=2,JE1
       DO I=(I_start+1),IE1
          UKO(I,J,K)=UKO(I,J,K)-UCOS(I,J)*DEUTIO(I,J)*DDUO(I,J,K)
       ENDDO
       ENDDO
!
       DO J=(J_start+1),JE1
       DO I=2,IE1
          VKE(I,J,K)=VKE(I,J,K)-VCOS(I,J)*DEUTIE(I,J)*DDUE(I,J,K)
       ENDDO
       ENDDO
!
2107   CONTINUE
!
!$OMP SINGLE
       CALL bounds_exch(UKO)
       CALL bounds_exch(VKE)
!$OMP END SINGLE
!
       IF(ITER.EQ.ITERMAX) GO TO 1000
!
!$OMP DO
       DO 9013 J=1,JE
       DO 9013 I=1,IE
         TLOW(I,J)=0.
         SLOW(I,J)=0.
         CONTRIJ(I,J) = 0
9013   CONTINUE
!
!$OMP DO
       DO 9015 K=KE,1,-1
!
       DO J=(J_start+1),JE
       DO I=(I_start+1),IE
         T1O(I,J,K)= SLOW(I,J) + DTI * WETO(I,J,K) * (                  &
     &      DTDXPO(I,J) * ( UKO(I-1,J,K)*DLYU(I-1,J)                    &
     &                     -UKO(I,J,K)  *DLYU(I,J)  )/DLYP(I,J)         &
     &    + DTDYO(I,J)  * ( VKE(I,J,K)  *DLXV(I,J)                      &
     &                     -VKE(I,J-1,K)*DLXV(I,J-1))/DLXP(I,J)  )
         SLOW(I,J)=T1O(I,J,K)
       ENDDO
       ENDDO
!
9015   CONTINUE
!
!$OMP SINGLE
       CALL bounds_exch(T1O)
!$OMP END SINGLE
!
!$OMP DO
       DO 30 K=2,KE
       DO 30 J=1,JE
       DO 30 I=1,IE
         OVE=WPO(I,J,K)
         TLOW(I,J)=TLOW(I,J)+G*STABIO(I,J,K)*DT*(                       &
     &       CONN*T1O(I,J,K)+WO(I,J,K))
         SPEED=2.*G*BRUVA(I,J)*(DTDXUO(I,J)**2+DTDYO(I,J)**2)
         BREMS=            CONN*STABN/(1.+SPEED)
         BREMS=UNDER*SPEED*CONN*STABN/(1.+SPEED)
         VORW=1.-BREMS
         WPO(I,J,K)=VORW*TLOW(I,J)+BREMS*WPO(I,J,K)
         CONTRIJ(I,J)=CONTRIJ(I,J)+(OVE-WPO(I,J,K))**2
30     CONTINUE
!
!$OMP SINGLE
       CALL bounds_exch(WPO)
!$OMP END SINGLE
!
       UNDER=1.
!
1000  CONTINUE
!
!$OMP DO
       DO 8881 K=1,KE
!
       DO J=1,JE
       DO I=I_start,IE
         uko(i,j,k)=uko(i,j,k)*amsuo(i,j,k)/(ALMZER+dduo(I,j,k))
         uk1o(i,j,k)=uk1o(i,j,k)*amsuo(i,j,k)/(ALMZER+dduo(i,j,k))
       ENDDO
       ENDDO
!
       DO J=J_start,JE
       DO I=1,IE
         vke(i,j,k)=vke(i,j,k)*amsue(i,j,k)/(ALMZER+ddue(i,j,k))
         vk1e(i,j,k)=vk1e(i,j,k)*amsue(i,j,k)/(ALMZER+ddue(i,j,k))
       ENDDO
       ENDDO
!
8881   CONTINUE

!$OMP END PARALLEL
!
      CALL gather_arr(contrij,contrij_g,p_io)

      IF(p_pe==p_io) THEN
        CONTRA=0.
        DO  J=1,JE_G
        DO  I=1,IE_G
          CONTRA=CONTRA+CONTRIJ_G(I,J)
        ENDDO
        ENDDO
!
        WRITE(IO_STDOUT,*)'CONTRA: ',CONTRA
!
        IF(CONTRA.GT.1.E-3) THEN
          WRITE(IO_STDOUT,*) 'CONTRI : ',SUM(CONTRIJ_G(1:IE_G,1:JE_G),1)
          WRITE(IO_STDOUT,*) ' '
          WRITE(IO_STDOUT,*) 'CONTRJ : ',SUM(CONTRIJ_G(1:IE_G,1:JE_G),2)
          CALL ABSTURZ
        ENDIF
      ENDIF
      ! test for NaN in occlit.f90
      !CALL ABSTURZ
     
      RETURN
      END
