      SUBROUTINE OCTIMF
      USE MO_PARAM1
      USE MO_COMMO1
      USE MO_COMMO2
      DO  K=1,KE
!
      DO  J=1,JE
      DO  I=I_start,IE
       UKO(I,J,K)=UOO(I,J,K)
      ENDDO
      ENDDO
!
      DO  J=J_start,JE
      DO  I=1,IE
       VKE(I,J,K)=VOE(I,J,K)
      ENDDO
      ENDDO
!
      ENDDO

      RETURN
      END
