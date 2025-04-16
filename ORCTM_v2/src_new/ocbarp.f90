      SUBROUTINE OCBARP
      USE MO_PARAM1
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
      IMPLICIT NONE
      INTEGER I,J
!
!
       DO J=1,JE
         DO I=I_start,IE
           UZO(I,J)=USO(I,J)
         ENDDO
       ENDDO

       DO J=J_start,JE
         DO I=1,IE
           VZE(I,J)=VSE(I,J)
         ENDDO
       ENDDO

!
!--------------------------------------------------------------------
!
!     B)
!
!     WINDSTRESS INPUT ON BAROCLINIC VELOCITIES
!
!=====================================================================
!
!     C)
!
       DO J=1,JE
       DO I=I_start,IE
       USO(I,J)=USO(I,J)*AMSUO(I,J,1)
       UZO(I,J)=UZO(I,J)*AMSUO(I,J,1)
       ENDDO
       ENDDO
  
       DO J=J_start,JE
       DO I=1,IE
       VSE(I,J)=VSE(I,J)*AMSUE(I,J,1)
       VZE(I,J)=VZE(I,J)*AMSUE(I,J,1)
       ENDDO
       ENDDO  
!	   
       RETURN
       END
