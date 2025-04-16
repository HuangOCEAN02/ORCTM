      SUBROUTINE ITPREP

! RJ: Rewritten for MPI-Parallelization Nov 2003
! Calculates some arrays needed for the iterative solver

      USE MO_PARAM1
      USE MO_PARA2
      USE MO_COMMO1
      USE MO_PARALLEL

      IMPLICIT NONE

      INTEGER I,J

      DO I=1,IE
      DO J=1,JE
        XX(I,J) = DLXP(I,J)*DLYP(I,J)
      ENDDO
      ENDDO

      DO I=I_start,IE
      DO J=1,JE
        UF(I,J) = DEUTO(I,J)*DLYU(I,J)/DLXU(I,J)*G*CONN*STABN*DT**2
      ENDDO
      ENDDO

      DO I=1,IE
      DO J=J_start,JE
        VF(I,J) = DEUTE(I,J)*DLXV(I,J)/DLYV(I,J)*G*CONN*STABN*DT**2
      ENDDO
      ENDDO

      FF(:,:) = 0
      DO J=(J_start+1),JE
      DO I=(I_start+1),IE
        FF(I,J) = 1./(XX(I,J)+UF(I,J)+UF(I-1,J)+VF(I,J)+VF(I,J-1))
      ENDDO
      ENDDO

      CALL bounds_exch(FF)

      END
