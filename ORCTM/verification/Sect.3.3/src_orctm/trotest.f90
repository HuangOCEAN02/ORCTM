      SUBROUTINE TROTEST
      USE MO_PARAM1
      USE MO_PARA2
      USE MO_COMMO1
      USE MO_UNITS

      USE MO_PARALLEL
!
!=======================================================================
!     PURPOSE :
!
!     
!UWE  PREPARATION FOR ITERATIVE SOLUTION OF ZETA FIELD
!     OPTIMISES SORPAR
!
! RJ: Rewritten for MPI-Parallelization Nov 2003
!

      IMPLICIT NONE
      INTEGER I, J, II, JJ, IS, ISOR, KITER
      REAL CHABAS, CONTRIC, CONVMIN, SORPAI, SORPMIN, SUMCHA, ZALT
      REAL SQRTRND, RND


      SORPAR=1.75
!
      CONVMIN=9.E99
      SORPMIN=SORPAR
!
      DO ISOR=1,200
      SORPAR=SORPAR+0.001
!
!UWE  INITIALIZE WITH RANDOM NUMBERS
!hh random_number doesn't work with fujitsu f90!
!hh    CALL RANDOM_NUMBER(XHH)
!hh    B1O(I,J)=(XHH-0.5)*WETO(I,J,1)
!OtB  Using a simple one

! RJ: We do that in the following way so that the result
!     does not depend on the decomposition
!     Why is that done for every step in the original version ????
!     We do it only once here

      IF (ISOR==1) THEN
        DO JJ=1,JE_G
        DO II=1,IE_G
          I = II - p_ioff
          J = JJ - p_joff
          RND = SQRTRND()
          IF(I>=1 .AND. I<=IE .AND. J>=1 .AND. J<=JE) THEN
            B1O(I,J)=(RND-0.5)*WETO(I,J,1)*XX(I,J)
          ENDIF
        ENDDO
        ENDDO

        CALL bounds_exch(B1O)
      ENDIF

      Z1O(:,:) = 0.

      SORPAI=1.-SORPAR

      DO KITER=1,200

        SUMCHA=0.

        ! odd/odd and even/even Z1O elements

        DO J=2,JE-1
          IS = MOD(p_ioff+p_joff+j,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch(Z1O)

        ! odd/even and even/odd Z1O elements

        DO J=2,JE-1
          IS = MOD(p_ioff+p_joff+j+1,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch(Z1O)

        CALL global_sum(SUMCHA)

        IF(KITER.EQ.1) CHABAS=SUMCHA

      ENDDO ! KITER Loop

      CONTRIC=SUMCHA/CHABAS
      IF(CONTRIC.LT.CONVMIN)THEN
        SORPMIN=SORPAR
        CONVMIN=CONTRIC
      ENDIF

      ENDDO ! ISOR Loop

      SORPAR=SORPMIN
      WRITE(IO_STDOUT,*) 'TROTEST SORPAR=',SORPAR
      RETURN
      END
