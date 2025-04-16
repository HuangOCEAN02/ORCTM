      SUBROUTINE TRONEU
!**********************************************************************
!
!
!     TTTTTT  RRRRR    OOO    NN   N  EEEEEE  U    U
!       TT    R    R  O   O   N N  N  E       U    U
!       TT    RRRRR   O   O   N  N N  EEEEE   U    U
!       TT    R RR    O   O   N   NN  E       U    U
!       TT    R   RR   OOO    N    N  EEEEEE  UUUUUU
!
!
!**********************************************************************
!
      USE MO_PARAM1
      USE MO_PARA2
      USE MO_COMMO1
      USE MO_UNITS
      USE MO_PARALLEL
      USE MO_OBCS

      IMPLICIT NONE
      INTEGER I, J, IS, KITER
      REAL SORPAI
!
!=======================================================================
!     SBR TROPIT
!
!     PURPOSE :
!
!     A) ITERATIVE SOLUTION OF ZETA FIELD
!
!
!   ITERATIVE SOLUTION OF THE SYSTEM
!   Z*(DX*DY+G*DT**2(HW*DYW/DXW+HO*DYO/DXO+HS*DXS/DYS+HN*DXN/DYN) )
!
! =            G*DT**2*(HW*ZWW*DYW/DXW+HO*ZOO*DYO/DXO
!                      +HS*ZSS*DXS/DYS+HN*ZNN*DXN/DYN  )
!+G*DT**3*(F*HW*(ZSW-ZNW)+F*HO*(ZNO-ZSO)+FS*HS*(ZSO-ZSW)+FN*HN*(ZNW-ZNO)
!                 +B
!
!   WHERE B CONTAINS THE WINSTRESS, DIVERGENCE OF OLD FLOW AND OLD Z
!
!
! RJ: Rewritten for MPI-Parallelization Nov 2003
!
!

!
      DO J=(J_start+1),JE
      DO I=(I_start+1),IE
      
      B1O(I,J) = WETO(I,J,1) * (                                       &
     & DTDXPO(I,J) * ((CONO*U1O(I-1,J) +  UZO(I-1,J)*CONN)*DLYU(I-1,J) &
     &       - (U1O(I,J)*CONO  +  UZO(I,J)*CONN)*DLYU(I,J) )/DLYP(I,J) &
     & + DTDYO(I,J) * (                                                &
     &  ( V1E(I,J) *CONO  +  CONN*VZE(I,J))*DLXV(I,J)                  &
     &- ( V1E(I,J-1) *CONO+CONN*VZE(I,J-1))*DLXV(I,J-1))/DLXP(I,J)  )

        
      B1O(I,J) = B1O(I,J)*XX(I,J)

      ENDDO
      ENDDO

      Z1O(:,:) = 0.0

      SORPAI=1.-SORPAR
! This routine is used to the pure velocity forcing at OB.
#ifdef OBC_UV_FORC
#ifndef OBC_ETA_FORC
      IF (WEST .and. have_g_is) THEN
        DO J=1,JE
          Z1O(1,J)  = B1O(1,J)/XX(1,J)
        ENDDO
     ENDIF
      IF (East .and. have_g_ie) THEN
        DO J=1,JE
          Z1O(IE,J)  = B1O(IE,J)/XX(IE,J)
        ENDDO
      ENDIF
      IF (North .and. have_g_js) THEN
        DO I=1,IE
          Z1O(I,1)  = B1O(I,1)/XX(I,1)        
        ENDDO
      ENDIF
      IF (South .and. have_g_je) THEN
        DO I=1,IE
          Z1O(I,JE)  = B1O(I,JE)/XX(I,JE)
        ENDDO
      ENDIF
#endif
#endif
      DO KITER=1,250

        ! odd/odd and even/even Z1O elements

        DO J=2,JE-1
          IS = MOD(p_ioff+p_joff+j,2)+2
          DO I=IS,IE-1,2
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
          ENDDO
        ENDDO

        CALL bounds_exch(Z1O)

        ! odd/even and even/odd Z1O elements

        DO J=2,JE-1
          IS = MOD(p_ioff+p_joff+j+1,2)+2
          DO I=IS,IE-1,2
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
          ENDDO
        ENDDO

        CALL bounds_exch(Z1O)

      ENDDO

     ! IF (have_g_is) THEN
     !   write(IO_STDOUT,*) 'TRONEU END West'
     !   write(IO_STDOUT,*) 'Z1O(1:3,:)',Z1O(1:3,:)
     ! ENDIF
     ! IF (have_g_ie) THEN
     !   write(IO_STDOUT,*) 'TRONEU END East'
     !   write(IO_STDOUT,*) 'Z1O(IE-2:IE,:)',Z1O(IE-2:IE,:)
     ! ENDIF

      END
