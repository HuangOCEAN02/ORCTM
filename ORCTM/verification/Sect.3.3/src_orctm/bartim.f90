      SUBROUTINE BARTIM
!*********************************************************
! UNTERPROGR. ZUR BERECHNUNG DES WASSERSTANDES
! BAROTROPES SYSTEM (IMPLIZITES VERFAHREN,
! BEWEGUNGSGL. UND KONTI. GL.)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_ELICOM
      USE MO_COMMO1
      USE MO_UNITS

      IMPLICIT NONE

      INTEGER i,j,l,ll
      REAL B1O_G(IE_G,JE_G), bbmm, xr
!
      Z1O(:,:) = 0
      B1O(:,:) = 0

      DO 201 J=(J_start+1),JE
      DO 201 I=(I_start+1),IE
      B1O(I,J) = WETO(I,J,1) * (                                        &
     & DTDXPO(I,J) * ((CONO*U1O(I-1,J) +  UZO(I-1,J)*CONN)*DLYU(I-1,J)  &
     &       - (U1O(I,J)*CONO  +  UZO(I,J)*CONN)*DLYU(I,J) )/DLYP(I,J)  &
     & + DTDYO(I,J) * (                                                 &
     &  ( V1E(I,J) *CONO  +  CONN*VZE(I,J))*DLXV(I,J)                   &
     &- ( V1E(I,J-1) *CONO+CONN*VZE(I,J-1))*DLXV(I,J-1))/DLXP(I,J)  )
201    CONTINUE
!
      CALL gather_arr(B1O,B1O_G,p_io)
!
      IF(p_pe==0) THEN

      DO 2031 J=1,JE_G
      DO 2031 I=1,IE_G
      IF(NUM_G(I,J).EQ.0) GO TO 2031
      L=NUM_G(I,J)
      B(L)=B1O_G(I,J) *SKAL(L)
2031  CONTINUE

      DO 21 J=1,MATR-1
      BBMM=B(J)
      DO 22 I=1,KB
      B(I+J)=B(I+J)-PGL(I,J)*BBMM
22    CONTINUE
21    CONTINUE

      X(MATR)=B(MATR)/PGL(KM,MATR)

      DO 23 LL=1,MATR-1
      L=MATR-LL
      XR=B(L)
      DO 24 J=1,KB
      XR=XR-PGL(KM+J,L)*X(L+J)
24    CONTINUE
      X(L)=XR/PGL(KM,L)
23    CONTINUE

      ENDIF

      CALL p_bcast(X,p_io)

      DO 25 J=1,JE
      DO 25 I=1,IE
      IF(NUM(I,J).EQ.0) GO TO 25
      L=NUM(I,J)
      Z1O(I,J)=X(L)
25    CONTINUE

#ifndef OBC_ETA_FORC
      IF (have_g_is) THEN
        DO J=1,JE
          Z1O(1,J)  = B1O(1,J)
        ENDDO
      ENDIF
      IF (have_g_ie) THEN
        DO J=1,JE
          Z1O(IE,J)  = B1O(IE,J)
        ENDDO
      ENDIF
      IF (have_g_js) THEN
        DO I=1,IE
          Z1O(I,1)  = B1O(I,1)
        ENDDO
      ENDIF
      IF (have_g_je) THEN
        DO I=1,IE
          Z1O(I,JE)  = B1O(I,JE)
        ENDDO
      ENDIF
#endif

      CALL bounds_exch(Z1O)

      RETURN
      END
