       SUBROUTINE TRIAN
!**********************************************************
! UNTERPROGRAMM ZUR TRIANGULARISIERUNG DER MATRIX
! LOESUNG LINEARER GL., GAUSS-SCHER ALGOR.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_ELICOM
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
      IMPLICIT NONE

      INTEGER I, IN, IO, IS, IW, J, JA, JJ, JO, K, K1, KITER, L, LL
      REAL FAKE, FAKN, FAKS, FAKW, RE, RR1, WILKIN, ZZZ
      REAL FW(IE,JE), FE(IE,JE), FN(IE,JE), FS(IE,JE)
      REAL FW_G(IE_G,JE_G), FE_G(IE_G,JE_G)
      REAL FN_G(IE_G,JE_G), FS_G(IE_G,JE_G)

      DO 4 J=1,ILL
      DO 4 I=1,KBM
      PGL(I,J)=0.
4     CONTINUE

      DO J=2,JE-1
      DO I=2,IE-1
         FAKW=DEUTO(I-1,J)*DTI*SIFTFO(J)*DTDXPO(I,J)*DTDXUO(I-1,J)      &
     &   *DLYU(I-1,J)/DLYP(I,J)
         FAKE=DEUTO(I,J)*DTI*SIFTFO(J)*DTDXPO(I,J)*DTDXUO(I,J)          &
     &   *DLYU(I,J)/DLYP(I,J)
         FAKN=DTI*SIFTFE(J-1)*DEUTE(I,J-1)*DPYO(I,J)*DPYE(I,J-1)        &
     &   *DTDXPO(I,J)/DTDXUE(I,J-1)
         FAKS=DTI*SIFTFE(J)*DEUTE(I,J)*DPYO(I,J)*DPYE(I,J)              &
     &   *DTDXPO(I,J)/DTDXUE(I,J)
!
         ZZZ=-G*CONN*STABN
         FW(I,J)=FAKW*ZZZ
         FE(I,J)=FAKE*ZZZ
         FN(I,J)=FAKN*ZZZ
         FS(I,J)=FAKS*ZZZ
      ENDDO
      ENDDO

      if (have_g_is) then
       I=1
        DO J=2,JE-1
         FAKW=DEUTO(I-1,J)*DTI*SIFTFO(J)*DTDXPO(I,J)*DTDXUO(I-1,J)      &
     &   *DLYU(I-1,J)/DLYP(I,J)
         FAKE=DEUTO(I,J)*DTI*SIFTFO(J)*DTDXPO(I,J)*DTDXUO(I,J)          &
     &   *DLYU(I,J)/DLYP(I,J)
         FAKN=DTI*SIFTFE(J-1)*DEUTE(I,J-1)*DPYO(I,J)*DPYE(I,J-1)        &
     &   *DTDXPO(I,J)/DTDXUE(I,J-1)
         FAKS=DTI*SIFTFE(J)*DEUTE(I,J)*DPYO(I,J)*DPYE(I,J)              &
     &   *DTDXPO(I,J)/DTDXUE(I,J)
         ZZZ=-G*CONN*STABN
         FW(I,J)=FAKW*ZZZ
         FE(I,J)=FAKE*ZZZ
         FN(I,J)=FAKN*ZZZ
         FS(I,J)=FAKS*ZZZ
        ENDDO
      endif

      if (have_g_ie) then
       I=IE
        DO J=2,JE-1
         FAKW=DEUTO(I-1,J)*DTI*SIFTFO(J)*DTDXPO(I,J)*DTDXUO(I-1,J)      &
     &   *DLYU(I-1,J)/DLYP(I,J)
         FAKE=DEUTO(I,J)*DTI*SIFTFO(J)*DTDXPO(I,J)*DTDXUO(I,J)          &
     &   *DLYU(I,J)/DLYP(I,J)
         FAKN=DTI*SIFTFE(J-1)*DEUTE(I,J-1)*DPYO(I,J)*DPYE(I,J-1)        &
     &   *DTDXPO(I,J)/DTDXUE(I,J-1)
         FAKS=DTI*SIFTFE(J)*DEUTE(I,J)*DPYO(I,J)*DPYE(I,J)              &
     &   *DTDXPO(I,J)/DTDXUE(I,J)
         ZZZ=-G*CONN*STABN
         FW(I,J)=FAKW*ZZZ
         FE(I,J)=FAKE*ZZZ
         FN(I,J)=FAKN*ZZZ
         FS(I,J)=FAKS*ZZZ
        ENDDO
      endif

      if (have_g_js) then
       J=1
        DO I=2,IE-1
         FAKW=DEUTO(I-1,J)*DTI*SIFTFO(J)*DTDXPO(I,J)*DTDXUO(I-1,J)      &
     &   *DLYU(I-1,J)/DLYP(I,J)
         FAKE=DEUTO(I,J)*DTI*SIFTFO(J)*DTDXPO(I,J)*DTDXUO(I,J)          &
     &   *DLYU(I,J)/DLYP(I,J)
         FAKN=DTI*SIFTFE(J-1)*DEUTE(I,J-1)*DPYO(I,J)*DPYE(I,J-1)        &
     &   *DTDXPO(I,J)/DTDXUE(I,J-1)
         FAKS=DTI*SIFTFE(J)*DEUTE(I,J)*DPYO(I,J)*DPYE(I,J)              &
     &   *DTDXPO(I,J)/DTDXUE(I,J)
         ZZZ=-G*CONN*STABN
         FW(I,J)=FAKW*ZZZ
         FE(I,J)=FAKE*ZZZ
         FN(I,J)=FAKN*ZZZ
         FS(I,J)=FAKS*ZZZ
        ENDDO
      endif

      if (have_g_je) then
       J=JE
        DO I=2,IE-1
         FAKW=DEUTO(I-1,J)*DTI*SIFTFO(J)*DTDXPO(I,J)*DTDXUO(I-1,J)      &
     &   *DLYU(I-1,J)/DLYP(I,J)
         FAKE=DEUTO(I,J)*DTI*SIFTFO(J)*DTDXPO(I,J)*DTDXUO(I,J)          &
     &   *DLYU(I,J)/DLYP(I,J)
         FAKN=DTI*SIFTFE(J-1)*DEUTE(I,J-1)*DPYO(I,J)*DPYE(I,J-1)        &
     &   *DTDXPO(I,J)/DTDXUE(I,J-1)
         FAKS=DTI*SIFTFE(J)*DEUTE(I,J)*DPYO(I,J)*DPYE(I,J)              &
     &   *DTDXPO(I,J)/DTDXUE(I,J)
         ZZZ=-G*CONN*STABN
         FW(I,J)=FAKW*ZZZ
         FE(I,J)=FAKE*ZZZ
         FN(I,J)=FAKN*ZZZ
         FS(I,J)=FAKS*ZZZ
        ENDDO
      endif

      CALL gather_arr(FW,FW_G,p_io)
      CALL gather_arr(FE,FE_G,p_io)
      CALL gather_arr(FN,FN_G,p_io)
      CALL gather_arr(FS,FS_G,p_io)

      IF(p_pe==p_io) THEN
        DO J=2,JE_G-1
        DO I=2,IE_G-1
          IF(NUM_G(I,J).EQ.0) CYCLE
          L=NUM_G(I,J)
          IN=NUM_G(I,J-1)-L
          IS=NUM_G(I,J+1)-L
          IW=NUM_G(I-1,J)-L
          IO=NUM_G(I+1,J)-L
          IF(NUM_G(I,J-1).LE.0) IN=0
          IF(NUM_G(I,J+1).LE.0) IS=0
          IF(NUM_G(I-1,J).LE.0) IW=0
          IF(NUM_G(I+1,J).LE.0) IO=0
!
          PGL(KM,L)=1.-FW_G(I,J)-FE_G(I,J)-FN_G(I,J)-FS_G(I,J)
          PGL(KM+IW,L)=FW_G(I,J)+PGL(KM+IW,L)
          PGL(KM+IO,L)=FE_G(I,J)+PGL(KM+IO,L)
          PGL(KM+IN,L)=PGL(KM+IN,L)+FN_G(I,J)
          PGL(KM+IS,L)=PGL(KM+IS,L)+FS_G(I,J)
!     
          SKAL(L)=1./PGL(KM,L)
          DO LL=1,KBM
            PGL(LL,L)=PGL(LL,L)*SKAL(L)
          ENDDO
        ENDDO
        ENDDO

        I=1
        DO J=1,JE_G
          IF(NUM_G(I,J).EQ.0) CYCLE
          L=NUM_G(I,J)
          PGL(KM,L)=1.0
        ENDDO

        I=IE_G
        DO J=1,JE_G
          IF(NUM_G(I,J).EQ.0) CYCLE
          L=NUM_G(I,J)
          PGL(KM,L)=1.0
        ENDDO

        J=1
        DO I=1,IE_G
          IF(NUM_G(I,J).EQ.0) CYCLE
          L=NUM_G(I,J)
          PGL(KM,L)=1.0
        ENDDO

        J=JE_G
        DO I=1,IE_G
          IF(NUM_G(I,J).EQ.0) CYCLE
          L=NUM_G(I,J)
          PGL(KM,L)=1.0
        ENDDO

        WILKIN=1.
        WRITE(IO_STDOUT,*)' LAUF ',L,MATR
!
        DO KITER=2,MATR
          K1=KITER-1
          RE=PGL(KM,K1)
          JA=KM
          DO K=1,KB
            PGL(K,K1)=0.
            JA=JA-1
            IF(K+K1.GT.MATR) CYCLE
            JO=JA+KB
            RR1=PGL(JA,K+K1)/RE
            PGL(K,K1)=RR1
            IF(ABS(RR1).GT.WILKIN) WILKIN=ABS(RR1)
            DO JJ=JA,JO
              PGL(JJ,K+K1)=PGL(JJ,K+K1)-RR1*PGL(JJ+K,K1)
            ENDDO
          ENDDO
        ENDDO
602     FORMAT(' MAX.ELIMI',E14.5)
        WRITE(IO_STDOUT,602) WILKIN
      ENDIF

      CALL p_bcast(PGL,p_io)
!
      RETURN
      END
