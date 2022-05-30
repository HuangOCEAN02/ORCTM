      SUBROUTINE FINDBOT
      USE MO_PARAM1
      USE MO_COMMO1
!:: SBR FINDS DEEPEST WET LEVEL ON PRESSURE POINT
!:: JHJ SEPT. 2, 1999
      DO 100 I=1,IE
       DO 100 J=1,JE
        IF (WETO(I,J,1).LT.0.5) THEN
         KBOT(I,J)=0
         GOTO 100
        ENDIF
        KBOT(I,J)=KE
        DO 200 K=2,KE
         IF (WETO(I,J,K) .LT. 0.5) THEN
           KBOT(I,J)=K-1
           GOTO 100
         ENDIF
 200    CONTINUE
 100  CONTINUE
      RETURN
      END
