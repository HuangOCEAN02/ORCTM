      SUBROUTINE BODEN
!*************************************************************
!
!     BBBBB    OOO   DDDDD   EEEEEE  NN   N
!     B    B  O   O  D    D  E       NNN  N
!     BBBBB   O   O  D    D  EEEEE   N NN N
!     B    B  O   O  D    D  E       N  NNN
!     BBBBB    OOO   DDDDD   EEEEEE  N   NN
!
!
!--------------------------------------------------------------
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
      IMPLICIT NONE
!
      CHARACTER*6 BLORD
      REAL TOPHILF(IE,JE)
      CHARACTER*1 CDEPTO(IE_G,JE_G)
      REAL DEPTO_G(IE_G,JE_G)
      REAL CMUE,CMUO,CMWE,CMWO,CUE,CUO,CWE,CWO,DEPALT,DEPMAX,DEPMIN
      REAL SUGGO,SUM,SUPPO,TMAX,dddwo
      INTEGER I,I1,I2,IANF,IEND,III,IPR,IPRANZ,J,JEV,JJ,JJJ,K
      INTEGER LAUF,NDIF1
      REAL minvalue
!
!------------------------------------------------------------------------
!
!   HERE   :  READ TOPOGRAPHY FROM EXTERNAL FILE ON ARRAY DEPTH
!
!      (BE CAREFUL WITH THE SIGN OF THE TOPOGR. DATA ! , THIS VERSION
!       USES POSITIVE DEPTH DATA BECAUSE SYNBAPS COMES WITH POSITIVE
!       VALUES)
!
      minvalue=MIN(1.0,DZW(1)/10.0)
!
      SUM=0.
      DO K=1,KE
         SUM=SUM+DZW(K)
      ENDDO
!
      DEPMAX=SUM
      DEPMIN=DZW(1)+DZW(2)+minvalue
!
      WRITE(IO_STDOUT,*)' --------------------------------------------- &
     &                    --------'
      WRITE(IO_STDOUT,*)' MINIMUM WATER DEPTH IS SET TO ',DEPMIN
      WRITE(IO_STDOUT,*)' MAXIMUM WATER DEPTH IS SET TO ',DEPMAX
      WRITE(IO_STDOUT,*)' --------------------------------------------- &
     &                    --------'
!
      depto_g(:,:)=0.0

      IF (ISTART .EQ. 0) THEN
      DO 1296 J=1,JE
      DO 1296 I=1,IE
        IF(DEPTO(I,J).LE.minvalue) DEPTO(I,J)=0.
        IF(DEPTO(I,J).GT.minvalue.AND.DEPTO(I,J).LT.DEPMIN) DEPTO(I,J)=DEPMIN
1296  CONTINUE
      ENDIF
!
      TMAX=0.
      DO J=1,JE
       DO I=1,IE
        TOPHILF(I,J)=DEPTO(I,J)
        TMAX=MAX(TMAX,DEPTO(I,J))
       ENDDO
      ENDDO
      CALL global_max(TMAX)
      WRITE(IO_STDOUT,*)'TOPOGRAPHIE CHECK: ',TMAX
!
      IF (ISTART .EQ. 0) CALL gather_arr(DEPTO,DEPTO_G,p_io)
!
      IF (p_pe==p_io) THEN
        OPEN(IO_IN_GIGI,FILE='topo'                                     &
     &              ,ACCESS='SEQUENTIAL',FORM='FORMATTED')
        DO I1=1,IE_G,20
          I2=MIN(I1+19,IE_G)
          WRITE(IO_STDOUT,*) 'LESSTREIFEN ',I1,I2
          IF (ISTART .EQ. 0) THEN
            WRITE(IO_IN_GIGI,*)'STREIFEN ',I1,I2
            DO J=1,JE_G
              WRITE(IO_IN_GIGI,63885)J,(DEPTO_G(I,J),I=I1,I2)
            ENDDO
          ELSE
            READ(IO_IN_GIGI,*) BLORD
            DO J=1,JE_G
              READ(IO_IN_GIGI,63885)JJJ,(DEPTO_G(I,J),I=I1,I2)
            ENDDO
          ENDIF
63885     FORMAT(I5,20F8.2)
        ENDDO
        CLOSE(IO_IN_GIGI)
      ENDIF
!
      IF (ISTART .NE. 0) CALL scatter_arr(DEPTO_G,DEPTO,p_io)
!
      DO J=1,JE
       DO I=1,IE
        DO K=1,KE
         IF(ABS(DEPTO(I,J)-TIESTU(K)).LT.minvalue/2) DEPTO(I,J)=TIESTU(K)+minvalue
        ENDDO
       ENDDO
      ENDDO

      CALL bounds_exch(DEPTO)
!
! Now the calculations on DEPTO are finished, we gather and distribute global DEPTO_G 
!
      CALL gather_arr(DEPTO,DEPTO_G,p_io)
      CALL p_bcast(DEPTO_G,p_io)

      DO 296 J=1,JE
      DO 296 I=1,IE1
        DEUTO(I,J)=MIN(DEPTO(I,J),DEPTO(I+1,J))
296   CONTINUE

      DEUTO(I_start,:)=DEUTO(1,:)
      DEUTO(IE,:)=DEUTO(IE1,:)

      DO 297 J=1,JE1
      DO 297 I=1,IE
        DEUTE(I,J)=MIN(DEPTO(I,J),DEPTO(I,J+1))
297   CONTINUE

      DEUTE(:,J_start)=DEUTE(:,1)
      DEUTE(:,JE)=DEUTE(:,JE1)

      CALL bounds_exch(DEUTO)
      CALL bounds_exch(DEUTE)
!
!----------------------------------------------------------------------
!     COMPUTE BOUNDARY TOPOGRAPHY , I.E. 3 ROWS
!     FROM THE NORTHERN AND FROM THE SOUTHERN BOUNDARY ARE SET
!     TO ZERO DEPTH
!
!     CALCULATION OF DEPTH AT SCALAR POINTS (HP : DEPTH AT PRESSURE PT.)
!        HP IS THE MAXIMUM DEPTH OF THE 4 SURROUNDING U-POINT DEPTHS (H)

      LAUF=0
      JJ=0
      NUM_G(:,:) = 0
      DO III=0,IE_G-1
         JJ=JJ+1
         IF(MOD(III,2).EQ.0)THEN
            I=III/2+1
         ELSE
            I=IE_G-(III/2)
         ENDIF
         DO J=1,JE_G
            IF(DEPTO_G(I,J).GE.1.) THEN
               LAUF=LAUF+1
               NUM_G(I,J)=LAUF
            ENDIF
         ENDDO
      ENDDO

!     Set local NUM

      NUM(1:ie,1:je) = NUM_G(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)
!
      WRITE(IO_STDOUT,*)' FEUCHTPUNKTE',LAUF
!
!============WERTE FUER TRIANGULARISIERUNG
      MATR=LAUF
      NMAX=0
!
      DO J=1,JE_G
      DO I=2,IE_G
        IF(NUM_G(I,J).EQ.0 .OR. NUM_G(I-1,J).EQ.0) CYCLE
        NDIF1=ABS(NUM_G(I,J)-NUM_G(I-1,J))
        IF(NDIF1.GT.NMAX) NMAX=NDIF1
      ENDDO
      ENDDO
      WRITE(IO_STDOUT,*)' MAXDIF ', NMAX
!
!======LAND/OZEAN-STEUERFELDER UND SCHICHTDICKEN
!
      DO 444 K=1,KE
#ifndef PARTIALCELL      
!
!----------------------------------p----------------------------
      DO J=1,JE
      DO I=1,IE
        IF(TIESTU(K).LT.DEPTO(I,J)) WETO(I,J,K)=1.
        IF(TIESTU(K+1).LT.DEPTO(I,J)) DDPO(I,J,K)=DZW(K)
        IF(TIESTU(K).LT.DEPTO(I,J) .AND. TIESTU(K+1).GE. DEPTO(I,J))    &
     &    DDPO(I,J,K)=DEPTO(I,J)-TIESTW(K)
        IF(DDPO(I,J,K).GT.ZERO) DPIO(I,J,K)=1./DDPO(I,J,K)
      ENDDO
      ENDDO
!----------------------------------u----------------------------
      DO J=1,JE
      DO I=I_start,IE
        IF(TIESTU(K+1).LT.DEUTO(I,J)) DDUO(I,J,K)=DZW(K)
        IF(TIESTU(K).LT.DEUTO(I,J) .AND. TIESTU(K+1).GE. DEUTO(I,J))    &
     &    DDUO(I,J,K)=DEUTO(I,J)-TIESTW(K)
        IF(DDUO(I,J,K).GT.ZERO) AMSUO(I,J,K)=1.
      ENDDO
      ENDDO
!----------------------------------v----------------------------
      DO J=J_start,JE
      DO I=1,IE
        IF(TIESTU(K+1).LT.DEUTE(I,J)) DDUE(I,J,K)=DZW(K)
        IF(TIESTU(K).LT.DEUTE(I,J) .AND. TIESTU(K+1).GE. DEUTE(I,J))    &
     &    DDUE(I,J,K)=DEUTE(I,J)-TIESTW(K)
        IF(DDUE(I,J,K).GT.ZERO) AMSUE(I,J,K)=1.
      ENDDO
      ENDDO
!
#else
!
!  DO 444 K=1,KE
!----------------------------------p----------------------------
      DO J=1,JE
      DO I=1,IE
        IF(TIESTU(K) .LT. DEPTO(I,J)) WETO(I,J,K)=1.
        IF(TIESTW(K+1) .LT. DEPTO(I,J)) DDPO(I,J,K)=DZW(K)
        ! bottom
        IF(TIESTW(K) .LT. DEPTO(I,J) .AND. TIESTW(K+1) .GE. DEPTO(I,J)) then
        
           IF (((DEPTO(I,J) - TIESTW(K))/DZW(K)) .GE. 0.2 ) then
               WETO(I,J,K) = 1.0
               DDPO(I,J,K) = DEPTO(i,j)-TIESTW(K)
           ELSE
               DEPTO(I,J) = TIESTW(k)
           ENDIF
           
        ENDIF
        IF(DDPO(I,J,K).GT.ZERO) DPIO(I,J,K)=1./DDPO(I,J,K)
      ENDDO
      ENDDO
!----------------------------------u----------------------------
      DO J=1,JE
      DO I=I_start,IE
        IF(TIESTW(K+1) .LT. DEUTO(I,J)) DDUO(I,J,K)=DZW(K)
         ! bottom
        IF(TIESTW(K) .LT. DEUTO(I,J) .AND. TIESTW(K+1) .GE. DEUTO(I,J)) then        
           IF (((DEUTO(I,J) - TIESTW(K))/DZW(K)) .GE. 0.2 ) then
               DDUO(I,J,K) = DEUTO(i,j)-TIESTW(K)
           ELSE
               DEUTO(I,J) = TIESTW(K)
           ENDIF 
        ENDIF   
        IF(DDUO(I,J,K).GT.ZERO) AMSUO(I,J,K)=1.
      ENDDO
      ENDDO
!----------------------------------v----------------------------
      DO J=J_start,JE
      DO I=1,IE
        IF(TIESTW(K+1) .LT. DEUTE(I,J)) DDUE(I,J,K)=DZW(K)
         ! bottom
        IF(TIESTW(K) .LT. DEUTE(I,J) .AND. TIESTW(K+1) .GE. DEUTE(I,J)) then  
           IF (((DEUTE(I,J) - TIESTW(K))/DZW(K)) .GE. 0.2 ) then
               DDUE(I,J,K) = DEUTE(i,j)-TIESTW(K)
           ELSE
               DEUTE(I,J) = TIESTW(k)
           ENDIF 
        ENDIF
        IF(DDUE(I,J,K).GT.ZERO) AMSUE(I,J,K)=1.
      ENDDO
      ENDDO

#endif
444   CONTINUE

!
      DO 1457 K=1,KE
      DO 1457 J=1,JE1
      DO 1457 I=1,IE1
      IF(DDUE(I,J,K).GT.0..AND.DDPO(I,J,K)*DDPO(I,J+1,K).LT.minvalue**2) THEN
        WRITE(IO_STDOUT,*) ' DDUE ',I+p_ioff,J+p_joff,K                 &
     &                  ,DDUE(I,J,K),DDPO(I,J,K),DDPO(I,J+1,K)
      ENDIF
      IF(DDUO(I,J,K).GT.0..AND.DDPO(I,J,K)*DDPO(I+1,J,K).LT.minvalue**2) THEN
        WRITE(IO_STDOUT,*) ' DDUO ',I+p_ioff,J+p_joff,K                 &
     &                  ,DDUO(I,J,K),DDPO(I,J,K),DDPO(I+1,J,K)
      ENDIF
      SUPPO=0.
      SUGGO=0.
      IF(DDPO(I,J,K).GT.minvalue) THEN
      SUPPO=SUPPO+DDPO(I,J,K)
      SUGGO=SUGGO+1.
      ENDIF
      IF(DDPO(I+1,J,K).GT.minvalue) THEN
      SUPPO=SUPPO+DDPO(I+1,J,K)
      SUGGO=SUGGO+1.
      ENDIF
      IF(DDPO(I,J+1,K).GT.minvalue) THEN
      SUPPO=SUPPO+DDPO(I,J+1,K)
      SUGGO=SUGGO+1.
      ENDIF
      IF(DDPO(I+1,J+1,K).GT.minvalue) THEN
      SUPPO=SUPPO+DDPO(I+1,J+1,K)
      SUGGO=SUGGO+1.
      ENDIF
      IF(SUGGO.GE.1.) DDPSIO(I,J,K)=SUPPO/SUGGO
1457  CONTINUE

      DO K=1,KE
        DO I=1,IE1
          DDPSIO(I,JE,K)=DDPSIO(I,JE1,K)
        ENDDO
        DO J=1,JE
          DDPSIO(IE,J,K)=DDPSIO(IE1,J,K)
        ENDDO
      ENDDO

      CALL bounds_exch(DDPSIO)
!
      DO J=1,JE
      DO I=1,IE
        KCONDEP(I,J)=NINT(WETO(I,J,1))-99*NINT(1.-WETO(I,J,1))
      ENDDO
      ENDDO
!
      DO J=1,JE
      DO I=I_start,IE
        IF(DEUTO(I,J).GT.minvalue) DEUTIO(I,J)=1./DEUTO(I,J)
      ENDDO
      ENDDO
!
      DO J=J_start,JE
      DO I=1,IE
        IF(DEUTE(I,J).GT.minvalue) DEUTIE(I,J)=1./DEUTE(I,J)
      ENDDO
      ENDDO
!
!====================================================================
!     
      CALL bounds_exch(DEPTO) 
      CALL bounds_exch(DEUTE,DEUTIE)
      CALL bounds_exch(DEUTO,DEUTIO)
      CALL bounds_exch(DDUO)
      CALL bounds_exch(DDUE)
      CALL bounds_exch(AMSUE)
      CALL bounds_exch(AMSUO)
      CALL bounds_exch(WETO)
      CALL bounds_exch(DDPO)
   
#ifdef NON_HYDROSTATIC    
      DO 1458 K=1,KEP
      DO 1458 J=1,JE
      DO 1458 I=1,IE
        if (k .eq. 1) then
          if (DDPO(I,J,k) .eq. 0) then
            dddwo=dzw(k)/2.0
          else
            dddwo=DDPO(I,J,k)/2.0
          endif
          
        elseif (k .eq. kep) then
          if (DDPO(I,J,KE) .eq. 0) then
            dddwo=dzw(ke)
          else
            dddwo=DDPO(I,J,KE)
          endif
          
        else
          if (DDPO(I,J,K-1) .eq. 0) then
            dddwo=(dzw(k)+dzw(k-1))/2.0
          elseif (DDPO(I,J,K) .eq. 0) then
            dddwo=DDPO(I,J,K-1)
          else
            dddwo=(DDPO(I,J,K-1)+DDPO(I,J,K))/2.0
          endif
          
        endif
        
        DZPNH(I,J,K)=dddwo
1458  CONTINUE
       CALL bounds_exch(DZPNH)  
#endif /*NON_HYDROSTATIC*/      
      
!
!----------------------------------------------------------------------
!  PRINT DEPTHS AT VECTOR CELLS
!
      CALL gather_arr(DEPTO,DEPTO_G,p_io)
      CALL p_bcast(DEPTO_G,p_io)
      DO I=1,IE_G
         DO J=1,JE_G
            CDEPTO(I,J)='7'
            IF(DEPTO_G(I,J).LE.6000.) CDEPTO(I,J)='6'
            IF(DEPTO_G(I,J).LE.5000.) CDEPTO(I,J)='5'
            IF(DEPTO_G(I,J).LE.4000.) CDEPTO(I,J)='4'
            IF(DEPTO_G(I,J).LE.3000.) CDEPTO(I,J)='3'
            IF(DEPTO_G(I,J).LE.2000.) CDEPTO(I,J)='2'
            IF(DEPTO_G(I,J).LE.1000.) CDEPTO(I,J)='1'
            IF(DEPTO_G(I,J).EQ.   0.) CDEPTO(I,J)='*'
         ENDDO
      ENDDO
      IPRANZ=(2*IE_G-1)/120 + 1
      DO IPR=1,IPRANZ
       IANF =    (IPR-1)*60    + 1
       IEND =     IPR   *60
       WRITE(IO_STDOUT,*) '    1M - 1000M : 1'
       WRITE(IO_STDOUT,*) ' 1001M - 2000M : 2'
       WRITE(IO_STDOUT,*) ' 2001M - 3000M : 3'
       WRITE(IO_STDOUT,*) ' 3001M - 4000M : 4'
       WRITE(IO_STDOUT,*) ' 4001M - 5000M : 5'
       WRITE(IO_STDOUT,6337) (I,I=IANF+1,IEND,2) 
       DO J=1,JE_G
        IF(MOD(J,10).NE.0) THEN
         WRITE(IO_STDOUT,6334) J,(CDEPTO(MOD(I-1,IE_G)+1,J),I=IANF,IEND)
        ELSE
         WRITE(IO_STDOUT,6335) J,(CDEPTO(MOD(I-1,IE_G)+1,J),I=IANF,IEND)
        ENDIF
       ENDDO
      ENDDO
!
 6334 FORMAT(6X,I3,3X,6(9(A1,1X),(A1,1H.)))
 6335 FORMAT(6X,I3,3X,6(9(A1,1H-),(A1,1H.)))
 6337 FORMAT(6X,'DEPTHS AT VECTOR CELLS:',/,11X,30I4,/)
!
!=====================================================================
      WRITE(IO_STDOUT,*)'SR BODEN UEBERLEBT!!!'
!
      RETURN
      END
