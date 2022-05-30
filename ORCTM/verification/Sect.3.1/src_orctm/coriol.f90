      SUBROUTINE CORIOL
!********************************************************************
!
!
!      CCCCC   OOO   RRRRR   II   OOO   L
!     C       O   O  R    R  II  O   O  L
!     C       O   O  RRRRR   II  O   O  L
!     C       O   O  R  RR   II  O   O  L
!      CCCCC   OOO   R   RR  II   OOO   LLLLLL
!
!
!********************************************************************
!
! CORIOLIS-PARAMETER
! TIEFE*FAKTOR(CORIOLIS, ZEITSCHRITT)
!------------------------------------------------------------
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
      IMPLICIT NONE
      INTEGER i,j,jj,ioff,joff
      REAL rad
!
      WRITE(IO_STDOUT,*)'IN CORIOL'
!
      RAD=PI/180.
!
      IF(p_pe==p_io) THEN
         OPEN(IO_IN_ANTA,FILE='anta'                                    &
     &                  ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
         READ(IO_IN_ANTA) IBLA
         READ(IO_IN_ANTA) GILA_G
         READ(IO_IN_ANTA) IBLA
         READ(IO_IN_ANTA) GIPH_G
         CLOSE(IO_IN_ANTA)
      ENDIF

      CALL p_bcast(GILA_G,p_io)
      CALL p_bcast(GIPH_G,p_io)

      ioff = 2*p_ioff
      joff = 2*p_joff
      GILA(:,:) = GILA_G(ioff+1:ioff+2*ie,joff+1:joff+2*je)
      GIPH(:,:) = GIPH_G(ioff+1:ioff+2*ie,joff+1:joff+2*je)
!
      DO 1478 J=1,JE
      DO 1478 I=1,IE
      ALAT(I,J)=GIPH(2*I,2*J)/RAD
1478  CONTINUE

      DO 1 JJ=1,JE
      SIFTFE(JJ)=DT
      COFTFO(JJ)=0.
      COFTFE(JJ)=0.    
      SIFTFO(JJ)=DT
1     CONTINUE
!
      RETURN
      END
