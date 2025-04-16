      SUBROUTINE AUFR
!****************************************************************
!
!     AAAAAA  U    U  FFFFFF  RRRRRR
!     A    A  U    U  F       R    R
!     AAAAAA  U    U  FFFFF   RRRRRR
!     A    A  U    U  F       R  RR
!     A    A  UUUUUU  F       R   RR
!
!
!*****************************************************************
!
!-----------------------------------------------------------------
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
#ifdef OBCS_UV_FLATHER
      USE MO_OBCS_NEST
#endif
!
      IMPLICIT NONE
      INTEGER(KIND=i4) I4DATE(4)
      REAL(KIND=sp) RDT27, RDT28, RDT(2)
      INTEGER i, j, k, IDATE
!
      IF (p_pe==p_io) THEN
        OPEN(IO_IN_Z370,FILE='Z37000',STATUS='UNKNOWN'                  &
     &               ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

        OPEN(IO_IN_Z380,FILE='Z38000',STATUS='UNKNOWN'                  &
     &               ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

        REWIND(IO_IN_Z370)
        REWIND(IO_IN_Z380)
!
        READ(IO_IN_Z370,ERR=99,END=99) I4DATE
        READ(IO_IN_Z370,ERR=99,END=99) RDT27
        GOTO 1
   99   RDT27=0.0
        WRITE(IO_STDOUT,*)'NO DATA ON RESTART FILE Z37000 !!!'
    1   READ(IO_IN_Z380,ERR=88,END=88) I4DATE
        READ(IO_IN_Z380,ERR=88,END=88) RDT28
        GOTO 2
   88   RDT28=0.0
        WRITE(IO_STDOUT,*)'NO DATA ON RESTART FILE Z38000 !!!'
    2   WRITE(IO_STDOUT,6000)RDT27,RDT28
 6000   FORMAT(' Z37000 TIMESTEPS : ',F10.0                             &
     &        ,' Z38000 TIMESTEPS : ',F10.0)
!
        IF(RDT27.GT.RDT28)THEN
          IUNIT=IO_IN_Z370
          IFLAG=1
        ELSE
          IUNIT=IO_IN_Z380
          IFLAG=-1
        ENDIF
!
        REWIND(IUNIT)
!
        READ(IUNIT) I4DATE
        READ(IUNIT) RDT

        LDT=NINT(RDT(1))
        IDATE=I4DATE(1)
      ENDIF

      CALL p_bcast(IFLAG,p_io)
      CALL p_bcast(IUNIT,p_io)
      CALL p_bcast(IDATE,p_io)
      CALL p_bcast(LDT,p_io)

      LYEARS=IDATE/10000
      LMONTS=(IDATE-LYEARS*10000)/100
      LDAYS=IDATE-(LYEARS*10000)-(LMONTS*100)
      IF (ISTART .EQ. 0) LYEARS=-1
      LYEARS=-1
!
!   INTENTIONAL CRASH IN YEAR LY_END
      IF(LYEARS.EQ.LY_END) THEN
         WRITE(IO_STDOUT,*)' INTENTIONAL CRASH AFTER ',LYEARS,' YEARS'
         CALL ABSTURZ
      ENDIF
!::
!:: SET YEAR TO LY_START
      IF(LY_START .GE. 0) THEN
       LYEARS=LY_START
       LMONTS=0
       LDAY=0
       WRITE(IO_STDOUT,*)'ATTN !! SET OFFSET --> LYEARS TO: ',LYEARS
      ENDIF
!:: SET MONTH TO LM_START
      IF(LM_START .GE. 0) THEN
       LMONTS=LM_START-1
       LDAY=0
       WRITE(IO_STDOUT,*)'ATTN !! SET OFFSET --> LMONTH TO: ',LMONTH
      ENDIF
!::
      WRITE(IO_STDOUT,6001) IUNIT,LYEARS,LMONTS,LDAYS,LDT
 6001 FORMAT('READ FROM UNIT : ',I5,' START YEAR,MONTH,DAY,LDT : ',4I10)
!CCC
!  READ ZONAL VELOCITY COMPONENT
      DO K=1,KE
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,UOO(:,:,K),1)
      ENDDO
!  READ MERIDIONAL  VELOCITY COMPONENT
      DO K=1,KE
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,VOE(:,:,K),2)
      ENDDO
!:: UPDATE VELOCITY FIELDS
      CALL OCTIMF
!  READ TEMPERATURE
      DO K=1,KE
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,THO(:,:,K))
      ENDDO
!  READ SALINITY
      DO K=1,KE
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,SAO(:,:,K))
      ENDDO
!:: READ 2-D FIELDS
!  READ SEA SURFACE ELEVATION
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,ZO)
!  READ SEA SURFACE ELEVATION CHANGE
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,Z1O)
! READ BAROTROPIC ZONAL VELOCITY
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,USO)
! READ BAROTROPIC MERIDIONAL VELOCITY
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,VSE)
      
#ifdef OBCS_UV_FLATHER
       ZOOLD  =  ZO
       USOOLD = USO
       VSEOLD = VSE
#endif

!  READ SEA ICE THICKNESS
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,SICTHO)
!  READ SEA ICE CONCENTRATION         
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,SICOMO)
!  READ ZONAL SEA ICE VELOCITY                
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,SICUO)
!  READ MERIDIONAL SEA ICE VELOCITY                
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,SICVE)
!  READ SNOW                                       
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,SICSNO)
!  READ HIBLER ETA/ZETA FIELDS                     
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,HIBETE)
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,HIBETO)
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,HIBZETE)
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,HIBZETO)
!  READ VERTICAL DIFFUSIVITY DVO                            
      DO K=1,KEP
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,DVO(:,:,K))
      ENDDO
!  READ VERTICAL FRICTION AVO
      DO K=1,KEP
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,AVO(:,:,K))
      ENDDO
!  READ VERTICAL VELOCITY (DIAGNOSTIC)
      DO K=1,KEP
      IF(p_pe==p_io) READ(IUNIT) I4DATE
      CALL read_slice(IUNIT,WO(:,:,K))
      ENDDO
!
      DO J=1,JE
       DO I=1,IE
        DVO(I,J,1)=0.
        DVO(I,J,KEP)=0.
        AVO(I,J,1)=0.
        AVO(I,J,KEP)=0.
       ENDDO
      ENDDO
!   
      DO 1467 K=1,KE
      DO 1467 J=1,JE
      DO 1467 I=1,IE
      SAO(I,J,K)=MAX(SAO(I,J,K),20.)
      THO(I,J,K)=MAX(THO(I,J,K),-2.)
      SAO(I,J,K)=MIN(SAO(I,J,K),70.)
      THO(I,J,K)=MIN(THO(I,J,K),70.)
1467  CONTINUE
!::
      IF(p_pe==p_io) THEN
        REWIND(IUNIT)
        CLOSE(IO_IN_Z370)
        CLOSE(IO_IN_Z380)
      ENDIF
!
      RETURN
      END
