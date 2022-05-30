      SUBROUTINE AUFR_FRAGMENT
#ifdef FRAGMENT
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
!
      IMPLICIT NONE
      INTEGER(KIND=i4) I4DATE(4)
      REAL(KIND=sp) RDT27, RDT28, RDT(2)
      INTEGER i, j, k, IDATE
      character(len=20) str
      
        write(str,'("_",I3.3)') p_pe
        OPEN(IO_IN_Z370,FILE='../restart_file/Z37000'//trim(adjustl(str)),    &
     &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

        OPEN(IO_IN_Z380,FILE='../restart_file/Z38000'//trim(adjustl(str)),    &
     &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
     
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


!      CALL p_bcast(IFLAG,p_io)
!      CALL p_bcast(IUNIT,p_io)
!      CALL p_bcast(IDATE,p_io)
!      CALL p_bcast(LDT,p_io)

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

!     CALL Sleep(5)
!     CALL ABSTURZ
!CCC
!  READ ZONAL VELOCITY COMPONENT
      DO K=1,KE
        READ(IUNIT) I4DATE
        READ(IUNIT) DUMMY_U
        UOO(:,:,K) = DUMMY_U
      ENDDO

!  READ MERIDIONAL  VELOCITY COMPONENT
      DO K=1,KE
        READ(IUNIT) I4DATE
        READ(IUNIT) DUMMY_V
        VOE(:,:,K) = DUMMY_V
      ENDDO
!:: UPDATE VELOCITY FIELDS
      CALL OCTIMF
!  READ TEMPERATURE
      DO K=1,KE
      	READ(IUNIT) I4DATE
        READ(IUNIT) DUMMY_P
        THO(:,:,K) = DUMMY_P
      ENDDO
!  READ SALINITY
      DO K=1,KE
        READ(IUNIT) I4DATE
        READ(IUNIT) DUMMY_P
        SAO(:,:,K) = DUMMY_P
      ENDDO

!      CALL Sleep(5)
!      CALL ABSTURZ
!:: READ 2-D FIELDS
!  READ SEA SURFACE ELEVATION
      READ(IUNIT) I4DATE
      READ(IUNIT) DUMMY_P
      ZO = DUMMY_P      
       
!  READ SEA SURFACE ELEVATION CHANGE
      READ(IUNIT) I4DATE
      READ(IUNIT) DUMMY_P
      Z1O = DUMMY_P         

!  READ SEA ICE THICKNESS
      READ(IUNIT) I4DATE
      READ(IUNIT) DUMMY_P
      SICTHO = DUMMY_P      

!  READ SEA ICE CONCENTRATION         
      READ(IUNIT) I4DATE
      READ(IUNIT) DUMMY_P
      SICOMO = DUMMY_P         

!  READ ZONAL SEA ICE VELOCITY                
      READ(IUNIT) I4DATE
      READ(IUNIT) DUMMY_P 
      SICUO = DUMMY_P       

!  READ MERIDIONAL SEA ICE VELOCITY                
      READ(IUNIT) I4DATE
      READ(IUNIT) DUMMY_P
      SICVE = DUMMY_P

!  READ SNOW                                       
      READ(IUNIT) I4DATE
      READ(IUNIT) DUMMY_P 
      SICSNO = DUMMY_P        

!  READ HIBLER ETA/ZETA FIELDS                     
      READ(IUNIT) I4DATE
      READ(IUNIT) DUMMY_P
      HIBETE = DUMMY_P          

      READ(IUNIT) I4DATE
      READ(IUNIT) DUMMY_P
      HIBETO = DUMMY_P          

      READ(IUNIT) I4DATE
      READ(IUNIT) DUMMY_P
      HIBZETE = DUMMY_P          

      READ(IUNIT) I4DATE
      READ(IUNIT) DUMMY_P
      HIBZETO = DUMMY_P  

!  READ VERTICAL DIFFUSIVITY DVO                            
      DO K=1,KEP
        READ(IUNIT) I4DATE
        READ(IUNIT) DUMMY_P
        DVO(:,:,K) = DUMMY_P  
      ENDDO
!  READ VERTICAL FRICTION AVO
      DO K=1,KEP
        READ(IUNIT) I4DATE
        READ(IUNIT) DUMMY_P
        AVO(:,:,K) = DUMMY_P         
      ENDDO
!  READ VERTICAL VELOCITY (DIAGNOSTIC)
      DO K=1,KEP
        READ(IUNIT) I4DATE
        READ(IUNIT) DUMMY_P
        WO(:,:,K) = DUMMY_P       
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

        REWIND(IUNIT)
        CLOSE(IO_IN_Z370)
        CLOSE(IO_IN_Z380)

!
#endif
      RETURN
      END SUBROUTINE AUFR_FRAGMENT
