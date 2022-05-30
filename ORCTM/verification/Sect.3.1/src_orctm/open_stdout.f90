      SUBROUTINE OPEN_STDOUT(kout,yfilename)

      USE MO_PARALLEL

!----------------------------------------------------------------------
!
!**** *OPEN_OCEOUT*    -  assign ocean model output file.
!
!     PURPOSE.
!     -------
!     
!
!
!     INTERFACE.
!     ----------
!
!----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) yfilename
      CHARACTER*256 outfile
      INTEGER       kout, iost
!
!
!*      0.    Assign output file.
!             -------------------

!  For the root processor: Write to yfilename
!  For all others: Write to yfilename with processor number appended

      IF(p_pe==p_io) THEN
         outfile = yfilename
      ELSE
         WRITE(outfile,'(A,"_",I3.3)')                                  &
     &                  yfilename(1:LEN_TRIM(yfilename)),p_pe
      ENDIF

      iost = 0

      OPEN (kout,FILE = outfile                                         &
     &     ,STATUS='UNKNOWN',FORM ='FORMATTED',ERR = 110,IOSTAT = iost)

 110  CONTINUE
      IF (iost .ne. 0) THEN
          PRINT*, '        ***WARNING***'
          PRINT*,' ===>>> : ERROR opening output FILE'
          PRINT*,' ======   =====                ===='
          PRINT*,' Logical unit ',kout,' error number = ',iost
          PRINT*, ' '
          PRINT*,' We STOP!!! Verify the file ',outfile
          PRINT*, ' '
          CALL STOP_ALL('STOP')
      ENDIF

      RETURN
      END
