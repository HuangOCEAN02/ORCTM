      SUBROUTINE ABSTURZ
!
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
      IMPLICIT NONE
!
      WRITE(IO_STDOUT,*)'IN ABSTURZ: INTENDED ERROR EXIT!!!!'
      CALL STOP_ALL('IN ABSTURZ: INTENDED ERROR EXIT!!!!')

      RETURN
      END