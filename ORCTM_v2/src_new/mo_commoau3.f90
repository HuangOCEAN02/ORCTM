      MODULE MO_COMMOAU3
!==>COMMOAU3
!     S. Legutke               *MPI-MaD, HH*          03.08.01
!
!
!*    *COMMON* ** -  control variables.
!
!*    VARIABLE    TYPE     PURPOSE.
!     --------    ----     --------
!     *brine*    *REAL*    - fr.flux due to th.dyn. ice growth [m/sec].
!
! ---------------------------------------------------------------------
!
      USE MO_PARAM1

      IMPLICIT NONE

      REAL, POINTER :: brine(:,:)
!<==END COMMOAU3

      CONTAINS

      SUBROUTINE alloc_mem_commoau3

      ALLOCATE(brine(ie,je))

      END SUBROUTINE alloc_mem_commoau3
      END MODULE MO_COMMOAU3
