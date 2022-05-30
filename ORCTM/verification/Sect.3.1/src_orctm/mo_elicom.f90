      MODULE MO_ELICOM
      USE MO_PARAM1

      IMPLICIT NONE

      REAL, POINTER :: PGL(:,:)

      CONTAINS

      SUBROUTINE alloc_mem_elicom

      ALLOCATE (PGL(IMM,ILL))

      END SUBROUTINE alloc_mem_elicom
      END MODULE MO_ELICOM
