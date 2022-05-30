      MODULE MO_COMMOBBl
      USE MO_PARAM1
      IMPLICIT NONE

      REAL,POINTER::alpx(:,:),alpy(:,:),bblflx(:,:),ubbl(:,:),vbbl(:,:)
      INTEGER,POINTER::kupubbl(:,:),kdwubbl(:,:),kupvbbl(:,:),          &
     &      kdwvbbl(:,:)

      CONTAINS

      SUBROUTINE alloc_mem_commobbl

      ALLOCATE(alpx(IE,JE),alpy(IE,JE),bblflx(IE,JE),ubbl(IE,JE),       &
     &      vbbl(IE,JE),kupubbl(IE,JE),kdwubbl(IE,JE),kupvbbl(IE,JE),   &
     &      kdwvbbl(IE,JE) )

      END SUBROUTINE alloc_mem_commobbl
      END MODULE MO_COMMOBBl
