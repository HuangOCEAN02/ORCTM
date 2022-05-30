      MODULE MO_COMMO2
      USE MO_PARAM1

      IMPLICIT NONE

!==> COMMO2
!UWE***FORCING FIELDS
      REAL, POINTER :: TXO1(:,:),TYE1(:,:),                             &
     &   FCLOU1(:,:),FSWR1(:,:),FU101(:,:),FPREC1(:,:),                 &
     &   FTDEW1(:,:),FSLP1(:,:),TAFO1(:,:)
!<==END COMMO2

      CONTAINS

      SUBROUTINE alloc_mem_commo2

      ALLOCATE(TXO1(I_start:IE,1:JE),TYE1(1:IE,J_start:JE),             &
     &   FCLOU1(IE,JE),FSWR1(IE,JE),FU101(IE,JE),FPREC1(IE,JE),         &
     &   FTDEW1(IE,JE),FSLP1(IE,JE),TAFO1(IE,JE))

      END SUBROUTINE alloc_mem_commo2

      END MODULE MO_COMMO2
