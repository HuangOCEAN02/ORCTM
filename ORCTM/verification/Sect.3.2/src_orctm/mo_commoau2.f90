      MODULE MO_COMMOAU2
      USE MO_PARAM1
      IMPLICIT NONE
!==>COMMOAU2
      REAL, POINTER ::                                                  &
     &   TICEO(:,:),ACLO(:,:),PAO(:,:),RPRECO(:,:),FRSE(:,:),           &
     &   PRECO(:,:),QSWO(:,:),QLWO(:,:),QSEO(:,:),                      &
     &   QLAO(:,:),TAIRO(:,:),TDO(:,:),                                 &
     &   SICUDO(:,:),SICVDE(:,:),PRECH(:,:)
!<==END COMMOAU2

      CONTAINS

      SUBROUTINE alloc_mem_commoau2

      ALLOCATE(                                                         &
     &  TICEO(IE,JE),ACLO(IE,JE),PAO(IE,JE),RPRECO(IE,JE),FRSE(IE,JE)   &
     & ,PRECO(IE,JE),QSWO(IE,JE),QLWO(IE,JE),QSEO(IE,JE)                &
     & ,QLAO(IE,JE),TAIRO(IE,JE),TDO(IE,JE)                             &
     & ,SICUDO(IE,JE),SICVDE(IE,JE),PRECH(IE,JE))

      END SUBROUTINE alloc_mem_commoau2
      END MODULE MO_COMMOAU2
