      MODULE MO_OCTDIFF

      USE mo_param1

      IMPLICIT NONE

!:: THREE-DIMENSIONAL FIELDS
       REAL,POINTER ::                                                  &
     &   vol_term(:,:,:),                                               &
     &   SLOPLOy(:,:,:), SLOPLUy(:,:,:), SLOPROy(:,:,:),                &
     &   SLOPRUy(:,:,:), SLOPLOx(:,:,:), SLOPLUx(:,:,:),                &
     &   SLOPROx(:,:,:), SLOPRUx(:,:,:),                                &
     &   XCOEFF_LO(:,:,:),XCOEFF_LU(:,:,:),XCOEFF_RO(:,:,:),            &
     &   XCOEFF_RU(:,:,:),ZCOEFF_LOy(:,:,:),ZCOEFF_LUy(:,:,:),          &
     &   ZCOEFF_ROy(:,:,:),ZCOEFF_RUy(:,:,:),ZCOEFF_LOx(:,:,:),         &
     &   ZCOEFF_LUx(:,:,:),ZCOEFF_ROx(:,:,:),ZCOEFF_RUx(:,:,:),         &
     &   YCOEFF_LO(:,:,:),YCOEFF_LU(:,:,:),YCOEFF_RO(:,:,:),            &
     &   YCOEFF_RU(:,:,:),XFLUX(:,:,:),YFLUX(:,:,:),TRIDSY(:,:,:,:)
!
!:: TWO-DIMENSIONAL FIELDS
      REAL,POINTER :: AULX(:,:),AULY(:,:),ZBOTT(:,:),ZSURF(:,:)

      INTEGER,POINTER :: KO(:,:),KU(:,:) 
!
      CONTAINS

      SUBROUTINE alloc_mem_octdiff

!:: THREE-DIMENSIONAL FIELDS
       ALLOCATE(                                                        &
     &   vol_term(IE,JE,KE),                                            &
     &   SLOPLOy(IE,JE,KE), SLOPLUy(IE,JE,KE), SLOPROy(IE,JE,KE),       &
     &   SLOPRUy(IE,JE,KE), SLOPLOx(IE,JE,KE), SLOPLUx(IE,JE,KE),       &
     &   SLOPROx(IE,JE,KE), SLOPRUx(IE,JE,KE),                          &
     &   XCOEFF_LO(IE,JE,KE),XCOEFF_LU(IE,JE,KE),XCOEFF_RO(IE,JE,KE),   &
     &   XCOEFF_RU(IE,JE,KE),ZCOEFF_LOy(IE,JE,KE),ZCOEFF_LUy(IE,JE,KE), &
     &   ZCOEFF_ROy(IE,JE,KE),ZCOEFF_RUy(IE,JE,KE),ZCOEFF_LOx(IE,JE,KE),&
     &   ZCOEFF_LUx(IE,JE,KE),ZCOEFF_ROx(IE,JE,KE),ZCOEFF_RUx(IE,JE,KE),&
     &   YCOEFF_LO(IE,JE,KE),YCOEFF_LU(IE,JE,KE),YCOEFF_RO(IE,JE,KE),   &
     &   YCOEFF_RU(IE,JE,KE),XFLUX(IE,JE,KE),YFLUX(IE,JE,KE),           &
     &   TRIDSY(IE,JE,KE,3)                                             &
     & )
!
!:: TWO-DIMENSIONAL FIELDS
       ALLOCATE(AULX(IE,JE),AULY(IE,JE),ZBOTT(JE,KE),ZSURF(JE,KE),      &
     &             KO(JE,KE),KU(JE,KE) )

      END SUBROUTINE alloc_mem_octdiff

      END MODULE MO_OCTDIFF
