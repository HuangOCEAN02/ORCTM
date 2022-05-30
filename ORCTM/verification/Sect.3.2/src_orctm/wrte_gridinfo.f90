      SUBROUTINE WRTE_GRIDINFO
!
! Modified by Peterspy,referred to aufw.f90
!
!**************************************************************************

      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS
      USE MO_KIND

      INTEGER(KIND=i4) IDATE,ICODE,KCODE,IEDIM,IEDIMU,IEDIMV
      REAL   rad,    LONP_G(IE_G , JE_G),   LATP_G(IE_G , JE_G),        &
     &               LONU_G(0:IE_G,JE_G),   LATU_G(0:IE_G,JE_G),        &
     &               LONV_G(IE_G,0:JE_G),   LATV_G(IE_G,0:JE_G)
!
! Write to disk (EXTRA format)
!
      IDATE=0
      KCODE=-100
      IEDIM=(IE_G*JE_G)
      IEDIMU=(IE_G+1)*JE_G
      IEDIMV=IE_G*(JE_G+1)

      RAD=PI/180.
      DO I=1,IE_G
        DO J=1,JE_G
          LONP_G(I,J)=GILA_G(2*I,2*J)/RAD
          LATP_G(I,J)=GIPH_G(2*I,2*J)/RAD
        ENDDO
      ENDDO
      DO I=0,IE_G
        DO J=1,JE_G
          LONU_G(I,J)=GILA_G(2*I+1,2*J)/RAD
          LATU_G(I,J)=GIPH_G(2*I+1,2*J)/RAD
        ENDDO
      ENDDO
      DO I=1,IE_G
        DO J=0,JE_G
          LONV_G(I,J)=GILA_G(2*I,2*J+1)/RAD
          LATV_G(I,J)=GIPH_G(2*I,2*J+1)/RAD
        ENDDO
      ENDDO

!:: WRITE GRID INFORMATION AND WET/DRY FIELDS:
      IF(p_pe==p_io) THEN
        OPEN(IO_OU_GRID,FILE='GRID_INFO.ext',STATUS='UNKNOWN'           &
     &                 ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      ENDIF
!
!:: PETERSPY: ICODE (ABC)
!     A: 1-> lon/lat 2 -> water depth; 3 -> sea mask; 4 -> grid interval
!     B: 0 -> no direction; 1/2/3 -> X/Y/Z direction
!     C: 0/1/2 -> variable locates at P/U/V point
!
!::   WRITE LONP
        ICODE=110
        IF (p_pe==p_io) THEN
          WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIM
          WRITE(IO_OU_GRID) real(LONP_G,sp)
        ENDIF
!::   WRITE LATP
        ICODE=120
        IF(p_pe==p_io) THEN
          WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIM
          WRITE(IO_OU_GRID) real(LATP_G,sp)
        ENDIF
!::   WRITE LONU
        ICODE=111
        IF (p_pe==p_io) THEN
          WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMU
          WRITE(IO_OU_GRID) real(LONU_G,sp)
        ENDIF
!::   WRITE LATU
        ICODE=121
        IF(p_pe==p_io) THEN
          WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMU
          WRITE(IO_OU_GRID) real(LATU_G,sp)
        ENDIF
!::   WRITE LONV
        ICODE=112
        IF (p_pe==p_io) THEN
          WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMV
          WRITE(IO_OU_GRID) real(LONV_G,sp)
        ENDIF
!::   WRITE LATV
        ICODE=122
        IF(p_pe==p_io) THEN
          WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMV
          WRITE(IO_OU_GRID) real(LATV_G,sp)
        ENDIF
!::   WRITE DEPTO
        ICODE=200
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIM
          CALL write_slice(IO_OU_GRID,DEPTO)
!::   WRITE DEUTO
        ICODE=201
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMU
          CALL write_slice(IO_OU_GRID,DEUTO,1)
!::   WRITE DEUTE
        ICODE=202
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMV
          CALL write_slice(IO_OU_GRID,DEUTE,2)
!::   WRITE DLXP
        ICODE=410
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIM
          CALL write_slice(IO_OU_GRID,DLXP)
!::   WRITE DLXU
        ICODE=411
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMU
          CALL write_slice(IO_OU_GRID,DLXU,1)
!::   WRITE DLXV
        ICODE=412
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMV
          CALL write_slice(IO_OU_GRID,DLXV,2)
!::   WRITE DLYP
        ICODE=420
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIM
          CALL write_slice(IO_OU_GRID,DLYP)
!::   WRITE DLYU
        ICODE=421
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMU
          CALL write_slice(IO_OU_GRID,DLYU,1)
!::   WRITE DLYV
        ICODE=422
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMV
          CALL write_slice(IO_OU_GRID,DLYV,2)
!::   WRITE WETO
        ICODE=300
        DO K=1,KE
          KCODE=K
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIM
          CALL write_slice(IO_OU_GRID,WETO(:,:,K))
        ENDDO
!::   WRITE AMSUO
        ICODE=301
        DO K=1,KE
          KCODE=K
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMU
          CALL write_slice(IO_OU_GRID,AMSUO(:,:,K),1)
        ENDDO
!::   WRITE AMSUE
        ICODE=302
        DO K=1,KE
          KCODE=K
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMV
          CALL write_slice(IO_OU_GRID,AMSUE(:,:,K),2)
        ENDDO
!::   WRITE DDPO
        ICODE=430
        DO K=1,KE
          KCODE=K
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIM
          CALL write_slice(IO_OU_GRID,DDPO(:,:,K))
        ENDDO
!::   WRITE DDUO
        ICODE=431
        DO K=1,KE
          KCODE=K
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMU
          CALL write_slice(IO_OU_GRID,DDUO(:,:,K),1)
        ENDDO
!::   WRITE DDUE
        ICODE=432
        DO K=1,KE
          KCODE=K
          IF(p_pe==p_io) WRITE(IO_OU_GRID) IDATE,ICODE,KCODE,IEDIMV
          CALL write_slice(IO_OU_GRID,DDUE(:,:,K),2)
        ENDDO

      IF(p_pe==p_io) CLOSE(IO_OU_GRID)

      RETURN
      END
