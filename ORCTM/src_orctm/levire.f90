      SUBROUTINE LEVIRE(LLL)
!
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_UNITS
      USE MO_LEVITUS

      IMPLICIT NONE
      INTEGER lll,i,j,k,ierr
      REAL SAOMAX0, SAOMAX1, SAOMIN0, SAOMIN1
      REAL THOMAX0, THOMAX1, THOMIN0, THOMIN1

!
      WRITE(IO_STDOUT,*) 'CALLING LEVIRE MIT ',LLL,' ...'
      IF (LLL .EQ.  0 ) THEN
         WRITE(IO_STDOUT,*)                                             &
     &       '==> HORIZONTAL INITIAL STRATIFICATION'
      ELSEIF (LLL .EQ. 13 ) THEN
         WRITE(IO_STDOUT,*)                                             &
     &       '==> READING INITIAL 3D STRATIFICATION'
      ELSEIF (LLL .EQ. -1 ) THEN
         WRITE(IO_STDOUT,*)                                             &
     &       '==> READING LEVITUS STRATIFICATION FOR 3D RESTORING.'
      ELSEIF (LLL .EQ. -2 ) THEN
         WRITE(IO_STDOUT,*)                                             &
     &       '==> READING SURFACE SALINITY'
      ELSEIF (LLL .EQ. -3 ) THEN
         WRITE(IO_STDOUT,*)                                             &
     &       '==> READING SURFACE TEMPERATURE'
      ELSE
         WRITE(IO_STDOUT,*)                                             &
     &       '==> INVALID PARAMETER IN CALL OF LEVIRE'
         WRITE(IO_STDOUT,*) 'STOP IN LEVIRE '
         CALL ABSTURZ
      ENDIF
!
!
! 3D RESTORING
!
      IF( LLL .EQ. -1 ) THEN
      IF (I3DREST .GT. 0) THEN
      DO K=1,KE
       IF(p_pe==p_io) THEN
        READ(IO_IN_INIS,IOSTAT=IERR) IBLA
        IF (IERR.GT.0) THEN
         WRITE(IO_STDOUT,*)'NO MONTHLY LEVITUS FIELDS FOUND IN MONTH= ' &
     &                     ,LMONTS
         WRITE(IO_STDOUT,*)                                             &
     &      'CHECK OPTION DREIDREST_MON AND FILES INISAL/TEM !!! '
         CALL ABSTURZ
        ENDIF
       ENDIF
        CALL read_slice(IO_IN_INIS,SLEVI(:,:,K))
        IF(p_pe==p_io) READ(IO_IN_INIT) ibla
        CALL read_slice(IO_IN_INIT,TLEVI(:,:,K))
       ENDDO
         DO K=1,KE
            DO J=1,JE
               DO I=1,IE
                  LSLEV(I,J,K)=.TRUE.
                  LTLEV(I,J,K)=.TRUE.
                  IF(ABS(SLEVI(I,J,K)).GT.40..OR.WETO(I,J,K).LT.0.5)    &
     &                LSLEV(I,J,K)=.FALSE.
                  IF(ABS(TLEVI(I,J,K)).GT.40..OR.WETO(I,J,K).LT.0.5)    &
     &                LTLEV(I,J,K)=.FALSE.
               ENDDO
            ENDDO
         ENDDO
         WRITE(IO_STDOUT,*)                                             &
     &   'LEVITUS DATA READ TO SLEVI/TLEVI FOR 3D RESTORING IN MONTH',  &
     &    LMONTS
!
         ELSE !I3DREST
!
         WRITE(IO_STDOUT,*)                                             &
     &   'LEVIRE IS CALLED WITH LLL = ',LLL,                            &
     &   ' BUT CODE IS NOT ACTIVATED. THE PROGRAM IS STOPPED.'
         STOP 'STOP IN LEVIRE: CODE NOT ACTIVATED'
!
      ENDIF !I3DREST
      ENDIF
!
! HORIZONTAL STRATIFICATION
!
      IF( LLL .EQ. 0 ) THEN
         DO K=1,KE
            DO J=1,JE
               DO I=1,IE
                  SAO(I,J,K)=SAF(K)
                  THO(I,J,K)=TAF(K)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
!
! 3D STRATIFICATION
!
      IF( LLL .EQ. 13 ) THEN
         IF(p_pe==p_io) REWIND(IO_IN_INIT)
         IF(p_pe==p_io) REWIND(IO_IN_INIS)
         DO K=1,KE
           IF(p_pe==p_io) READ(IO_IN_INIS)
           CALL read_slice(IO_IN_INIS,SAO(:,:,K))
           IF(p_pe==p_io) READ(IO_IN_INIT)
           CALL read_slice(IO_IN_INIT,THO(:,:,K))
         ENDDO
!
         SAOMAX0 = -1.0
         SAOMIN0 = 99.0
         SAOMAX1 = -1.0
         SAOMIN1 = 99.0
         THOMAX0 = -1.0
         THOMIN0 = 99.0
         THOMAX1 = -1.0
         THOMIN1 = 99.0
         DO K=1,KE
            DO J=1,JE
               DO I=1,IE
                  IF( WETO(I,J,K) .LT. 0.5) THEN
                     SAOMAX0 = MAX( SAOMAX0, SAO(I,J,K))
                     SAOMIN0 = MIN( SAOMIN0, SAO(I,J,K))
                     THOMAX0 = MAX( THOMAX0, THO(I,J,K))
                     THOMIN0 = MIN( THOMIN0, THO(I,J,K))
!                    SAO(I,J,K) = 999.99
!                    THO(I,J,K) = 999.99
                  ELSE
                     SAOMAX1 = MAX( SAOMAX1, SAO(I,J,K))
                     SAOMIN1 = MIN( SAOMIN1, SAO(I,J,K))
                     THOMAX1 = MAX( THOMAX1, THO(I,J,K))
                     THOMIN1 = MIN( THOMIN1, THO(I,J,K))
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         CALL global_max(SAOMAX0,SAOMAX1,THOMAX0,THOMAX1)
         CALL global_min(SAOMIN0,SAOMIN1,THOMIN0,THOMIN1)
         WRITE(IO_STDOUT,*) 'MIN/MAX AT DRY CELLS (S): ',SAOMIN0,SAOMAX0
         WRITE(IO_STDOUT,*) 'MIN/MAX AT WET CELLS (S): ',SAOMIN1,SAOMAX1
         WRITE(IO_STDOUT,*) 'MIN/MAX AT DRY CELLS (T): ',THOMIN0,THOMAX0
         WRITE(IO_STDOUT,*) 'MIN/MAX AT WET CELLS (T): ',THOMIN1,THOMAX1
!
      ENDIF
!
!
! SURFACE SALINITY
!
      IF( LLL .EQ. -2 ) THEN
          IF(p_pe==p_io) THEN
           READ(IO_IN_SURS,IOSTAT=IERR)IBLA
           IF (IERR.GT.0) THEN
            WRITE(IO_STDOUT,*)'NO MONTHLY SSS FIELDS FOUND IN MONTH= '  &
     &                        ,LMONTS
            WRITE(IO_STDOUT,*)                                          &
     &            'CHECK OPTION RESTORE_MON AND FILE SURSAL!!! '
            CALL ABSTURZ
           ENDIF
          ENDIF
          CALL read_slice(IO_IN_SURS,RELSAO)
      ENDIF

! SURFACE TEMPERATURE
      IF( LLL .EQ. -3 ) THEN
          IF(p_pe==p_io) THEN
           READ(IO_IN_SURT,IOSTAT=IERR)IBLA
           IF (IERR.GT.0) THEN
            WRITE(IO_STDOUT,*)'NO MONTHLY SST FIELDS FOUND IN MONTH= '  &
     &                        ,LMONTS
            WRITE(IO_STDOUT,*)                                          &
     &            'CHECK OPTION RESTORE_MON AND FILE SURTEM!!! '
            CALL ABSTURZ
           ENDIF
          ENDIF
          CALL read_slice(IO_IN_SURT,RELTHO)
      ENDIF
!
      RETURN
      END

