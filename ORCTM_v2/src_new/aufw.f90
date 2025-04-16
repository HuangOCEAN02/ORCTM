      SUBROUTINE AUFW
!****************************************************************
!
!     AAAAAA  U    U  FFFFFF  W    W
!     A    A  U    U  F       W    W
!     AAAAAA  U    U  FFFFF   W WW W
!     A    A  U    U  F       WWWWWW
!     A    A  UUUUUU  F       WW  WW
!
!
!*****************************************************************
! SUBROUTINE AUFW
!
!     THIS SBR WRITES ROUTINELY AT CERTAIN TIME INTERVALS
!     A RESTART FILE ON Z37000 OR Z38000
!
!-----------------------------------------------------------------
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMO2
      USE MO_COMMOAU2
!
#ifdef SLOPECON_ADPO
      USE MO_COMMOBBL
#endif /*SLOPECON_ADPO*/
!::
      USE MO_UNITS

#ifdef OBCS_UV_FLATHER
      USE MO_OBCS_NEST
#endif

      IMPLICIT NONE
      INTEGER(KIND=i4) IDATE,ICODE,KCODE,IEDIM00,IEDIM,IEDIMU,IEDIMV
      INTEGER k, l
      REAL(KIND=sp) rldt, rldays
!
      IF(p_pe==p_io) THEN
        OPEN(IO_IN_Z370,FILE='Z37000',STATUS='UNKNOWN'                  &
     &                 ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

        OPEN(IO_IN_Z380,FILE='Z38000',STATUS='UNKNOWN'                  &
     &                 ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

        REWIND(IO_IN_Z370)
        REWIND(IO_IN_Z380)
      ENDIF
!
      IUNIT=IUNIT+IFLAG
      WRITE(IO_STDOUT,*) 'AUFW: IUNIT IFLAG',IUNIT,IFLAG
      IF(p_pe==p_io) REWIND(IUNIT)
      WRITE(IO_STDOUT,*)' ++++++ WRITING RESTART FILE ON ',IUNIT        &
     &,' AFTER ',LDT                                                    &
     &,' TIME STEPS . CALCULATED YEAR : ',LYEARS,' MONTH : ',LMONTS
!:: WRITE ALL DATA IN EXTRA FORMAT
      IDATE= (LYEARS*10000)+(LMONTS*100)+LDAYS
      RLDT=REAL(LDT)
      RLDAYS=REAL(LDAYS)
      IEDIM00=2
      IEDIM=IE_G*JE_G
      IEDIMU=(IE_G+1)*JE_G
      IEDIMV=IE_G*(JE_G+1)
!:: WRITE TIME STEP INFORMATION
      KCODE=-100
      ICODE=999
      IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM00
      IF(p_pe==p_io) WRITE(IUNIT) RLDT,RLDAYS
!  WRITE ZONAL VELOCITY COMPONENT
      ICODE=3
      DO K=1,KE
       KCODE=TIESTU(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIMU
       CALL write_slice(IUNIT,UKO(:,:,K),1)
      ENDDO
!  WRITE MERIDIONAL  VELOCITY COMPONENT
      ICODE=4
      DO K=1,KE
       KCODE=TIESTU(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIMV
       CALL write_slice(IUNIT,VKE(:,:,K),2)
      ENDDO
!  WRITE TEMPERATURE                     
      ICODE=2
      DO K=1,KE
       KCODE=TIESTU(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,THO(:,:,K))
      ENDDO
!  WRITE SALINITY                        
      ICODE=5
      DO K=1,KE
       KCODE=TIESTU(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,SAO(:,:,K))
      ENDDO
!:: WRITE 2-D FIELDS
      KCODE=-100
!  WRITE SEA SURFACE ELEVATION           
      ICODE=1
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,ZO)
!  WRITE SEA SURFACE ELEVATION CHANGE    
      ICODE=82
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,Z1O)
!  WRITE BAROTROPIC ZONAL VELOCITY
      ICODE=63
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIMU
       CALL write_slice(IUNIT,USO,1)
!  WRITE BAROTROPIC MERIDIONAL VELOCITY
      ICODE=64
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIMV
       CALL write_slice(IUNIT,VSE,2)
!  WRITE SEA ICE THICKNESS               
      ICODE=13
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,SICTHO)
!  WRITE SEA ICE CONCENTRATION           
      ICODE=15
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,SICOMO)
!  WRITE ZONAL SEA ICE VELOCITY                
      ICODE=35
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,SICUO)
!  WRITE MERIDIONAL SEA ICE VELOCITY                
      ICODE=36
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,SICVE)
!  WRITE SNOW                                       
      ICODE=141
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,SICSNO)
!  WRITE HIBLER ETA/ZETA FIELDS                     
      ICODE=501 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,HIBETE)
      ICODE=502 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,HIBETO)
      ICODE=503 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,HIBZETE)
      ICODE=504 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,HIBZETO)
!  WRITE VERTICAL DIFFUSIVITY DVO                            
      ICODE=111
      DO K=1,KEP
       KCODE=TIESTW(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,DVO(:,:,K))
      ENDDO
!  WRITE VERTICAL FRICTION AVO
      ICODE=110
      DO K=1,KEP
       KCODE=TIESTW(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,AVO(:,:,K))
      ENDDO
!  WRITE VERTICAL VELOCITY (DIAGNOSTIC)
      ICODE=7
      DO K=1,KEP
       KCODE=TIESTW(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,WO(:,:,K))
      ENDDO
!:: WRITE 2-D DIAGNOSTIC FIELDS
      KCODE=-100
!:: WRITE MAX. CONVECTION DEPTH    
      ICODE=69 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,FLOAT(KCONDEP(:,:)))
!:: WRITE ZONAL WIND STRESS              
      ICODE=52 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIMU
       CALL write_slice(IUNIT,TXO,1)
!:: WRITE MERIDIONAL WIND STRESS              
      ICODE=53 
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIMV
       CALL write_slice(IUNIT,TYE,2)
!:: WRITE SHORT WAVE RAD FLUX      
      ICODE=176
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,QSWO)
!:: WRITE LONG WAVE RAD FLUX      
      ICODE=177
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,QLWO)
!:: WRITE SENS HEAT FLUX          
      ICODE=146
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,QSEO)
!:: WRITE LATENT HEAT FLUX          
      ICODE=147
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,QLAO)
#ifdef SLOPECON_ADPO
      ICODE=199
      KCODE=-100
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,BBLFLX)
#endif
!:: WRITE GRID INFORMATION AND WET/DRY FIELDS:
!:: WRITE WETO
      ICODE=84
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,DEPTO)
!:: WRITE DLXP
      ICODE=85
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,DLXP)
!:: WRITE DLYP
      ICODE=86
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,DLYP)
!:: WRITE WETO
      ICODE=506
      DO K=1,KE
       KCODE=TIESTU(K)
       IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       CALL write_slice(IUNIT,WETO(:,:,K))
      ENDDO
!::
      IFLAG=-IFLAG
      IF(p_pe==p_io) THEN
        REWIND(IUNIT)
        CLOSE(IO_IN_Z370)
        CLOSE(IO_IN_Z380)
      ENDIF
!
      RETURN
      END
