      SUBROUTINE AUFW_FRAGMENT
#ifdef FRAGMENT
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
      character(len=20) str
!

      write(str,'("_",I3.3)') p_pe
      OPEN(IO_IN_Z370,FILE='./restart_file/Z37000'//trim(adjustl(str)),    &
     &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      OPEN(IO_IN_Z380,FILE='./restart_file/Z38000'//trim(adjustl(str)),    &
     &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      REWIND(IO_IN_Z370)
      REWIND(IO_IN_Z380)

!
      IUNIT=IUNIT+IFLAG
      WRITE(IO_STDOUT,*) 'AUFW: IUNIT IFLAG',IUNIT,IFLAG
      REWIND(IUNIT)
      WRITE(IO_STDOUT,*)' ++++++ WRITING RESTART FILE ON ',IUNIT        &
     &,' AFTER ',LDT                                                    &
     &,' TIME STEPS . CALCULATED YEAR : ',LYEARS,' MONTH : ',LMONTS
!:: WRITE ALL DATA IN EXTRA FORMAT
      IDATE= (LYEARS*10000)+(LMONTS*100)+LDAYS
      RLDT=REAL(LDT)
      RLDAYS=REAL(LDAYS)
      IEDIM00=2
      IEDIM=IE*JE
      
      if (have_g_is) then
         IEDIMU=(IE+1)*JE
      else
         IEDIMU= IE*JE
      ENDIF
      
      if (have_g_js) then
         IEDIMV=IE*(JE+1)
      else
         IEDIMV= IE*JE
      ENDIF      

!:: WRITE TIME STEP INFORMATION
      KCODE=-100
      ICODE=999
      WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM00
      WRITE(IUNIT) RLDT,RLDAYS
!  WRITE ZONAL VELOCITY COMPONENT
      ICODE=3
      DO K=1,KE
       KCODE=TIESTU(K)
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIMU
       WRITE(IUNIT) real(UKO(:,:,k),sp)
      ENDDO
!  WRITE MERIDIONAL  VELOCITY COMPONENT
      ICODE=4
      DO K=1,KE
       KCODE=TIESTU(K)
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIMV
       WRITE(IUNIT) real(VKE(:,:,k),sp)       
      ENDDO
!  WRITE TEMPERATURE                     
      ICODE=2
      DO K=1,KE
       KCODE=TIESTU(K)
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(THO(:,:,k),sp)        
      ENDDO
!  WRITE SALINITY                        
      ICODE=5
      DO K=1,KE
       KCODE=TIESTU(K)
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(SAO(:,:,k),sp)        
      ENDDO
!:: WRITE 2-D FIELDS
      KCODE=-100
!  WRITE SEA SURFACE ELEVATION           
      ICODE=1
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(ZO,sp)        
!  WRITE SEA SURFACE ELEVATION CHANGE    
      ICODE=82
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(Z1O,sp)           
!  WRITE BAROTROPIC ZONAL VELOCITY
      ICODE=63
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIMU
       WRITE(IUNIT) real(USO,sp)           
!  WRITE BAROTROPIC MERIDIONAL VELOCITY
      ICODE=64
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIMV
       WRITE(IUNIT) real(VSE,sp)           
!  WRITE SEA ICE THICKNESS               
      ICODE=13
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(SICTHO,sp)        
!  WRITE SEA ICE CONCENTRATION           
      ICODE=15
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(SICOMO,sp)        
!  WRITE ZONAL SEA ICE VELOCITY                
      ICODE=35
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(SICUO,sp)        
!  WRITE MERIDIONAL SEA ICE VELOCITY                
      ICODE=36
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(SICVE,sp)        
!  WRITE SNOW                                       
      ICODE=141
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(SICSNO,sp)        
!  WRITE HIBLER ETA/ZETA FIELDS                     
      ICODE=501 
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(HIBETE,sp)        
      ICODE=502 
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(HIBETO,sp)        
      ICODE=503 
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(HIBZETE,sp)        
      ICODE=504 
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(HIBZETO,sp)        
!  WRITE VERTICAL DIFFUSIVITY DVO                            
      ICODE=111
      DO K=1,KEP
       KCODE=TIESTW(K)
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(DVO(:,:,K),sp)        
      ENDDO
!  WRITE VERTICAL FRICTION AVO
      ICODE=110
      DO K=1,KEP
       KCODE=TIESTW(K)
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(AVO(:,:,K),sp)        
      ENDDO
!  WRITE VERTICAL VELOCITY (DIAGNOSTIC)
      ICODE=7
      DO K=1,KEP
       KCODE=TIESTW(K)
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(WO(:,:,K),sp)        
      ENDDO
!:: WRITE 2-D DIAGNOSTIC FIELDS
      KCODE=-100
!:: WRITE MAX. CONVECTION DEPTH    
      ICODE=69 
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(FLOAT(KCONDEP(:,:)),sp)        
!:: WRITE ZONAL WIND STRESS              
      ICODE=52 
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIMU
       WRITE(IUNIT) real(TXO,sp)        
!:: WRITE MERIDIONAL WIND STRESS              
      ICODE=53 
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIMV
       WRITE(IUNIT) real(TYE,sp)        
!:: WRITE SHORT WAVE RAD FLUX      
      ICODE=176
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(QSWO,sp)        
!:: WRITE LONG WAVE RAD FLUX      
      ICODE=177
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(QLWO,sp)        
!:: WRITE SENS HEAT FLUX          
      ICODE=146
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(QSEO,sp)        
!:: WRITE LATENT HEAT FLUX          
      ICODE=147
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(QLAO,sp)        

#ifdef SLOPECON_ADPO
      ICODE=199
      KCODE=-100
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(BBLFLX,sp)        

#endif
!:: WRITE GRID INFORMATION AND WET/DRY FIELDS:
!:: WRITE WETO
      ICODE=84
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(DEPTO,sp)         
!:: WRITE DLXP
      ICODE=85
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(DLXP,sp)         
!:: WRITE DLYP
      ICODE=86
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(DLYP,sp)         
!:: WRITE WETO
      ICODE=506
      DO K=1,KE
       KCODE=TIESTU(K)
       WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
       WRITE(IUNIT) real(WETO(:,:,K),sp)         
      ENDDO
!::
      IFLAG=-IFLAG

      REWIND(IUNIT)
      CLOSE(IO_IN_Z370)
      CLOSE(IO_IN_Z380)

!
#endif
      RETURN
      END SUBROUTINE AUFW_FRAGMENT
