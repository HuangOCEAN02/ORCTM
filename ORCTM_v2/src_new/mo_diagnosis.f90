MODULE MO_DIAGNOSIS

#ifdef DIAGNOSIS    
  USE MO_PARAM1
  USE MO_MPI
  USE MO_PARALLEL
  USE MO_COMMO1
  USE MO_KIND
  
  IMPLICIT NONE      
! OCECTL Namelist variables in OCEDIAG      
  INTEGER::DMEAN
  INTEGER,PARAMETER::POINT=10
  INTEGER,POINTER :: PE_POINT(:)  ! PROCESSOR NUMBER
  INTEGER,POINTER :: i_point(:),j_point(:)
  INTEGER,POINTER :: i_point_G(:),j_point_G(:)
! Diagnostic output IO-NUMBER in diagnosis subdirectory
  INTEGER,PARAMETER::                                   &
         &   IO_DIAG_Uu=73,IO_DIAG_Vv=74,IO_DIAG_Ww=146 &  
         &  ,IO_DIAG_Tt=71,IO_DIAG_Ss=72,IO_DIAG_Po=76  &
#ifdef DIFFDIAG  
         &  ,IO_DIAG_Av=144,IO_DIAG_Dv=145              &
#endif 
#ifdef BAROTROPIC
         &  ,IO_DIAG_US=63,IO_DIAG_VS=64                &
#endif
#ifdef NON_HYDROSTATIC       
         &  ,IO_DIAG_Pn=77,IO_DIAG_Div=78               &
#endif       
         &  ,IO_DIAG_Zz=82              
  REAL,POINTER::Uu_point(:,:),Vv_point(:,:),Ww_point(:,:), &
     &          Tt_point(:,:),Ss_point(:,:),Po_point(:,:), &
#ifdef DIFFDIAG      
     &          Av_point(:,:),dv_point(:,:),               &
#endif 
#ifdef BAROTROPIC
     &          Us_point(:),Vz_point(:),                   &
#endif
#ifdef NON_HYDROSTATIC       
     &          Pn_point(:,:),div_point(:,:),              & 
#endif               
     &          Zz_point(:)
     
     


  CONTAINS
         
    SUBROUTINE alloc_mem_diagnosis   

      ALLOCATE(PE_POINT(POINT))
      ALLOCATE(i_point(POINT),j_point(POINT))      
      allocate(i_point_G(POINT),j_point_G(POINT))
      allocate(Uu_point(POINT,KE),Vv_point(POINT,KE),Ww_point(POINT,KEP) &
      &       ,Tt_point(POINT,KE),Ss_point(POINT,KE),Po_point(POINT,KE) &
      
#ifdef DIFFDIAG      
      &       ,Av_point(POINT,KEP),Dv_point(POINT,KEP) &
#endif   
#ifdef BAROTROPIC
      &       ,Us_point(POINT),Vz_point(POINT)         &
#endif
#ifdef NON_HYDROSTATIC      
      &       ,Pn_point(POINT,KE),div_point(POINT,KE)  &
#endif      

      &       ,Zz_point(POINT))
      
        
      Uu_point(:,:)=0.0
      Vv_point(:,:)=0.0
      Zz_point(:)=0.0
#ifdef BAROTROPIC      
      Us_point(:)=0.0
      Vz_point(:)=0.0
#endif      
      Tt_point(:,:)=0.0
      Ss_point(:,:)=0.0
      Po_point(:,:)=0.0               
      Ww_point(:,:)=0.0
      
#ifdef NON_HYDROSTATIC       
      Pn_point(:,:)=0.0
      Div_point(:,:)=0.0      
#endif 
#ifdef DIFFDIAG     
      Av_point(:,:)=0.0
      Dv_point(:,:)=0.0 
#endif                             


  
    END SUBROUTINE alloc_mem_diagnosis     
            
 
    SUBROUTINE INI_diagnosis
    
    INTEGER::P,N
! find the processor number and processor's index 
    if (p_pe == p_io) THEN
     DO N = 1, POINT
      DO P = 0, nprocxy-1
          IF (i_point_G(N) >= p_Ioff_GG(P)+1 .AND. i_point_G(N) < p_Ioff_GG(P)+ie ) THEN
            IF (j_point_G(N) >= p_Joff_GG(P)+1 .AND. j_point_G(N) < p_Joff_GG(P)+je ) THEN
              PE_POINT(N) = P
              i_point(N) = i_point_G(N) - p_ioff_gg(P)
              j_point(N) = j_point_G(N) - p_joff_gg(P)             
            !ELSE
            ! PRINT*, '        ***WARNING***'
            ! PRINT*,' ===>>> : ERROR finding the processor'
            ! PRINT*,' ======   =====                ===='
            ! PRINT*,' j_point_G ',j_point_G(I),' I_point_G = ',I_point_G(I)
            ! PRINT*, ' '
            ! PRINT*,' We STOP!!! Verify the OCECTL '
            ! PRINT*, ' '
            ! CALL STOP_ALL('STOP')        
           ENDIF
         ENDIF
     ENDDO
    ENDDO
   ENDIF
       
      CALL p_bcast(PE_POINT,p_io)
      CALL p_bcast(i_point,p_io)
      CALL p_bcast(j_point,p_io) 
       
   !  IF (p_pe == p_io) then
   !    DO N = 1, POINT
   !       WRITE(IO_STDOUT,*) 'p_pe number for diagnostic ouput: '
   !       WRITE(IO_STDOUT,*) PE_POINT(N)
   !       WRITE(IO_STDOUT,*) 'Index of p_pe number : '
   !       WRITE(IO_STDOUT,*) 'i_point: ', i_point(N)
   !       WRITE(IO_STDOUT,*) 'j_point: ', j_point(N)                                
   !    ENDDO
   !  ENDIF

 
   END SUBROUTINE INI_diagnosis
 
 
 
 
 
 
                                 
    SUBROUTINE MMEAN_DIAG2D                                           &
     &(KDAYS,KMONTS,KYEARS,MNNNDT,MNNN720,IUN2D,                      &
     &                   FIELD2D,SUM_FIELD2D)
      
!    *KDAYS*       - actual day.
!    *KMONTS*      - actual month.
!    *KYEARS*      - actual year.
!    *MNNNDT*      - actual number of day loop(begin).
!    *MNNN720*     - average number(all = begin + MNNN720)
!    *IUN2D*       - IO-name
!    *FIELD2D*     - field 2D variable
!    *SUM_FIELD2D* - global field 2D variable

      INTEGER,INTENT(IN)::KDAYS,KMONTS,KYEARS,MNNNDT,MNNN720,IUN2D
      REAL,INTENT(IN)::FIELD2D(:,:)
      REAL,INTENT(INOUT)::SUM_FIELD2D(:)
      INTEGER(KIND=i4) IDATE,ICODE,KLEV,IEDIM
      INTEGER N
      character(len=180) str
      REAL::FIELD21,SUM_FIELD21
      REAL(KIND=sp)::FF_G
           
      DO N = 1, POINT       
        SUM_FIELD21=SUM_FIELD2D(N)
        FIELD21 = 0.

        IF (p_pe .eq. PE_POINT(N)) THEN
           FIELD21=FIELD2D(i_point(N),j_point(N))

           ICODE=N
           IEDIM=1

           IF (MNNNDT.EQ.1) THEN
              SUM_FIELD21=0.
           ENDIF

           IF (MNNNDT.NE.MNNN720) THEN
              SUM_FIELD21=SUM_FIELD21+FIELD21
           ELSE
              SUM_FIELD21=SUM_FIELD21+FIELD21
              SUM_FIELD21=SUM_FIELD21/FLOAT(MNNN720)
              IF (ABS(SUM_FIELD21).LT.(1.E-35)) THEN
                SUM_FIELD21=0.
              ENDIF                                                  
              FF_G = real(SUM_FIELD21,sp)

           !   WRITE(IO_STDOUT,*)'OPEN DIAG UNIT: ',IUN2D
              write(str,'(I4,"_",I3.3)') IUN2D,p_pe
           !   write(IO_STDOUT,*)'diagnosis.'//trim(adjustl(str))
              OPEN(IUN2D,FILE='./diagnosis/diagnosis.'//trim(adjustl(str)), & 
     &           STATUS='UNKNOWN', ACCESS='SEQUENTIAL', & 
     &           POSITION='APPEND', FORM='UNFORMATTED')

              KLEV=-100
              IDATE= (KYEARS*10000)+(KMONTS*100)+KDAYS 

              WRITE(IO_STDOUT,*)'DIAG2D: Unit=',IUN2D
              WRITE(IO_STDOUT,'(1x,A11,I8)') 'YYYYMMDD=',IDATE
              WRITE(IO_STDOUT,*)'         No.=',icode
              WRITE(IO_STDOUT,*)'    nanf/nend=',MNNNDT,MNNN720

              WRITE(IUN2D) IDATE,ICODE,KLEV,IEDIM
              WRITE(IUN2D) FF_G
              CLOSE(IUN2D)
           ENDIF 
        ENDIF ! p_pe == PE_POINT(I) 
           
        SUM_FIELD2D(N)=SUM_FIELD21 
      ENDDO

      RETURN
    END SUBROUTINE MMEAN_DIAG2D
    
    
    
    
    

    SUBROUTINE MMEAN_DIAG3D                                            &
     & (KDAYS,KMONTS,KYEARS,MNNNDT,MNNN720,IUN3D,                       &
     &   FIELD,SUM_FIELD,IVEC)
     
!    *KDAYS*     - actual day.
!    *KMONTS*    - actual month.
!    *KYEARS*    - actual year.
!    *MNNNDT*    - actual number of day loop(begin).
!    *MNNN720*   - average number(all = begin + MNNN720)
!    *iun3d*     - IO-name
!    *field*     - field variable
!    *sum_field* - field variable sum
!    *ievc*      - check number for u/v/w/p-variable
          
      INTEGER,INTENT(IN)::KDAYS,KMONTS,KYEARS,MNNNDT,MNNN720,IUN3D,  &
     &   IVEC
      REAL,INTENT(IN)::FIELD(:,:,:)
      REAL,INTENT(INOUT)::SUM_FIELD(:,:)
      INTEGER(KIND=i4) IDATE,KLEV,ICODE,IEDIM      
      INTEGER K,N,IEVC1,KK
      character(len=180) str
      REAL,allocatable::FIELD1(:),SUM_FIELD1(:)
      REAL(KIND=sp),allocatable::FF_G(:)
      
      ievc1=IVEC
     
      IF (ievc1 .eq. -1) then
        kk = kep
      ELSE
        kk = ke
      ENDIF
     
      allocate(FIELD1(1:KK)) 
      allocate(SUM_FIELD1(1:KK))

     FIELD1(:)=0.

      DO N = 1, POINT       
           SUM_FIELD1(:)=SUM_FIELD(N,:) 
        IF (p_pe .eq. PE_POINT(N)) THEN   ! The PE_point processor can proceed the entire DO-ENDDO                    
           FIELD1(:)=FIELD(i_point(N),j_point(N),:)

           ICODE=N                    
           IEDIM=1
 
          IF (MNNNDT.EQ.1) THEN
              SUM_FIELD1(:)=0.      
          ENDIF

          IF (MNNNDT.NE.MNNN720) THEN
              SUM_FIELD1(:)=SUM_FIELD1(:)+FIELD1(:)
          ELSE
              ALLOCATE(FF_G(KK))

              DO K=1,KE
                SUM_FIELD1(k)=SUM_FIELD1(k)+FIELD1(k)
                SUM_FIELD1(k)=SUM_FIELD1(k)/FLOAT(MNNN720) 
                IF (ABS(SUM_FIELD1(k)).LT.(1.E-35)) THEN
                    SUM_FIELD1(k)=0.
                ENDIF
              ENDDO

              IDATE= (KYEARS*10000)+(KMONTS*100)+KDAYS 
              WRITE(IO_STDOUT,*)'DIAG3D: Unit=',IUN3D
              WRITE(IO_STDOUT,'(1X,A11,I8)') '  YYYYMMDD=',IDATE
              WRITE(IO_STDOUT,*)'         No.=',icode
              WRITE(IO_STDOUT,*)'    nanf/nend=',MNNNDT,MNNN720

              FF_G = REAL(SUM_FIELD1,sp)    
              write(str,'(I4,"_",I3.3)') IUN3D,p_pe 
            !  write(IO_STDOUT,*)'diagnosis.'//trim(adjustl(str))
              OPEN(IUN3D,FILE='./diagnosis/diagnosis.'//trim(adjustl(str)), & 
     &         STATUS='UNKNOWN', ACCESS='SEQUENTIAL', &
     &         POSITION='APPEND', FORM='UNFORMATTED')

              DO K=1,KE
                IF (ievc1 .eq. -1) THEN
                  KLEV=NINT(TIESTW(K))
                ELSE
                  KLEV=NINT(TIESTU(K))
                ENDIF

                IF (DZW(1) .LT. 1) KLEV=K ! spy

                WRITE(IUN3D) IDATE,ICODE,KLEV,IEDIM
                WRITE(IUN3D) FF_G(k)           
              ENDDO      
              CLOSE(IUN3D)

             DEALLOCATE(FF_G)
          ENDIF
        ENDIF ! p_pe==PE_POINT(N)
        SUM_FIELD(N,:)=SUM_FIELD1(:)
      ENDDO
     
      DEALLOCATE(FIELD1,SUM_FIELD1)


      RETURN     
    END SUBROUTINE MMEAN_DIAG3D    
#endif
END MODULE MO_DIAGNOSIS                  
