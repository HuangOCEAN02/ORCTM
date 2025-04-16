       MODULE MO_MEAN
       USE MO_PARAM1
#ifdef GLOBAL_NETCDF
       USE MO_GLUE
#endif
       IMPLICIT NONE

       REAL, POINTER :: FLUM(:,:),PEM(:,:),SICTRU(:,:)                  &
     &  ,SICTRV(:,:),SUM_EMINPO(:,:),SUM_ZO(:,:)                        &
     &  ,SUM_FLUM(:,:),SUM_PEM(:,:),SUM_SICTHO(:,:)                     &
     &  ,SUM_SICOMO(:,:),SUM_SICUO(:,:),SUM_SICVE(:,:)                  &
     &  ,SUM_SICTRU(:,:),SUM_SICTRV(:,:),SUM_THO(:,:,:)                 &
     &  ,SUM_SAO(:,:,:),SUM_WO(:,:,:),SUM_PO(:,:,:)                     &
     &  ,SUM_SICSNO(:,:),SUM_PRECH(:,:),SUM_RIVRUN(:,:)                 &
     &  ,SUM_QSWO(:,:),SUM_QLWO(:,:),SUM_QLAO(:,:),SUM_QSEO(:,:)        &
     &  ,SUM_UKO(:,:,:),SUM_VKE(:,:,:)                                  &
! PO added by peterspy
#ifdef NON_HYDROSTATIC
     &  ,SUM_PNH(:,:,:),SUM_DIVG(:,:,:)                                 &
#ifdef CG2D4NONHY      
     &  ,SUM_PNHBOT(:,:),SUM_DIVGBOT(:,:)                               & 
#endif /*CG2D4NONHY*/      
#endif
#ifdef DENSITYIO
   &  ,ROO(:,:,:), SUM_ROO(:,:,:)                                                       &
#endif
#ifdef DIFFDIAG
     &  ,SUM_AVO(:,:,:),SUM_DVO(:,:,:) &
     &  ,WTMIX(:,:,:),SUM_WTMIX(:,:,:) &
#endif /*DIFFDIAG*/
#ifdef BAROTROPIC
     &  ,SUM_USO(:,:),SUM_VSE(:,:)                                      & 
#endif /*BAROTROPIC*/
#ifdef GMBOLUS      
     &  ,SUM_WGO(:,:,:),SUM_BOLX(:,:),SUM_BOLY(:,:)                     &
#endif /*GMBOLUS*/
     &  

      CONTAINS

      SUBROUTINE alloc_mem_mean

       ALLOCATE( FLUM(IE,JE),PEM(IE,JE),SICTRU(IE,JE)                   &
     &  ,SICTRV(IE,JE),SUM_EMINPO(IE,JE),SUM_ZO(IE,JE)                  &
     &  ,SUM_FLUM(IE,JE),SUM_PEM(IE,JE),SUM_SICTHO(IE,JE)               &
     &  ,SUM_SICOMO(IE,JE),SUM_SICUO(IE,JE),SUM_SICVE(IE,JE)            &
     &  ,SUM_SICTRU(IE,JE),SUM_SICTRV(IE,JE),SUM_THO(IE,JE,KE)          &
     &  ,SUM_SAO(IE,JE,KE),SUM_WO(IE,JE,KE),SUM_PO(IE,JE,KE)            &
     &  ,SUM_SICSNO(IE,JE),SUM_PRECH(IE,JE),SUM_RIVRUN(IE,JE)           &
     &  ,SUM_QSWO(IE,JE),SUM_QLWO(IE,JE),SUM_QLAO(IE,JE),SUM_QSEO(IE,JE)&
     &  ,SUM_UKO(I_start:IE,JE,KE),SUM_VKE(IE,J_start:JE,KE)            &
! PO added by peterspy
#ifdef NON_HYDROSTATIC
     &  ,SUM_PNH(IE,JE,KE),SUM_DIVG(IE,JE,KE)                           &
#ifdef CG2D4NONHY          
     &  ,SUM_PNHBOT(IE,JE),SUM_DIVGBOT(IE,JE)                           &
#endif /*CG2D4NONHY*/                                                   
#endif /*NON_HYDROSTATIC*/
#ifdef DENSITYIO
     &  ,ROO(IE,JE,KE),SUM_ROO(IE,JE,KE)                                                &
#endif
#ifdef DIFFDIAG
     &  ,SUM_AVO(IE,JE,KE),SUM_DVO(IE,JE,KE)                            &
     &  ,WTMIX(IE,JE,KE),SUM_WTMIX(IE,JE,KE)                            &
#endif /*DIFFDIAG*/
#ifdef BAROTROPIC
     &  ,SUM_USO(I_start:IE,JE),SUM_VSE(IE,J_start:JE)                  &
#endif
#ifdef GMBOLUS      
     &  ,SUM_WGO(IE,JE,KE),SUM_BOLX(IE,JE),SUM_BOLY(IE,JE)              &
#endif /*GMBOLUS*/
     &  )
      END SUBROUTINE alloc_mem_mean

      SUBROUTINE MMEAN2D                                                &
     &(kdays,kmonts,kyears,mnnndt,mnnn720,iun2d,                        &
     &                   field2,sum_field2,                             &
     &                   i4code,ivec)
!****************************************************************
!
!**** *MMEAN2D* - average data and save.
!
!     CH,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    01.10.01
!     - netCDF version (with cond.comp. PNETCDFO)
!
!     Purpose
!     -------
!     
!
!     Method
!     -------
!     
! 
!**   Interface.
!     ----------
!
!     *CALL*       *MMEAN2D(kdays,kmonts,kyears,mnnndt,mnnn720,
!                           iun2d,field2,sum_field2,i4code)
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean
!     model.
!     *COMMON*     *COMMO1.h*     - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS.h*      - std I/O logical units.
!
!**   Interface to calling routine (parameter list):
!     ----------------------------------------------
!
!     *INTEGER* *kdays*    - .
!     *INTEGER* *kmonts*    - .
!     *INTEGER* *kyears*    - .
!     *INTEGER* *mnnndt*    - actual number of day loop.
!     *INTEGER* *mnnn720*    - .
!     *INTEGER* *iun2d*    - .
!     *REAL*    *field2*   - .
!     *REAL*    *sum_field2*   - .
!     *INTEGER* *i4code*    - .
!
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS
      USE MO_KIND

      INTEGER,INTENT(IN)::kdays,kmonts,kyears,mnnndt,mnnn720,           &
     &        iun2d,i4code,ivec
      REAL,INTENT(IN)::FIELD2(:,:)
      REAL,INTENT(INOUT)::SUM_FIELD2(:,:)

      INTEGER I,J,K,I_start_G,J_start_G,ivec1
      REAL AWET
      INTEGER(KIND=i4) IDATE,ICODE,KLEV,IEDIM      
      REAL,allocatable::FIELD21(:,:),SUM_FIELD21(:,:)  
#ifndef FRAGMENT       
      REAL,allocatable::SUM_FIELD2_G(:,:)
      REAL(KIND=sp),allocatable::FF_G(:,:)
#else
      character(len=90) str,str1     
#endif /*FRAGMENT*/  

      I_start_G=1
      J_start_G=1
      if (ivec .eq. 1) then
        I_start_G=0
        allocate(FIELD21(I_start:IE,1:JE),SUM_FIELD21(I_start:IE,1:JE))
      elseif (ivec .eq. 2) then
        J_start_G=0
        allocate(FIELD21(1:IE,J_start:JE),SUM_FIELD21(1:IE,J_start:JE))
      elseif (ivec .eq. 0) then
        allocate(FIELD21(1:IE,1:JE),SUM_FIELD21(1:IE,1:JE))
      else
        CALL STOP_ALL('mo_parallel, wrong ivec number')
      endif

      ivec1=ivec

      FIELD21(:,:)=FIELD2(:,:)
      SUM_FIELD21(:,:)=SUM_FIELD2(:,:)

      ICODE=I4CODE
      IEDIM=(IE_G+1-I_start_G)*(JE_G+1-J_start_G)
!
!H At beginning of averaging period : initialize the field
!
      IF (MNNNDT.EQ.1) THEN
        if (icode.eq.1) WRITE(IO_STDOUT,*)'MMEAN2D: reset at ',mnnndt
        if (ivec .eq. 1) then
          DO I=I_start,IE
            DO J=1,JE
              SUM_FIELD21(I,J)=0.      
            ENDDO
          ENDDO
        elseif (ivec .eq. 2) then
          DO I=1,IE
            DO J=J_start,JE
              SUM_FIELD21(I,J)=0.      
            ENDDO
          ENDDO
        else
          DO I=1,IE
            DO J=1,JE
              SUM_FIELD21(I,J)=0.     
            ENDDO
          ENDDO
        endif
      ENDIF
!
!H Summation and averaging of the field; Check out of range for IEEE
!
      IF (MNNNDT.NE.MNNN720) THEN
        if (icode.eq.1) WRITE(IO_STDOUT,*)'MMEAN2D: accu at ',mnnndt
        if (ivec .eq. 1) then
          DO I=I_start,IE
            DO J=1,JE
              SUM_FIELD21(I,J)=SUM_FIELD21(I,J)+FIELD21(I,J)
            ENDDO
          ENDDO
        elseif (ivec .eq. 2) then
          DO I=1,IE
            DO J=J_start,JE
              SUM_FIELD21(I,J)=SUM_FIELD21(I,J)+FIELD21(I,J)
            ENDDO
          ENDDO
        else
          DO I=1,IE
            DO J=1,JE
              SUM_FIELD21(I,J)=SUM_FIELD21(I,J)+FIELD21(I,J)
            ENDDO
          ENDDO
        endif
      ELSE
        if (icode.eq.1) WRITE(IO_STDOUT,*)'MMEAN2D: accu at ',mnnndt
        if (icode.eq.1) WRITE(IO_STDOUT,*)'MMEAN2D: norm at ',mnnndt
        

        if (ivec .eq. 1) then
          DO I=I_start,IE
            DO J=1,JE
              SUM_FIELD21(I,J)=SUM_FIELD21(I,J)+FIELD21(I,J)
              SUM_FIELD21(I,J)=SUM_FIELD21(I,J)/FLOAT(MNNN720)
              AWET= AMSUO(I,J,1)
!              IF (AWET.LT.0.5) THEN SUM_FIELD21(I,J)=9.e9
              IF (ABS(SUM_FIELD21(I,J)).LT.(1.E-35)) THEN
                SUM_FIELD21(I,J)=0.
              ENDIF              
            ENDDO
          ENDDO
        elseif (ivec .eq. 2) then
          DO I=1,IE
            DO J=J_start,JE
              SUM_FIELD21(I,J)=SUM_FIELD21(I,J)+FIELD21(I,J)
              SUM_FIELD21(I,J)=SUM_FIELD21(I,J)/FLOAT(MNNN720)
              AWET= AMSUE(I,J,1)
!              IF (AWET.LT.0.5) THEN SUM_FIELD21(I,J)=9.e9
              IF (ABS(SUM_FIELD21(I,J)).LT.(1.E-35)) THEN
                SUM_FIELD21(I,J)=0.
              ENDIF              
            ENDDO
          ENDDO
        else
          DO I=1,IE
            DO J=1,JE
              SUM_FIELD21(I,J)=SUM_FIELD21(I,J)+FIELD21(I,J)
              SUM_FIELD21(I,J)=SUM_FIELD21(I,J)/FLOAT(MNNN720)
              AWET= WETO(I,J,1)
!              IF (AWET.LT.0.5) THEN SUM_FIELD21(I,J)=9.e9
              IF (ABS(SUM_FIELD21(I,J)).LT.(1.E-35)) THEN
                SUM_FIELD21(I,J)=0.
              ENDIF              
            ENDDO
          ENDDO
        endif
        
#ifndef FRAGMENT       
        ALLOCATE(SUM_FIELD2_G(I_start_G:IE_G,J_start_G:JE_G))
        ALLOCATE(FF_G(I_start_G:IE_G,J_start_G:JE_G))
        CALL gather_arr(SUM_FIELD21,SUM_FIELD2_G,p_io,ivec1)

!H Write in EXTRA format
        IF (p_pe == p_io) THEN
          WRITE(IO_STDOUT,*)'OPEN UNIT: ',IUN2D
          OPEN(IUN2D, STATUS='UNKNOWN',ACCESS='SEQUENTIAL',             &
     &                POSITION='APPEND',FORM='UNFORMATTED')

          KLEV=-100
          IDATE= (KYEARS*10000)+(KMONTS*100)+KDAYS 

          WRITE(IO_STDOUT,*)'MMEAN2D: Unit=',IUN2D
          WRITE(IO_STDOUT,'(1x,A11,I8)') 'YYYYMMDD=',IDATE
          WRITE(IO_STDOUT,*)'         Code=',icode
          WRITE(IO_STDOUT,*)'    nanf/nend=',mnnndt,mnnn720

          DO J=J_start_G,JE_G
            DO i=I_start_G,IE_G
              FF_G(I,J)=real(SUM_FIELD2_G(I,J),sp)
            ENDDO
          ENDDO

          WRITE(IUN2D) IDATE,ICODE,KLEV,IEDIM
          WRITE(IUN2D)((FF_G(I,J),I=I_start_G,IE_G),J=J_start_G,JE_G)
          CLOSE(IUN2D)
        ENDIF  ! p_pe==p_io
        DEALLOCATE(SUM_FIELD2_G,FF_G)
#elif defined GLOBAL_NETCDF
        CALL MEAN2D_GLOBAL(KDAYS, KMONTS, KYEARS, & 
        &        SUM_FIELD21, ICODE,ivec1)
#else

        write(str,'(I4,"/fort.")') IUN2D
        write(str1,'(I4,"_",I3.3)') IUN2D,p_pe
                 
      !  WRITE(IO_STDOUT,*) './fort.'//trim(adjustl(str))        &
      !  &  //trim(adjustl(str1))
        
        OPEN(IUN2D,    &
     &  FILE='./fort.'//trim(adjustl(str))//trim(adjustl(str1)),&
     &  STATUS='UNKNOWN',ACCESS='SEQUENTIAL',    & 
     &  POSITION='APPEND',FORM='UNFORMATTED')

        KLEV=-100
        IDATE= (KYEARS*10000)+(KMONTS*100)+KDAYS 

        WRITE(IO_STDOUT,*)'MMEAN2D: Unit=',IUN2D
        WRITE(IO_STDOUT,'(1x,A11,I8)') 'YYYYMMDD=',IDATE
        WRITE(IO_STDOUT,*)'         Code=',icode
        WRITE(IO_STDOUT,*)'    nanf/nend=',mnnndt,mnnn720
        
        WRITE(IUN2D) IDATE,ICODE,KLEV,INT(SIZE(SUM_FIELD21),i4)
        WRITE(IUN2D) real(SUM_FIELD21,sp)
        CLOSE(IUN2D)

#endif /*FRAGMENT*/          
      ENDIF 

      SUM_FIELD2(:,:)=SUM_FIELD21(:,:)
      DEALLOCATE(FIELD21,SUM_FIELD21)

      RETURN
      END SUBROUTINE MMEAN2D

      SUBROUTINE MMEAN3D                                                &
     &(kdays,kmonts,kyears,mnnndt,mnnn720,iun3d,                        &
     &   field,sum_field,i4code,ivec)
!
!**** *MMEAN3D* - save mean diagnostic output.
!
!     CHH,    *MPI-Met, HH*   14.01.99
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    01.10.01
!     - separate routine extracted from OLLIE (MAIN)
!     - NetCDF version (cond. comp. PNETCDF)
!
!     Purpose
!     -------
!     Accumulate 3d fields, average, and save.
!
!     Method
!     -------
!     Field is set to 0 at the first time step of each month,
!     accumulated each step, and normalized at the end of the
!     averaging time period. The filed is written to disk in
!     EXTRA format (default, or in NetCDF).
!
!**   Interface.
!     ----------
!
!     *CALL*       *MMEAN3D(kdays,kmonts,kyears,mnnndt,mnnn720,iun3d,
!                        field,sum_field,i4code)*
!
!     *PARAMETER*  *MO_PARAM1*     - grid size parameters for ocean
!     model.
!     *COMMON*     *MO_COMMO1*     - ocean/sediment tracer arrays.
!     *COMMON*     *MO_UNITS*      - std I/O logical units.
!
!**   Interface to calling routine (parameter list):
!     ----------------------------------------------
!
!     *INTEGER* *KYEARS   - actual year.
!     *INTEGER* *KMONTS*   - actual month.
!     *INTEGER* *KDAYS    - actual day.
!
!
!     Externals
!     ---------
!     none.
!
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS
      USE MO_KIND

      INTEGER,INTENT(IN)::kdays,kmonts,kyears,mnnndt,mnnn720,           &
     &        iun3d,i4code,ivec
      REAL,INTENT(IN)::FIELD(:,:,:)
      REAL,INTENT(INOUT)::SUM_FIELD(:,:,:)

      INTEGER I,J,K,I_start_G,J_start_G,ivec1,Kvec
      REAL AWET
      REAL,allocatable::FIELD1(:,:,:),SUM_FIELD1(:,:,:)

      INTEGER(KIND=i4) IDATE,KLEV,ICODE,IEDIM
#ifndef FRAGMENT      
      REAL,allocatable::SUM_FIELD_G(:,:,:)    
      REAL(KIND=sp),allocatable::FF_G(:,:)
#else
      character(len=90) str,str1 
#endif /*FRAGMENT*/ 

      I_start_G=1
      J_start_G=1
      if (ivec .eq. 1) then
        I_start_G=0
        allocate(FIELD1(I_start:IE,1:JE,1:KE))
        allocate(SUM_FIELD1(I_start:IE,1:JE,1:KE))
      elseif (ivec .eq. 2) then
        J_start_G=0
        allocate(FIELD1(1:IE,J_start:JE,1:KE))
        allocate(SUM_FIELD1(1:IE,J_start:JE,1:KE))
      elseif (ivec .eq. -1) then
        allocate(FIELD1(1:IE,1:JE,1:KEP))
        allocate(SUM_FIELD1(1:IE,1:JE,1:KEP))
      elseif (ivec .eq. 0) then
        allocate(FIELD1(1:IE,1:JE,1:KE))
        allocate(SUM_FIELD1(1:IE,1:JE,1:KE))
      else
        CALL STOP_ALL('mo_parallel, wrong ivec number')
      endif

      ivec1=ivec
      Kvec = KE
      if (ivec1 .lt. 0) then
         ivec1=0
         Kvec = KEP
      endif

      FIELD1(:,:,:)=FIELD(:,:,:)
      SUM_FIELD1(:,:,:)=SUM_FIELD(:,:,:)

      IEDIM=(IE_G+1-I_start_G)*(JE_G+1-J_start_G)
      ICODE=I4CODE

!H    INITIALIZE THE FIELD
  
      IF (MNNNDT.EQ.1) THEN
        if (icode.eq.2) WRITE(IO_STDOUT,*)'MMEAN3D: reset at ',mnnndt
        DO K=1,KE
          if (ivec .eq. 1) then
            DO I=I_start,IE
              DO J=1,JE
                SUM_FIELD1(I,J,K)=0.      
              ENDDO
            ENDDO
          elseif (ivec .eq. 2) then
            DO I=1,IE
              DO J=J_start,JE
                SUM_FIELD1(I,J,K)=0.      
              ENDDO
            ENDDO
          else
            DO I=1,IE
              DO J=1,JE
                SUM_FIELD1(I,J,K)=0.      
              ENDDO
            ENDDO
          endif
        ENDDO
      ENDIF

!H    SUMMATION AND MONTHLY MEAN OF THE FIELD

      IF (MNNNDT.NE.MNNN720) THEN
        if (icode.eq.2) WRITE(IO_STDOUT,*)'MMEAN3D: accu at ',mnnndt
        DO K=1,KE
          if (ivec .eq. 1) then
            DO I=I_start,IE
              DO J=1,JE
                SUM_FIELD1(I,J,K)=SUM_FIELD1(I,J,K)+FIELD1(I,J,K)
              ENDDO
            ENDDO
          elseif (ivec .eq. 2) then
            DO I=1,IE
              DO J=J_start,JE
                SUM_FIELD1(I,J,K)=SUM_FIELD1(I,J,K)+FIELD1(I,J,K)
              ENDDO
            ENDDO
          else
            DO I=1,IE
              DO J=1,JE
                SUM_FIELD1(I,J,K)=SUM_FIELD1(I,J,K)+FIELD1(I,J,K)
              ENDDO
            ENDDO
          endif
        ENDDO
      ELSE
        if (icode.eq.2) WRITE(IO_STDOUT,*)'MMEAN3D: accu at ',mnnndt
        if (icode.eq.2) WRITE(IO_STDOUT,*)'MMEAN3D: norm at ',mnnndt


        DO K=1,KE
          if (ivec .eq. 1) then
            DO I=I_start,IE
              DO J=1,JE
                SUM_FIELD1(I,J,K)=SUM_FIELD1(I,J,K)+FIELD1(I,J,K)
                SUM_FIELD1(I,J,K)=SUM_FIELD1(I,J,K)/FLOAT(MNNN720)
                AWET= AMSUO(I,J,K)
!                IF (AWET.LT.0.5) THEN SUM_FIELD1(I,J,K)=9.e9
                IF (ABS(SUM_FIELD1(I,J,K)).LT.(1.E-35)) THEN
                  SUM_FIELD1(I,J,K)=0.
                ENDIF
              ENDDO
            ENDDO
          elseif (ivec .eq. 2) then
            DO I=1,IE
              DO J=J_start,JE
                SUM_FIELD1(I,J,K)=SUM_FIELD1(I,J,K)+FIELD1(I,J,K)
                SUM_FIELD1(I,J,K)=SUM_FIELD1(I,J,K)/FLOAT(MNNN720)
                AWET= AMSUE(I,J,K)
!                IF (AWET.LT.0.5) THEN SUM_FIELD1(I,J,K)=9.e9
                IF (ABS(SUM_FIELD1(I,J,K)).LT.(1.E-35)) THEN
                  SUM_FIELD1(I,J,K)=0.
                ENDIF
              ENDDO
            ENDDO
          else
            DO I=1,IE
              DO J=1,JE
                SUM_FIELD1(I,J,K)=SUM_FIELD1(I,J,K)+FIELD1(I,J,K)
                SUM_FIELD1(I,J,K)=SUM_FIELD1(I,J,K)/FLOAT(MNNN720)
                AWET= WETO(I,J,K)
!                IF (AWET.LT.0.5) THEN SUM_FIELD1(I,J,K)=9.e9
                IF (ABS(SUM_FIELD1(I,J,K)).LT.(1.E-35)) THEN
                  SUM_FIELD1(I,J,K)=0.
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO

        IDATE= (KYEARS*10000)+(KMONTS*100)+KDAYS 
        WRITE(IO_STDOUT,*)'MMEAN3D: Unit=',IUN3D
        WRITE(IO_STDOUT,'(1X,A11,I8)') '  YYYYMMDD=',IDATE
        WRITE(IO_STDOUT,*)'         Code=',icode
        WRITE(IO_STDOUT,*)'    nanf/nend=',mnnndt,mnnn720
        
#ifndef FRAGMENT 

        ALLOCATE(SUM_FIELD_G(I_start_G:IE_G,J_start_G:JE_G,KE))
        ALLOCATE(FF_G(I_start_G:IE_G,J_start_G:JE_G))

        DO K=1,KE
        CALL gather_arr(SUM_FIELD1(:,:,K),SUM_FIELD_G(:,:,K),p_io,ivec1)
        ENDDO

        IF (p_pe==p_io) THEN
          OPEN(IUN3D, STATUS='UNKNOWN', ACCESS='SEQUENTIAL',            &
     &                POSITION='APPEND', FORM='UNFORMATTED')

          DO K=1,KE
            DO j=J_start_G,JE_G
              DO i=I_start_G,IE_G
                FF_G(i,j)=REAL(SUM_FIELD_G(I,J,K),sp)
              ENDDO
            ENDDO  

            IF (ICODE.EQ.7.OR.ICODE.EQ.207.OR.                          &
     &          ICODE.EQ.110.OR.ICODE.EQ.111) THEN
              KLEV=NINT(TIESTW(K))
            ELSE
              KLEV=NINT(TIESTU(K))
            ENDIF

            IF (DZW(1) .LT. 1) KLEV=K ! spy

            WRITE(IUN3D) IDATE,ICODE,KLEV,IEDIM
            WRITE(IUN3D)((FF_G(I,J),I=I_start_G,IE_G),J=J_start_G,JE_G)
          ENDDO

          CLOSE(IUN3D)
        ENDIF ! p_pe==p_io
        DEALLOCATE(SUM_FIELD_G,FF_G)
#elif defined GLOBAL_NETCDF
        CALL MEAN3D_GLOBAL(KDAYS, KMONTS, KYEARS, &
        &        SUM_FIELD1, ICODE,ivec1)
#else

        write(str,'(I4,"/fort.")') IUN3D
        write(str1,'(I4,"_",I3.3)') IUN3D,p_pe
        
        IEDIM = INT(SIZE(SUM_FIELD1),i4)/Kvec         
      !  WRITE(IO_STDOUT,*) './fort.'//trim(adjustl(str))          &
     ! &    //trim(adjustl(str1))
        
        OPEN(IUN3D, &
     &  FILE='./fort.'//trim(adjustl(str))//trim(adjustl(str1)),  &
     &  STATUS='UNKNOWN',ACCESS='SEQUENTIAL',    & 
     &  POSITION='APPEND',FORM='UNFORMATTED')
     
        DO K=1,KE
          IF (ICODE.EQ.7.OR.ICODE.EQ.207.OR.                          &
     &        ICODE.EQ.110.OR.ICODE.EQ.111) THEN
              KLEV=NINT(TIESTW(K))
          ELSE
              KLEV=NINT(TIESTU(K))
          ENDIF

          IF (DZW(1) .LT. 1) KLEV=K ! spy
          WRITE(IUN3D) IDATE,ICODE,KLEV,IEDIM
          WRITE(IUN3D) real(SUM_FIELD1(:,:,K),sp)
        ENDDO
          
          CLOSE(IUN3D)
#endif /*FRAGMENT*/
        
      ENDIF

      SUM_FIELD(:,:,:)=SUM_FIELD1(:,:,:)
      DEALLOCATE(FIELD1,SUM_FIELD1)

      RETURN
      END SUBROUTINE MMEAN3D


      END MODULE MO_MEAN
