! GLUE: This is f77 code used petcdf(parallel netcdf) library to   
! merge the fields of each processor into a global field.
! Attendion !!! 
! If adopted Periodic boundary condition(2.5D/Globe) or purely
! 2D domian, array indexs need to be modified for saving time.
! JULY 24 2021  H.H  - OLD VERSION: Single processor(long time)
! JULY 25 2021  H.H  - NEW VERSION: Pnetcdf (relatively short time)
MODULE MO_GLUE
  INCLUDE 'mpif.h'
  INCLUDE "pnetcdf.inc"

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)  
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
  INTEGER, PARAMETER :: wp = dp   ! working precision
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)

!|-----------------BASIC--INFO----NEED-TO-BE-Modified-------------|     
  INTEGER,PARAMETER::IN_ETA_FILE=82,IN_TEM_FILE=71,IN_SAL_FILE=72, &
   & IN_UKO_FILE=73,IN_VKE_FILE=74,IN_WWW_FILE=46,IN_POO_FILE=76,  &
   & IN_PNH_FILE=77,IN_UZO_FILE=63,IN_VSE_FILE=64 
  INTEGER,PARAMETER::STDOUT=6  
  INTEGER,PARAMETER::nprocx=120,nprocy=1
  INTEGER,PARAMETER::ie_g=4000,je_g=5,KE=300,TE=1440
  LOGICAL::DENSITY=.TRUE.,BAROTROPIC=.TRUE.,VERTICAL=.TRUE.,     &
   &  Nonhydrostatic=.TRUE.,Pressure=.TRUE.,Salinity=.TRUE.,     &
   &  Tradition=.TRUE.,EastWest=.TRUE.,NorthSouth=.FALSE.
!|-----------------BASIC--INFO--END-------------------------------|  


  INTEGER::p_pe,p_io,ierr,nprocxy,nodes,KEP
  INTEGER::ie,je,p_ioff,p_joff
  INTEGER,ALLOCATABLE::p_lim_x(:), p_lim_y(:), p_num_x(:), p_num_y(:) &
                & ,p_ioff_g(:), p_joff_g(:), p_size_x(:), p_size_y(:)                 
  LOGICAL::have_g_is, have_g_ie, have_g_js, have_g_je   
  INTEGER::I_START,J_START,iss,jss,iee,jee
  REAL(sp),ALLOCATABLE::UUU(:,:,:),VVV(:,:,:),WWW(:,:,:),TTT(:,:,:),  &
   &    SSS(:,:,:),POO(:,:,:),ROU(:,:,:),ETA(:,:)
  REAL(sp),ALLOCATABLE::dummyE(:,:),dummyE_3d(:,:,:), &
   & dummyE_3d_U(:,:,:),dummyE_3d_V(:,:,:),dummyE_2d_U(:,:), &
   & dummyE_2d_V(:,:)
   

!|-----------------BASIC--INFO--OF--Pnetcdf-----------------------|   
  INTEGER*8::global_nx, global_ny, global_nz   
  INTEGER::cmode,varid,ncid,dimidxy(2),dimidxyz(3), &
    &  varidT,varidS,varidP,varidR    
  INTEGER*8::starts_2d(2), counts_2d(2)   
  INTEGER*8::starts_3d(3), counts_3d(3)    
  INTEGER::subarray_p_2d,subarray_p_3d,subarray_u,subarray_v,nTypes 
  INTEGER::subarray_u_2d,subarray_v_2d
  INTEGER::array_of_sizes(2), array_of_subsizes(2)
  INTEGER::array_of_starts(2)   
  INTEGER::array_of_sizes_3d(3), array_of_subsizes_3d(3)
  INTEGER::array_of_starts_3d(3)   
  
!|-----------------BASIC--INFO--OF--Pnetcdf--END------------------|  

  CONTAINS
  
  
  SUBROUTINE MY_MPI
  
   CALL MPI_Init(ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,p_pe,ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nodes,ierr)

   p_io = 0
   ! CALL MPI_BCAST(nprocx,1,MPI_INTEGER4,p_io,MPI_COMM_WORLD,ierr)
   ! CALL MPI_BCAST(nprocy,1,MPI_INTEGER4,p_io,MPI_COMM_WORLD,ierr)
  
   nprocxy = nprocx*nprocy

   if (p_pe .eq. p_io) then
    if (nprocxy .ne. nodes) then
       WRITE(STDOUT,*) '  Error about the proceseors prescribed !!!'
       CALL sleep(5)
       stop
    else 
        WRITE(STDOUT,*) '   MPI SUCCESS !!! ', p_pe
    endif
   endif  
      
  END SUBROUTINE MY_MPI   
  
  SUBROUTINE p_deco

   INTEGER i, nx, ny
  
   ALLOCATE(p_lim_x(0:nprocx), p_lim_y(0:nprocy))
   DO i=0,nprocx
        p_lim_x(i) = 2 + i*(ie_g-2)/nprocx
   ENDDO
   DO i=0,nprocy
        p_lim_y(i) = 2 + i*(je_g-2)/nprocy
   ENDDO

! Set number of processors in x and y direction
   ALLOCATE(p_num_x(0:nprocxy-1), p_num_y(0:nprocxy-1))
   DO i=0,nprocx-1
         p_num_x(i:nprocxy-1:nprocx) = i
   ENDDO
   DO i=0,nprocy-1
        p_num_y(i*nprocx:(i+1)*nprocx-1) = i
   ENDDO

! Offsets and sizes
   ALLOCATE(p_ioff_g(0:nprocxy-1), p_joff_g(0:nprocxy-1))
   ALLOCATE(p_size_x(0:nprocxy-1), p_size_y(0:nprocxy-1))   
   DO i=0,nprocxy-1
      nx = p_num_x(i)
      ny = p_num_y(i)

      p_ioff_g(i) = p_lim_x(nx)-2
      p_joff_g(i) = p_lim_y(ny)-2
      p_size_x(i) = p_lim_x(nx+1) - p_lim_x(nx) + 2
      p_size_y(i) = p_lim_y(ny+1) - p_lim_y(ny) + 2
   ENDDO

      ! Get our own values

   ie = p_size_x(p_pe)
   je = p_size_y(p_pe)
   p_ioff = p_ioff_g(p_pe)
   p_joff = p_joff_g(p_pe)
   have_g_is = p_num_x(p_pe) == 0
   have_g_ie = p_num_x(p_pe) == nprocx-1
   have_g_js = p_num_y(p_pe) == 0
   have_g_je = p_num_y(p_pe) == nprocy-1

    if (have_g_is) then
      I_start = 0
    else
      I_start = 1
    endif
    if (have_g_js) then
      J_start = 0
    else
      J_start = 1
    endif
    
  END SUBROUTINE p_deco  

  SUBROUTINE ADISITJ(TH,SH,PA)
!     TRANSFORMATION FROM POTENTIAL TO IN SITU DENSITY
!     from adisitj.f90 in MPIOM
!     changed by peterspy, used in post process

      PARAMETER(A1=3.6504E-4,A2=8.3198E-5,A3=5.4065E-7,A4=4.0274E-9,    &
     & B1=1.7439E-5,B2=2.9778E-7,C1=8.9309E-7,C2=3.1628E-8,             &
     & C3=2.1987E-10,D=4.1057E-9,E1=1.6056E-10,E2=5.0484E-12)


      REAL(kind=SELECTED_REAL_KIND(6,37))::TH(ie,je),SH(ie,je),PA(ie,je),     &
     &  PR(ie,je),TPO(ie,je),T(ie,je),QC(ie,je),QV(ie,je),DC(ie,je),DV(ie,je),&
     &  QNQ(ie,je),QN3(ie,je),QVS(ie,je),DVS(ie,je),FNE(ie,je),FST(ie,je)


      PR=PA
      QC = PR*(A1 + PR*(C1 - E1*PR))
      QV = PR*(B1 - D*PR)
      DC = 1. + PR*(-A2 + PR*(C2 -E2*PR))
      DV = B2*PR
      QNQ  = -PR*(-A3 + PR*C3)
      QN3  = -PR*A4
!
      QVS = QV*(SH - 35.) + QC
      DVS = DV*(SH - 35.) + DC
      TPO     = TH
      TH      = (TH + QVS)/DVS
      T       = TH
      FNE     = - QVS + T*(DVS + T*(QNQ + T*QN3)) - TPO
      FST     = DVS + T*(2*QNQ + 3*QN3*T)
      TH      =TH   -FNE/FST
!
      RETURN
   END SUBROUTINE ADISITJ



   SUBROUTINE RHO(T,S,P,RH)
!     calculate density by using tem sal and pres
!     from rho2.f90 in MPIOM
!     changed by peterspy, used in post process
      REAL(kind=SELECTED_REAL_KIND(6,37))::B0,B1,B2,B3,B4,C0,C1,C2,D0,  &
               A0,A1,A2,A3,A4,A5,F0,F1,F2,F3,G0,G1,G2,AI0,AI1,AI2,AJ0,  &
               AM0,AM1,AM2,E0,E1,E2,E3,E4,H0,H1,H2,H3,AK0,AK1,AK2
      REAL(kind=SELECTED_REAL_KIND(6,37))::S(ie,je),T(ie,je),P(ie,je),  &
     & RH(ie,je),S3H(ie,je)
      DATA B0,B1,B2,B3,B4/8.24493E-1,-4.0899E-3,7.6438E-5,              &
     &-8.2467E-7,5.3875E-9/
      DATA C0,C1,C2/-5.72466E-3,1.0227E-4,-1.6546E-6/
      DATA D0/4.8314E-4/
      DATA A0,A1,A2,A3,A4,A5/999.842594,6.793952E-2,                    &
     &-9.095290E-3,1.001685E-4,-1.120083E-6,6.536332E-9/
      DATA F0,F1,F2,F3/54.6746,-0.603459,1.09987E-2,-6.1670E-5/
      DATA G0,G1,G2/7.944E-2,1.6483E-2,-5.3009E-4/
      DATA AI0,AI1,AI2/2.2838E-3,-1.0981E-5,-1.6078E-6/
      DATA AJ0/1.91075E-4/
      DATA AM0,AM1,AM2/-9.9348E-7,2.0816E-8,9.1697E-10/
      DATA E0,E1,E2,E3,E4/19652.21,148.4206,-2.327105,                  &
     &1.360477E-2,-5.155288E-5/
      DATA H0,H1,H2,H3/3.239908,1.43713E-3,1.16092E-4,-5.77905E-7/
      DATA AK0,AK1,AK2/8.50935E-5,-6.12293E-6,5.2787E-8/
!
      S3H=SQRT(S **3)
      RH =(A0+T *(A1+T                                      &
     &       *(A2+T *(A3+T *(A4+T *A5))))                   &
     &       +S *(B0+T *(B1+T                               &
     &      *(B2+T *(B3+T *B4))))+D0*S **2                  &
     &      +S3H*(C0+T *(C1+C2*T )) )                       &
     &       /(1.-P /(P *(                                  &
     &   H0+T *(H1+T *(H2+T *H3))                           &
     &  +S *(AI0+T *(AI1+AI2*T ))+AJ0*S3H                   &
     &  +(AK0+T *(AK1+T *AK2)                               &
     &  +S *(AM0+T *(AM1+T *AM2)))*P )+                     &
     &    E0+T *(E1+T *(E2+T *(E3+T *E4)))                  &
     &      +S *(F0+T *(F1+T *(F2+T *F3)))                  &
     &      +S3H*(G0+T *(G1+G2*T ))))

      RETURN
   END SUBROUTINE RHO


   SUBROUTINE check(err, message)
      implicit none
      include "mpif.h"
      include "pnetcdf.inc"
      integer err
      character message*(*)

      ! It is a good idea to check returned value for possible error
      if (err .NE. NF_NOERR) then
          write(6,*) message//' '//nfmpi_strerror(err)
          call MPI_Abort(MPI_COMM_WORLD, -1, err)
      endif
    end SUBROUTINE check ! subroutine check


   SUBROUTINE find_index_P

     cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)  
     nTypes = 1

!    2-D Field  defaul because of the ghost points
     ! x-direction
     iss = 2
     iee = ie-1
     starts_2d(1) = p_ioff+2
     counts_2d(1) = ie-2  
     
     ! y-direction 
     jss = 2  
     jee = je-1               
     starts_2d(2) = p_joff+2
     counts_2d(2) = je-2  
                               
     ! West
     if (have_g_is) then
        iss = 1
        iee = ie-1
        starts_2d(1) = p_ioff+1
        counts_2d(1) = ie-1         
     endif    
      ! East    
     if (have_g_ie) then
        iss = 2     
        iee = ie    
        starts_2d(1) = p_ioff+2
        counts_2d(1) = ie-1             
     endif     
 
     IF(nprocy .ne. 1) then  
      ! North
      if (have_g_js) then
        jss = 1   
        jee = je-1   
        starts_2d(2) = p_joff+1
        counts_2d(2) = je-1             
      endif      
      ! South
      if (have_g_je) then
        jss = 2   
        jee = je        
        starts_2d(2) = p_joff+2
        counts_2d(2) = je-1              
      endif    
     
     ELSE ! non_parallel in y-direction for Purely 2D or 2.5D(Periodic Boundary Condition)
       jss = 3         ! Purely 2D simulation           
       jee = 3         ! Purely 2D simulation           
     ! jss = 1         ! Periodic Boundary Condition
     ! jee = je - 2    ! Periodic Boundary Condition  
       starts_2d(2) = p_joff+1
       counts_2d(2) = 1        ! Purely 2D simulation    
     ! counts_2d(2) = je -2    ! Periodic Boundary Condition    
     ENDIF  
     
     ALLOCATE(dummyE(counts_2d(1),counts_2d(2)))
          
! define an MPI datatype using MPI_Type_create_subarray()
     array_of_sizes(1)    = ie_g
     array_of_sizes(2)    = 1           ! Purely 2D simulation     
   ! array_of_sizes(2) = je_g  - 2      ! Periodic Boundary Condition
     
     array_of_subsizes(1) = counts_2d(1)
     array_of_subsizes(2) = counts_2d(2)
     
     array_of_starts(1) = 0  
     array_of_starts(2) = 0      
     
     call MPI_Type_create_subarray(2, array_of_sizes, &
      &   array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, &
      &   MPI_REAL, subarray_p_2d, err)
     call MPI_Type_commit(subarray_p_2d, err)
        
        
!    3-D Field          
     starts_3d(1) = p_ioff+2
     counts_3d(1) = ie-2       
     starts_3d(2) = p_joff+2
     counts_3d(2) = je-2         
     starts_3d(3) = 1
     counts_3d(3) = ke           
     ! West
     if (have_g_is) then
        starts_3d(1) = p_ioff+1
        counts_3d(1) = ie-1     
     endif    
      ! East    
     if (have_g_ie) then
        starts_3d(1) = p_ioff+2
        counts_3d(1) = ie-1     
     endif     
 
     IF(nprocy .ne. 1) then  
      ! North
      if (have_g_js) then
        starts_3d(2) = p_joff+1
        counts_3d(2) = je-1     
      endif      
      ! South
      if (have_g_je) then
        starts_3d(2) = p_joff+2
        counts_3d(2) = je-1     
      endif    
     
     ELSE ! non_parallel in y-direction for  Purely 2D or 2.5D(Periodic Boundary Condition)
       starts_3d(2) = p_joff+1
       counts_3d(2) =  1         ! Purely 2D simulation    
     ! counts_3d(2) = je -2      ! Periodic Boundary Condition
     ENDIF          

     ALLOCATE(dummyE_3d(counts_3d(1),counts_3d(2),counts_3d(3)))
      
   ! define an MPI datatype using MPI_Type_create_subarray()
   array_of_sizes_3d(1)    = ie_g
   array_of_sizes_3d(2)    = 1          ! Purely 2D simulation      
 ! array_of_sizes_3d(2) = je_g  - 2     ! Periodic Boundary Condition
   array_of_sizes_3d(3)    = ke  
   array_of_subsizes_3d(1) = counts_3d(1)
   array_of_subsizes_3d(2) = counts_3d(2)
   array_of_subsizes_3d(3) = counts_3d(3)  
   array_of_starts_3d(1)   = 0   ! MPI start index starts with 0
   array_of_starts_3d(2)   = 0
   array_of_starts_3d(3)   = 0   
   call MPI_Type_create_subarray(3, array_of_sizes_3d, &
     &   array_of_subsizes_3d, array_of_starts_3d, MPI_ORDER_FORTRAN, &
     &   MPI_REAL, subarray_p_3d, err)
   call MPI_Type_commit(subarray_p_3d, err)
                 
   END SUBROUTINE find_index_P

   SUBROUTINE find_index_U
 
   ! x-direction
     starts_3d(1) = p_ioff+2
     counts_3d(1) = ie-2     
   ! y-direction     
     starts_3d(2) = p_joff+2
     counts_3d(2) = je-2       
   ! z-direction        
     starts_3d(3) = 1
     counts_3d(3) = ke   
                
     ! West
     if (have_g_is) then
        iss = I_start      ! U-point  
        starts_3d(1) = p_ioff+1
        counts_3d(1) = ie   ! U-point
     endif    
      ! East    
     if (have_g_ie) then     
        starts_3d(1) = p_ioff+2
        counts_3d(1) = ie-1     
     endif     
 
     IF(nprocy .ne. 1) then  
      ! North
      if (have_g_js) then
        starts_3d(2) = p_joff+1
        counts_3d(2) = je-1     
      endif      
      ! South
      if (have_g_je) then
        starts_3d(2) = p_joff+2
        counts_3d(2) = je-1     
      endif    
     
     ELSE ! non_parallel in y-direction for  Purely 2D or 2.5D(Periodic Boundary Condition)
       starts_3d(2) = p_joff+1
       counts_3d(2) =  1     ! Purely 2D simulation    
     ! counts_3d(2) = je -2  ! Periodic Boundary Condition
     ENDIF         

     ALLOCATE(dummyE_3d_U(counts_3d(1),counts_3d(2),counts_3d(3)))
     
     
   ! define an MPI datatype using MPI_Type_create_subarray()
   array_of_sizes_3d(1)    = ie_g+1
   array_of_sizes_3d(2)    = 1         ! Purely 2D simulation      
 ! array_of_sizes_3d(2) = je_g  - 2    ! Periodic Boundary Condition
   array_of_sizes_3d(3)    = ke  
   array_of_subsizes_3d(1) = counts_3d(1)
   array_of_subsizes_3d(2) = counts_3d(2)
   array_of_subsizes_3d(3) = counts_3d(3)  
   array_of_starts_3d(1)   = 0   ! MPI start index starts with 0
   array_of_starts_3d(2)   = 0
   array_of_starts_3d(3)   = 0   
   call MPI_Type_create_subarray(3, array_of_sizes_3d, &
     &   array_of_subsizes_3d, array_of_starts_3d, MPI_ORDER_FORTRAN, &
     &   MPI_REAL, subarray_u, err)
   call MPI_Type_commit(subarray_u, err)     
     
   END SUBROUTINE find_index_U
   
   SUBROUTINE find_index_V  
   
     starts_3d(1) = p_ioff+2
     counts_3d(1) = ie-2       
     starts_3d(2) = p_joff+2
     counts_3d(2) = je-2         
     starts_3d(3) = 1
     counts_3d(3) = ke 
               
     ! West
     if (have_g_is) then
        starts_3d(1) = p_ioff+1
        counts_3d(1) = ie-1     
     endif    
      ! East    
     if (have_g_ie) then
        starts_3d(1) = p_ioff+2
        counts_3d(1) = ie-1     
     endif     
 
     IF(nprocy .ne. 1) then  
      ! North
      if (have_g_js) then
        jss = J_start      ! V-point 
        starts_3d(2) = p_joff+1
        counts_3d(2) = je  ! V-point  
      endif      
      ! South
      if (have_g_je) then      
        starts_3d(2) = p_joff+2
        counts_3d(2) = je-1     
      endif    
     
     ELSE ! non_parallel in y-direction for  Purely 2D or 2.5D(Periodic Boundary Condition)
       jss = J_start+2    !  Purely 2D simulation      
     ! jss = J_start      ! V-point & Periodic Boundary Condition
       starts_3d(2) = p_joff+1
       counts_3d(2) = 2         !  Purely 2D simulation   
     ! counts_3d(2) = je+1 - 2  ! V-point & Periodic Boundary Condition
     ENDIF         

     ALLOCATE(dummyE_3d_V(counts_3d(1),counts_3d(2),counts_3d(3)))

   ! define an MPI datatype using MPI_Type_create_subarray()
   array_of_sizes_3d(1)    = ie_g
   array_of_sizes_3d(2)    = 2            !  Purely 2D simulation   
 ! array_of_sizes_3d(2)    = je_g + 1 - 2 ! Periodic Boundary Condition
   array_of_sizes_3d(3)    = ke  
   array_of_subsizes_3d(1) = counts_3d(1)
   array_of_subsizes_3d(2) = counts_3d(2)
   array_of_subsizes_3d(3) = counts_3d(3)  
   array_of_starts_3d(1)   = 0   ! MPI start index starts with 0
   array_of_starts_3d(2)   = 0
   array_of_starts_3d(3)   = 0   
   call MPI_Type_create_subarray(3, array_of_sizes_3d, &
     &   array_of_subsizes_3d, array_of_starts_3d, MPI_ORDER_FORTRAN, &
     &   MPI_REAL, subarray_v, err)
   call MPI_Type_commit(subarray_v, err)        
     
   END SUBROUTINE find_index_V   


   SUBROUTINE find_index_U_2D
   
   ! x-direction
     starts_2d(1) = p_ioff+2
     counts_2d(1) = ie-2     
   ! y-direction     
     starts_2d(2) = p_joff+2
     counts_2d(2) = je-2       
     
     ! West
     if (have_g_is) then
        iss = I_start      ! U-point  
        starts_2d(1) = p_ioff+1
        counts_2d(1) = ie   ! U-point
     endif    
      ! East    
     if (have_g_ie) then     
        starts_2d(1) = p_ioff+2
        counts_2d(1) = ie-1     
     endif     
 
     IF(nprocy .ne. 1) then  
      ! North
      if (have_g_js) then
        starts_2d(2) = p_joff+1
        counts_2d(2) = je-1     
      endif      
      ! South
      if (have_g_je) then
        starts_2d(2) = p_joff+2
        counts_2d(2) = je-1     
      endif    
     
     ELSE ! non_parallel in y-direction for  Purely 2D or 2.5D(Periodic Boundary Condition)
       starts_2d(2) = p_joff+1
       counts_2d(2) = 1        !  Purely 2D simulation 
     ! counts_2d(2) = je  - 2  ! Periodic Boundary Condition
     ENDIF         


     ALLOCATE(dummyE_2d_U(counts_2d(1),counts_2d(2)))
  
     
   ! define an MPI datatype using MPI_Type_create_subarray()
     array_of_sizes(1)    = ie_g+1
     array_of_sizes(2)    = 1         !  Purely 2D simulation 
   ! array_of_sizes(2)    = je_g - 2  ! Periodic Boundary Condition

     array_of_subsizes(1) = counts_2d(1)
     array_of_subsizes(2) = counts_2d(2)
     
     array_of_starts(1) = 0  
     array_of_starts(2) = 0     
     
     call MPI_Type_create_subarray(2, array_of_sizes, &
      &   array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, &
      &   MPI_REAL, subarray_u_2d, err)
     call MPI_Type_commit(subarray_u_2d, err)     
   
   END SUBROUTINE find_index_U_2D


   SUBROUTINE find_index_V_2D
   ! x-direction   
     starts_2d(1) = p_ioff+2
     counts_2d(1) = ie-2  
   ! y-direction          
     starts_2d(2) = p_joff+2
     counts_2d(2) = je-2         
               
     ! West
     if (have_g_is) then
        starts_2d(1) = p_ioff+1
        counts_2d(1) = ie-1     
     endif    
      ! East    
     if (have_g_ie) then
        starts_2d(1) = p_ioff+2
        counts_2d(1) = ie-1     
     endif     
 
     IF(nprocy .ne. 1) then  
      ! North
      if (have_g_js) then
        jss = J_start      ! V-point 
        starts_2d(2) = p_joff+1
        counts_2d(2) = je  ! V-point  
      endif      
      ! South
      if (have_g_je) then      
        starts_2d(2) = p_joff+2
        counts_2d(2) = je-1     
      endif    
     
     ELSE ! non_parallel in y-direction for  Purely 2D or 2.5D(Periodic Boundary Condition)
       jss = J_start+2      !  Purely 2D simulation 
     ! jss = J_start        ! V-point & Periodic Boundary Condition
       starts_2d(2) = p_joff+1
       counts_2d(2) = 2   !  Purely 2D simulation 
     ! counts_2d(2) = je+1 - 2  ! V-point & Periodic Boundary Condition
     ENDIF         

     ALLOCATE(dummyE_2d_V(counts_2d(1),counts_2d(2)))


   ! define an MPI datatype using MPI_Type_create_subarray()
  	 array_of_sizes(1)    = ie_g
  	 array_of_sizes(2)    = 2 	  !  Purely 2D simulation 
   ! array_of_sizes(2)    = je_g + 1 - 2 ! Periodic Boundary Condition

  	 array_of_subsizes(1) = counts_2d(1)
  	 array_of_subsizes(2) = counts_2d(2)
 
     array_of_starts(1) = 0  
     array_of_starts(2) = 0       ! MPI start index starts with 0
 
     call MPI_Type_create_subarray(2, array_of_sizes, &
      &   array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, &
      &   MPI_REAL, subarray_v_2d, err)
     call MPI_Type_commit(subarray_v_2d, err)   
        
   END SUBROUTINE find_index_V_2D
END MODULE MO_GLUE




PROGRAM GLUE
   USE MO_GLUE
   INTEGER(KIND=i4)::I4DATE(4)   
   INTEGER :: i,j,k,t,err 
   character(len=120) :: str
   real(sp),allocatable::dummyT(:,:),dummyS(:,:),dummyP(:,:),&
    & dummyR(:,:)
   
   KEP=KE+1  
        
   CALL my_mpi 
   CALL p_deco 
      
   ALLOCATE(UUU(I_start:ie,je,ke),VVV(ie,J_start:je,ke),WWW(ie,je,kep), &
   &       ETA(ie,je),TTT(ie,je,ke),SSS(ie,je,ke),&
   &       POO(ie,je,ke),ROU(ie,je,ke))   
   ALLOCATE(dummyT(ie,je),dummyS(ie,je),dummyP(ie,je),dummyR(ie,je))
   UUU = 0.0;VVV = 0.0;WWW = 0.0
   ETA = 0.0;TTT = 0.0;SSS = 0.0
   POO = 0.0;ROU = 0.0
   dummyT = 0.0;dummyS = 0.0
   dummyP = 0.0;dummyR = 0.0
!|------------------------------------------------------------------|
   
   write(str,'("_",I3.3)') p_pe
  ! write(STDOUT,*) '../fort.82/fort.82'//trim(adjustl(str))
  ! write(STDOUT,*) '../fort.146/fort.146'//trim(adjustl(str))

   IF(Tradition) then
   
    OPEN(IN_ETA_FILE,FILE='../../fort.82/fort.82'//trim(adjustl(str)),    &
      &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

    OPEN(IN_TEM_FILE,FILE='../../fort.71/fort.71'//trim(adjustl(str)),    &
      &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

    IF (SALINITY) THEN
     OPEN(IN_SAL_FILE,FILE='../../fort.72/fort.72'//trim(adjustl(str)),    &
       &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')    
    ENDIF
 
    IF (Pressure) THEN  
     OPEN(IN_POO_FILE,FILE='../../fort.76/fort.76'//trim(adjustl(str)),    &
       &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    ENDIF
   
    IF (EastWest) THEN
     OPEN(IN_UKO_FILE,FILE='../../fort.73/fort.73'//trim(adjustl(str)),    &
       &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    ENDIF
    
    IF (NorthSouth) THEN
     OPEN(IN_VKE_FILE,FILE='../../fort.74/fort.74'//trim(adjustl(str)),    &
       &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
    ENDIF
    
   ENDIF
   
   IF (VERTICAL) THEN
    OPEN(IN_WWW_FILE,FILE='../../fort.146/fort.146'//trim(adjustl(str)), &
     &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
   ENDIF 

  
   if (BAROTROPIC) then
   
     IF (EastWest) THEN
       OPEN(IN_UZO_FILE,FILE='../../fort.63/fort.63'//trim(adjustl(str)),    &
         &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
     ENDIF
     
     IF (NorthSouth) THEN
       OPEN(IN_VSE_FILE,FILE='../../fort.64/fort.64'//trim(adjustl(str)),    &
         &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
     ENDIF
   ENDIF  

   IF (Nonhydrostatic) THEN
    OPEN(IN_PNH_FILE,FILE='../../fort.77/fort.77'//trim(adjustl(str)), &
     &    STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
   ENDIF 


   !  READ Eta  
   CALL find_index_P      
   global_nx = ie_g
   global_ny = 1           ! Purely 2D simulation
  ! global_ny = je_g  - 2  ! Periodic Boundary Condition
   global_nz = KE
   
   if (Tradition) then
   DO T=1,TE
      READ(IN_ETA_FILE) I4DATE
      READ(IN_ETA_FILE) ETA    
      
      dummyE(:,:) = ETA(iss:iee,jss:jee)
            
      write(str,'("_",I6.6)') t            
      err = nfmpi_create(MPI_COMM_WORLD,'./globalfiles/eta'//  &
       &  trim(adjustl(str))//'.nc',cmode,MPI_INFO_NULL,ncid)
       call check(err, 'In nfmpi_create: ')    
        
      ! define dimensions x and y
      err = nfmpi_def_dim(ncid, "lon", global_nx, dimidxy(1))
       call check(err, 'In nfmpi_def_dim x: ')
      err = nfmpi_def_dim(ncid, "lat", global_ny, dimidxy(2))
       call check(err, 'In nfmpi_def_dim y: ')    
         
      ! define a 2D variable of Real type
      err = nfmpi_def_var(ncid, "ETA", NF_DOUBLE, 2, dimidxy, varid)
       call check(err, 'In nfmpi_def_var: ')
       
      ! do not forget to exit define mode
      err = nfmpi_enddef(ncid)
       call check(err, 'In nfmpi_enddef: ')     
      ! now we are in data mode
      ! Note that in Fortran, array indices start with 1               
      err = nfmpi_put_vara_real_all(ncid, varid, starts_2d,  &
      &  counts_2d, dummyE, nType, subarray_p_2d)
       call check(err, 'In nfmpi_put_var_all: ')     

      err = nfmpi_close(ncid)
      call check(err, 'In nfmpi_close: ')
                
      IF(p_pe .eq. p_io) WRITE(STDOUT,*) 'Eta Timestep = ',t 
      
   ENDDO
   
   CALL MPI_Type_free(subarray_p_2d, err)      


   !  READ TS and Calculate Rou   
   DO t=1,TE 
      DO K=1,KE  
         ! tem     
         READ(IN_TEM_FILE) I4DATE
         READ(IN_TEM_FILE) dummyT
         TTT(:,:,k) = dummyT
         
         ! sal 
        if (Salinity) then         
         READ(IN_SAL_FILE) I4DATE
         READ(IN_SAL_FILE) dummyS
         SSS(:,:,k) = dummyS
        endif    
            
         ! po
        if (Pressure) then
         READ(IN_POO_FILE) I4DATE
         READ(IN_POO_FILE) dummyP 
         POO(:,:,k) = dummyP
        endif     
        
         ! rou
         if (Density .and. Salinity .and. Pressure) then          
           dummyP = dummyP/100._sp        
           CALL ADISITJ(dummyT,dummyS,dummyP)
           CALL RHO(dummyT,dummyS,dummyP,dummyR) 
           ROU(:,:,k) = dummyR
         endif           
      ENDDO
      
      
      write(str,'("_",I6.6)') t  
                
      err = nfmpi_create(MPI_COMM_WORLD,'./globalfiles/Tspr'//  &
       &  trim(adjustl(str))//'.nc',cmode,MPI_INFO_NULL,ncid)
      call check(err, 'In nfmpi_create: ')     
      ! define dimensions x ,y, z  
      err = nfmpi_def_dim(ncid, "lon", global_nx, dimidxyz(1))
      call check(err, 'In nfmpi_def_dim x: ')
      err = nfmpi_def_dim(ncid, "lat", global_ny, dimidxyz(2))
      call check(err, 'In nfmpi_def_dim y: ')  
      err = nfmpi_def_dim(ncid, "z", global_nz, dimidxyz(3))
      call check(err, 'In nfmpi_def_dim z: ')              
      ! define a 3D variable of Real type

       err = nfmpi_def_var(ncid, "THO", NF_DOUBLE, 3, &
        &  dimidxyz, varidT)
       call check(err, 'In nfmpi_def_var: ')
    
      if (Salinity) then
       err = nfmpi_def_var(ncid, "SAO", NF_DOUBLE, 3,    &
        &  dimidxyz, varidS)
       call check(err, 'In nfmpi_def_var: ')
      endif
    
      if (Pressure) then       
       err = nfmpi_def_var(ncid, "PO", NF_DOUBLE, 3,    &
        &  dimidxyz, varidP)
       call check(err, 'In nfmpi_def_var: ')   
      endif

      if (Density) then           
       err = nfmpi_def_var(ncid, "ROU", NF_DOUBLE, 3,     &
        &  dimidxyz, varidR)
       call check(err, 'In nfmpi_def_var: ')   
      endif    
                 
      ! do not forget to exit define mode
      err = nfmpi_enddef(ncid)
      call check(err, 'In nfmpi_enddef: ')     
  
      dummyE_3d(:,:,:)=TTT(iss:iee,jss:jee,1:ke)                               
      err = nfmpi_put_vara_real_all(ncid, varidT, starts_3d,   &
       &  counts_3d, dummyE_3d, nTypes, subarray_p_3d)
      call check(err, 'In nfmpi_put_var_all: ')     

      if (Salinity) then
       dummyE_3d(:,:,:)=SSS(iss:iee,jss:jee,1:ke)       
       err = nfmpi_put_vara_real_all(ncid, varidS, starts_3d,   &
        &  counts_3d, dummyE_3d, nTypes, subarray_p_3d)
       call check(err, 'In nfmpi_put_var_all: ')      
      endif
      
      if (Pressure) then 
       dummyE_3d(:,:,:)=POO(iss:iee,jss:jee,1:ke)          
       err = nfmpi_put_vara_real_all(ncid, varidP, starts_3d,   &
        &  counts_3d, dummyE_3d, nTypes, subarray_p_3d)
       call check(err, 'In nfmpi_put_var_all: ')
      endif
      
      if (Density) then 
       dummyE_3d(:,:,:)=ROU(iss:iee,jss:jee,1:ke)       	   
       err = nfmpi_put_vara_real_all(ncid, varidR, starts_3d,  &
        &  counts_3d, dummyE_3d, nTypes, subarray_p_3d)
       call check(err, 'In nfmpi_put_var_all: ')
      endif
                   
      err = nfmpi_close(ncid)
      call check(err, 'In nfmpi_close: ')
                
      IF(p_pe .eq. p_io) WRITE(STDOUT,*) 'Tspr Timestep = ',t    
    
   ENDDO      
   ENDIF
   
   ! Read Pnh
   if ( nonhydrostatic ) then
     DO t=1,TE
       ! Pnh   
        DO K=1,KE  
         READ(IN_PNH_FILE) I4DATE
         READ(IN_PNH_FILE) dummyP
         POO(:,:,k) = dummyP
       ENDDO  
   
       dummyE_3d(:,:,:)=POO(iss:iee,jss:jee,1:ke) 
       
       write(str,'("_",I6.6)') t                   
       err = nfmpi_create(MPI_COMM_WORLD,'./globalfiles/pnh'//  &
        &  trim(adjustl(str))//'.nc',cmode,MPI_INFO_NULL,ncid)
       call check(err, 'In nfmpi_create: ')      
       err = nfmpi_def_dim(ncid, "lon", global_nx, dimidxyz(1))
       call check(err, 'In nfmpi_def_dim x: ')
       err = nfmpi_def_dim(ncid, "lat", global_ny, dimidxyz(2))
       call check(err, 'In nfmpi_def_dim y: ')                
       err = nfmpi_def_dim(ncid, "z", global_nz, dimidxyz(3))
       call check(err, 'In nfmpi_def_dim z: ')                 
       err = nfmpi_def_var(ncid, "PNH", NF_DOUBLE,  &
        & 3, dimidxyz, varid)
       call check(err, 'In nfmpi_def_var: ')                        
       err = nfmpi_enddef(ncid)
       call check(err, 'In nfmpi_enddef: ')     
        
       err = nfmpi_put_vara_real_all(ncid, varid, starts_3d,   &
        &  counts_3d, dummyE_3d, nTypes, subarray_p_3d)
       call check(err, 'In nfmpi_put_var_all: ')  
          
       err = nfmpi_close(ncid)
       call check(err, 'In nfmpi_close: ')
                
       IF(p_pe .eq. p_io) WRITE(STDOUT,*) 'PNH Timestep = ',t       

     enddo
   endif
   
   
   !  READ W 
   if (VERTICAL) then 
    DO t=1,TE
       ! woo    	
       DO K=1,KE  
         READ(IN_WWW_FILE) I4DATE
         READ(IN_WWW_FILE) WWW(:,:,K)     
       ENDDO
       
       dummyE_3d = WWW(iss:iee,jss:jee,1:ke)
       
       write(str,'("_",I6.6)') t                   
       err = nfmpi_create(MPI_COMM_WORLD,'./globalfiles/woo'//  &
        &  trim(adjustl(str))//'.nc',cmode,MPI_INFO_NULL,ncid)
       call check(err, 'In nfmpi_create: ')      
       err = nfmpi_def_dim(ncid, "lon", global_nx, dimidxyz(1))
       call check(err, 'In nfmpi_def_dim x: ')
       err = nfmpi_def_dim(ncid, "lat", global_ny, dimidxyz(2))
       call check(err, 'In nfmpi_def_dim y: ')                
       err = nfmpi_def_dim(ncid, "z", global_nz, dimidxyz(3))
       call check(err, 'In nfmpi_def_dim z: ')                 
       err = nfmpi_def_var(ncid, "WO", NF_DOUBLE,  &
        & 3, dimidxyz, varid)
       call check(err, 'In nfmpi_def_var: ')                        
       err = nfmpi_enddef(ncid)
       call check(err, 'In nfmpi_enddef: ')     
        
       err = nfmpi_put_vara_real_all(ncid, varid, starts_3d,   &
        &  counts_3d, dummyE_3d, nTypes, subarray_p_3d)
       call check(err, 'In nfmpi_put_var_all: ')  
          
       err = nfmpi_close(ncid)
       call check(err, 'In nfmpi_close: ')
                
       IF(p_pe .eq. p_io) WRITE(STDOUT,*) 'WOO Timestep = ',t       
    
    ENDDO  
   ENDIF
         
   CALL MPI_Type_free(subarray_p_3d, err)   
   
  !  READ U
   IF (Tradition) THEN
      IF (EastWest) THEN
       CALL find_index_U
       global_nx = ie_g + 1
       global_ny = 1          ! Purely 2D simulation
     !  global_ny = je_g - 2  !  Periodic Boundary Condition
       DO T = 1, TE 
       ! uuu    	
       DO K=1,KE
         READ(IN_UKO_FILE) I4DATE
         READ(IN_UKO_FILE) UUU(:,:,K)     
       ENDDO
       
       dummyE_3d_U = UUU(iss:iee,jss:jee,1:ke)
       write(str,'("_",I6.6)') t                   
       err = nfmpi_create(MPI_COMM_WORLD,'./globalfiles/uko'//  &
        &  trim(adjustl(str))//'.nc',cmode,MPI_INFO_NULL,ncid)
       call check(err, 'In nfmpi_create: ')      
       err = nfmpi_def_dim(ncid, "lon", global_nx, dimidxyz(1))
       call check(err, 'In nfmpi_def_dim x: ')
       err = nfmpi_def_dim(ncid, "lat", global_ny, dimidxyz(2))
       call check(err, 'In nfmpi_def_dim y: ')                
       err = nfmpi_def_dim(ncid, "z", global_nz, dimidxyz(3))
       call check(err, 'In nfmpi_def_dim z: ')     
            
       err = nfmpi_def_var(ncid, "UKO", NF_DOUBLE,  &
        & 3, dimidxyz, varid)
       call check(err, 'In nfmpi_def_var: ')                        
       err = nfmpi_enddef(ncid)
       call check(err, 'In nfmpi_enddef: ')     
           
       err = nfmpi_put_vara_real_all(ncid, varid, starts_3d,   &
        &  counts_3d, dummyE_3d_U, ntypes, subarray_u)
       call check(err, 'In nfmpi_put_var_all: ')  
                
       err = nfmpi_close(ncid)
       call check(err, 'In nfmpi_close: ')
                
       IF(p_pe .eq. p_io) WRITE(STDOUT,*) 'UKO Timestep = ',t          
              
       ENDDO       
       CALL MPI_Type_free(subarray_u, err)     
      ENDIF
   ENDIF
   
   IF (BAROTROPIC) THEN
  ! Read the U
    IF (EastWest) THEN
     CALL find_index_U_2D
     global_nx = ie_g + 1
     global_ny = 1        ! Purely 2D simulation
   ! global_ny = je_g - 2 ! Periodic Boundary Condition
     DO T = 1, TE 
       READ(IN_UZO_FILE) I4DATE
       READ(IN_UZO_FILE) UUU(:,:,1)     
      
       dummyE_2d_U(:,:) = UUU(iss:iee,jss:jee,1)
       write(str,'("_",I6.6)') t                   
       err = nfmpi_create(MPI_COMM_WORLD,'./globalfiles/uzo'//  &
        &  trim(adjustl(str))//'.nc',cmode,MPI_INFO_NULL,ncid)
       call check(err, 'In nfmpi_create: ')      
       err = nfmpi_def_dim(ncid, "lon", global_nx, dimidxy(1))
       call check(err, 'In nfmpi_def_dim x: ')
       err = nfmpi_def_dim(ncid, "lat", global_ny, dimidxy(2))
       call check(err, 'In nfmpi_def_dim y: ')                     
       err = nfmpi_def_var(ncid, "UZO", NF_DOUBLE,  &
        & 2, dimidxy, varid)
       call check(err, 'In nfmpi_def_var: ')                        
       err = nfmpi_enddef(ncid)
       call check(err, 'In nfmpi_enddef: ')     
           
      err = nfmpi_put_vara_real_all(ncid, varid, starts_2d,   &
       &  counts_2d, dummyE_2d_U, nTypes, subarray_u_2d)
      call check(err, 'In nfmpi_put_var_all: ')  
                
      err = nfmpi_close(ncid)
      call check(err, 'In nfmpi_close: ')
                
      IF(p_pe .eq. p_io) WRITE(STDOUT,*) 'UZO Timestep = ',t          
              
     ENDDO       
     CALL MPI_Type_free(subarray_u_2d, err)     
    ENDIF
   ENDIF
   

!  READ V 
   IF (Tradition) THEN 
     IF (NorthSouth) THEN
      CALL find_index_V   
      global_nx = ie_g
      global_ny = 2                ! Purely 2D simulation
     ! global_ny = je_g + 1  - 2   ! Periodic Boundary Condition  
      DO T = 1, TE 		
       DO K = 1, KE
         READ(IN_VKE_FILE) I4DATE
         READ(IN_VKE_FILE) VVV(:,:,K)     	
       ENDDO      
       dummyE_3d_V = VVV(iss:iee,jss:jee,1:KE)
       
       write(str,'("_",I6.6)') t                  
       err = nfmpi_create(MPI_COMM_WORLD,'./globalfiles/vke'//  &
        &  trim(adjustl(str))//'.nc',cmode,MPI_INFO_NULL,ncid)
       call check(err, 'In nfmpi_create: ')      
       err = nfmpi_def_dim(ncid, "lon", global_nx, dimidxyz(1))
       call check(err, 'In nfmpi_def_dim x: ')
       err = nfmpi_def_dim(ncid, "lat", global_ny, dimidxyz(2))
       call check(err, 'In nfmpi_def_dim y: ')                
       err = nfmpi_def_dim(ncid, "z", global_nz, dimidxyz(3))
       call check(err, 'In nfmpi_def_dim z: ')                 
       err = nfmpi_def_var(ncid, "VKE", NF_DOUBLE,  &
        & 3, dimidxyz, varid)
       call check(err, 'In nfmpi_def_var: ')                        
       err = nfmpi_enddef(ncid)
       call check(err, 'In nfmpi_enddef: ')     
           
       err = nfmpi_put_vara_real_all(ncid, varid, starts_3d,   &
        &  counts_3d, dummyE_3d_V, nTypes, subarray_v)
       call check(err, 'In nfmpi_put_var_all: ')  
          
       err = nfmpi_close(ncid)
       call check(err, 'In nfmpi_close: ')
                
       IF(p_pe .eq. p_io) WRITE(STDOUT,*) 'VKE Timestep = ',t          
              
      ENDDO       
      CALL MPI_Type_free(subarray_v, err)   
     ENDIF     
   ENDIF 
 
 
   IF ( BAROTROPIC ) THEN
  ! Read the V
    IF (NorthSouth) THEN
     CALL find_index_V_2D
     global_nx = ie_g 
     global_ny = 2  ! Purely 2D simulation
   ! global_ny = je_g +1  - 2 ! Periodic Boundary Condition

      DO T = 1, TE 
       READ(IN_VSE_FILE) I4DATE
       READ(IN_VSE_FILE) VVV(:,:,1)     
      
       dummyE_2d_V(:,:) = VVV(iss:iee,jss:jee,1)
       write(str,'("_",I6.6)') t                   
       err = nfmpi_create(MPI_COMM_WORLD,'./globalfiles/vse'//  &
        &  trim(adjustl(str))//'.nc',cmode,MPI_INFO_NULL,ncid)
       call check(err, 'In nfmpi_create: ')      
       err = nfmpi_def_dim(ncid, "lon", global_nx, dimidxy(1))
       call check(err, 'In nfmpi_def_dim x: ')
       err = nfmpi_def_dim(ncid, "lat", global_ny, dimidxy(2))
       call check(err, 'In nfmpi_def_dim y: ')                
       err = nfmpi_def_var(ncid, "VSE", NF_DOUBLE,  &
        & 2, dimidxy, varid)
       call check(err, 'In nfmpi_def_var: ')                        
       err = nfmpi_enddef(ncid)
       call check(err, 'In nfmpi_enddef: ')     
           
       err = nfmpi_put_vara_real_all(ncid, varid, starts_2d,   &
        &  counts_2d, dummyE_2d_V, nTypes, subarray_v_2d)
       call check(err, 'In nfmpi_put_var_all: ')  
                
       err = nfmpi_close(ncid)
       call check(err, 'In nfmpi_close: ')
                
       IF(p_pe .eq. p_io) WRITE(STDOUT,*) 'VSE Timestep = ',t          
              
      ENDDO    
      CALL MPI_Type_free(subarray_v_2d, err)  
     ENDIF
   ENDIF
           
   CALL MPI_FINALIZE(IERR)   
   stop  
   
 END PROGRAM GLUE   
