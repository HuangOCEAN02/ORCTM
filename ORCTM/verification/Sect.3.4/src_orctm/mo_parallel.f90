MODULE MO_PARALLEL
  USE mo_param1, ONLY: ke, ie, je, ie_g, je_g, ICYCLI_X, ICYCLI_Y, I_start, J_start, Z3DOC
  USE mo_mpi
  USE mo_kind
  use mo_units
  IMPLICIT NONE

  INTEGER, PARAMETER :: maxproc = 1024

! Number of subdivisions in x/y direction

  INTEGER nprocx, nprocy
  INTEGER nprocxy

! Our own offset

  INTEGER p_ioff, p_joff

! Flag if we have the boundaries

  LOGICAL have_g_is, have_g_ie, have_g_js, have_g_je

! Start of the inner domains in x/y direction
! p_lim_x(0) = 2
! p_lim_x(i) = start of inner domain i
! p_lim_x(nprocx) = ie_g

  INTEGER p_lim_x(0:maxproc), p_lim_y(0:maxproc)
  PRIVATE p_lim_x, p_lim_y

! For every processor: number in x/y direction (0-based)

  INTEGER p_num_x(0:maxproc), p_num_y(0:maxproc)
  PRIVATE p_num_x, p_num_y

! Global offsets and sizes for each processor (both for outer domains)

  INTEGER p_ioff_g(0:maxproc), p_joff_g(0:maxproc)
  INTEGER p_size_x(0:maxproc), p_size_y(0:maxproc)
  INTEGER,allocatable:: p_Ioff_GG(:),p_Joff_GG(:)    ! H.H for diagnosis output
  PRIVATE p_ioff_g, p_joff_g, p_size_x, p_size_y

  REAL*8  :: t2d=0, t3d=0, ts, te
  INTEGER :: n2d=0, n3d=0
  PRIVATE t2d, t3d, ts, te, n2d, n3d

  INTERFACE bounds_exch
     MODULE PROCEDURE bounds_exch_2d
     MODULE PROCEDURE bounds_exch_3d
  END INTERFACE

  INTERFACE global_sum
     MODULE PROCEDURE global_sum_i
     MODULE PROCEDURE global_sum_r
     MODULE PROCEDURE global_sum_1d
     MODULE PROCEDURE global_sum_2d
  END INTERFACE

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE p_deco

!   domain decomposition

    IMPLICIT NONE
    INTEGER i, nx, ny

    IF(p_nprocs > 1) THEN

      WRITE(nerr,*) 'Process ',p_pe,' of ',p_nprocs,' is alive'

      ! set some variables

      IF(p_pe==0) THEN

         ! nprocx and nprocy must be set by the calling process

         IF(nprocx==0 .OR. nprocy==0) THEN
           WRITE(nerr,*) 'ERROR: nprocx or nprocy not set'
           CALL p_abort
         ENDIF

         nprocxy = nprocx*nprocy

         IF(nprocxy /= p_nprocs .AND. nprocxy /= p_nprocs-1) THEN
            WRITE(nerr,*)'Number of processors = ',p_nprocs
            WRITE(nerr,*)'nprocx = ',nprocx,' nprocy = ',nprocy
            WRITE(nerr,*)'Number of processors doesnt fit!'
            CALL p_abort
         ENDIF

         IF(((ie_g-2)/nprocx)<3 .or. ((je_g-2)/nprocy)<3) THEN
            WRITE(nerr,*)'Decomposed domain gets too small'
            WRITE(nerr,*)'We need at least 3 rows in every direction'
            CALL p_abort
         ENDIF

      ENDIF

      ! broadcast nprocx and nprocy

      CALL p_bcast(nprocx,0)
      CALL p_bcast(nprocy,0)
      nprocxy = nprocx*nprocy

      ! Decomposition - compute domain limits

      DO i=0,nprocx
        p_lim_x(i) = 2 + i*(ie_g-2)/nprocx
      ENDDO
      DO i=0,nprocy
        p_lim_y(i) = 2 + i*(je_g-2)/nprocy
      ENDDO

      ! Set number of processors in x and y direction

      DO i=0,nprocx-1
         p_num_x(i:nprocxy-1:nprocx) = i
      ENDDO
      DO i=0,nprocy-1
        p_num_y(i*nprocx:(i+1)*nprocx-1) = i
      ENDDO

      ! Offsets and sizes

      DO i=0,nprocxy-1
        nx = p_num_x(i)
        ny = p_num_y(i)

        p_ioff_g(i) = p_lim_x(nx)-2
        p_joff_g(i) = p_lim_y(ny)-2
        p_size_x(i) = p_lim_x(nx+1) - p_lim_x(nx) + 2
        p_size_y(i) = p_lim_y(ny+1) - p_lim_y(ny) + 2
      ENDDO

       ALLOCATE(p_Ioff_GG(0:nprocxy-1),p_Joff_GG(0:nprocxy-1))
      ! Get our own values
      
      p_Ioff_GG(:) = p_ioff_g(0:nprocxy-1)
      p_Joff_GG(:) = p_joff_g(0:nprocxy-1)     
      
      CALL p_bcast(p_Ioff_GG,0)
      CALL p_bcast(p_Joff_GG,0)
      
      IF(p_pe<nprocxy) THEN
        ie = p_size_x(p_pe)
        je = p_size_y(p_pe)
        p_ioff = p_ioff_g(p_pe)
        p_joff = p_joff_g(p_pe)
        have_g_is = p_num_x(p_pe) == 0
        have_g_ie = p_num_x(p_pe) == nprocx-1
        have_g_js = p_num_y(p_pe) == 0
        have_g_je = p_num_y(p_pe) == nprocy-1
      ELSE
        ie = ie_g
        je = je_g
        p_ioff = 0
        p_joff = 0
        have_g_is = .TRUE.
        have_g_ie = .TRUE.
        have_g_js = .TRUE.
        have_g_je = .TRUE.
      ENDIF

      WRITE(nerr,'(a,i4,a,2i4,a,2i4)')'Proc ',p_pe,' offset: ',p_ioff,p_joff, &
        ' Size: ',ie,je

    ELSE

      nprocx = 1
      nprocy = 1
      nprocxy = 1
      p_lim_x(0) = 2
      p_lim_x(1) = ie_g
      p_lim_y(0) = 2
      p_lim_y(1) = je_g
      p_num_x(0) = 0
      p_num_y(0) = 0
      p_ioff_g(0) = 0
      p_joff_g(0) = 0
      p_size_x(0) = ie_g
      p_size_y(0) = je_g
      ie = ie_g
      je = je_g
      p_ioff = 0
      p_joff = 0
      have_g_is = .TRUE.
      have_g_ie = .TRUE.
      have_g_js = .TRUE.
      have_g_je = .TRUE.

      WRITE(nerr,*) 'Running on a single processor'

    ENDIF
    if (have_g_is) then
      I_start=0
    else
      I_start=1
    endif
    if (have_g_js) then
      J_start=0
    else
      J_start=1
    endif
   END SUBROUTINE p_deco

!-----------------------------------------------------------------------

   SUBROUTINE gather_arr(arrl, arrg, pe, ivec)

      ! Gathers all local array parts into a global array on PE pe

      IMPLICIT NONE

      REAL, INTENT(IN)  :: arrl(:,:)
      REAL, INTENT(INOUT) :: arrg(:,:)
      INTEGER, INTENT(IN) :: pe
      INTEGER, INTENT(IN), OPTIONAL :: ivec

      INTEGER n, iis, iie, jjs, jje, ioff, joff, I_start_G, J_start_G
      REAL, ALLOCATABLE :: aux(:,:), auxg(:,:)

      I_start_G=1
      J_start_G=1
      if (present(ivec)) then
        if (ivec .eq. 1) then
          I_start_G=0
        elseif (ivec .eq. 2) then
          J_start_G=0
        elseif (ivec .eq. 0) then
        else
          CALL STOP_ALL('mo_parallel, wrong ivec number')
        endif
      endif

      IF(p_pe/=pe) THEN
         IF(p_pe<nprocxy) CALL p_send(arrl,pe,1111)
      ELSE
         ALLOCATE(auxg(I_start_G:IE_G,J_start_G:JE_G))
         DO n=0,nprocxy-1
            iis = 1
            iie = p_size_x(n)
            IF (p_num_x(n).eq.0) iis=I_start_G

            jjs = 1
            jje = p_size_y(n)
            IF (p_num_y(n).eq.0) jjs=J_start_G

            ALLOCATE(aux(iis:iie,jjs:jje))

            IF (n==pe) THEN
               aux(:,:) = arrl(:,:)
            ELSE
               CALL p_recv(aux,n,1111)
            ENDIF
            ! Copy only outer limits into arrg
            IF(p_num_x(n) /= 0        ) iis = iis+1
            IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
            IF(p_num_y(n) /= 0        ) jjs = jjs+1
            IF(p_num_y(n) /= nprocy-1 ) jje = jje-1
            ioff = p_ioff_g(n)
            joff = p_joff_g(n)
            auxg(ioff+iis:ioff+iie,joff+jjs:joff+jje) = aux(iis:iie,jjs:jje)

            DEALLOCATE(aux)
         ENDDO
         arrg(:,:)=auxg(:,:)
         DEALLOCATE(auxg)
      ENDIF

      RETURN
   END SUBROUTINE gather_arr

!-----------------------------------------------------------------------

   SUBROUTINE scatter_arr(arrg, arrl, pe, ivec)

      ! Scatters global array (arrg) on PE pe to the local array (arrl) on all PEs
      ! Since this routine is not used in perfomance critical parts
      ! we use the simple broadcast version

      IMPLICIT NONE

      REAL, INTENT(INOUT) :: arrg(:,:)
      REAL, INTENT(INOUT) :: arrl(:,:)
      INTEGER, INTENT(IN) :: pe
      INTEGER, INTENT(IN), OPTIONAL :: ivec

      INTEGER :: I_start_G, J_start_G
      REAL, ALLOCATABLE :: aux(:,:), auxg(:,:)

      I_start_G=1
      J_start_G=1
      if (present(ivec)) then
        if (ivec .eq. 1) then
          I_start_G=0
          ALLOCATE(aux(I_start:IE,1:JE))
        elseif (ivec .eq. 2) then
          J_start_G=0
          ALLOCATE(aux(1:IE,J_start:JE))
        elseif (ivec .eq. 0) then
          ALLOCATE(aux(1:IE,1:JE))
        else
          CALL STOP_ALL('mo_parallel, wrong ivec number')
        endif
      else
        ALLOCATE(aux(1:IE,1:JE))
      endif
      ALLOCATE(auxg(I_start_G:IE_G,J_start_G:JE_G))

      auxg(:,:)=arrg(:,:)
      CALL p_bcast(auxg,pe)

      if (present(ivec)) then
        if (ivec .eq. 1) then
          aux(:,:) = auxg(p_ioff+I_start:p_ioff+ie,p_joff+1:p_joff+je)
        elseif (ivec .eq. 2) then
          aux(:,:) = auxg(p_ioff+1:p_ioff+ie,p_joff+J_start:p_joff+je)
        elseif (ivec .eq. 0) then
          aux(:,:) = auxg(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)
        else
          CALL STOP_ALL('mo_parallel, wrong ivec number')
        endif
      else
        aux(:,:) = auxg(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)
      endif

      arrl(:,:)=aux(:,:)

      DEALLOCATE(aux,auxg)

   END SUBROUTINE scatter_arr

!-----------------------------------------------------------------------

   SUBROUTINE bounds_exch_2d(a0,a1,a2,a3,a4,a5)

      ! Exchanges boundaries of 2D arrays
      ! Since 2D boundary exchange has an impact on performance,
      ! we allow that up to 6 2D arrays are exchanged at the same time

! peterspy: Note that a0 to a9 should be located at the same point (one of p/u/v)

      REAL, INTENT(INOUT) :: a0(:,:)
      REAL, INTENT(INOUT), OPTIONAL :: a1(:,:),a2(:,:),a3(:,:),a4(:,:),a5(:,:)
      REAL, DIMENSION(:,:), ALLOCATABLE :: a01,a11,a21,a31,a41,a51
      REAL, DIMENSION(:,:), ALLOCATABLE :: xr1,xr2,xs1,xs2,yr1,yr2,ys1,ys2
      INTEGER :: nm,np,n,i_x1,i_x2,i_y1,i_y2,J_INDEX

      IF (UBOUND(a0,1) .eq. (ie+1)) THEN
        ALLOCATE(a01(0:ie,je))
        ALLOCATE(xr1(je,6), xr2(je,6), xs1(je,6), xs2(je,6))
        ALLOCATE(yr1(0:ie,6), yr2(0:ie,6), ys1(0:ie,6), ys2(0:ie,6))
      ELSEIF (UBOUND(a0,2) .eq. (je+1)) THEN
        ALLOCATE(a01(ie,0:je))
        ALLOCATE(xr1(0:je,6), xr2(0:je,6), xs1(0:je,6), xs2(0:je,6))
        ALLOCATE(yr1(ie,6), yr2(ie,6), ys1(ie,6), ys2(ie,6))
      ELSE
        ALLOCATE(a01(ie,je))
        ALLOCATE(xr1(je,6), xr2(je,6), xs1(je,6), xs2(je,6))
        ALLOCATE(yr1(ie,6), yr2(ie,6), ys1(ie,6), ys2(ie,6))
      ENDIF

      a01(:,:)=a0(:,:)
      i_x1=LBOUND(a01,1);i_x2=UBOUND(a01,1);i_y1=LBOUND(a01,2);i_y2=UBOUND(a01,2);
      IF(PRESENT(a1)) THEN; ALLOCATE(a11(i_x1:i_x2,i_y1:i_y2)); a11(:,:)=a1(:,:); ENDIF
      IF(PRESENT(a2)) THEN; ALLOCATE(a21(i_x1:i_x2,i_y1:i_y2)); a21(:,:)=a2(:,:); ENDIF
      IF(PRESENT(a3)) THEN; ALLOCATE(a31(i_x1:i_x2,i_y1:i_y2)); a31(:,:)=a3(:,:); ENDIF
      IF(PRESENT(a4)) THEN; ALLOCATE(a41(i_x1:i_x2,i_y1:i_y2)); a41(:,:)=a4(:,:); ENDIF
      IF(PRESENT(a5)) THEN; ALLOCATE(a51(i_x1:i_x2,i_y1:i_y2)); a51(:,:)=a5(:,:); ENDIF

#ifdef TIMECHECK
      ts = MPI_Wtime()
#endif

      ! x-direction

      xs1(:,1) = a01(2,:)
      xs2(:,1) = a01(ie-1,:)
      n = 1
      IF(PRESENT(a1)) THEN; n=n+1; xs1(:,n) = a11(2,:); xs2(:,n) = a11(ie-1,:); ENDIF
      IF(PRESENT(a2)) THEN; n=n+1; xs1(:,n) = a21(2,:); xs2(:,n) = a21(ie-1,:); ENDIF
      IF(PRESENT(a3)) THEN; n=n+1; xs1(:,n) = a31(2,:); xs2(:,n) = a31(ie-1,:); ENDIF
      IF(PRESENT(a4)) THEN; n=n+1; xs1(:,n) = a41(2,:); xs2(:,n) = a41(ie-1,:); ENDIF
      IF(PRESENT(a5)) THEN; n=n+1; xs1(:,n) = a51(2,:); xs2(:,n) = a51(ie-1,:); ENDIF

      IF (nprocx>1 .AND. p_pe<nprocxy ) THEN
        ! Get processor numbers of neighbors
        nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
        np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)

#ifdef bounds_exch_isend
        CALL p_isend(xs1(:,1:n),nm,1)
        CALL p_isend(xs2(:,1:n),np,2)
        CALL p_recv(xr1(:,1:n),np,1)
        CALL p_recv(xr2(:,1:n),nm,2)
        CALL p_wait
#else
        CALL p_sendrecv(xs1(:,1:n),nm,xr1(:,1:n),np,1)
        CALL p_sendrecv(xs2(:,1:n),np,xr2(:,1:n),nm,2)
#endif
      ELSE
        xr1(:,1:n) = xs1(:,1:n)
        xr2(:,1:n) = xs2(:,1:n)
      ENDIF

      IF(ICYCLI_X/=0 .OR. .NOT. have_g_is) THEN
         a01( 1,:) = xr2(:,1)
         n = 1
         IF(PRESENT(a1)) THEN; n=n+1; a11( 1,:) = xr2(:,n); ENDIF
         IF(PRESENT(a2)) THEN; n=n+1; a21( 1,:) = xr2(:,n); ENDIF
         IF(PRESENT(a3)) THEN; n=n+1; a31( 1,:) = xr2(:,n); ENDIF
         IF(PRESENT(a4)) THEN; n=n+1; a41( 1,:) = xr2(:,n); ENDIF
         IF(PRESENT(a5)) THEN; n=n+1; a51( 1,:) = xr2(:,n); ENDIF
      ENDIF
      IF(ICYCLI_X/=0 .OR. .NOT. have_g_ie) THEN
        a01(ie,:) = xr1(:,1)
        n = 1
        IF(PRESENT(a1)) THEN; n=n+1; a11(ie,:) = xr1(:,n); ENDIF
        IF(PRESENT(a2)) THEN; n=n+1; a21(ie,:) = xr1(:,n); ENDIF
        IF(PRESENT(a3)) THEN; n=n+1; a31(ie,:) = xr1(:,n); ENDIF
        IF(PRESENT(a4)) THEN; n=n+1; a41(ie,:) = xr1(:,n); ENDIF
        IF(PRESENT(a5)) THEN; n=n+1; a51(ie,:) = xr1(:,n); ENDIF
      ENDIF

      ! y-direction
      ! Note that there is no action required if nprocy==1
      ! IF nprocy>1 changed by peterspy, same as x-direction
      IF (nprocy>1 .AND. p_pe<nprocxy) THEN

        ys1(:,1) = a01(:,2)
        ys2(:,1) = a01(:,je-1)
        n = 1

        IF(PRESENT(a1)) THEN; n=n+1; ys1(:,n) = a11(:,2); ys2(:,n) = a11(:,je-1); ENDIF
        IF(PRESENT(a2)) THEN; n=n+1; ys1(:,n) = a21(:,2); ys2(:,n) = a21(:,je-1); ENDIF
        IF(PRESENT(a3)) THEN; n=n+1; ys1(:,n) = a31(:,2); ys2(:,n) = a31(:,je-1); ENDIF
        IF(PRESENT(a4)) THEN; n=n+1; ys1(:,n) = a41(:,2); ys2(:,n) = a41(:,je-1); ENDIF
        IF(PRESENT(a5)) THEN; n=n+1; ys1(:,n) = a51(:,2); ys2(:,n) = a51(:,je-1); ENDIF

! H.H - Get processor numbers of neighbors
! H.H - For the sake of simplicity  the p_sendrecv call is as
! H.H - for periodic boundary conditions
        nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
        np = MOD(p_pe+nprocxy+nprocx,nprocxy)

#ifdef bounds_exch_isend
        CALL p_isend(ys1(:,1:n),nm,1)
        CALL p_isend(ys2(:,1:n),np,2)
        CALL p_recv(yr1(:,1:n),np,1)
        CALL p_recv(yr2(:,1:n),nm,2)
        CALL p_wait
#else
        CALL p_sendrecv(ys1(:,1:n),nm,yr1(:,1:n),np,1)
        CALL p_sendrecv(ys2(:,1:n),np,yr2(:,1:n),nm,2)
#endif
! H.H - N-S Periodic Boundary Condition for single processor in y-direction
! H.H - changed by H.H
      ELSEIF(nprocy == 1 .AND. p_pe<nprocxy .AND. ICYCLI_Y .EQ. 1) then 
       
        IF (UBOUND(a0,2) .eq. (je+1)) THEN      	
            J_INDEX = 1      
        ELSE
            J_INDEX = 2  
        ENDIF   
              
        ys1(:,1) = a01(:,J_INDEX)
        ys2(:,1) = a01(:,je-1)
        n = 1

        IF(PRESENT(a1)) THEN; n=n+1; ys1(:,n) = a11(:,J_INDEX); ys2(:,n) = a11(:,je-1); ENDIF
        IF(PRESENT(a2)) THEN; n=n+1; ys1(:,n) = a21(:,J_INDEX); ys2(:,n) = a21(:,je-1); ENDIF
        IF(PRESENT(a3)) THEN; n=n+1; ys1(:,n) = a31(:,J_INDEX); ys2(:,n) = a31(:,je-1); ENDIF
        IF(PRESENT(a4)) THEN; n=n+1; ys1(:,n) = a41(:,J_INDEX); ys2(:,n) = a41(:,je-1); ENDIF
        IF(PRESENT(a5)) THEN; n=n+1; ys1(:,n) = a51(:,J_INDEX); ys2(:,n) = a51(:,je-1); ENDIF
                         
        yr1(:,1:n) = ys1(:,1:n)
        yr2(:,1:n) = ys2(:,1:n)
! H.H - Record the j=2&4, make j=2&4,j=5&1 exchanged message      
      ELSE
        yr1(:,1:n) = ys1(:,1:n)
        yr2(:,1:n) = ys2(:,1:n)
      ENDIF

! H.H - ICYCLI_Y is assigned by one when using N-S Periodic Boundary Condition 
    IF (nprocy == 1 .AND. ICYCLI_Y .EQ. 1) THEN  
     	
           IF (UBOUND(a0,2) .eq. (je+1)) THEN      	
              J_INDEX = 0     
           ELSE
              J_INDEX = 1 
           ENDIF   
          
           a01(:, J_INDEX) = yr2(:,1)
           n = 1
           IF(PRESENT(a1)) THEN; n=n+1; a11(:, J_INDEX) = yr2(:,n); ENDIF
           IF(PRESENT(a2)) THEN; n=n+1; a21(:, J_INDEX) = yr2(:,n); ENDIF
           IF(PRESENT(a3)) THEN; n=n+1; a31(:, J_INDEX) = yr2(:,n); ENDIF
           IF(PRESENT(a4)) THEN; n=n+1; a41(:, J_INDEX) = yr2(:,n); ENDIF
           IF(PRESENT(a5)) THEN; n=n+1; a51(:, J_INDEX) = yr2(:,n); ENDIF

           a01(:,je) = yr1(:,1)
           n = 1
           IF(PRESENT(a1)) THEN; n=n+1; a11(:,je) = yr1(:,n); ENDIF
           IF(PRESENT(a2)) THEN; n=n+1; a21(:,je) = yr1(:,n); ENDIF
           IF(PRESENT(a3)) THEN; n=n+1; a31(:,je) = yr1(:,n); ENDIF
           IF(PRESENT(a4)) THEN; n=n+1; a41(:,je) = yr1(:,n); ENDIF
           IF(PRESENT(a5)) THEN; n=n+1; a51(:,je) = yr1(:,n); ENDIF
    ELSE       
! original
        IF(.NOT. have_g_js) THEN
           a01(:, 1) = yr2(:,1)
           n = 1
           IF(PRESENT(a1)) THEN; n=n+1; a11(:, 1) = yr2(:,n); ENDIF
           IF(PRESENT(a2)) THEN; n=n+1; a21(:, 1) = yr2(:,n); ENDIF
           IF(PRESENT(a3)) THEN; n=n+1; a31(:, 1) = yr2(:,n); ENDIF
           IF(PRESENT(a4)) THEN; n=n+1; a41(:, 1) = yr2(:,n); ENDIF
           IF(PRESENT(a5)) THEN; n=n+1; a51(:, 1) = yr2(:,n); ENDIF
        ENDIF
        IF(.NOT. have_g_je) THEN
          a01(:,je) = yr1(:,1)
          n = 1
          IF(PRESENT(a1)) THEN; n=n+1; a11(:,je) = yr1(:,n); ENDIF
          IF(PRESENT(a2)) THEN; n=n+1; a21(:,je) = yr1(:,n); ENDIF
          IF(PRESENT(a3)) THEN; n=n+1; a31(:,je) = yr1(:,n); ENDIF
          IF(PRESENT(a4)) THEN; n=n+1; a41(:,je) = yr1(:,n); ENDIF
          IF(PRESENT(a5)) THEN; n=n+1; a51(:,je) = yr1(:,n); ENDIF
        ENDIF
    ENDIF


#ifdef TIMECHECK
      te = MPI_Wtime()
#endif
      t2d = t2d + te-ts
      n2d = n2d + 1

      a0(:,:)=a01(:,:)
      IF(PRESENT(a1)) a1(:,:)=a11(:,:)
      IF(PRESENT(a2)) a2(:,:)=a21(:,:)
      IF(PRESENT(a3)) a3(:,:)=a31(:,:)
      IF(PRESENT(a4)) a4(:,:)=a41(:,:)
      IF(PRESENT(a5)) a5(:,:)=a51(:,:)

   END SUBROUTINE bounds_exch_2d

!-----------------------------------------------------------------------

   SUBROUTINE bounds_exch_3d(arr)

      ! Exchanges boundaries of 3D arrays

      REAL, INTENT(INOUT) :: arr(:,:,:)
      REAL, ALLOCATABLE :: arr1(:,:,:)
      REAL, DIMENSION(:,:), ALLOCATABLE :: x3r1,x3r2,xs1,xs2,y3r1,y3r2,ys1,ys2
      INTEGER nm, np, kk

      kk = UBOUND(arr,3)
      IF (UBOUND(arr,1) .eq. (ie+1)) THEN
        ALLOCATE(arr1(0:ie,je,kk))
        ALLOCATE(x3r1(je,kk), x3r2(je,kk), xs1(je,kk), xs2(je,kk))
        ALLOCATE(y3r1(0:ie,kk), y3r2(0:ie,kk), ys1(0:ie,kk), ys2(0:ie,kk))
      ELSEIF (UBOUND(arr,2) .eq. (je+1)) THEN
        ALLOCATE(arr1(ie,0:je,kk))
        ALLOCATE(x3r1(0:je,kk), x3r2(0:je,kk), xs1(0:je,kk), xs2(0:je,kk))
        ALLOCATE(y3r1(ie,kk), y3r2(ie,kk), ys1(ie,kk), ys2(ie,kk))
      ELSE
        ALLOCATE(arr1(ie,je,kk))
        ALLOCATE(x3r1(je,kk), x3r2(je,kk), xs1(je,kk), xs2(je,kk))
        ALLOCATE(y3r1(ie,kk), y3r2(ie,kk), ys1(ie,kk), ys2(ie,kk))
      ENDIF
      arr1(:,:,:)=arr(:,:,:)

#ifdef TIMECHECK
      ts = MPI_Wtime()
#endif

      ! x-direction

      xs1(:,:) = arr1(2,:,:)
      xs2(:,:) = arr1(ie-1,:,:)

      IF (nprocx>1 .AND. p_pe<nprocxy) THEN
        ! Get processor numbers of neighbors
        nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
        np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)
#ifdef bounds_exch_isend
        CALL p_isend(xs1,nm,1)
        CALL p_isend(xs2,np,2)
        CALL p_recv(x3r1,np,1)
        CALL p_recv(x3r2,nm,2)
        CALL p_wait
#else
        CALL p_sendrecv(xs1,nm,x3r1,np,1)
        CALL p_sendrecv(xs2,np,x3r2,nm,2)
#endif
      ELSE
        x3r1 = xs1
        x3r2 = xs2
      ENDIF

      IF(ICYCLI_X/=0 .OR. p_num_x(p_pe)/=0)        arr1( 1,:,:) = x3r2(:,:)
      IF(ICYCLI_X/=0 .OR. p_num_x(p_pe)/=nprocx-1) arr1(ie,:,:) = x3r1(:,:)

      ! y-direction
      ! Note that there is no action required if nprocy==1

      IF (nprocy>1 .AND. p_pe<nprocxy) THEN
        ys1(:,:) = arr1(:,2,:)
        ys2(:,:) = arr1(:,je-1,:)
! H.H - Get processor numbers of neighbors
! H.H - For the sake of simplicity  the p_sendrecv call is as
! H.H - for periodic boundary conditions
        nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
        np = MOD(p_pe+nprocxy+nprocx,nprocxy)
#ifdef bounds_exch_isend
        CALL p_isend(ys1,nm,1)
        CALL p_isend(ys2,np,2)
        CALL p_recv(y3r1,np,1)
        CALL p_recv(y3r2,nm,2)
        CALL p_wait
#else
        CALL p_sendrecv(ys1,nm,y3r1,np,1)
        CALL p_sendrecv(ys2,np,y3r2,nm,2)
#endif
! H.H&QS - N-S Periodic Boundary Condition for single processor in y-direction
! H.H - Record the j=2&4, make j=2&4,j=5&1 exchanged message  
      ELSEIF(nprocy == 1 .AND. p_pe<nprocxy .AND. ICYCLI_Y .EQ. 1) then 
      	
        IF (UBOUND(arr,2) .eq. (je+1)) THEN      	
            y3r1(:,:) = arr1(:,1,:)        
        ELSE
            y3r1(:,:) = arr1(:,2,:)
        ENDIF         
            y3r2(:,:) = arr1(:,je-1,:)     
      ELSE
        y3r1 = ys1
        y3r2 = ys2
      ENDIF

! H.H - ICYCLI_Y is assigned by one when using N-S Periodic Boundary Condition 
     IF (nprocy == 1 .AND. ICYCLI_Y .EQ. 1) THEN    	
      	  IF (UBOUND(arr,2) .eq. (je+1)) THEN          
              arr1(:,0,:) = y3r2(:,:)
          ELSE
              arr1(:,1,:) = y3r2(:,:)
          ENDIF                
      	  arr1(:,je,:) = y3r1(:,:)
      	     
     ELSE

! H.H - ICYCLI_Y is assigned by zero when not using N-S Periodic Boundary Condition
      IF(p_num_y(p_pe)/=0)        arr1(:, 1,:) = y3r2(:,:)
      IF(p_num_y(p_pe)/=nprocy-1) arr1(:,je,:) = y3r1(:,:)
      
     ENDIF

#ifdef TIMECHECK
      te = MPI_Wtime()
#endif
      t3d = t3d + te-ts
      n3d = n3d + 1

      arr(:,:,:)=arr1(:,:,:)

   END SUBROUTINE bounds_exch_3d

!-----------------------------------------------------------------------

   SUBROUTINE read_slice(iunit,arr,ivec)

      ! Reads a 2D array and scatters it to all

      INTEGER,INTENT(IN) :: iunit
      REAL,INTENT(INOUT) :: arr(:,:)
      INTEGER,INTENT(IN),OPTIONAL :: ivec

      INTEGER ::I_start_G,J_start_G
      REAL,ALLOCATABLE :: arr0(:,:),arr_g(:,:)
      
! H.H --for reading the restart file Z37000&Z38000 perfectly    
      REAL(Kind=sp),ALLOCATABLE :: dummy_g(:,:)

      I_start_G=1
      J_start_G=1
      if (present(ivec)) then
        if (ivec .eq. 1) then
          I_start_G=0
          ALLOCATE(arr0(I_start:IE,1:JE))
        elseif (ivec .eq. 2) then
          J_start_G=0
          ALLOCATE(arr0(1:IE,J_start:JE))
        elseif (ivec .eq. 0) then
          ALLOCATE(arr0(1:IE,1:JE))
        else
          CALL STOP_ALL('mo_parallel, wrong ivec number')
        endif
      else
        ALLOCATE(arr0(1:IE,1:JE))
      endif

      ALLOCATE(arr_g(I_start_G:ie_g,J_start_G:je_g))   
      ALLOCATE(dummy_g(I_start_G:ie_g,J_start_G:je_g))

      IF(p_pe==p_io) THEN
        ! for restart
        IF (Z3DOC .eq. 3) then
            READ(iunit) dummy_g
            arr_g = real(dummy_g,8)
        ENDIF
        ! for no restart
        IF (Z3DOC .eq. 1) then
            READ(iunit) arr_g
        ENDIF
      ENDIF

      if (present(ivec)) then
        CALL scatter_arr(arr_g,arr0,p_io,ivec)
      else
        CALL scatter_arr(arr_g,arr0,p_io)
      endif

      arr(:,:)=arr0(:,:)

      DEALLOCATE(arr0,arr_g,dummy_g)

   END SUBROUTINE read_slice

!-----------------------------------------------------------------------

   SUBROUTINE write_slice(iunit,arr,ivec)

      ! Gathers a 2D array and writes it to iunit

      INTEGER,INTENT(IN) :: iunit
      REAL,INTENT(IN) :: arr(:,:)
      INTEGER,INTENT(IN),OPTIONAL :: ivec

      INTEGER::I_start_G,J_start_G

      REAL,ALLOCATABLE :: arr0(:,:),arr_g(:,:)

      I_start_G=1
      J_start_G=1
      if (present(ivec)) then
        if (ivec .eq. 1) then
          I_start_G=0
          ALLOCATE(arr0(I_start:IE,1:JE))
        elseif (ivec .eq. 2) then
          J_start_G=0
          ALLOCATE(arr0(1:IE,J_start:JE))
        elseif (ivec .eq. 0) then
          ALLOCATE(arr0(1:IE,1:JE))
        else
          CALL STOP_ALL('mo_parallel, wrong ivec number')
        endif
      else
        ALLOCATE(arr0(1:IE,1:JE))
      endif

      arr0(:,:)=arr(:,:)

      ALLOCATE(arr_g(I_start_G:ie_g,J_start_G:je_g))

      if (present(ivec)) then
        CALL gather_arr(arr,arr_g,p_io,ivec)
      else
        CALL gather_arr(arr,arr_g,p_io)
      endif

      IF(p_pe==p_io) THEN
        WRITE(iunit) real(arr_g,sp)
      ENDIF

      DEALLOCATE(arr0,arr_g)

   END SUBROUTINE write_slice

!-----------------------------------------------------------------------

   SUBROUTINE global_sum_r(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9)

   ! Build global sum for real scalars
   ! For performance reasons we permit up to 10 arguments in a single call

   REAL, INTENT(INOUT):: s0
   REAL, INTENT(INOUT), OPTIONAL:: s1,s2,s3,s4,s5,s6,s7,s8,s9

   REAL s(10)
   INTEGER n

   s(1) = s0
   n = 1
   IF (PRESENT(s1)) THEN; n = n+1; s(n) = s1; ENDIF
   IF (PRESENT(s2)) THEN; n = n+1; s(n) = s2; ENDIF
   IF (PRESENT(s3)) THEN; n = n+1; s(n) = s3; ENDIF
   IF (PRESENT(s4)) THEN; n = n+1; s(n) = s4; ENDIF
   IF (PRESENT(s5)) THEN; n = n+1; s(n) = s5; ENDIF
   IF (PRESENT(s6)) THEN; n = n+1; s(n) = s6; ENDIF
   IF (PRESENT(s7)) THEN; n = n+1; s(n) = s7; ENDIF
   IF (PRESENT(s8)) THEN; n = n+1; s(n) = s8; ENDIF
   IF (PRESENT(s9)) THEN; n = n+1; s(n) = s9; ENDIF

   CALL global_sum_1d(s(1:n))

   s0 = s(1)
   n = 1
   IF (PRESENT(s1)) THEN; n = n+1; s1 = s(n); ENDIF
   IF (PRESENT(s2)) THEN; n = n+1; s2 = s(n); ENDIF
   IF (PRESENT(s3)) THEN; n = n+1; s3 = s(n); ENDIF
   IF (PRESENT(s4)) THEN; n = n+1; s4 = s(n); ENDIF
   IF (PRESENT(s5)) THEN; n = n+1; s5 = s(n); ENDIF
   IF (PRESENT(s6)) THEN; n = n+1; s6 = s(n); ENDIF
   IF (PRESENT(s7)) THEN; n = n+1; s7 = s(n); ENDIF
   IF (PRESENT(s8)) THEN; n = n+1; s8 = s(n); ENDIF
   IF (PRESENT(s9)) THEN; n = n+1; s9 = s(n); ENDIF

   END SUBROUTINE global_sum_r

!-----------------------------------------------------------------------

   SUBROUTINE global_sum_i(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9)

   ! Build global sum for integer scalars
   ! For performance reasons we permit up to 10 arguments in a single call

   INTEGER, INTENT(INOUT):: s0
   INTEGER, INTENT(INOUT), OPTIONAL:: s1,s2,s3,s4,s5,s6,s7,s8,s9

   REAL s(10)
   INTEGER n

   s(1) = s0
   n = 1
   IF (PRESENT(s1)) THEN; n = n+1; s(n) = s1; ENDIF
   IF (PRESENT(s2)) THEN; n = n+1; s(n) = s2; ENDIF
   IF (PRESENT(s3)) THEN; n = n+1; s(n) = s3; ENDIF
   IF (PRESENT(s4)) THEN; n = n+1; s(n) = s4; ENDIF
   IF (PRESENT(s5)) THEN; n = n+1; s(n) = s5; ENDIF
   IF (PRESENT(s6)) THEN; n = n+1; s(n) = s6; ENDIF
   IF (PRESENT(s7)) THEN; n = n+1; s(n) = s7; ENDIF
   IF (PRESENT(s8)) THEN; n = n+1; s(n) = s8; ENDIF
   IF (PRESENT(s9)) THEN; n = n+1; s(n) = s9; ENDIF

   CALL global_sum_1d(s(1:n))

   s0 = NINT(s(1))
   n = 1
   IF (PRESENT(s1)) THEN; n = n+1; s1 = NINT(s(n)); ENDIF
   IF (PRESENT(s2)) THEN; n = n+1; s2 = NINT(s(n)); ENDIF
   IF (PRESENT(s3)) THEN; n = n+1; s3 = NINT(s(n)); ENDIF
   IF (PRESENT(s4)) THEN; n = n+1; s4 = NINT(s(n)); ENDIF
   IF (PRESENT(s5)) THEN; n = n+1; s5 = NINT(s(n)); ENDIF
   IF (PRESENT(s6)) THEN; n = n+1; s6 = NINT(s(n)); ENDIF
   IF (PRESENT(s7)) THEN; n = n+1; s7 = NINT(s(n)); ENDIF
   IF (PRESENT(s8)) THEN; n = n+1; s8 = NINT(s(n)); ENDIF
   IF (PRESENT(s9)) THEN; n = n+1; s9 = NINT(s(n)); ENDIF

   END SUBROUTINE global_sum_i

!-----------------------------------------------------------------------

   SUBROUTINE global_sum_1d(s)

   ! Build global sum for real 1D array

   REAL, INTENT(INOUT):: s(:)

   REAL r(SIZE(s)), q(SIZE(s)), errmax, err
   INTEGER n, i, ierr

   ! If we are running in test mode, use result from last PE and
   ! check if the others are approximatly right, else do real summation

   n = SIZE(s)

   IF(p_nprocs>nprocxy) THEN

      ! Test mode
      IF(p_pe>=nprocxy) THEN
         q(:) = s(:) ! Save s
         s(:) = 0    ! Remove our contribution to global sum
      ENDIF

      ! Get sum on working PEs
      r = p_sum(s)

      ! Check against saved value on test PE
      IF(p_pe==nprocxy) THEN
         errmax = 0;
         DO i=1,n
           IF(q(i)/=r(i)) THEN
             err = ABS(q(i)-r(i))/MAX(ABS(q(i)),ABS(r(i)))
             errmax = MAX(errmax,err)
           ENDIF
         ENDDO
         WRITE(nerr,*) 'Global Sum Max Err = ',errmax
         IF(errmax>1.e-5) CALL p_abort
         s(:) = q(:) ! Restore s
      ENDIF

      ! We use the value from test PE in order to get always
      ! identical results during test
      CALL p_bcast(s,nprocxy)

   ELSE

      r = p_sum(s)
      s(:) = r(:)

   ENDIF

   END SUBROUTINE global_sum_1d

!-----------------------------------------------------------------------

   SUBROUTINE global_sum_2d(s)

   ! Build global sum for real 2D array

   REAL, INTENT(INOUT):: s(:,:)

   REAL r(SIZE(s))
   INTEGER ierr

   IF(p_nprocs>nprocxy) THEN
      ! Test mode - use global_sum_1d
      r = RESHAPE(s,(/SIZE(s)/))
      CALL global_sum_1d(r)
      s = RESHAPE(r,SHAPE(s))
   ELSE
      r = p_sum(RESHAPE(s,(/SIZE(s)/)))
      s = RESHAPE(r,SHAPE(s))
   ENDIF

   END SUBROUTINE global_sum_2d

!-----------------------------------------------------------------------

   SUBROUTINE global_max(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9)

   ! Build global max for real scalars
   ! For performance reasons we permit up to 10 arguments in a single call

   REAL, INTENT(INOUT):: s0
   REAL, INTENT(INOUT), OPTIONAL:: s1,s2,s3,s4,s5,s6,s7,s8,s9

   REAL s(10), r(10)
   INTEGER n, ierr

   s(1) = s0
   n = 1
   IF (PRESENT(s1)) THEN; n = n+1; s(n) = s1; ENDIF
   IF (PRESENT(s2)) THEN; n = n+1; s(n) = s2; ENDIF
   IF (PRESENT(s3)) THEN; n = n+1; s(n) = s3; ENDIF
   IF (PRESENT(s4)) THEN; n = n+1; s(n) = s4; ENDIF
   IF (PRESENT(s5)) THEN; n = n+1; s(n) = s5; ENDIF
   IF (PRESENT(s6)) THEN; n = n+1; s(n) = s6; ENDIF
   IF (PRESENT(s7)) THEN; n = n+1; s(n) = s7; ENDIF
   IF (PRESENT(s8)) THEN; n = n+1; s(n) = s8; ENDIF
   IF (PRESENT(s9)) THEN; n = n+1; s(n) = s9; ENDIF

   r = p_max(s)

   s0 = r(1)
   n = 1
   IF (PRESENT(s1)) THEN; n = n+1; s1 = r(n); ENDIF
   IF (PRESENT(s2)) THEN; n = n+1; s2 = r(n); ENDIF
   IF (PRESENT(s3)) THEN; n = n+1; s3 = r(n); ENDIF
   IF (PRESENT(s4)) THEN; n = n+1; s4 = r(n); ENDIF
   IF (PRESENT(s5)) THEN; n = n+1; s5 = r(n); ENDIF
   IF (PRESENT(s6)) THEN; n = n+1; s6 = r(n); ENDIF
   IF (PRESENT(s7)) THEN; n = n+1; s7 = r(n); ENDIF
   IF (PRESENT(s8)) THEN; n = n+1; s8 = r(n); ENDIF
   IF (PRESENT(s9)) THEN; n = n+1; s9 = r(n); ENDIF

   END SUBROUTINE global_max

!-----------------------------------------------------------------------

   SUBROUTINE global_min(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9)

   ! Build global min for real scalars
   ! For performance reasons we permit up to 10 arguments in a single call

   REAL, INTENT(INOUT):: s0
   REAL, INTENT(INOUT), OPTIONAL:: s1,s2,s3,s4,s5,s6,s7,s8,s9

   REAL s(10), r(10)
   INTEGER n, ierr

   s(1) = s0
   n = 1
   IF (PRESENT(s1)) THEN; n = n+1; s(n) = s1; ENDIF
   IF (PRESENT(s2)) THEN; n = n+1; s(n) = s2; ENDIF
   IF (PRESENT(s3)) THEN; n = n+1; s(n) = s3; ENDIF
   IF (PRESENT(s4)) THEN; n = n+1; s(n) = s4; ENDIF
   IF (PRESENT(s5)) THEN; n = n+1; s(n) = s5; ENDIF
   IF (PRESENT(s6)) THEN; n = n+1; s(n) = s6; ENDIF
   IF (PRESENT(s7)) THEN; n = n+1; s(n) = s7; ENDIF
   IF (PRESENT(s8)) THEN; n = n+1; s(n) = s8; ENDIF
   IF (PRESENT(s9)) THEN; n = n+1; s(n) = s9; ENDIF

   r = p_min(s)

   s0 = r(1)
   n = 1
   IF (PRESENT(s1)) THEN; n = n+1; s1 = r(n); ENDIF
   IF (PRESENT(s2)) THEN; n = n+1; s2 = r(n); ENDIF
   IF (PRESENT(s3)) THEN; n = n+1; s3 = r(n); ENDIF
   IF (PRESENT(s4)) THEN; n = n+1; s4 = r(n); ENDIF
   IF (PRESENT(s5)) THEN; n = n+1; s5 = r(n); ENDIF
   IF (PRESENT(s6)) THEN; n = n+1; s6 = r(n); ENDIF
   IF (PRESENT(s7)) THEN; n = n+1; s7 = r(n); ENDIF
   IF (PRESENT(s8)) THEN; n = n+1; s8 = r(n); ENDIF
   IF (PRESENT(s9)) THEN; n = n+1; s9 = r(n); ENDIF

   END SUBROUTINE global_min

!-----------------------------------------------------------------------

   REAL FUNCTION global_array_sum(arr)

   ! Builds the global sum of a 2D array by gathering it on one PE,
   ! calculating th sum on this PE and broadcasting the result
   ! As opposed to calculating the local sum and the using a global_sum call,
   ! this method should give identical results indepentently of the number
   ! of processors, at least when optimization is switched off.

   REAL, INTENT(IN) :: arr(:,:)

   INTEGER i,j
   REAL arr_g(ie_g,je_g), s, sum

   CALL gather_arr(arr,arr_g,0)

   IF(p_pe==0) THEN
      sum = 0
      DO j=2,je_g-1
        DO i=2,ie_g-1
          sum = sum + arr_g(i,j)
        ENDDO
      ENDDO
   ENDIF

   CALL p_bcast(sum,0)

   global_array_sum = sum

   IF(p_pe==nprocxy) THEN
      s = 0
      DO j=2,je-1
        DO i=2,ie-1
          s = s + arr(i,j)
        ENDDO
      ENDDO

      IF(s/=sum) THEN
         WRITE(nerr,*) 'global_array_sum: ',s,sum,abs(s-sum)
         call p_abort
      ENDIF
   ENDIF

   END FUNCTION global_array_sum

!-----------------------------------------------------------------------

   SUBROUTINE stop_all(text)

   ! Does an emergency stop by calling p_abort

   CHARACTER (LEN=*), INTENT(IN) :: text

   write(nerr,*) text

   CALL p_abort

   END SUBROUTINE stop_all

!-----------------------------------------------------------------------

   SUBROUTINE print_stats

   REAL*8 t2, t3
   INTEGER n

   IF(p_pe==p_io) THEN
     WRITE(nerr,*) 'Number of 2D boundary exchanges: ',n2d
     WRITE(nerr,*) 'Number of 3D boundary exchanges: ',n3d
   ENDIF

   DO n=0,nprocxy-1
     t2 = t2d
     t3 = t3d
     CALL p_bcast(t2,n)
     CALL p_bcast(t3,n)
     IF(p_pe==p_io) THEN
       WRITE(nerr,'(a,i4,a,2f10.3)') 'PE: ',n, &
         ' Times for 2D/3D boundary exchanges: ',t2,t3
     ENDIF
   ENDDO

   END SUBROUTINE print_stats

!-----------------------------------------------------------------------
END MODULE MO_PARALLEL
