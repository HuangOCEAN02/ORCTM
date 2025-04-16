      MODULE MO_GLUE
      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_KIND
#ifdef NON_HYDROSTATIC
      USE MO_PETSC
#endif
#ifdef GLOBAL_NETCDF
      INCLUDE "pnetcdf.inc"
      REAL(sp),ALLOCATABLE :: dummyE(:,:),dummyE_3d(:,:,:), &
      & dummyE_3d_U(:,:,:),dummyE_3d_V(:,:,:)
      CHARACTER(len=120) :: nter_str
      INTEGER::NITER

      !|-----------------BASIC--INFO--OF--Pnetcdf-----------------------|
      INTEGER::err,cmode,varid,ncid,dimidxy(2),dimidxyz(3)
      INTEGER*8::starts_2d(2), counts_2d(2)
      INTEGER*8::starts_3d_p(3), counts_3d_p(3)
      INTEGER*8::starts_3d_u(3), counts_3d_u(3)
      INTEGER*8::starts_3d_v(3), counts_3d_v(3)
      INTEGER::iss,iee,jss,jee,ius,iue,jus,jue,ivs,ive,jvs,jve
      INTEGER::subarray_p_2d,subarray_p_3d,subarray_u,subarray_v,nTypes
      INTEGER::array_of_sizes(2), array_of_subsizes(2)
      INTEGER::array_of_starts(2)
      INTEGER::array_of_sizes_3d(3), array_of_subsizes_3d(3)
      INTEGER::array_of_starts_3d(3)
      !|-----------------BASIC--INFO--OF--Pnetcdf--END------------------|

! IN-SITU Density calculation variables
#ifdef DENSITYIO
      REAL(sp),ALLOCATABLE ::TTO(:,:),SSO(:,:) &
      ,PPO(:,:),RRO(:,:)
#endif /*DENSITYIO*/

      CONTAINS

      SUBROUTINE alloc_mem_glue
      INCLUDE 'mpif.h'
      
! IN-SITU Density calculation variables
#ifdef DENSITYIO
     ALLOCATE(TTO(IE,JE),SSO(IE,JE) &
      ,PPO(IE,JE),RRO(IE,JE))
#endif /*DENSITYIO*/

      cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
      nTypes = 1
!-------------------------------------  p-fields default --------------------------------------------------
      				!    2-D Field  default because of the ghost points
      				! x-direction
      				iss = 2
      				iee = ie-1

      				starts_2d(1) = p_ioff+2
      				counts_2d(1) = ie-2

      				starts_3d_p(1) = p_ioff+2
      				counts_3d_p(1) = ie-2

      				! y-direction
      				jss = 2
      				jee = je-1

      				starts_2d(2) = p_joff+2
      				counts_2d(2) = je-2

      				starts_3d_p(2) = p_joff+2
      				counts_3d_p(2) = je-2

      				! z-direction
      				starts_3d_p(3) = 1
      				counts_3d_p(3) = ke


      				! West
      				if (have_g_is) then
      					iss = 1
      					iee = ie-1
      					starts_2d(1) = p_ioff+1
      					counts_2d(1) = ie-1
      					starts_3d_p(1) = p_ioff+1
      					counts_3d_p(1) = ie-1
      				endif
      				! East
      				if (have_g_ie) then
      					iss = 2
      					iee = ie
      					starts_2d(1) = p_ioff+2
      					counts_2d(1) = ie-1
      					starts_3d_p(1) = p_ioff+2
      					counts_3d_p(1) = ie-1
      				endif

      			! North
      			 if (have_g_js) then
      						jss = 1
      						jee = je-1
      						starts_2d(2) = p_joff+1
      						counts_2d(2) = je-1
      						starts_3d_p(2) = p_joff+1
      						counts_3d_p(2) = je-1
      			endif

      			! South
      			if (have_g_je) then
      						jss = 2
      						jee = je
      						starts_2d(2) = p_joff+2
      						counts_2d(2) = je-1
      						starts_3d_p(2) = p_joff+2
      						counts_3d_p(2) = je-1
      			endif


      				ALLOCATE(dummyE(counts_2d(1),counts_2d(2)))


      				! define an MPI datatype using MPI_Type_create_subarray()
      				array_of_sizes(1)    = ie_g
      				array_of_sizes(2)    = je_g ! - 2 ! Periodic Boundary Condition

      				array_of_subsizes(1) = counts_2d(1)
      				array_of_subsizes(2) = counts_2d(2)

      				array_of_starts(1) = 0
      				array_of_starts(2) = 0

      			call MPI_Type_create_subarray(2, array_of_sizes, &
      			&   array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, &
      			&   MPI_REAL, subarray_p_2d, err)
      			call MPI_Type_commit(subarray_p_2d, err)

      			ALLOCATE(dummyE_3d(counts_3d_p(1),counts_3d_p(2),counts_3d_p(3)))
	          ! define an MPI datatype using MPI_Type_create_subarray()
      				array_of_sizes_3d(1)    = ie_g
      				array_of_sizes_3d(2)    = je_g ! - 2 ! Periodic Boundary Condition
      				array_of_sizes_3d(3)    = ke
      				array_of_subsizes_3d(1) = counts_3d_p(1)
      				array_of_subsizes_3d(2) = counts_3d_p(2)
      				array_of_subsizes_3d(3) = counts_3d_p(3)
      				array_of_starts_3d(1)   = 0   ! MPI start index starts with 0
      				array_of_starts_3d(2)   = 0
      				array_of_starts_3d(3)   = 0
      		   call MPI_Type_create_subarray(3, array_of_sizes_3d, &
      		& array_of_subsizes_3d, array_of_starts_3d, MPI_ORDER_FORTRAN, &
      		& MPI_REAL, subarray_p_3d, err)
      		   call MPI_Type_commit(subarray_p_3d, err)

!-------------------------------------  u-fields default --------------------------------------------------
!         u-fields excluded ghost points
      				! x-direction
      				ius = 2
      				iue = ie-1
      				starts_3d_u(1) = p_ioff+3 ! Global matrix from 1
      				counts_3d_u(1) = ie-2

      				! y-direction
              jus = 2
              jue = je-1
      				starts_3d_u(2) = p_joff+2
      				counts_3d_u(2) = je-2

      				! z-direction
      				starts_3d_u(3) = 1
      				counts_3d_u(3) = ke

      				! West
      				if (have_g_is) then
      					ius = 0      ! U-point
      				  iue = ie-1
      					starts_3d_u(1) = p_ioff+1  ! Global matrix from 1
      					counts_3d_u(1) = ie   ! U-point
      				endif
      				! East
      				if (have_g_ie) then
                 ius =2
                 iue = ie
      					starts_3d_u(1) = p_ioff+3
      					counts_3d_u(1) = ie-1
      				endif

      				! IF(nprocy .ne. 1) then

      					! North
      					if (have_g_js) then
                   jus = 1
                   jue = je-1
      						starts_3d_u(2) = p_joff+1
      						counts_3d_u(2) = je-1
      					endif
      					! South
      					if (have_g_je) then
                   jus = 2
                   jue = je
      						starts_3d_u(2) = p_joff+2
      						counts_3d_u(2) = je-1
      					endif

      				! ELSE ! non_parallel in y-direction
      				! starts_3d_u(2) = p_joff+1
      				! 	counts_3d_u(2) = je !  - 2  ! Periodic Boundary Condition
      				! ENDIF

      				ALLOCATE(dummyE_3d_U(counts_3d_u(1),counts_3d_u(2),counts_3d_u(3)))

      				! define an MPI datatype using MPI_Type_create_subarray()
      				array_of_sizes_3d(1)    = ie_g+1
      				array_of_sizes_3d(2)    = je_g ! - 2  ! Periodic Boundary Condition
      				array_of_sizes_3d(3)    = ke
      				array_of_subsizes_3d(1) = counts_3d_u(1)
      				array_of_subsizes_3d(2) = counts_3d_u(2)
      				array_of_subsizes_3d(3) = counts_3d_u(3)
      				array_of_starts_3d(1)   = 0   ! MPI start index starts with 0
      				array_of_starts_3d(2)   = 0
      				array_of_starts_3d(3)   = 0
      			  call MPI_Type_create_subarray(3, array_of_sizes_3d, &
      				&   array_of_subsizes_3d, array_of_starts_3d, MPI_ORDER_FORTRAN, &
      				&   MPI_REAL, subarray_u, err)
      				call MPI_Type_commit(subarray_u, err)

!-------------------------------------  v-fields default --------------------------------------------------
      				! x-direction
      				ivs = 2
      				ive = ie-1
              starts_3d_v(1) = p_ioff+2
      				counts_3d_v(1) = ie-2

      				! y-direction
      				jvs = 2
      				jve = je-1
      				starts_3d_v(2) = p_joff+3 ! Global matrix from 1
      				counts_3d_v(2) = je-2

      				! z-direction
      				starts_3d_v(3) = 1
      				counts_3d_v(3) = ke

      				! West
      				if (have_g_is) then
                ivs = 1
                ive = ie-1
      					starts_3d_v(1) = p_ioff+1
      					counts_3d_v(1) = ie-1
      				endif
      				! East
      				if (have_g_ie) then
                 ivs = 2
                 ive = ie
      					starts_3d_v(1) = p_ioff+2
      					counts_3d_v(1) = ie-1
      				endif

      				! IF(nprocy .ne. 1) then
      					! North
      					if (have_g_js) then
      						jvs = 0      ! V-point
                  jve = je-1
      						starts_3d_v(2) = p_joff+1 ! Global matrix from 1
      						counts_3d_v(2) = je  ! V-point
      					endif

      					! South
      					if (have_g_je) then
                  jvs = 2
                  jve = je
      						starts_3d_v(2) = p_joff+3 ! Global matrix from 1
      						counts_3d_v(2) = je-1
      					endif

      				! ELSE ! non_parallel in y-direction
      				! 	jss = J_start        ! V-point
      				! 	starts_3d_v(2) = p_joff+1
      				! 	counts_3d(2) = je+1 ! - 2  ! V-point & Periodic Boundary Condition
      				! ENDIF

      				ALLOCATE(dummyE_3d_V(counts_3d_v(1),counts_3d_v(2),counts_3d_v(3)))

      				! define an MPI datatype using MPI_Type_create_subarray()
      				array_of_sizes_3d(1)    = ie_g
      				array_of_sizes_3d(2)    = je_g + 1 ! - 2  Periodic Boundary Condition
      				array_of_sizes_3d(3)    = ke
      				array_of_subsizes_3d(1) = counts_3d_v(1)
      				array_of_subsizes_3d(2) = counts_3d_v(2)
      				array_of_subsizes_3d(3) = counts_3d_v(3)
      				array_of_starts_3d(1)   = 0   ! MPI start index starts with 0
      				array_of_starts_3d(2)   = 0
      				array_of_starts_3d(3)   = 0
      			  call MPI_Type_create_subarray(3, array_of_sizes_3d, &
      				&   array_of_subsizes_3d, array_of_starts_3d, MPI_ORDER_FORTRAN, &
      				&   MPI_REAL, subarray_v, err)
      				call MPI_Type_commit(subarray_v, err)

      END SUBROUTINE alloc_mem_glue


      SUBROUTINE GLUE_END

      CALL MPI_Type_free(subarray_p_2d, err)
      CALL MPI_Type_free(subarray_p_3d, err)
      CALL MPI_Type_free(subarray_u, err)
      CALL MPI_Type_free(subarray_v, err)

      END SUBROUTINE  GLUE_END

      SUBROUTINE check(ierr, message)
      INCLUDE 'mpif.h'
      INTEGER ierr
      character message*(*)

      ! It is a good idea to check returned value for possible error
      if (err .NE. NF_NOERR) then
        write(IO_STDOUT,*)  message//' '//nfmpi_strerror(ierr)
      	CALL sleep(120)
      	CALL p_abort
#ifdef NON_HYDROSTATIC
      	CALL petsc_end
#endif
      	CALL ABSTURZ
      endif
     END SUBROUTINE check ! subroutine check

      SUBROUTINE MEAN2D_GLOBAL(kdays,kmonts,kyears, FIELD2D, I4CODE,&
               & ivec)
      USE MO_PARAM1
      include 'mpif.h'
      INTEGER,INTENT(IN)::kdays,kmonts,kyears,I4CODE,ivec
      REAL,INTENT(IN)::FIELD2D(:,:)

      REAL(sp),ALLOCATABLE::FIELD21(:,:)
      INTEGER I_start_G,J_start_G,ICODE,ITER
      INTEGER*8::global_nx, global_ny,global_nz
      CHARACTER(LEN=3) :: VAR

      ICODE=I4CODE
      ITER = NITER

      I_start_G=1
      J_start_G=1
      if (ivec .eq. 1) then
      	I_start_G=0
        allocate(FIELD21(I_start:IE,1:JE))
        FIELD21(:,:) = FIELD2D(:,:)
      	dummyE_3d_U(:,:,1) = FIELD21(ius:iue,jus:jue)
      	elseif (ivec .eq. 2) then
      		J_start_G=0
          allocate(FIELD21(1:IE,J_start:JE))
          FIELD21(:,:) = FIELD2D(:,:)
      		dummyE_3d_V(:,:,1) = FIELD21(ivs:ive,jvs:jve)
      		elseif (ivec .eq. 0) then
            allocate(FIELD21(1:IE,1:JE))
            FIELD21(:,:) = FIELD2D(:,:)
      			dummyE(:,:) = FIELD21(iss:iee,jss:jee)
      		else
      			CALL STOP_ALL('mo_parallel, wrong ivec number')
      		endif
         
         ! This is for IO testing of u and v fields 
         ! IF (have_g_is .AND. ICODE .EQ.  203 ) THEN
         !   WRITE(IO_STDOUT,*) 'dummyE_3d_U(ie,:,1)', dummyE_3d_U(ie,:,1)
         !   WRITE(IO_STDOUT,*) 'FIELD21(iue,:)', FIELD21(iue,:)
         !   WRITE(IO_STDOUT,*) 'ius,iue,jus,jue', ius,iue, jus,jue
         ! ENDIF   
          
         ! IF (have_g_js .AND. ICODE .EQ.  203 .and. p_pe ==1 ) THEN
         !   WRITE(IO_STDOUT,*) 'dummyE_3d_U(1,:,1)', dummyE_3d_U(1,:,1)
         !   WRITE(IO_STDOUT,*) 'FIELD21(ius,:)', FIELD21(ius,:)
         !   WRITE(IO_STDOUT,*) 'ius,iue,jus,jue', ius,iue, jus,jue
         ! ENDIF   
                    
          global_nx =  INT(IE_G+1-I_start_G)
          global_ny =  INT(JE_G+1-J_start_G)
          global_nz =  1

      		WRITE(nter_str,'("_",I8.8 )') ITER

      		if (ICODE .EQ.  1 ) then
      			VAR = 'ZOO'
      			err = nfmpi_create(MPI_COMM_WORLD,'./zoo/ZOO'//  &
      			&  trim(adjustl(nter_str))//'.nc',cmode,MPI_INFO_NULL,ncid)
      		endif

      		if (ICODE .EQ.  203 ) then
      			VAR = 'USO'
      			err = nfmpi_create(MPI_COMM_WORLD,'./uso/USO'//  &
      			&  trim(adjustl(nter_str))//'.nc',cmode,MPI_INFO_NULL,ncid)

      		endif

      		if (ICODE .EQ.  204 ) then
      			VAR = 'VSE'
      			err = nfmpi_create(MPI_COMM_WORLD,'./vse/VSE'//  &
      			&  trim(adjustl(nter_str))//'.nc',cmode,MPI_INFO_NULL,ncid)
      		endif
      		CALL check(err, 'In nfmpi_create: ')

      		IF (ICODE .EQ.  203 .OR. ICODE .EQ.  204 ) THEN
      		  ! define dimensions x and y and z
      		  err = nfmpi_def_dim(ncid, "lon", global_nx, dimidxyz(1))
      		   call check(err, 'In nfmpi_def_dim x: ')
      		  err = nfmpi_def_dim(ncid, "lat", global_ny, dimidxyz(2))
      		   call check(err, 'In nfmpi_def_dim y: ')
            err = nfmpi_def_dim(ncid, "z", global_nz, dimidxyz(3))
             call check(err, 'In nfmpi_def_dim z: ')
       		  ! define a 3D variable of Real type
      		  err = nfmpi_def_var(ncid, VAR, NF_DOUBLE, 3, dimidxyz, varid)
      		   call check(err, 'In nfmpi_def_var: ')
      		ENDIF

      		if (ICODE .EQ.  1 ) THEN
      		  ! define dimensions x and y
      		  err = nfmpi_def_dim(ncid, "lon", global_nx, dimidxy(1))
      		   call check(err, 'In nfmpi_def_dim x: ')
      		  err = nfmpi_def_dim(ncid, "lat", global_ny, dimidxy(2))
      		   call check(err, 'In nfmpi_def_dim y: ')
      		  ! define a 2D variable of Real type
      		  err = nfmpi_def_var(ncid, VAR, NF_DOUBLE, 2, dimidxy, varid)
      		   call check(err, 'In nfmpi_def_var: ')
      		ENDIF

      		! do not forget to exit define mode
      		err = nfmpi_enddef(ncid)
      		call check(err, 'In nfmpi_enddef: ')

      		! Now we are in data mode
      		! Note that in Fortran, array indices start with 1
      		if (ICODE .EQ.  1 ) then
      			err = nfmpi_put_vara_real_all(ncid, varid, starts_2d,  &
      			&  counts_2d, dummyE, nType, subarray_p_2d)
      		ENDIF
      		if (ICODE .EQ.  203 ) then
      			counts_3d_u(3) = 1
      			err = nfmpi_put_vara_real_all(ncid, varid, starts_3d_u,  &
      			&  counts_3d_u, dummyE_3d_U(:,:,1))
      		ENDIF
      		if (ICODE .EQ.  204 ) then
      			counts_3d_v(3) = 1
      			err = nfmpi_put_vara_real_all(ncid, varid, starts_3d_v,  &
      			&  counts_3d_v, dummyE_3d_V(:,:,1))
      		ENDIF
      		call check(err, 'In nfmpi_put_var_all: ')

      		err = nfmpi_close(ncid)
      		call check(err, 'In nfmpi_close: ')

      		IF(p_pe .eq. p_io) then
      			WRITE(IO_STDOUT,*) ' MEAN  2D ', kyears, kmonts, kdays, ITER
      			WRITE(IO_STDOUT,*) ' GLUE NITER  ', VAR
      		ENDIF
          DEALLOCATE(FIELD21)

      		END SUBROUTINE MEAN2D_GLOBAL


      SUBROUTINE MEAN3D_GLOBAL(KDAYS,KMONTS,KYEARS, FIELD3D, I4CODE,&
               & ivec)
      USE MO_PARAM1
      INCLUDE 'mpif.h'
      		INTEGER,INTENT(IN)::KDAYS,KMONTS,KYEARS,I4CODE,ivec
      		REAL,INTENT(IN)::FIELD3D(:,:,:)

          REAL(sp),allocatable::FIELD1(:,:,:)
      		INTEGER I_start_G,J_start_G,ICODE,ITER
      		INTEGER*8::global_nx, global_ny, global_nz
      		CHARACTER(LEN=3) :: VAR

      		ICODE=I4CODE
      		ITER = NITER

      		I_start_G = 1
      		J_start_G = 1

      		if (ivec .eq. 1) then
      			I_start_G=0
      			allocate(FIELD1(I_start:IE,1:JE,KE))
      			FIELD1(:,:,:) = FIELD3D(:,:,:)
      			dummyE_3d_U(:,:,:) = FIELD1(ius:iue,jus:jue,1:ke)
      			elseif (ivec .eq. 2) then
      				J_start_G=0
      				allocate(FIELD1(1:IE,J_start:JE,KE))
      			  FIELD1(:,:,:) = FIELD3D(:,:,:)
      				dummyE_3d_V(:,:,:) = FIELD1(ivs:ive,jvs:jve,1:ke)
        			elseif (ivec .eq. -1) then
       				allocate(FIELD1(1:IE,1:JE,KEP))
      			   FIELD1(:,:,:) = FIELD3D(:,:,:)
      				dummyE_3d(:,:,:) = FIELD1(iss:iee,jss:jee,1:ke)
      				elseif (ivec .eq. 0) then
      					 allocate(FIELD1(1:IE,1:JE,KE))
      			     FIELD1(:,:,:) = FIELD3D(:,:,:)
      				  dummyE_3d(:,:,:) = FIELD1(iss:iee,jss:jee,1:ke)
      				else
      					CALL STOP_ALL('mo_parallel, wrong ivec number')
      				endif

      				global_nx =  INT(IE_G+1-I_start_G)
      				global_ny = INT(JE_G+1-J_start_G)
      				global_nz = INT(KE)

              counts_3d_u(3) = KE ! Re-assign due to  corrections of 2D variables
              counts_3d_v(3) = KE
              WRITE(nter_str,'("_",I8.8 )') ITER
                                if (ICODE .EQ.  2 ) then
      					VAR = 'THO'
      					err = nfmpi_create(MPI_COMM_WORLD,'./tho/THO'//  &
      					&  trim(adjustl(nter_str))//'.nc',cmode,MPI_INFO_NULL,ncid)
      				endif

      				if (ICODE .EQ.  5 ) then
      					VAR = 'SAO'
      					err = nfmpi_create(MPI_COMM_WORLD,'./sao/SAO'//  &
      					&  trim(adjustl(nter_str))//'.nc',cmode,MPI_INFO_NULL,ncid)
      				endif

      				if (ICODE .EQ.  6 ) then
      					VAR = 'POO'
      					err = nfmpi_create(MPI_COMM_WORLD,'./poo/POO'//  &
      					&  trim(adjustl(nter_str))//'.nc',cmode,MPI_INFO_NULL,ncid)
      				endif

      				if (ICODE .EQ.  8 ) then
      				    VAR = 'ROO'
      			     err = nfmpi_create(MPI_COMM_WORLD,'./roo/ROO'//  &
      		      &  trim(adjustl(nter_str))//'.nc',cmode,MPI_INFO_NULL,ncid)
      				endif

      				if (ICODE .EQ.  3 ) then
      					VAR = 'UOO'
      					err = nfmpi_create(MPI_COMM_WORLD,'./uoo/UOO'//  &
      					&  trim(adjustl(nter_str))//'.nc',cmode,MPI_INFO_NULL,ncid)
      				endif

      				if (ICODE .EQ.  4 ) then
      					VAR = 'VOE'
      					err = nfmpi_create(MPI_COMM_WORLD,'./voe/VOE'//  &
      					&  trim(adjustl(nter_str))//'.nc',cmode,MPI_INFO_NULL,ncid)
      				endif

      				if (ICODE .EQ.  7 ) then
      					VAR = 'WOO'
      					err = nfmpi_create(MPI_COMM_WORLD,'./woo/WOO'//  &
      					&  trim(adjustl(nter_str))//'.nc',cmode,MPI_INFO_NULL,ncid)
      				endif

#ifdef NON_HYDROSTATIC
      				if (ICODE .EQ.  16 ) then
      					VAR = 'PNH'
      					err = nfmpi_create(MPI_COMM_WORLD,'./pnh/PNH'//  &
      					&  trim(adjustl(nter_str))//'.nc',cmode,MPI_INFO_NULL,ncid)
      				endif
#endif

#ifdef DIFFDIAG
      				if (ICODE .EQ.  110 ) then
      					VAR = 'AVO'
      					err = nfmpi_create(MPI_COMM_WORLD,'./avo/AVO'//  &
      					&  trim(adjustl(nter_str))//'.nc',cmode,MPI_INFO_NULL,ncid)
      				endif
      				if (ICODE .EQ.  111 ) then
      					VAR = 'DVO'
      					err = nfmpi_create(MPI_COMM_WORLD,'./dvo/DVO'//  &
      					&  trim(adjustl(nter_str))//'.nc',cmode,MPI_INFO_NULL,ncid)
      				endif
#endif

      				CALL check(err, 'In nfmpi_create: ')

      				err = nfmpi_def_dim(ncid, "lon", global_nx, dimidxyz(1))
      				call check(err, 'In nfmpi_def_dim x: ')
      				err = nfmpi_def_dim(ncid, "lat", global_ny, dimidxyz(2))
      				call check(err, 'In nfmpi_def_dim y: ')
      				err = nfmpi_def_dim(ncid, "z", global_nz, dimidxyz(3))
      				call check(err, 'In nfmpi_def_dim z: ')

      				err = nfmpi_def_var(ncid, VAR, NF_DOUBLE,  &
      				& 3, dimidxyz, varid)
      				call check(err, 'In nfmpi_def_var: ')

      				err = nfmpi_enddef(ncid)
      				call check(err, 'In nfmpi_enddef: ')

      		if (ICODE .EQ.  2 .OR. ICODE .EQ.  5 .OR. ICODE .EQ.  6 .OR. &
      		& ICODE .EQ.  7 .OR. ICODE .EQ.  16 .OR. ICODE .EQ.  110  .OR. &
      			& ICODE .EQ.  111 .OR.  ICODE .EQ.  8 ) then
      					err = nfmpi_put_vara_real_all(ncid, varid, starts_3d_p,   &
      					&  counts_3d_p, dummyE_3d, nTypes, subarray_p_3d)
      				endif

      				if (ICODE .EQ.  3 ) then
      					err = nfmpi_put_vara_real_all(ncid, varid, starts_3d_u,   &
      					&  counts_3d_u, dummyE_3d_U, ntypes, subarray_u)
      				endif

      				if (ICODE .EQ.  4 ) then
      					err = nfmpi_put_vara_real_all(ncid, varid, starts_3d_v,   &
      					&  counts_3d_v, dummyE_3d_V, ntypes, subarray_v)
      				endif

      				call check(err, 'In nfmpi_put_var_all: ')

      				err = nfmpi_close(ncid)
      				call check(err, 'In nfmpi_close: ')

      				IF(p_pe .eq. p_io) then
      					WRITE(IO_STDOUT,*) ' MEAN  3D ', kyears, kmonts, kdays, ITER
      					WRITE(IO_STDOUT,*) ' GLUE NITER  ', VAR
      				ENDIF
              DEALLOCATE(FIELD1)

      				END SUBROUTINE MEAN3D_GLOBAL

              SUBROUTINE INSITU_THO(TH,SH,PA)
        !     TRANSFORMATION FROM POTENTIAL TO IN SITU DENSITY
        !     from adisitj.f90 in MPIOM
        !     changed by peterspy, used in post process
        		  REAL(kind=SELECTED_REAL_KIND(6,37)):: A1,A2,A3,A4,B1,B2,C1,C2,C3,D,E1,E2
        !      PARAMETER(A1=3.6504E-4,A2=8.3198E-5,A3=5.4065E-7,A4=4.0274E-9,    &
        !     & B1=1.7439E-5,B2=2.9778E-7,C1=8.9309E-7,C2=3.1628E-8,             &
        !     & C3=2.1987E-10,D=4.1057E-9,E1=1.6056E-10,E2=5.0484E-12)
              REAL(kind=SELECTED_REAL_KIND(6,37))::TH(ie,je),SH(ie,je),PA(ie,je),     &
             &  PR(ie,je),TPO(ie,je),T(ie,je),QC(ie,je),QV(ie,je),DC(ie,je),DV(ie,je),&
             &  QNQ(ie,je),QN3(ie,je),QVS(ie,je),DVS(ie,je),FNE(ie,je),FST(ie,je)
        
             DATA A1,A2,A3,A4/3.6504E-4,8.3198E-5,5.4065E-7,4.0274E-9/
             DATA B1,B2/1.7439E-5,2.9778E-7/
             DATA C1,C2,C3/8.9309E-7,3.1628E-8,2.1987E-10/
             DATA D/4.1057E-9/
             DATA E1,E2/1.6056E-10,5.0484E-12/
             
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
              END SUBROUTINE INSITU_THO
        
              SUBROUTINE INSITU_ROO(T,S,P,RH)
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
              END SUBROUTINE INSITU_ROO
#endif /*GLOBAL_NETCDF*/

      				END MODULE MO_GLUE
