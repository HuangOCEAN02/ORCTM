MODULE mo_mpi

  ! Comment: Please use basic WRITE to nerr for messaging in the whole
  !          MPI package to achieve proper output. 

  USE mo_kind

  IMPLICIT NONE

  PRIVATE                          ! all declarations are private

  PUBLIC :: nerr
  INTEGER, PARAMETER :: nerr = 0

  ! subroutines defined, overloaded depending on argument type 

  PUBLIC :: p_start, p_stop, p_abort
  PUBLIC :: p_send, p_recv, p_sendrecv, p_bcast, p_barrier 
  PUBLIC :: p_isend, p_irecv, p_wait
#ifdef __ONESIDED
  PUBLIC :: p_win_create, p_win_fence, p_win_free, p_put
#endif
  PUBLIC :: p_gather, p_max, p_min, p_sum
  PUBLIC :: p_set_communicator
  PUBLIC :: p_probe

#ifndef NOMPI
  PUBLIC :: MPI_INTEGER, MPI_STATUS_SIZE, MPI_SUCCESS
  PUBLIC :: MPI_SUM, MPI_MIN, MPI_MAX, MPI_Wtime
#endif

  ! real data type matching real type of MPI implementation

  PUBLIC :: p_real_dp

  ! logical switches

  PUBLIC :: p_parallel, p_parallel_io

  ! PE identifier

  PUBLIC :: p_pe, p_io, p_nprocs

  ! communicator

  PUBLIC :: p_all_comm
  PUBLIC :: p_communicator_a, p_communicator_b, p_communicator_d

  ! old fashioned method (MPI-1)

#ifndef NOMPI
  INCLUDE 'mpif.h'
#endif

  ! general run time information

#ifndef NOMPI
  INTEGER :: version, subversion   ! MPI version
#endif

  ! MPI call inherent variables

#ifndef NOMPI
  INTEGER :: p_error                     ! MPI error number

  INTEGER :: p_status(MPI_STATUS_SIZE)   ! standard information of MPI_RECV
#endif

  ! public parallel run information

  LOGICAL :: p_parallel, p_parallel_io

  INTEGER :: p_pe                  ! this is the PE number of this task
  INTEGER :: p_io                  ! PE number of PE handling IO
  INTEGER :: p_nprocs              ! number of available PEs (processors)

  ! communicator sets 

  INTEGER :: p_all_comm       ! replaces MPI_COMM_WORLD in one application
  INTEGER :: p_communicator_a ! for Set A       
  INTEGER :: p_communicator_b ! for Set B       
  INTEGER :: p_communicator_d ! for debug node
            
  ! non blocking calls

#ifndef NOMPI
  INTEGER, ALLOCATABLE :: p_request(:) ! request values for non blocking calls
  INTEGER :: p_irequest ! the first p_irequest values of p_request are in use
  INTEGER :: p_mrequest ! actual size of p_request
  INTEGER, PARAMETER :: p_request_alloc_size = 1024
#endif

  ! module intrinsic names

  INTEGER :: mype                  ! this is the PE number of this task
#ifndef NOMPI
  INTEGER :: iope                  ! PE able to do IO
#endif
  INTEGER :: npes                  ! number of available PEs

  INTEGER :: nbcast                ! counter for broadcasts for debugging 

  ! MPI transfer types

  INTEGER :: p_real_dp
#ifndef NOMPI
  INTEGER :: p_real_sp
  INTEGER :: p_int_i4
  INTEGER :: p_int_i8

  ! native types

  INTEGER :: p_int       ! maybe switched by compiler options therefor reset 
  INTEGER :: p_real      ! maybe switched by compiler options therefor reset 

  PUBLIC :: p_int, p_real

  INTEGER :: p_bool
  INTEGER :: p_char

  INTEGER :: p_ig, p_rg
  INTEGER :: p_i4, p_i8
  INTEGER :: p_sp, p_dp
#endif

  ! for checking out integer and real variables separat. KIND values usually
  ! overlap and give the byte size or are defined as sequence separate for
  ! both groups.

  INTEGER, PARAMETER :: real_type    = 1
  INTEGER, PARAMETER :: integer_type = 2

  ! define generic interfaces to allow proper compiling
  ! with picky compilers like NAG f95 for clean argument checking and 
  ! shortening the call sequence.

  INTERFACE p_send
     MODULE PROCEDURE p_send_real
     MODULE PROCEDURE p_send_int
     MODULE PROCEDURE p_send_bool
     MODULE PROCEDURE p_send_real_1d
     MODULE PROCEDURE p_send_int_1d
     MODULE PROCEDURE p_send_bool_1d
     MODULE PROCEDURE p_send_real_2d
     MODULE PROCEDURE p_send_int_2d
     MODULE PROCEDURE p_send_bool_2d
     MODULE PROCEDURE p_send_real_3d
     MODULE PROCEDURE p_send_int_3d
     MODULE PROCEDURE p_send_bool_3d
     MODULE PROCEDURE p_send_real_4d
     MODULE PROCEDURE p_send_int_4d
     MODULE PROCEDURE p_send_bool_4d
     MODULE PROCEDURE p_send_char
     MODULE PROCEDURE p_send_real_5d
  END INTERFACE

  INTERFACE p_isend
     MODULE PROCEDURE p_isend_real
     MODULE PROCEDURE p_isend_int
     MODULE PROCEDURE p_isend_bool
     MODULE PROCEDURE p_isend_real_1d
     MODULE PROCEDURE p_isend_int_1d
     MODULE PROCEDURE p_isend_bool_1d
     MODULE PROCEDURE p_isend_real_2d
     MODULE PROCEDURE p_isend_int_2d
     MODULE PROCEDURE p_isend_bool_2d
     MODULE PROCEDURE p_isend_real_3d
     MODULE PROCEDURE p_isend_int_3d
     MODULE PROCEDURE p_isend_bool_3d
     MODULE PROCEDURE p_isend_real_4d
     MODULE PROCEDURE p_isend_int_4d
     MODULE PROCEDURE p_isend_bool_4d
     MODULE PROCEDURE p_isend_char
     MODULE PROCEDURE p_isend_real_5d
  END INTERFACE

#ifdef __ONESIDED
  INTERFACE p_win_create
     MODULE PROCEDURE p_win_create_real_2d
     MODULE PROCEDURE p_win_create_real_3d
  END INTERFACE

  INTERFACE p_win_fence
     MODULE PROCEDURE p_win_fence_int
     MODULE PROCEDURE p_win_fence_int_1d
  END INTERFACE

  INTERFACE p_win_free
     MODULE PROCEDURE p_win_free_int
     MODULE PROCEDURE p_win_free_int_1d
  END INTERFACE

  INTERFACE p_put
     MODULE PROCEDURE p_put_real_2d
     MODULE PROCEDURE p_put_real_3d
  END INTERFACE
#endif
  INTERFACE p_recv
     MODULE PROCEDURE p_recv_real
     MODULE PROCEDURE p_recv_int
     MODULE PROCEDURE p_recv_bool
     MODULE PROCEDURE p_recv_real_1d
     MODULE PROCEDURE p_recv_int_1d
     MODULE PROCEDURE p_recv_bool_1d
     MODULE PROCEDURE p_recv_real_2d
     MODULE PROCEDURE p_recv_int_2d
     MODULE PROCEDURE p_recv_bool_2d
     MODULE PROCEDURE p_recv_real_3d
     MODULE PROCEDURE p_recv_int_3d
     MODULE PROCEDURE p_recv_bool_3d
     MODULE PROCEDURE p_recv_real_4d
     MODULE PROCEDURE p_recv_int_4d
     MODULE PROCEDURE p_recv_bool_4d
     MODULE PROCEDURE p_recv_char
     MODULE PROCEDURE p_recv_real_5d
  END INTERFACE

  INTERFACE p_irecv
     MODULE PROCEDURE p_irecv_real
     MODULE PROCEDURE p_irecv_real_1d
     MODULE PROCEDURE p_irecv_real_2d
     MODULE PROCEDURE p_irecv_real_3d
     MODULE PROCEDURE p_irecv_real_4d
  END INTERFACE

  INTERFACE p_sendrecv
     MODULE PROCEDURE p_sendrecv_real_1d
     MODULE PROCEDURE p_sendrecv_real_2d
     MODULE PROCEDURE p_sendrecv_real_3d
     MODULE PROCEDURE p_sendrecv_real_4d
  END INTERFACE

  INTERFACE p_bcast
     MODULE PROCEDURE p_bcast_real
     MODULE PROCEDURE p_bcast_int_i4
     MODULE PROCEDURE p_bcast_int_i8
     MODULE PROCEDURE p_bcast_bool
     MODULE PROCEDURE p_bcast_real_1d
     MODULE PROCEDURE p_bcast_int_1d
     MODULE PROCEDURE p_bcast_bool_1d
     MODULE PROCEDURE p_bcast_real_2d
     MODULE PROCEDURE p_bcast_int_2d
     MODULE PROCEDURE p_bcast_bool_2d
     MODULE PROCEDURE p_bcast_real_3d
     MODULE PROCEDURE p_bcast_int_3d
     MODULE PROCEDURE p_bcast_bool_3d
     MODULE PROCEDURE p_bcast_real_4d
     MODULE PROCEDURE p_bcast_int_4d
     MODULE PROCEDURE p_bcast_bool_4d
     MODULE PROCEDURE p_bcast_char
  END INTERFACE

  INTERFACE p_probe
     MODULE PROCEDURE p_probe_real
     MODULE PROCEDURE p_probe_int
     MODULE PROCEDURE p_probe_bool
     MODULE PROCEDURE p_probe_char
  END INTERFACE

  INTERFACE p_gather
     MODULE PROCEDURE p_gather_real_1d2d
  END INTERFACE

  INTERFACE p_max
     MODULE PROCEDURE p_max_0d
     MODULE PROCEDURE p_max_1d
     MODULE PROCEDURE p_max_2d
     MODULE PROCEDURE p_max_3d
  END INTERFACE

  INTERFACE p_min
     MODULE PROCEDURE p_min_0d
     MODULE PROCEDURE p_min_1d
     MODULE PROCEDURE p_min_2d
     MODULE PROCEDURE p_min_3d
  END INTERFACE

  INTERFACE p_sum
     MODULE PROCEDURE p_sum_0d
     MODULE PROCEDURE p_sum_1d
  END INTERFACE

CONTAINS

  SUBROUTINE p_start(model_name)

#ifndef NOMPI
#if defined (__prism) && defined (use_comm_MPI1)
    USE mod_prism_proto, ONLY: prism_ok
#endif
#endif

    CHARACTER(len=*), INTENT(in), OPTIONAL :: model_name

    ! variables are required for determing I/O size in bytes of the defined
    ! KIND types for assigning the right MPI data types with the used kinds

    INTEGER :: io_size, integer_io_size, integer_byte_size

    INTEGER      :: iig = 0  
    INTEGER (i4) :: ii4 = 0_i4
    INTEGER (i8) :: ii8 = 0_i8

    REAL         :: rrg = 0.0
    REAL (sp)    :: rsp = 0.0_sp
    REAL (dp)    :: rdp = 0.0_dp

    ! temporary array to distibute the determined MPI types

    INTEGER :: p_send(7)

    ! variables used for determing the I/O PE

    LOGICAL :: liope
    INTEGER, ALLOCATABLE :: iope_table(:)
    CHARACTER(len=132) :: io_pe_message 

    CHARACTER(len=132) :: yname

#ifndef NOMPI
#if defined (__prism) && defined (use_comm_MPI1)
    INTEGER :: prism_model_number
    CHARACTER(len=132) :: prism_model_name

    EXTERNAL :: prism_abort_proto
    EXTERNAL :: prism_init_comp_proto
    EXTERNAL :: prism_get_localcomm_proto
#endif
#endif

    ! Executable statements:

    IF (PRESENT(model_name)) THEN
      yname = TRIM(model_name)
    ELSE
      yname = 'mpi-om'
    END IF

    nbcast = 0                 

    io_pe_message(:) = ' '

    ! start MPI

#ifndef NOMPI

    WRITE (nerr,*) 'orctm: call mpi_init ...'
    CALL flush (nerr)

    CALL MPI_INIT (p_error)

    WRITE (nerr,*) 'orctm: call mpi_init finished ...'
    CALL flush (nerr)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_INIT failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF
#endif

    ! create communicator for this process alone before 
    ! potentially joining MPI2

#ifndef NOMPI
#if defined (__prism) && defined (use_comm_MPI1)

    prism_model_name = TRIM(yname)

    WRITE (nerr,*) 'orctm: call prism_init_comp_proto ...'
    CALL flush (nerr)

    CALL prism_init_comp_proto (prism_model_number, TRIM(prism_model_name), &
         p_error)

    WRITE (nerr,*) 'orctm: call prism_init_comp_proto finished ...'
    CALL flush (nerr)

    IF (p_error /= prism_ok) THEN
      WRITE (nerr,*) ' prism_init_comp_proto failed'
      CALL prism_abort_proto(prism_model_number, TRIM(yname),'abort1')
    ENDIF

    CALL prism_get_localcomm_proto(p_all_comm, p_error)

    IF (p_error /= prism_ok) THEN
      WRITE (nerr,*) ' prism_get_localcomm_proto failed'
      CALL prism_abort_proto(prism_model_number, TRIM(yname),'abort2')
    ENDIF

#else

    CALL MPI_COMM_DUP (MPI_COMM_WORLD, p_all_comm, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_DUP failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF
#endif
#endif

    ! get local PE identification

#ifndef NOMPI
    CALL MPI_COMM_RANK (p_all_comm, mype, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_RANK failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ELSE
#ifdef DEBUG       
       WRITE (nerr,'(a,i4,a)') ' PE ', mype, ' started.'
#endif
    END IF
#else
    mype = 0
#endif

    ! get number of available PEs

#ifndef NOMPI
    CALL MPI_COMM_SIZE (p_all_comm, npes, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF
#else
    npes = 1
#endif

    ! for non blocking calls

#ifndef NOMPI
    p_mrequest = p_request_alloc_size
    ALLOCATE(p_request(p_mrequest))
    p_irequest = 0
#endif

    ! look for a dedicated IO PE

#ifndef NOMPI
    CALL MPI_ATTR_GET (p_all_comm, MPI_IO, iope, liope, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_ATTR_GET failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    IF (iope == MPI_ANY_SOURCE) THEN

       ! all nodes can do IO

       IF (mype == 0) THEN
          WRITE (io_pe_message,'(a)') &
               '  All nodes can do I/O, select PE 0 for I/O handling.'
       END IF
       p_io = 0

    ELSE

       ALLOCATE (iope_table(0:npes-1))

       IF (liope) THEN
          iope_table(mype) = iope
          CALL MPI_GATHER (iope_table(mype), 1, MPI_INTEGER, &
               iope_table, 1, MPI_INTEGER,       &
               0, p_all_comm, p_error)
          
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GATHER failed.'
             WRITE (nerr,'(a,i4)') ' Error =  ', p_error
             CALL p_abort
          END IF

          IF (mype == 0) THEN
             ! Now select the first given iope from table as IO PE.
             WRITE (io_pe_message,'(a,i3,a)') &
                  '  Select PE ', iope_table(0), ' for I/O handling.'
             p_io = iope_table(0)
          END IF
          CALL MPI_BCAST (p_io, 1, MPI_INTEGER, 0, p_all_comm, p_error)

       ELSE
          ! if no dedicated IO PE is given, use PE 0
          p_io = 0

          WRITE (io_pe_message,'(a)') &
               '  No dedicated I/O PE, select PE 0 for I/O handling.'
       END IF

       DEALLOCATE (iope_table)

    END IF
#else    
    p_io = 0
#endif

    ! Information ...

    IF (mype == 0) THEN    
       WRITE (nerr,'(/,a,a,a)') ' ', &
            TRIM(yname), ' MPI interface runtime information:'
    END IF

     IF (npes < 2) THEN
       p_parallel = .FALSE.
       p_parallel_io = .TRUE.   ! can always do I/O
       IF (mype == 0) THEN
          WRITE (nerr,'(a)') '  Single processor run.'
       END IF
       p_pe = 0
       p_nprocs = 1
    ELSE 
       p_parallel = .TRUE.
       IF (mype == p_io) THEN
          p_parallel_io = .TRUE.
       ELSE
          p_parallel_io = .FALSE.
       END IF
       IF (mype == 0) THEN       
          WRITE (nerr,'(a,i4,a)') '  Run on ', npes, ' processors.'
       END IF
       p_pe = mype
       p_nprocs = npes
    END IF

    ! inform on I/O PE situation

    IF (mype == 0) THEN  
       WRITE (nerr, '(a)') io_pe_message(1:LEN_TRIM(io_pe_message))
    END IF

#ifndef NOMPI

    IF (p_parallel) THEN

       ! lets check the available MPI version

       CALL MPI_GET_VERSION (version, subversion, p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a)') ' MPI_GET_VERSION failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       IF (mype == 0) THEN 
          WRITE (nerr,'(a,i1,a1,i1)') &
               '  Used MPI version: ', version, '.', subversion
       END IF

       ! due to a possible circular dependency with mo_machine and other 
       ! modules, we determine here locally the I/O size of the different  
       ! kind types (assume 8 bit/byte. This is than used for determing 
       ! the right MPI send/receive type parameters. 

       ! first get the native INTEGER size

       integer_byte_size = BIT_SIZE(iig)/8

       ! and inquire for the I/O size (is independent of byte or word 
       ! values ...)

       INQUIRE (iolength=io_size) iig
       integer_io_size = io_size
       p_ig = io_size/integer_io_size*integer_byte_size

       ! and the native REAL size 

       INQUIRE (iolength=io_size) rrg
       p_rg = io_size/integer_io_size*integer_byte_size

       ! find now the size of usual 4 byte and 8 byte INTEGER 
       ! (might be 8 byte both, or only the 4 byte available ...

       INQUIRE (iolength=io_size) ii4
       p_i4 = io_size/integer_io_size*integer_byte_size
       INQUIRE (iolength=io_size) ii8
       p_i8 = io_size/integer_io_size*integer_byte_size

       ! find now the size of usual 4 byte and 8 byte REAL 
       ! (might be 8 byte both)

       INQUIRE (iolength=io_size) rsp
       p_sp = io_size/integer_io_size*integer_byte_size
       INQUIRE (iolength=io_size) rdp
       p_dp = io_size/integer_io_size*integer_byte_size

       ! testing this variables

       p_int_i4  = p_type (i4, integer_type) 
       p_int_i8  = p_type (i8, integer_type) 
       p_real_sp = p_type (sp, real_type) 
       p_real_dp = p_type (dp, real_type) 

       IF (mype == 0) THEN

          IF (p_ig == p_i4) THEN
             p_int = p_int_i4
          ELSE IF (p_ig == p_i8) THEN
             p_int = p_int_i8
          END IF

          IF (p_rg == p_sp) THEN
             p_real = p_real_sp
          ELSE IF (p_rg == p_dp) THEN
             p_real = p_real_dp
          END IF

       END IF

       p_send(1) = p_int
       p_send(2) = p_int_i4
       p_send(3) = p_int_i8
       p_send(4) = p_real
       p_send(5) = p_real_sp
       p_send(6) = p_real_dp

       CALL MPI_BCAST (p_send, 6, MPI_INTEGER, 0, p_all_comm, p_error)

#ifdef DEBUG
       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a)') ' MPI_BCAST for send/receive types failed.'
          WRITE (nerr,'(a,i4)') ' Error = ', p_error
       END IF
#endif

       p_int      = p_send(1) 
       p_int_i4   = p_send(2) 
       p_int_i8   = p_send(3) 
       p_real     = p_send(4) 
       p_real_sp  = p_send(5) 
       p_real_dp  = p_send(6) 

       ! set logical and character types to native types

       p_bool = MPI_LOGICAL
       p_char = MPI_CHARACTER

       IF (mype == 0) THEN
          WRITE (nerr,'(/)')
          IF (p_real_sp == MPI_REAL) THEN
             WRITE (nerr,'(a)') ' Selected type: MPI_REAL for KIND sp'
          ELSE IF (p_real_sp == MPI_DOUBLE_PRECISION) THEN
             WRITE (nerr,'(a)') &
                  ' Selected type: MPI_DOUBLE_PRECISION for KIND sp'
          END IF

          IF (p_real_dp == MPI_DOUBLE_PRECISION) THEN
             WRITE (nerr,'(a)') &
                  ' Selected type: MPI_DOUBLE_PRECISION for KIND dp'
          END IF

          IF (p_int_i4 == MPI_INTEGER4) THEN
             WRITE (nerr,'(a)') ' Selected type: MPI_INTEGER4 for KIND i4'
          ELSE IF (p_int_i4 == MPI_INTEGER8) THEN
             WRITE (nerr,'(a)') ' Selected type: MPI_INTEGER8 for KIND i4'
          END IF

          IF (p_int_i8 == MPI_INTEGER8) THEN
             WRITE (nerr,'(a)') ' Selected type: MPI_INTEGER8 for KIND i8'
          END IF
       END IF
    END IF
    WRITE (nerr,'(/)')
#endif



#ifdef DEBUG
    WRITE (nerr,'(a)')    ' Transfer types:'
    WRITE (nerr,'(a,i4)') '  INTEGER generic:', p_int
    WRITE (nerr,'(a,i4)') '  INTEGER 4 byte :', p_int_i4
    WRITE (nerr,'(a,i4)') '  INTEGER 8 byte :', p_int_i8
    WRITE (nerr,'(a,i4)') '  REAL generic   :', p_real
    WRITE (nerr,'(a,i4)') '  REAL single    :', p_real_sp
    WRITE (nerr,'(a,i4)') '  REAL double    :', p_real_dp
#endif

  END SUBROUTINE p_start

  SUBROUTINE p_stop

    ! finish MPI and clean up all PEs

#ifndef NOMPI
    ! to prevent abort due to unfinished communication
    CALL p_barrier(p_all_comm) 

    CALL MPI_FINALIZE (p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_FINALIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
    p_parallel = .FALSE.
    DEALLOCATE(p_request)
#endif

  STOP 'in p_stop'

  END SUBROUTINE p_stop

  SUBROUTINE p_abort

    ! this routine should be used instead of abort, util_abort() or STOP 
    ! in all routines for proper clean up of all PEs

#ifndef NOMPI
    CALL MPI_ABORT (MPI_COMM_WORLD, 0, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_ABORT failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF
#else
    ! p_abort should also abort if we are not running with MPI

    STOP 'p_abort called'
#endif

  END SUBROUTINE p_abort

  FUNCTION p_type (kind_type, var_type) RESULT (p_message_type)

    USE mo_kind 

    INTEGER              :: p_message_type
    INTEGER, INTENT(in)  :: kind_type, var_type

    IF (var_type == integer_type) THEN
       IF (kind_type == i8) THEN 
          p_message_type = check_type_i8 ()
       ELSE IF (kind_type == i4) THEN    
          p_message_type = check_type_i4 ()
       END IF
    ELSE IF (var_type == real_type) THEN
       IF (kind_type == dp) THEN    
          p_message_type = check_type_dp ()
       ELSE IF (kind_type == sp) THEN    
          p_message_type = check_type_sp ()
       END IF
    END IF

  END FUNCTION p_type

  FUNCTION check_type_i4 () RESULT (p_type)

    INTEGER :: p_type

#ifndef NOMPI
    INTEGER :: datatype

    INTEGER (i4) :: buf_int_i4(1,8)
    INTEGER (i4) :: a, b, c

    p_type = MPI_DATATYPE_NULL

    IF (mype == 0) THEN        

       a = HUGE(a)
       buf_int_i4(1,1) = a

       CALL MPI_SEND (buf_int_i4(1,1), 1, MPI_INTEGER4, 1, 1, &
            p_all_comm, p_error)
       CALL MPI_RECV (buf_int_i4(1,5), 1, MPI_INTEGER4, 1, 2, &
            p_all_comm, p_status, p_error) 
       b = buf_int_i4(1,5)

       IF (a == b) THEN
          CALL MPI_SEND (MPI_BYTE, 1, MPI_INTEGER, 1, 0, &
               p_all_comm, p_error) 
          CALL MPI_SEND (buf_int_i4(1,1), p_i4, MPI_BYTE, 1, 3, &
               p_all_comm, p_error) 
          CALL MPI_RECV (buf_int_i4(1,5), p_i4, MPI_BYTE, 1, 4, &
               p_all_comm, p_status, p_error) 
          c =  buf_int_i4(1,5)

          IF (a /= c) THEN
             WRITE (nerr,'(a)') &
                  ' Warning: MPI_INTEGER4 and MPI_BYTE not equivalent'
             WRITE (nerr,'(a)') ' Using MPI_INTEGER4 anyway'
          END IF

          p_type = MPI_INTEGER4

       ELSE
          CALL MPI_SEND (MPI_INTEGER8, 1, MPI_INTEGER, 1, 0, &
               p_all_comm, p_error)            
          CALL MPI_SEND (buf_int_i4(1,1), 1, MPI_INTEGER8, 1, 5, &
               p_all_comm, p_error) 
          CALL MPI_RECV (buf_int_i4(1,5), 1, MPI_INTEGER8, 1, 6, &
               p_all_comm, p_status, p_error) 
          b = buf_int_i4(1,5)

          IF (a == b) THEN
             CALL MPI_SEND (buf_int_i4(1,1), p_i4, MPI_BYTE, 1, 7, &
                  p_all_comm, p_error) 
             CALL MPI_RECV (buf_int_i4(1,5), p_i4, MPI_BYTE, 1, 8, &
                  p_all_comm, p_status, p_error) 
             c =  buf_int_i4(1,5)

             IF (a /= c) THEN
                WRITE (nerr,'(a)') &
                     ' Warning: MPI_INTEGER8 and MPI_BYTE not equivalent'
                WRITE (nerr,'(a)') &
                     ' Using MPI_INTEGER8 anyway'
             END IF

             p_type = MPI_INTEGER8

          END IF
       END IF
    END IF

    IF (mype == 1) THEN
       CALL MPI_RECV (buf_int_i4(1,7), 1, MPI_INTEGER4, 0, 1, &
            p_all_comm, p_status, p_error)
       buf_int_i4(1,3) = buf_int_i4(1,7)
       CALL MPI_SEND (buf_int_i4(1,3), 1, MPI_INTEGER4, 0, 2, &
            p_all_comm, p_error)

       CALL MPI_RECV (datatype, 1, MPI_INTEGER, 0, 0, &
            p_all_comm, p_status, p_error)

       IF (datatype == MPI_BYTE) THEN
          CALL MPI_RECV (buf_int_i4(1,7), p_i4, MPI_BYTE, 0, 3, &
               p_all_comm, p_status, p_error)
          buf_int_i4(1,3) = buf_int_i4(1,7)
          CALL MPI_SEND (buf_int_i4(1,3), p_i4, MPI_BYTE, 0, 4, &
               p_all_comm, p_error)
       ELSE IF (datatype == MPI_INTEGER8) THEN
          CALL MPI_RECV (buf_int_i4(1,7), 1, MPI_INTEGER8, 0, 5, &
               p_all_comm, p_status, p_error)
          buf_int_i4(1,3) = buf_int_i4(1,7)
          CALL MPI_SEND (buf_int_i4(1,3), 1, MPI_INTEGER8, 0, 6, &
               p_all_comm, p_error)
       END IF

       IF (datatype == MPI_INTEGER8) THEN
          CALL MPI_RECV (buf_int_i4(1,7), p_i4, MPI_BYTE, 0, 7, &
               p_all_comm, p_status, p_error)
          buf_int_i4(1,3) = buf_int_i4(1,7)
          CALL MPI_SEND (buf_int_i4(1,3), p_i4, MPI_BYTE, 0, 8, &
               p_all_comm, p_error)
       END IF
    END IF
#else
    p_type = 0
#endif

  END FUNCTION check_type_i4

  FUNCTION check_type_i8 () RESULT (p_type)

    INTEGER :: p_type

#ifndef NOMPI
    INTEGER (i8) :: buf_int_i8(1,8)
    INTEGER (i8) :: a, b, c

    p_type = MPI_DATATYPE_NULL

    IF (mype == 0) THEN        

       a = HUGE(a)
       buf_int_i8(1,1) = a

       CALL MPI_SEND (buf_int_i8(1,1), 1, MPI_INTEGER8, 1, 1, &
            p_all_comm, p_error)
       CALL MPI_RECV (buf_int_i8(1,5), 1, MPI_INTEGER8, 1, 2, &
            p_all_comm, p_status, p_error) 
       b = buf_int_i8(1,5)

       IF (a == b) THEN
          CALL MPI_SEND (buf_int_i8(1,1), p_i8, MPI_BYTE, 1, 3, &
               p_all_comm, p_error) 
          CALL MPI_RECV (buf_int_i8(1,5), p_i8, MPI_BYTE, 1, 4, &
               p_all_comm, p_status, p_error) 
          c =  buf_int_i8(1,5)

          IF (a /= c) THEN
             WRITE (nerr,'(a)') &
                  ' Warning: MPI_INTEGER8 and MPI_BYTE not equivalent'
             WRITE (nerr,'(a)') ' Using MPI_INTEGER8 anyway'
          END IF

          p_type = MPI_INTEGER8

       ELSE
          WRITE (nerr,'(a)') ' MPI_INTEGER8 not available.'
       END IF
    END IF

    IF (mype == 1) THEN
       CALL MPI_RECV (buf_int_i8(1,7), 1, MPI_INTEGER8, 0, 1, &
            p_all_comm, p_status, p_error)
       buf_int_i8(1,3) = buf_int_i8(1,7)
       CALL MPI_SEND (buf_int_i8(1,3), 1, MPI_INTEGER8, 0, 2, &
            p_all_comm, p_error)

       CALL MPI_RECV (buf_int_i8(1,7), p_i8, MPI_BYTE, 0, 3, &
            p_all_comm, p_status, p_error)
       buf_int_i8(1,3) = buf_int_i8(1,7)
       CALL MPI_SEND (buf_int_i8(1,3), p_i8, MPI_BYTE, 0, 4, &
            p_all_comm, p_error)
    END IF
#else
    p_type = 0
#endif

  END FUNCTION check_type_i8

  FUNCTION check_type_sp () RESULT (p_type)

    INTEGER :: p_type

#ifndef NOMPI
    INTEGER :: datatype

    REAL (sp) :: buf_real_sp(1,8)
    REAL (sp) :: a, b, c

    p_type = MPI_DATATYPE_NULL

    IF (mype == 0) THEN        

       a = HUGE(a)
       buf_real_sp(1,1) = a

       CALL MPI_SEND (buf_real_sp(1,1), 1, MPI_REAL, 1, 1, &
            p_all_comm, p_error)
       CALL MPI_RECV (buf_real_sp(1,5), 1, MPI_REAL, 1, 2, &
            p_all_comm, p_status, p_error) 
       b = buf_real_sp(1,5)

       IF (a == b) THEN
          CALL MPI_SEND (MPI_BYTE, 1, MPI_INTEGER, 1, 0, &
               p_all_comm, p_error) 
          CALL MPI_SEND (buf_real_sp(1,1), p_sp, MPI_BYTE, 1, 3, &
               p_all_comm, p_error) 
          CALL MPI_RECV (buf_real_sp(1,5), p_sp, MPI_BYTE, 1, 4, &
               p_all_comm, p_status, p_error) 
          c =  buf_real_sp(1,5)

          IF (a /= c) THEN
             WRITE (nerr,'(a)') ' Warning: MPI_REAL and MPI_BYTE not equivalent'
             WRITE (nerr,'(a)') ' Using MPI_REAL anyway'
          END IF

          p_type = MPI_REAL

       ELSE
          CALL MPI_SEND (MPI_DOUBLE_PRECISION, 1, MPI_INTEGER, 1, 0, &
               p_all_comm, p_error)            
          CALL MPI_SEND (buf_real_sp(1,1), 1, MPI_DOUBLE_PRECISION, 1, 5, &
               p_all_comm, p_error) 
          CALL MPI_RECV (buf_real_sp(1,5), 1, MPI_DOUBLE_PRECISION, 1, 6, &
               p_all_comm, p_status, p_error) 
          b = buf_real_sp(1,5)

          IF (a == b) THEN
             CALL MPI_SEND (buf_real_sp(1,1), p_sp, MPI_BYTE, 1, 7, &
                  p_all_comm, p_error) 
             CALL MPI_RECV (buf_real_sp(1,5), p_sp, MPI_BYTE, 1, 8, &
                  p_all_comm, p_status, p_error) 
             c =  buf_real_sp(1,5)

             IF (a /= c) THEN
                WRITE (nerr,'(a,a)') &
                     ' Warning: MPI_DOUBLE_PRECISION and MPI_BYTE ', &
                     'not equivalent'
                WRITE (nerr,'(a)') &
                     ' Using MPI_DOUBLE_PRECISION anyway'
             END IF

             p_type = MPI_DOUBLE_PRECISION

          END IF
       END IF
    END IF

    IF (mype == 1) THEN
       CALL MPI_RECV (buf_real_sp(1,7), 1, MPI_REAL, 0, 1, &
            p_all_comm, p_status, p_error)
       buf_real_sp(1,3) = buf_real_sp(1,7)
       CALL MPI_SEND (buf_real_sp(1,3), 1, MPI_REAL, 0, 2, &
            p_all_comm, p_error)

       CALL MPI_RECV (datatype, 1, MPI_INTEGER, 0, 0, &
            p_all_comm, p_status, p_error)

       IF (datatype == MPI_BYTE) THEN
          CALL MPI_RECV (buf_real_sp(1,7), p_sp, MPI_BYTE, 0, 3, &
               p_all_comm, p_status, p_error)
          buf_real_sp(1,3) = buf_real_sp(1,7)
          CALL MPI_SEND (buf_real_sp(1,3), p_sp, MPI_BYTE, 0, 4, &
               p_all_comm, p_error)
       ELSE IF (datatype == MPI_DOUBLE_PRECISION) THEN
          CALL MPI_RECV (buf_real_sp(1,7), 1, MPI_DOUBLE_PRECISION, 0, 5, &
               p_all_comm, p_status, p_error)
          buf_real_sp(1,3) = buf_real_sp(1,7)
          CALL MPI_SEND (buf_real_sp(1,3), 1, MPI_DOUBLE_PRECISION, 0, 6, &
               p_all_comm, p_error)
       END IF

       IF (datatype == MPI_DOUBLE_PRECISION) THEN
          CALL MPI_RECV (buf_real_sp(1,7), p_sp, MPI_BYTE, 0, 7, &
               p_all_comm, p_status, p_error)
          buf_real_sp(1,3) = buf_real_sp(1,7)
          CALL MPI_SEND (buf_real_sp(1,3), p_sp, MPI_BYTE, 0, 8, &
               p_all_comm, p_error)
       END IF
    END IF
#else
    p_type = 0
#endif

  END FUNCTION check_type_sp

  FUNCTION check_type_dp () RESULT (p_type)

    INTEGER :: p_type

#ifndef NOMPI
    REAL (dp) :: buf_real_dp(1,8)
    REAL (dp) :: a, b, c

    p_type = MPI_DATATYPE_NULL

    IF (mype == 0) THEN        

       a = HUGE(a)
       buf_real_dp(1,1) = a

       CALL MPI_SEND (buf_real_dp(1,1), 1, MPI_DOUBLE_PRECISION, 1, 1, &
            p_all_comm, p_error)
       CALL MPI_RECV (buf_real_dp(1,5), 1, MPI_DOUBLE_PRECISION, 1, 2, &
            p_all_comm, p_status, p_error) 
       b = buf_real_dp(1,5)

       IF (a == b) THEN
          CALL MPI_SEND (buf_real_dp(1,1), p_dp, MPI_BYTE, 1, 3, &
               p_all_comm, p_error) 
          CALL MPI_RECV (buf_real_dp(1,5), p_dp, MPI_BYTE, 1, 4, &
               p_all_comm, p_status, p_error) 
          c =  buf_real_dp(1,5)

          IF (a /= c) THEN
             WRITE (nerr,'(a)') &
                  ' Warning: MPI_DOUBLE_PRECISION and MPI_BYTE not equivalent'
             WRITE (nerr,'(a)') ' Using MPI_DOUBLE_PRECISION anyway'
          END IF

          p_type = MPI_DOUBLE_PRECISION

       ELSE
          WRITE (nerr,'(a)') ' MPI_DOUBLE_PRECISION not available.'
       END IF
    END IF

    IF (mype == 1) THEN
       CALL MPI_RECV (buf_real_dp(1,7), 1, MPI_DOUBLE_PRECISION, 0, 1, &
            p_all_comm, p_status, p_error)
       buf_real_dp(1,3) = buf_real_dp(1,7)
       CALL MPI_SEND (buf_real_dp(1,3), 1, MPI_DOUBLE_PRECISION, 0, 2, &
            p_all_comm, p_error)

       CALL MPI_RECV (buf_real_dp(1,7), p_dp, MPI_BYTE, 0, 3, &
            p_all_comm, p_status, p_error)
       buf_real_dp(1,3) = buf_real_dp(1,7)
       CALL MPI_SEND (buf_real_dp(1,3), p_dp, MPI_BYTE, 0, 4, &
            p_all_comm, p_error)
    END IF
#else
    p_type = 0
#endif

  END FUNCTION check_type_dp

  ! communicator set up

  SUBROUTINE p_set_communicator (nproca, nprocb, mapmesh, debug_parallel)

    INTEGER, INTENT(in) :: nproca, nprocb
    INTEGER, INTENT(in) :: mapmesh(0:,0:)
    INTEGER, INTENT(in) :: debug_parallel

#ifndef NOMPI
    INTEGER :: all_debug_pes(SIZE(mapmesh))

    INTEGER :: group_world, group_a, group_b, group_d
    INTEGER :: p_communicator_tmp
    
    INTEGER :: n, members

    ! first set global group

    CALL MPI_COMM_GROUP (p_all_comm, group_world, p_error)
       
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_GROUP failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF
       
    ! communicator is p_all_comm

    IF (debug_parallel >= 0 ) THEN

       CALL MPI_GROUP_INCL (group_world, 1, 0, group_d, p_error)
          
       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF
          
       CALL MPI_COMM_CREATE (p_all_comm, group_d, p_communicator_tmp, &
            p_error)
          
       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       IF (mype == 0) p_communicator_d = p_communicator_tmp

       DO n = 1, SIZE(mapmesh)
          all_debug_pes(n) = n
       END DO

       CALL MPI_GROUP_INCL (group_world, SIZE(mapmesh), all_debug_pes, &
            group_d, p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       CALL MPI_COMM_CREATE (p_all_comm, group_d, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       IF (mype /= 0) p_communicator_d = p_communicator_tmp

    ELSE
       p_communicator_d = p_all_comm
    END IF

    DO n = 0, nproca-1
       members = nprocb
       CALL MPI_GROUP_INCL (group_world, members, mapmesh(:,n), group_a, &
            p_error)
          
       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF
       
       CALL MPI_COMM_CREATE (p_all_comm, group_a, p_communicator_tmp, &
            p_error)
          
       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF
       IF(p_communicator_tmp/=MPI_COMM_NULL) &
         p_communicator_a = p_communicator_tmp
       
    END DO
       
    ! create groups for set Bs
       
    DO n = 0, nprocb-1
       members = nproca
       CALL MPI_GROUP_INCL (group_world, members, mapmesh(n,:), group_b, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF
       
       CALL MPI_COMM_CREATE (p_all_comm, group_b, p_communicator_tmp, &
            p_error)
       
       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF
       IF(p_communicator_tmp/=MPI_COMM_NULL) &
         p_communicator_b = p_communicator_tmp
       
    END DO

    CALL MPI_BARRIER (p_all_comm, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_BARRIER failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    IF (debug_parallel >= 0 .AND. mype == 0) THEN
      p_communicator_a = p_communicator_d
      p_communicator_b = p_communicator_d
    ENDIF
    
#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,3i8)') &
         'p_set_communicator on PE ', mype, ': ', &
         p_communicator_d, &
         p_communicator_a, &
         p_communicator_b
#endif
#endif   
  END SUBROUTINE p_set_communicator

!=========================================================================

  ! send implementation

  SUBROUTINE p_send_real (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, 1, p_real, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real

  SUBROUTINE p_send_real_1d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_1d

  SUBROUTINE p_send_real_2d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real, p_destination, p_tag, &
            p_comm, p_error)
    END IF
       
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_2d

  SUBROUTINE p_send_real_3d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_3d

  SUBROUTINE p_send_real_4d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_4d

  SUBROUTINE p_send_real_5d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm 
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_5d

  SUBROUTINE p_send_int (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, 1, p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int

  SUBROUTINE p_send_int_1d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_1d

  SUBROUTINE p_send_int_2d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_2d

  SUBROUTINE p_send_int_3d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_3d

  SUBROUTINE p_send_int_4d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_4d


  SUBROUTINE p_send_bool (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, 1, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool

  SUBROUTINE p_send_bool_1d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_1d

  SUBROUTINE p_send_bool_2d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_2d

  SUBROUTINE p_send_bool_3d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_3d

  SUBROUTINE p_send_bool_4d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_4d

  SUBROUTINE p_send_char (buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(in) :: buffer
    INTEGER,           INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_char, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, LEN(buffer), p_char, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_char

! non-blocking sends

  SUBROUTINE p_inc_request
    INTEGER, ALLOCATABLE :: tmp(:)

#ifndef NOMPI
    p_irequest = p_irequest + 1
    IF (p_irequest > p_mrequest) THEN
      ALLOCATE(tmp(p_mrequest))
      tmp(:) = p_request(:)
      DEALLOCATE(p_request)
      ALLOCATE(p_request(p_mrequest+p_request_alloc_size))
      p_request(1:p_mrequest) = tmp(:)
      p_mrequest = p_mrequest+p_request_alloc_size
      DEALLOCATE(tmp)
    ENDIF
#endif

  END SUBROUTINE p_inc_request

  SUBROUTINE p_isend_real (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, 1, p_real, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real

  SUBROUTINE p_isend_real_1d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_1d

  SUBROUTINE p_isend_real_2d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF
       
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_2d

  SUBROUTINE p_isend_real_3d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_3d

  SUBROUTINE p_isend_real_4d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_4d

  SUBROUTINE p_isend_real_5d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm 

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_5d

  SUBROUTINE p_isend_int (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, 1, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int

  SUBROUTINE p_isend_int_1d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_1d

  SUBROUTINE p_isend_int_2d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_2d

  SUBROUTINE p_isend_int_3d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_3d

  SUBROUTINE p_isend_int_4d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_4d


  SUBROUTINE p_isend_bool (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, 1, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool

  SUBROUTINE p_isend_bool_1d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_1d

  SUBROUTINE p_isend_bool_2d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_2d

  SUBROUTINE p_isend_bool_3d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_3d

  SUBROUTINE p_isend_bool_4d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_4d

  SUBROUTINE p_isend_char (buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(in) :: buffer
    INTEGER,           INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_char, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, LEN(buffer), p_char, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_char

#ifdef __ONESIDED
  ! some MPI2 things

  SUBROUTINE p_win_create_real_2d(buffer,p_win,comm)

    REAL (dp), INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(out) :: p_win
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_WIN_CREATE(buffer,SIZE(buffer)*p_rg,p_rg,MPI_INFO_NULL,p_comm,p_win,p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_WIN_CREATE on ', mype, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif

#endif
   END SUBROUTINE p_win_create_real_2d

  SUBROUTINE p_win_create_real_3d(buffer,p_win,comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(out) :: p_win
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_WIN_CREATE(buffer,SIZE(buffer)*p_rg,p_rg,MPI_INFO_NULL,p_comm,p_win,p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_WIN_CREATE on ', mype, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif

#endif

  END SUBROUTINE p_win_create_real_3d

  SUBROUTINE p_win_fence_int(p_win)

    INTEGER, INTENT(in) :: p_win

#ifndef NOMPI

    CALL MPI_WIN_FENCE(0,p_win,p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_WIN_FENCE on ', mype, ' for win', p_win, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif

#endif
  END SUBROUTINE p_win_fence_int

  SUBROUTINE p_win_fence_int_1d(p_win)

    INTEGER, INTENT(in) :: p_win(:)
    INTEGER :: i

#ifndef NOMPI

    do i=1,SIZE(p_win)
    CALL MPI_WIN_FENCE(0,p_win(i),p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_WIN_FENCE on ', mype, ' for win', p_win(i), ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
    end do

#endif
  END SUBROUTINE p_win_fence_int_1d

  SUBROUTINE p_win_free_int(p_win)

    INTEGER, INTENT(in) :: p_win

#ifndef NOMPI

    CALL MPI_WIN_FREE(p_win,p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_WIN_FREE on ', mype, ' for win', p_win, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif

#endif
  END SUBROUTINE p_win_free_int

  SUBROUTINE p_win_free_int_1d(p_win)

    INTEGER, INTENT(in) :: p_win(:)
    INTEGER :: i

#ifndef NOMPI

    do i=1,SIZE(p_win)
    CALL MPI_WIN_FREE(p_win(i),p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_WIN_FREE on ', mype, ' for win', p_win(i), ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
    end do

#endif
  END SUBROUTINE p_win_free_int_1d

  SUBROUTINE p_put_real_2d(buffer,p_destination,p_win,comm)

    REAL (dp), INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_win
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_PUT(buffer,SIZE(buffer),p_real,p_destination,0,SIZE(buffer),p_real,p_win,p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a)') ' MPI_PUT on ', mype,'to',p_destination',&
                               'through win', p_win',' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif

#endif
   END SUBROUTINE p_put_real_2d

  SUBROUTINE p_put_real_3d(buffer,p_destination,p_win,comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_win
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_PUT(buffer,SIZE(buffer),p_real,p_destination,0,SIZE(buffer),p_real,p_win,p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a)') ' MPI_PUT on ', mype,'to',p_destination',&
                               'through win', p_win',' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

   END SUBROUTINE p_put_real_3d
#endif
  ! recv implementation

  SUBROUTINE p_recv_real (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, 1, p_real, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real

  SUBROUTINE p_recv_real_1d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_1d

  SUBROUTINE p_recv_real_2d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_2d

  SUBROUTINE p_recv_real_3d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_3d

  SUBROUTINE p_recv_real_4d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_4d

  SUBROUTINE p_recv_real_5d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_5d

  SUBROUTINE p_recv_int (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, 1, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int

  SUBROUTINE p_recv_int_1d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_1d

  SUBROUTINE p_recv_int_2d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_2d

  SUBROUTINE p_recv_int_3d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_3d

  SUBROUTINE p_recv_int_4d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_4d


  SUBROUTINE p_recv_bool (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, 1, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool

  SUBROUTINE p_recv_bool_1d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_1d

  SUBROUTINE p_recv_bool_2d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_2d

  SUBROUTINE p_recv_bool_3d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_3d

  SUBROUTINE p_recv_bool_4d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_4d

  SUBROUTINE p_recv_char (buffer, p_source, p_tag, p_count, comm) 

    CHARACTER(len=*),  INTENT(out) :: buffer
    INTEGER,           INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, LEN(buffer), p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_char

! non-blocking receives

  SUBROUTINE p_irecv_real (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, 1, p_real, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real

  SUBROUTINE p_irecv_real_1d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_real, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_1d

  SUBROUTINE p_irecv_real_2d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_real, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF
       
#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_2d

  SUBROUTINE p_irecv_real_3d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_real, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_3d

  SUBROUTINE p_irecv_real_4d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_real, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_4d

  ! sendrecv implementation

  SUBROUTINE p_sendrecv_real_1d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm
    
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_1d

  SUBROUTINE p_sendrecv_real_2d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm
    
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_2d

  SUBROUTINE p_sendrecv_real_3d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm
    
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_3d

  SUBROUTINE p_sendrecv_real_4d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_4d

  ! bcast implementation

  SUBROUTINE p_bcast_real (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, 1, p_real, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real

  SUBROUTINE p_bcast_real_1d (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_1d

  SUBROUTINE p_bcast_real_2d (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_2d

  SUBROUTINE p_bcast_real_3d (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_3d

  SUBROUTINE p_bcast_real_4d (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_4d

  SUBROUTINE p_bcast_int_i4 (buffer, p_source, comm)

    INTEGER (i4), INTENT(inout) :: buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, 1, p_int_i4, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i4

  SUBROUTINE p_bcast_int_i8 (buffer, p_source, comm)

    INTEGER (i8), INTENT(inout) :: buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, 1, p_int_i8, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i8

  SUBROUTINE p_bcast_int_1d (buffer, p_source, comm)

    INTEGER, INTENT(inout) :: buffer(:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_1d

  SUBROUTINE p_bcast_int_2d (buffer, p_source, comm)

    INTEGER, INTENT(inout) :: buffer(:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_2d

  SUBROUTINE p_bcast_int_3d (buffer, p_source, comm)

    INTEGER, INTENT(inout) :: buffer(:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_3d

  SUBROUTINE p_bcast_int_4d (buffer, p_source, comm)

    INTEGER, INTENT(inout) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_4d


  SUBROUTINE p_bcast_bool (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, 1, p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool

  SUBROUTINE p_bcast_bool_1d (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer(:)
    INTEGER, INTENT(in) :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_1d

  SUBROUTINE p_bcast_bool_2d (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_2d

  SUBROUTINE p_bcast_bool_3d (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_3d

  SUBROUTINE p_bcast_bool_4d (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_4d

  SUBROUTINE p_bcast_char (buffer, p_source, comm)

    CHARACTER(len=*),  INTENT(inout) :: buffer
    INTEGER,           INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, LEN(buffer), p_char, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_char



  ! probe implementation

  SUBROUTINE p_probe_real (buffer, p_tagcount, p_tagtable, p_source, &
&                          p_tag, p_count, comm)

    REAL (dp), INTENT(in)  :: buffer
    INTEGER,   INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,   INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i = 1, p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
&                          flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', mype, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_real, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', mype, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO

#endif

  END SUBROUTINE p_probe_real

  SUBROUTINE p_probe_int (buffer, p_tagcount, p_tagtable, p_source, &
&                          p_tag, p_count, comm)

    INTEGER,   INTENT(in)  :: buffer
    INTEGER,   INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,   INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i=1,p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
&                          flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', mype, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_int, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', mype, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO


#endif

  END SUBROUTINE p_probe_int

  SUBROUTINE p_probe_bool (buffer, p_tagcount, p_tagtable, p_source, &
&                          p_tag, p_count, comm)

    LOGICAL,   INTENT(in)  :: buffer
    INTEGER,   INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,   INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i=1,p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
&                          flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', mype, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_bool, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', mype, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO


#endif

  END SUBROUTINE p_probe_bool

  SUBROUTINE p_probe_char (buffer, p_tagcount, p_tagtable, p_source, &
&                          p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(in)  :: buffer
    INTEGER,           INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,           INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i=1,p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
&                          flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', mype, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_char, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', mype, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO


#endif

  END SUBROUTINE p_probe_char

  SUBROUTINE p_wait
#ifndef NOMPI
    INTEGER :: i

    DO i = 1, p_irequest
        CALL MPI_WAIT(p_request(i), p_status, p_error)
    END DO 
    p_irequest = 0
#endif
  END SUBROUTINE p_wait

  SUBROUTINE p_barrier (comm)
  INTEGER ,INTENT(IN) ,OPTIONAL :: comm
#ifndef NOMPI
    INTEGER :: com
    com = MPI_COMM_WORLD; IF(PRESENT(comm)) com = comm
    CALL MPI_BARRIER (com, p_error)

!#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BARRIER on ', mype, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
!#endif
#endif

  END SUBROUTINE p_barrier

  FUNCTION p_sum_0d (zfield, comm) RESULT (p_sum)

    REAL(dp)                      :: p_sum  
    REAL(dp),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN 
       CALL MPI_ALLREDUCE (zfield, p_sum, 1, p_real, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_0d

  FUNCTION p_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN 
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_real, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_1d

  FUNCTION p_max_0d (zfield, comm) RESULT (p_max)

    REAL(dp)                      :: p_max  
    REAL(dp),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN 
       CALL MPI_ALLREDUCE (zfield, p_max, 1, p_real, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_0d

  FUNCTION p_max_1d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN 
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_1d

  FUNCTION p_max_2d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN 
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_2d

  FUNCTION p_max_3d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))  

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN 
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_3d

  FUNCTION p_min_0d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min  
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN 
       CALL MPI_ALLREDUCE (zfield, p_min, 1, p_real, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_0d

  FUNCTION p_min_1d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN 
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_1d

  FUNCTION p_min_2d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN 
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_2d

  FUNCTION p_min_3d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))  
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN 
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_3d

  SUBROUTINE p_gather_real_1d2d (sendbuf, recvbuf, p_dest, comm)
    REAL(dp),          INTENT(inout) :: sendbuf(:), recvbuf(:,:)
    INTEGER,           INTENT(in) :: p_dest
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

     CALL MPI_GATHER(sendbuf, SIZE(sendbuf), p_real, &
                     recvbuf, SIZE(recvbuf), p_real, &
                     p_dest, p_comm, p_error)
#else
     recvbuf(:,LBOUND(recvbuf,2)) = sendbuf(:)
#endif
   END SUBROUTINE p_gather_real_1d2d

END MODULE mo_mpi
