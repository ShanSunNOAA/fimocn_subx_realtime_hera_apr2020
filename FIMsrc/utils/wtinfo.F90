module module_wtinfo

  use fimnamelist,only:        &
    abort_on_bad_task_distrib, &
    check_omp_consistency,     &
    cpn,                       &
    debugmsg_on,               &
    max_write_tasks_per_node,  &
    mpipn,                     &
    nthreads,                  &
    num_write_tasks,           &
    readnl_serial,             &
    root_own_node

  implicit none

contains

  subroutine wtinfo (comm_in)

! Read namelist and return cores per node (cpn), number of write tasks
! (num_write_tasks), max write tasks per node (mwtpn), and whether root has node
! to himself (root_own_node).
!
! If this module is built in an SMS-parallelized context, insert and compile code
! to read the namelist only on the root task, package the write-task info in a
! user-defined type, and broadcast these values to the other tasks via an MPI
! derived type.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IMPORTANT: Much of the complication in this routine is due to the fact that it is called from
! fimcore.F90, which happens before SMS has been completely set up. Thus, many of the things SMS
! can normally do are not available here, and we have to revert to manual MPI calls. The likely
! failure mode if one tries to replace the MPI stuff in here with SMS is an SMS ASSERTION ERROR
! due to SMS not having been correctly initialized.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USEMPIMOD
!sms$insert use mpi
#endif

    implicit none

#ifndef USEMPIMOD
!sms$insert include 'mpif.h'
#endif

    type wtconf
      sequence
      integer::cpn
      integer::mpipn
      integer::mwtpn
      integer::nwt
      integer::nthreads ! number of OMP threads
      logical::eomp     ! check OpenMP consistency
      logical::aobtd
      logical::do
      logical::ron
    end type wtconf


    integer,intent(in),optional :: comm_in      ! an MPI intracommunicator

    integer       :: comm,ierrcode,istatus,me=0,comm_wtconf,size_i
    integer       :: count(2),offset(2),type(2)
    logical, save :: alreadyReadWriteTaskInfo = .false. ! initial value
    type(wtconf)  :: conf

! Prefer the passed-in MPI communicator, if provided.

!sms$insert if (present(comm_in)) then
!sms$insert comm=comm_in
!sms$insert else
!sms$insert comm=mpi_comm_world
!sms$insert endif

    if (.not.alreadyReadWriteTaskInfo) then

!sms$insert call mpi_comm_rank(comm,me,istatus)

!sms$insert if (me.eq.0) then

      call readnl_serial

      if (cpn == -1) then
!sms$ignore begin
        write(6,*) 'wtinfo ERROR: No value for cpn found in namelist. Stopping'
!sms$ignore end
        stop
      end if

      if (mpipn == -1) then
!sms$ignore begin
        write(6,*) 'wtinfo NOTE: No value for mpipn found in namelist. Setting to cpn=', cpn
!sms$ignore end
        mpipn = cpn
      end if

! Populate conf with individual values.

!sms$insert conf%aobtd    = abort_on_bad_task_distrib
!sms$insert conf%cpn      = cpn
!sms$insert conf%do       = debugmsg_on
!sms$insert conf%mwtpn    = max_write_tasks_per_node
!sms$insert conf%nwt      = num_write_tasks
!sms$insert conf%ron      = root_own_node
!sms$insert conf%eomp     = check_omp_consistency
!sms$insert conf%mpipn    = mpipn
!sms$insert conf%nthreads = nthreads

!sms$insert endif ! me.eq.0

! Determine the size in bytes of the MPI integer type

!sms$insert call mpi_type_size(mpi_integer,size_i,istatus)
!sms$insert if (istatus.ne.0) then
!sms$ignore begin
!sms$insert write (*,'(a,i0)') 'wtinfo: MPI_Type_extent returned ',istatus
!sms$insert call mpi_abort(comm,ierrcode,istatus)
!sms$insert stop
!sms$ignore end
!sms$insert endif

! Populate the arrays defining the MPI derived type
! COUNT THE NUMBER OF THINGS OF THAT TYPE IN WTCONF!!!!

!sms$insert count(1)=5
!sms$insert offset(1)=0
!sms$insert type(1)=mpi_integer
!sms$insert count(2)=4
!sms$insert offset(2)=offset(1)+(count(1)*size_i)
!sms$insert type(2)=mpi_logical

! Create the MPI derived-type struct

!sms$insert call mpi_type_struct(2,count,offset,type,comm_wtconf,istatus)
!sms$insert if (istatus.ne.0) then
!sms$ignore begin
!sms$insert write (*,'(a,i0)') 'wtinfo: MPI_Type_struct returned ',istatus
!sms$insert call mpi_abort(comm,ierrcode,istatus)
!sms$insert stop
!sms$ignore end
!sms$insert endif

! Commit the MPI derived type

!sms$insert call mpi_type_commit(comm_wtconf,istatus)
!sms$insert if (istatus.ne.0) then
!sms$ignore begin
!sms$insert write (*,'(a,i0)') 'wtinfo: MPI_Type_commit returned ',istatus
!sms$insert call mpi_abort(comm,ierrcode,istatus)
!sms$insert stop
!sms$ignore end
!sms$insert endif

! Broadcast the writetask configuration info

!sms$insert call mpi_bcast(conf,1,comm_wtconf,0,comm,istatus)
!sms$insert if (istatus.ne.0) then
!sms$ignore begin
!sms$insert write (*,'(a,i0)') 'wtinfo: MPI_Bcast returned ',istatus
!sms$insert call mpi_abort(comm,ierrcode,istatus)
!sms$insert stop
!sms$ignore end
!sms$insert endif

! Non-root tasks unpack conf to individual values.

!sms$insert if (me.ne.0) then
!sms$insert abort_on_bad_task_distrib = conf%aobtd
!sms$insert cpn                       = conf%cpn
!sms$insert debugmsg_on               = conf%do
!sms$insert max_write_tasks_per_node  = conf%mwtpn
!sms$insert num_write_tasks           = conf%nwt
!sms$insert root_own_node             = conf%ron
!sms$insert check_omp_consistency     = conf%eomp
!sms$insert mpipn                     = conf%mpipn
!sms$insert nthreads                  = conf%nthreads
!sms$insert endif ! me.ne.0

      alreadyReadWriteTaskInfo=.true.

    endif ! .not.alreadyReadWriteTaskInfo

    return

  end subroutine wtinfo

end module module_wtinfo
