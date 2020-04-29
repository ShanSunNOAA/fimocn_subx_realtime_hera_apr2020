! TODO Replace "write" nomenclature with "io", as "write" tasks now read, too.

module icosio

! NOTE Two cpp tokens, SERIAL and NOGRIB, guard portions of this code. When
! NOTE SERIAL is set, code sections that introduce a dependence on MPI are
! NOTE removed by cpp. Therefore, if SERIAL is *not* set, icosio must be linked
! NOTE with MPI. Likewise, when NOGRIB is set, code sections relating to grib
! NOTE output file production are removed by cpp. For SERIAL, the following
! NOTE convention has been followed: If a routine is unneeded in serial mode but
! NOTE introduces no MPI dependence, it is left alone; if it is unneeded but
! NOTE does introduce an MPI dependence, it is completely removed; and if it
! NOTE is needed in serial mode, only the lines introducing MPI dependencies
! NOTE are removed. If NOGRIB is not set, the three routines post_init_file,
! NOTE post_write_field, and post_finalize_file must be available at link time.

! TODO Ideally, the necessary grib-production code would be moved into and built
! TODO as part of icosio.

  use stringlist

#ifndef SERIAL
#ifdef USEMPIMOD
  use mpi
#endif /* USEMPIMOD */
#endif /* SERIAL */

  implicit none

  private

#ifndef SERIAL

#ifndef USEMPIMOD
  include 'mpif.h'
#endif /* USEMPIMOD */

! Public subroutines defining the icosio API: for parallel builds only

  public icosio_run

#endif /* SERIAL */

! Public subroutines defining the icosio API: for parallel and serial builds

  public icosio_in_receive,icosio_in_request,icosio_in_start,icosio_in_stop,&
    icosio_out,icosio_out_stop,icosio_setup,icosio_stop

! Module parameters (for serial and parallel builds)

#ifndef SERIAL
  integer,parameter::cmdsize      = 4    ! Size of frame command
  integer,parameter::num_i_params = 2    ! # of integer config parameters
  integer,parameter::num_l_params = 4    ! # of logical config parameters
#endif /* SERIAL */

! Module variables whose values will be computed, with default/dummy values
! where possible.

  character*2::tasktype        = 'ct'    ! May be set to 'wt' later
  integer::ierr                = -1      ! For error checks
  integer::me                  =  0      ! 0..(n-1) MPI rank
  integer::nct                 =  1      ! Number of compute tasks
  integer::nwt                 =  0      ! Number of write tasks
  integer::outfiles            =  0      ! Number of assigned write tasks
  integer::wtindex             =  0      ! Which write task receives data next
  logical::i_am_compute_root   = .true.  ! Am I serial / compute intracomm root?
  logical::i_am_compute_task   = .true.  ! Am I a compute task?
  logical::i_am_write_root     = .false. ! Am I root of a write intracomm?
  logical::icosio_setup_called = .false. ! Was icosio_setup() called?
  logical::in_receive_called   = .false. ! Was icosio_in_receive() called?
  logical::in_request_called   = .false. ! Was icosio_in_request() called?
  logical::in_start_called     = .false. ! Was icosio_in_start() called?
  logical::serial              = .true.  ! Am I running serial?
  logical::single              = .true.  ! Am I serial or alone in my intracomm?

  character*1000::msg
  integer,pointer::inv_perm_global(:)=>null()
  integer,pointer::perm_global(:)=>null()
  type(stringlist_node),pointer::current_filename=>null()
  type(stringlist_node),pointer::input_filenames=>null()
  type(stringlist_node),pointer::output_filenames=>null()

#ifndef SERIAL
  integer::first_reader        = -1      ! First round-robin assigned reader
  integer::intercomm           = -1      ! Compute/write intercomm
  integer::intracomm           = -1      ! Communicator for my group
  integer::reader              = -1      ! Currently-selected input task
  integer::size_c              = -1      ! Size of a character, per mpi
  integer::size_i              = -1      ! Size of an integer, per mpi
  integer::size_l              = -1      ! Size of a logical, per mpi
  integer::size_r              = -1      ! Size of a real, per mpi

  integer,allocatable::bounds(:,:)
  integer,allocatable::interior_sizes(:)
  logical,allocatable::readers(:)
#endif /* SERIAL */

! Caller-set variables (some sent from compute to write tasks) with defaults

  integer::ipe               = -1      ! Bounds: upper inner
  integer::ips               = -1      ! Bounds: lower inner
  integer::lun               =  1      ! Use this lun for disk writes
  integer::nip               = -1      ! Number of icosahedral points
  logical::binout            = .true.  ! Write non-grib history files?
  logical::client_server_io  = .true.  ! Use client/server io?
  logical::debugmsg_on       = .false. ! Print verbose status msgs?
  logical::gribout           = .false. ! Write grib files?
  logical::i_am_write_task   = .false. ! Am I a write task?
  logical::permute           = .false. ! Should output arrays be permuted?
  logical::print_diags       = .false. ! Only print diagnostics?
  logical::using_write_tasks = .false. ! Are we using write tasks?

  integer,pointer::inv_perm(:)=>null()
  integer,pointer::perm(:)=>null()

#ifndef SERIAL
  integer::comm              = -1      ! An MPI communicator
#endif /* SERIAL */

! MPI comm tags

#ifndef SERIAL
  integer,parameter::tag_cmd                      = 100
  integer,parameter::tag_collect_bounds           = 101
  integer,parameter::tag_collect_inv_perm_segment = 102
  integer,parameter::tag_collect_var_segment      = 103
  integer,parameter::tag_data                     = 104
  integer,parameter::tag_in_filename              = 105
  integer,parameter::tag_in_receive               = 106
  integer,parameter::tag_integer_parameters       = 107
  integer,parameter::tag_interior_sizes           = 108
  integer,parameter::tag_inv_perm_data            = 109
  integer,parameter::tag_logical_parameters       = 110
  integer,parameter::tag_metadata                 = 111
  integer,parameter::tag_perm_data                = 112
  integer,parameter::tag_permute_control          = 113
  integer,parameter::tag_scatter_input_var        = 114
#endif /* SERIAL */

! Command values

#ifndef SERIAL
  integer,parameter::cmd_in_receive = 1
  integer,parameter::cmd_in_request = 2
  integer,parameter::cmd_in_start   = 3
  integer,parameter::cmd_in_stop    = 4
  integer,parameter::cmd_out        = 5
  integer,parameter::cmd_stop       = 6
#endif /* SERIAL */

! Input buffers

  type in_buffer_node
    real,pointer::data(:)              => null() ! Variable data
    integer::data_size                 =  0      ! Size of data record
    character,pointer::filename(:)     => null() ! Input filename
    integer::filename_len              =  0      ! Filename length
    type(in_buffer_node),pointer::next => null() ! Pointer to next list node
  end type in_buffer_node

  type in_buffer
    type(in_buffer_node),pointer::tail=>null()
    type(in_buffer_node),pointer::head=>null()
  end type in_buffer

#ifndef SERIAL
  type(in_buffer),pointer::in_buffers ! Input buffers
#endif /* SERIAL */

! Output buffers

  type out_buffer_node
    real::scalefactor               =  1.     ! GRIB scale factor
    integer::accum_start            = -1      ! GRIB time-accumulation factor
    integer::filename_len           = -1      ! Actual filename length
    integer::header_len             = -1      ! Actual header length
    integer::levels                 = -1      ! Vertical levels in variable
    integer::segment_size           = -1      ! Bytes in variable segment
    integer::time                   = -1      ! Time variable represents
    integer::varname_len            = -1      ! Actual varname length
    character,pointer::filename(:)  => null() ! Output filename
    character,pointer::header(:)    => null() ! Header
    character,pointer::varname(:)   => null() ! Variable name
    real,pointer::segment(:)        => null() ! Variable data
    type(out_buffer_node),pointer::next => null() ! Pointer to next list node
  end type out_buffer_node

  type out_buffer
    integer::segments=0
    type(out_buffer_node),pointer::tail=>null()
    type(out_buffer_node),pointer::head=>null()
  end type out_buffer

  type(out_buffer),pointer::out_buffers(:) ! Output buffers

! Pointers to globally-sized data buffers

  type glbvar_ptrs
    real,pointer::ptr(:,:) => null()
  end type glbvar_ptrs

! Generic interfaces

  interface icosio_in_receive
    module procedure icosio_in_receive_2dr4
    module procedure icosio_in_receive_3dr4
  end interface icosio_in_receive

  interface icosio_out
    module procedure icosio_out_2dr4
    module procedure icosio_out_2dr8
    module procedure icosio_out_3dr4
    module procedure icosio_out_3dr8
  end interface icosio_out

contains

!-------------------------------------------------------------------------------
  subroutine append_to_list(its,task,varname,levels,segment_size,segment,&
    filename,header,time,scalefactor,accum_start)
!-------------------------------------------------------------------------------

! Attach a new node to the linked list of output buffer nodes for the specified
! write task. If the optional variable name is not specified an empty node is
! appended; otherwise, it is assumed that the rest of the optional arguments are
! also present, and the metadata and field data are filled in with the supplied
! values. This is a private routine -- not part of the icosio API. It is called
! by both compute and write tasks.

    character(len=*),intent(in),optional::filename,header,varname
    integer,intent(in),optional::levels,segment_size,time,accum_start
    integer,intent(in)::its,task
    real,intent(in),optional::scalefactor
    real,pointer,optional::segment(:)

    character(len=100)::debughdr
    character(len=2)::recipient_type
    character(len=14)::this='append_to_list'
    type(out_buffer),pointer::bp
    type(out_buffer_node),pointer::newnode

    call setdebughdr(debughdr,this,its)

! Allocate a new node, set its metadata and point it at the passed-in variable
! data segment.

    allocate(newnode,stat=ierr)
    write (msg,'(2a)') trim(debughdr),&
      ' Failed to allocate new output buffer node.'
    call die_if(ierr.ne.0,ierr,msg)

! Note that, if "varname" is present, we assume that all the optional arguments
! are present. We deallocate "segment" in clear_output_list().

    if (present(varname)) then
      newnode%scalefactor=scalefactor
      newnode%accum_start=accum_start
      newnode%filename_len=len(filename)
      newnode%header_len=len(header)
      newnode%levels=levels
      newnode%segment_size=segment_size
      newnode%time=time
      newnode%varname_len=len(varname)
      allocate(newnode%filename(len(filename)),stat=ierr)
      write (msg,'(2a)') trim(debughdr),' Failed to allocate filename member.'
      call die_if(ierr.ne.0,ierr,msg)
      newnode%filename=str_to_arr(filename)
      allocate(newnode%header(len(header)),stat=ierr)
      write (msg,'(2a)') trim(debughdr),' Failed to allocate header member.'
      call die_if(ierr.ne.0,ierr,msg)
      newnode%header=str_to_arr(header)
      allocate(newnode%varname(len(varname)),stat=ierr)
      write (msg,'(2a)') trim(debughdr),' Failed to allocate varname member.'
      call die_if(ierr.ne.0,ierr,msg)
      newnode%varname=str_to_arr(varname)
      newnode%segment=>segment
      write (msg,'(4a)') trim(debughdr),' Set variable metadata for ',&
        trim(varname),'.'
      call debugmsg
    end if

! Append the new node to the appropriate linked list.

    bp=>out_buffers(task)
    if (associated(bp%head)) then
      bp%tail%next=>newnode
    else
      bp%head=>newnode
    end if
    bp%tail=>newnode
    bp%segments=bp%segments+1

! Compute tasks buffer data for write tasks, and vice versa. Get messages right.

    if (tasktype=='ct') then
      recipient_type='wt'
    else
      recipient_type='ct'
    end if
    if (present(varname)) then
      write (msg,'(6a,i0,a)') trim(debughdr),' Appended ',trim(varname),&
        ' to list for ',recipient_type,' ',task,'.'
      call debugmsg
    else
      write (msg,'(4a,i0,a)') trim(debughdr),&
        ' Appended blank node to list for ',recipient_type,' ',task,'.'
      call debugmsg
    end if

  end subroutine append_to_list

!-------------------------------------------------------------------------------
  function arr_eq_str(arr,str,count)
!-------------------------------------------------------------------------------

! Returns true if characters 1..count of the supplied character array and string
! match, and false otherwise. This is a private routine -- not part of the
! icosio API. It may be called by either compute or write tasks.

    character,intent(in)::arr(:)
    character(len=*),intent(in)::str
    integer,intent(in)::count

    logical::arr_eq_str
    integer::i

    arr_eq_str=.true.
    do i=1,count
      if (arr(i).ne.str(i:i)) arr_eq_str=.false.
    end do

  end function arr_eq_str

!-------------------------------------------------------------------------------
  function arr_to_str(arr,count)
!-------------------------------------------------------------------------------

! Return a character string corresponding to the supplied character array. This
! is a private routine -- not part of the icosio API. It may be called by either
! compute or write tasks.

    character,intent(in)::arr(:)
    integer,intent(in)::count

    character(len=count)::arr_to_str
    integer::i

    do i=1,count
      arr_to_str(i:i)=arr(i)
    end do

  end function arr_to_str

!-------------------------------------------------------------------------------
  subroutine buffer_var(its,varname,var,filename,header,time,scalefactor,&
    accum_start)
!-------------------------------------------------------------------------------

! Buffer data for eventual transmission to write task(s), which is triggered by
! a flush_all() call. This is a private routine -- not part of the icosio API.
! It is called only by compute tasks.

    character(len=*),intent(in)::filename,header,varname
    integer,intent(in)::accum_start,its,time
    real,intent(in)::scalefactor
    real,intent(in)::var(:,:)

    character(len=100)::debughdr
    character(len=10)::this='buffer_var'
    integer::i,j,levels,offset,segment_size,wt
    integer,parameter::notfound=-1
    real,pointer::segment(:)
    type(out_buffer_node),pointer::node

    call setdebughdr(debughdr,this,its)

    write (msg,'(4a)') trim(debughdr),' varname=',trim(varname),' entry.'
    call debugmsg

! Lay out variable data in a 1D buffer.

    levels=size(var,1)
    segment_size=levels*interior_size()
    allocate(segment(segment_size),stat=ierr)
    write (msg,'(2a)') trim(debughdr),' Failed to allocate segment.'
    call die_if(ierr.ne.0,ierr,msg)
    offset=1
    do i=1,interior_size()
      do j=1,levels
        segment(offset)=var(j,i)
        offset=offset+1
      end do
    end do
    write (msg,'(4a)') trim(debughdr),' Laid out ',trim(varname),&
      ' in 1D buffer.'
    call debugmsg

! Determine which write task will receive this variable. mod() ensures round-
! robin assignment of responsibility for output files to write tasks. When there
! are fewer write tasks than output files, multiple variables will be assigned
! to some write task(s). When there are more write tasks than output files,
! write tasks not used in one output frame may be used in the next. In a future
! enhancement, it may be possible to split a single variable across two or more
! write tasks for output, for scalability.

! Search the output buffers for a matching filename. If we find one, the same
! write task should receive this variable segment as well.

    wt=notfound
    do i=0,nwt-1
      node=>out_buffers(i)%head
      do while (associated(node).and.wt.eq.notfound)
        if (arr_eq_str(node%filename,filename,node%filename_len)) wt=i
        node=>node%next
      end do
      if (wt.ne.notfound) exit
    end do

! If no matching filename was found, do round-robin assignment.

    if (wt.eq.notfound) then
      outfiles=outfiles+1
      wtindex=mod(wtindex,nwt)
      wt=wtindex
      wtindex=wtindex+1
      write (msg,'(4a,i0,a)') trim(debughdr),' ',trim(varname),&
        ' assigned to wt ',wt,' (round-robin).'
      call debugmsg
    else
      write (msg,'(4a,i0,a)') trim(debughdr),' ',trim(varname),&
        ' assigned to wt ',wt,' (filename match).'
      call debugmsg
    end if

! Buffer data to send later.

    call append_to_list(its,wt,varname,levels,segment_size,segment,filename,&
      header,time,scalefactor,accum_start)

    write (msg,'(4a)') trim(debughdr),' varname=',trim(varname),' exit.'
    call debugmsg

  end subroutine buffer_var

!-------------------------------------------------------------------------------
  subroutine check_setup_called(caller,its)
!-------------------------------------------------------------------------------

! Return with error if icosio_setup() has not been called. This is a private
! routine -- not part of the icosio API. It is called by both compute and write
! tasks.

    character(len=*),intent(in)::caller
    integer,intent(in)::its

    character(len=100)::debughdr
    character(len=18)::this='check_setup_called'

    call setdebughdr(debughdr,this,its)

    write (msg,'(4a)') trim(debughdr),' icosio_setup has not been called in ',&
      trim(caller),'.'
    call die_if(.not.icosio_setup_called,ierr,msg)

  end subroutine check_setup_called

!-------------------------------------------------------------------------------
  subroutine clear_output_list(its)
!-------------------------------------------------------------------------------

! Walk the list of output data buffers for the specified write task and de-
! allocate dynamic buffers. Reset the buffer-head pointers & counter. This is a
! private routine -- not part of the icosio API. It is called by both compute
! and write tasks.

    integer,intent(in)::its

    character(len=100)::debughdr
    character(len=17)::this='clear_output_list'
    integer::i
    type(out_buffer_node),pointer::next,node

    call setdebughdr(debughdr,this,its)

    write (msg,'(2a)') trim(debughdr),' Entry.'
    call debugmsg

    if (associated(out_buffers)) then
      do i=0,size(out_buffers)-1
        if (associated(out_buffers(i)%head)) then
          node=>out_buffers(i)%head
          do while (associated(node))
            next=>node%next
            deallocate(node%filename)
            deallocate(node%header)
            deallocate(node%varname)
            if (associated(node%segment)) deallocate(node%segment)
            nullify(node%segment)
            deallocate(node)
            node=>next
          end do
        end if
        out_buffers(i)%segments=0
        nullify(out_buffers(i)%tail)
        nullify(out_buffers(i)%head)
        nullify(node)
      end do
    else
      write (msg,'(2a)') trim(debughdr),&
        ' Called but output buffers not allocated.'
      call die(ierr,msg)
    end if

    write (msg,'(2a)') trim(debughdr),' Exit.'
    call debugmsg

  end subroutine clear_output_list

#ifndef SERIAL

!-------------------------------------------------------------------------------
  subroutine collect_bounds(its)
!-------------------------------------------------------------------------------

! This is a private routine -- not part of the icosio API. It is called only by
! compute tasks, and only when write tasks are not in use.

    integer,intent(in)::its

    character(len=100)::debughdr
    character(len=14)::this='collect_bounds'
    integer::bound(2),ct,expected,mpistatus(mpi_status_size),mpitype,received,&
      tag
    logical,save::already_called=.false.

    call setdebughdr(debughdr,this,its)

    if ((.not.already_called).and.(.not.single)) then

      write (msg,'(2a)') trim(debughdr),' Entry.'
      call debugmsg

      if (i_am_compute_root) then

        allocate(bounds(2,nct-1),stat=ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to allocate bounds.'
        call die_if(ierr.ne.0,ierr,msg)

        do ct=1,nct-1 ! First non-root compute task has task id 1.

          write (msg,'(2a,i0,a)') trim(debughdr),&
            ' Receiving bounds from ct ',ct,'...'
          call debugmsg

          expected=2
          mpitype=mpi_integer
          tag=tag_collect_bounds

          call mpi_recv(bound,expected,mpitype,ct,tag,intracomm,mpistatus,ierr)
          write (msg,'(2a,i0,a)') trim(debughdr),&
            ' Failed to receive bounds from ',ct,'.'
          call die_if(ierr.ne.mpi_success,ierr,msg)

          call mpi_get_count(mpistatus,mpitype,received,ierr)
          write (msg,'(2a,2(i0,a))') trim(debughdr),' bound expected ',&
            expected,' elements, received ',received,'.'
          call die_if(received.ne.expected,ierr,msg)

          write (msg,'(2a,3(i0,a))') trim(debughdr),' Received  bounds (',&
            bound(1),':',bound(2),') from ct ',ct,'.'
          call debugmsg

          bounds(1,ct)=bound(1)
          bounds(2,ct)=bound(2)

        end do

      else ! not the compute root...

        write (msg,'(2a,2(i0,a))') trim(debughdr),' Sending bounds (',&
          bound(1),':',bound(2),') to root.'
        call debugmsg

        bound(1)=ips
        bound(2)=ipe

        call mpi_ssend(bound,2,mpi_integer,0,tag_collect_bounds,intracomm,ierr)
        write (msg,'(2a)') trim(debughdr),&
          ' Failed to send bounds to compute root.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        write (msg,'(2a,2(i0,a))') trim(debughdr),' Sent    bounds (',&
          bound(1),':',bound(2),') to root.'
        call debugmsg

      end if

      already_called=.true.

      write (msg,'(2a)') trim(debughdr),' Exit.'
      call debugmsg

    end if

  end subroutine collect_bounds

#endif /* SERIAL */

!-------------------------------------------------------------------------------
  subroutine collect_inv_perm(its)
!-------------------------------------------------------------------------------

! Collect the various segments of the inverse grid permutation array on the
! compute root. This is a private routine -- not part of the icosio API. It is
! called only by compute tasks, and only when write tasks are not in use.

! TODO Consider using MPI 2 intercomm collective operations (probably MPI_Gather
! TODO and MPI_Gatherv) to get interior sizes and then global inv_perm.

    character(len=100)::debughdr
    character(len=16)::this='collect_inv_perm'
    integer,intent(in)::its
#ifndef SERIAL
    integer::ct,expected,inv_perm_global_size,mpistatus(mpi_status_size),&
      mpitype,received,segment_size,tag
#endif /* SERIAL */

    call setdebughdr(debughdr,this,its)

#ifdef SERIAL

    inv_perm_global => inv_perm

#else

    write (msg,'(2a)') trim(debughdr),' Entry.'
    call debugmsg

! Do some initial error checking.

    write (msg,'(4a)') trim(debughdr),' Calling ',this,' twice is a bug.'
    call die_if(associated(inv_perm_global),ierr,msg)

    write (msg,'(4a)') trim(debughdr),' A write task calling ',this,' is a bug.'
    call die_if(associated(inv_perm_global),ierr,msg)

! The association status of inv_perm_global is used to signal that the compute
! root has already collected the necessary global-size inv_perm, so both root
! and non-root compute tasks need to allocate it. Allocate the inv_perm_global
! of the compute root at global size, and allocate those of other compute task
! at unit size.

    if (i_am_compute_root) then
      inv_perm_global_size=nip
    else
      inv_perm_global_size=1
    end if
    allocate(inv_perm_global(inv_perm_global_size),stat=ierr)
    write (msg,'(2a)') trim(debughdr),' Failed to allocate inv_perm_global.'
    call die_if(ierr.ne.0,ierr,msg)

    call collect_bounds(its)

    compute_root: if (i_am_compute_root) then

! The compute root already has its own local segment: Copy it into the correct
! location in inv_perm_global.

      inv_perm_global(ips:ipe)=inv_perm(ips:ipe)

! If running in single mode (serial or one-task MPI), inv_perm_global=inv_perm
! on the one and only task, and the previous assignment statement is sufficient,
! so simply return here.

      if (single) then
        write (msg,'(2a)') trim(debughdr),' Single task returning...'
        call debugmsg
        return
      end if

! Otherwise, the compute root collects the local inv_perm segments from each
! compute task to assemble a global-size inv_perm.

      do ct=1,nct-1 ! First non-root compute task has task id 1.

! Receive a segment into the appropriate location in the global array.

        write (msg,'(2a,i0,a)') trim(debughdr),&
          ' Receiving inv_perm segment from ',ct,'.'
        call debugmsg

        segment_size=bounds(2,ct)-bounds(1,ct)+1
        expected=segment_size
        mpitype=mpi_integer
        tag=tag_collect_inv_perm_segment

        call mpi_recv(inv_perm_global(bounds(1,ct):bounds(2,ct)),expected,&
          mpitype,ct,tag,intracomm,mpistatus,ierr)
        write (msg,'(2a,i0,a)') trim(debughdr),&
          ' Failed to receive inv_perm segment from ',ct,'.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_get_count(mpistatus,mpitype,received,ierr)
        write (msg,'(2a,2(i0,a))') trim(debughdr),&
          ' inv_perm_global expected ',expected,' elements, received ',&
          received,'.'
        call die_if(received.ne.expected,ierr,msg)

        write (msg,'(2a,i0,a)') trim(debughdr),&
          ' Received  inv_perm segment from ',ct,'.'
        call debugmsg

      end do

    else ! I am a non-root compute task...

      inv_perm_global=-1 ! A default value.

! Send the segment of this task to the compute root.

      write (msg,'(2a)') trim(debughdr),' Sending inv_perm segment to root...'
      call debugmsg

      segment_size=interior_size()

      call mpi_ssend(inv_perm(ips:ipe),segment_size,mpi_integer,0,&
        tag_collect_inv_perm_segment,intracomm,ierr)
      write (msg,'(2a)') trim(debughdr),&
        ' Failed to send inv_perm segment to compute root.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a)') trim(debughdr),' Sent    inv_perm segment to root.'
      call debugmsg

    end if compute_root

    write (msg,'(2a)') trim(debughdr),' Exit.'
    call debugmsg

#endif /* SERIAL */

  end subroutine collect_inv_perm

#ifndef SERIAL

!-------------------------------------------------------------------------------
  subroutine collect_var(its,var,glbvar)
!-------------------------------------------------------------------------------

! Collect the various segments of a global model field array on the compute
! root. This is a private routine -- not part of the icosio API. It is called
! only by compute tasks, and only when write tasks are not in use.

! TODO Consider using MPI 2 intercomm collective operations (probably MPI_Gather
! TODO and MPI_Gatherv) to get interior sizes and then global fields.

    integer,intent(in)::its
    real,intent(in)::var(:,:)
    real,pointer::glbvar(:,:)

    character(len=100)::debughdr
    character(len=11)::this='collect_var'
    integer::ct,expected,levels,mpistatus(mpi_status_size),mpitype,received,&
      segment_size,tag

    call setdebughdr(debughdr,this,its)

    call collect_bounds(its)

    levels=size(var,1)

    compute_root: if (i_am_compute_root) then

! Allocate space for a global-size variable.

      allocate(glbvar(levels,nip),stat=ierr)
      write (msg,'(2a)') trim(debughdr),' Failed to allocate glbvar.'
      call die_if(ierr.ne.0,ierr,msg)

! Copy my local segment into the global array.

      glbvar(:,1:interior_size())=var(:,1:interior_size())

! If nct=1 (i.e. we are single) this loop will not execute.

      do ct=1,nct-1 ! First non-root compute task has task id 1.

        segment_size=levels*(bounds(2,ct)-bounds(1,ct)+1)

! Receive a segment into the appropriate location in the global array.

        write (msg,'(2a,2(i0,a))') trim(debughdr),' Receiving ',segment_size,&
          '-element variable segment from ',ct,'...'
        call debugmsg

        expected=segment_size
        mpitype=mpi_real
        tag=tag_collect_var_segment

        call mpi_recv(glbvar(:,bounds(1,ct):bounds(2,ct)),expected,mpitype,&
          ct,tag,intracomm,mpistatus,ierr)
        write (msg,'(2a,i0,a)') trim(debughdr),&
          ' Failed to receive var segment from ',ct,'.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_get_count(mpistatus,mpitype,received,ierr)
        write (msg,'(2a,2(i0,a))') trim(debughdr),' glbvar expected ',&
          expected,' elements, received ',received,'.'
        call die_if(received.ne.expected,ierr,msg)

        write (msg,'(2a,2(i0,a))') trim(debughdr),' Received  ',segment_size,&
          '-element variable segment from ',ct,'.'
        call debugmsg

      end do

    else ! I am a non-root compute task...

      segment_size=levels*interior_size()

! Send the segment of this task to the compute root.

      write (msg,'(2a,i0,a)') trim(debughdr),' Sending ',segment_size,&
        '-element variable segment to root...'
      call debugmsg

      call mpi_ssend(var(:,1:interior_size()),segment_size,mpi_real,0,&
        tag_collect_var_segment,intracomm,ierr)
      write (msg,'(2a)') trim(debughdr),&
        ' Failed to send var segment to compute root.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a,i0,a)') trim(debughdr),' Sent    ',segment_size,&
        '-element variable segment to root.'
      call debugmsg

    end if compute_root

  end subroutine collect_var

#endif /* SERIAL */

#ifndef SERIAL

!-------------------------------------------------------------------------------
  subroutine comm_type_sizes
!-------------------------------------------------------------------------------

! Set size variables for communication purposes. Query MPI for the sizes of
! intrinsic types. This is a private routine -- not part of the icosio API. It
! is called by both compute and write tasks.

    logical::set=.false.

    if (.not.set) then
      call mpi_type_size(mpi_character,size_c,ierr)
      call mpi_type_size(mpi_integer,  size_i,ierr)
      call mpi_type_size(mpi_logical,  size_l,ierr)
      call mpi_type_size(mpi_real,     size_r,ierr)
      set=.true.
    end if

  end subroutine comm_type_sizes

#endif /* SERIAL */

!-------------------------------------------------------------------------------
  subroutine debugmsg
!-------------------------------------------------------------------------------

! Print a message to stdout if verbose debugging is enabled. This is a private
! routine -- not part of the icosio API. It is called by both compute and write
! tasks.

    if (debugmsg_on) then
      write (6,'(2a)') 'ICOSIO DEBUG: ',trim(msg)
      call flush(6)
    end if

  end subroutine debugmsg

!-------------------------------------------------------------------------------
  subroutine die(i,message,oldstatus)
!-------------------------------------------------------------------------------

! End program after printing the specified message to stdout. This is a private
! routine -- not part of the icosio API. It is called by both compute and write
! tasks.

    integer,intent(inout)::i
    character(len=*),intent(in)::message
    integer,intent(in),optional::oldstatus

    if (present(oldstatus)) then
      write (*,'(2a,i0,a,i0)') trim(message),' i=',i,' oldstatus=',oldstatus
    else
      write (*,*) trim(message)
    end if
    call flush(6)

#ifndef SERIAL

    call mpi_abort(mpi_comm_world,i,ierr)

#endif /* SERIAL */

    stop

  end subroutine die

!-------------------------------------------------------------------------------
  subroutine die_if(condition,i,message,oldstatus)
!-------------------------------------------------------------------------------

! End program and print message if condition is true. This is a private
! routine -- not part of the icosio API. It is called by both compute and write
! tasks.

    logical,intent(in)::condition
    integer,intent(inout)::i
    character(len=*),intent(in)::message
    integer,intent(in),optional::oldstatus

    if (condition) then
      if (present(oldstatus)) then
        call die(i,message,oldstatus)
      else
        call die(i,message)
      end if
    end if

  end subroutine die_if

#ifndef SERIAL

!-------------------------------------------------------------------------------
  subroutine flush_all(its)
!-------------------------------------------------------------------------------

! Send accumulated metadata and field-data segments to all write tasks involed
! in the output of data for this frame. This is a private routine -- not part of
! the icosio API. If is called only by compute tasks.

    integer,intent(in)::its

    character(len=100)::debughdr
    character(len=9)::this='flush_all'
    integer::i,j,k
    logical::warning_given=.false.

    call setdebughdr(debughdr,this,its)

! Return if buffers have not accumulated any output data

    if (out_buffers(0)%segments.eq.0) return

! Print one-per-file write task layout instructions.

    if (i_am_compute_root.and..not.warning_given) then
      if (nwt.eq.1) then
        write (*,'(a,i0,a)') 'NOTE: Using one write task for ',outfiles,&
          ' output files.'
      else if (outfiles.ne.nwt) then
        write (*,'(a,i0,a)') 'NOTE: Set num_write_tasks=',outfiles,&
          ' in namelist file for a one-per-file write-task layout.'
      end if
      warning_given=.true.
    end if

! Flush the output buffers of all write tasks holding data from the current
! frame.

    if (outfiles.ge.nwt) then
      i=0
      j=nwt-1
    else
      i=wtindex-outfiles
      if (i.lt.0) i=i+nwt
      j=i+outfiles-1
    end if
    do k=i,j
      call flush_one(its,mod(k,nwt))
    end do

! Clean up and reset for next frame.

    call clear_output_list(its)
    outfiles=0

  end subroutine flush_all

#endif /* SERIAL */

#ifndef SERIAL

!-------------------------------------------------------------------------------
  subroutine flush_one(its,wt)
!-------------------------------------------------------------------------------

! Send the metadata and field-data segments accumulated in one write-task
! linked-list output buffer to the appropriate write task. This is a private
! routine -- not part of the icosio API. It is called only by compute tasks.

    integer,intent(in)::its,wt

    character(len=100)::debughdr
    character(len=9)::this='flush_one'
    character,allocatable::m(:)
    integer::cmd(cmdsize),count,datasize,ic,metasize,p,segments
    real,allocatable::d(:)
    type(out_buffer_node),pointer::head,node

    call setdebughdr(debughdr,this,its)

    write (msg,'(2a,i0,a)') trim(debughdr),' wt=',wt,' entry.'
    call debugmsg

    head=>out_buffers(wt)%head

    segments=out_buffers(wt)%segments

    associated_head: if (associated(head)) then

      if ((wt.lt.0).or.(wt.gt.(nwt-1))) then
        write (msg,'(2a,i0,a)') trim(debughdr),' wt=',wt,' out of range.'
        call die(ierr,msg)
      end if

! Calculate the size of the metadata buffer.

      metasize=0
      node=>head
      do while (associated(node))
        metasize=metasize           &
          +size_r*1                 & ! real      :: scalefactor
          +size_i*1                 & ! integer   :: accum_start
          +size_i*1                 & ! integer   :: filename_len
          +size_i*1                 & ! integer   :: header_len
          +size_i*1                 & ! integer   :: levels
          +size_i*1                 & ! integer   :: segment_size
          +size_i*1                 & ! integer   :: time
          +size_i*1                 & ! integer   :: varname_len
          +size_c*node%filename_len & ! character :: filename
          +size_c*node%header_len   & ! character :: header
          +size_c*node%varname_len    ! character :: varname
        node=>node%next
      end do

! Allocate metadata buffer to hold worth of metadata sets, as we will send all
! metadata to this task in one contiguous buffer.

      allocate(m(metasize),stat=ierr)
      write (msg,'(2a)') trim(debughdr),' Failed to allocate metadata buffer.'
      call die_if(ierr.ne.0,ierr,msg)

      datasize=0
      ic=intercomm
      p=0

      write (msg,'(2a,i0,a)') trim(debughdr),' Packing wt ',wt,' metadata...'
      call debugmsg

! Walk the list and pack the metadata of each node into the single metadata
! buffer.

      node=>head

      metadata_pack_loop: do while (associated(node))

        write (msg,'(2a,i0,3a)') trim(debughdr),' Packing wt ',wt,' ',&
          arr_to_str(node%varname,node%varname_len),' metadata...'
        call debugmsg

        call mpi_pack(node%scalefactor,1,mpi_real,m,metasize,p,ic,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to pack scalefactor.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_pack(node%accum_start,1,mpi_integer,m,metasize,p,ic,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to pack accum_start.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_pack(node%filename_len,1,mpi_integer,m,metasize,p,ic,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to pack filename_len.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_pack(node%header_len,1,mpi_integer,m,metasize,p,ic,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to pack header_len.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_pack(node%levels,1,mpi_integer,m,metasize,p,ic,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to pack levels.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_pack(node%segment_size,1,mpi_integer,m,metasize,p,ic,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to pack segment_size.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_pack(node%time,1,mpi_integer,m,metasize,p,ic,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to pack time.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_pack(node%varname_len,1,mpi_integer,m,metasize,p,ic,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to pack varname_len.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_pack(node%filename,node%filename_len,mpi_character,m,metasize,&
          p,ic,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to pack filename.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_pack(node%header,node%header_len,mpi_character,m,metasize,p,&
          ic,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to pack header.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_pack(node%varname,node%varname_len,mpi_character,m,metasize,p,&
          ic,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to pack varname.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        write (msg,'(2a,i0,3a)') trim(debughdr),' Packed  wt ',wt,' ',&
          arr_to_str(node%varname,node%varname_len),' metadata.'
        call debugmsg

! Keep track of the storage requirement for the variable data buffer.

        datasize=datasize+node%segment_size

        node=>node%next

      end do metadata_pack_loop

      write (msg,'(2a,i0,a)') trim(debughdr),' Packed  wt ',wt,' metadata.'
      call debugmsg

! If the data segment of more than one variable will be sent to this write task,
! pack the segments into a contiguous buffer.

      if (segments.gt.1) then

        allocate(d(datasize),stat=ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to allocate data buffer.'
        call die_if(ierr.ne.0,ierr,msg)

        count=0
        p=0

        write (msg,'(2a,i0,a)') trim(debughdr),' Packing wt ',wt,&
          ' variable data...'
        call debugmsg

! Walk the list and pack the variable data of each node into the single
! variable-data buffer.

        node=>head

        vardata_pack_loop: do while (associated(node))

          write (msg,'(2a,2(i0,a),2a)') trim(debughdr),' Packing segment ',&
            count,' wt ',wt,' ',arr_to_str(node%varname,node%varname_len),&
            ' variable data...'
          call debugmsg

          call mpi_pack(node%segment,node%segment_size,mpi_real,d,&
            datasize*size_r,p,ic,ierr)
          write (msg,'(4a)') trim(debughdr),&
            ' Failed to pack data segment for ',&
            arr_to_str(node%varname,node%varname_len),'.'
          call die_if(ierr.ne.mpi_success,ierr,msg)

          write (msg,'(2a,2(i0,a),2a)') trim(debughdr),' Packed  segment ',&
            count,' wt ',wt,' ',arr_to_str(node%varname,node%varname_len),&
            ' variable data.'
          call debugmsg

          count=count+1

          node=>node%next

        end do vardata_pack_loop

        write (msg,'(2a,i0,a)') trim(debughdr),' Packed  wt ',wt,&
          ' variable data.'
        call debugmsg

      end if

! Send frame command. A command specifying a number of segments > 0 tells the
! write tasks that a frame of output metadata and variable data will follow, and
! how many segments of different variables will be sent. The current its value
! is also included. Note that only the compute root sends the frame command.

      if (i_am_compute_root) then

        cmd(1)=cmd_out
        cmd(2)=segments
        cmd(3)=its
        cmd(4)=metasize

        write (msg,'(2a,2(i0,a))') trim(debughdr),' Sending wt ',wt,&
          ' frame command for ',cmd(2),' segments...'
        call debugmsg

        call mpi_ssend(cmd,cmdsize,mpi_integer,wt,tag_cmd,intercomm,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to send frame command.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        write (msg,'(2a,2(i0,a))') trim(debughdr),' Sent    wt ',wt,&
          ' frame command for ',cmd(2),' segments.'
        call debugmsg

      end if

! Send metadata.

      write (msg,'(2a,i0,a)') trim(debughdr),' Sending wt ',wt,' metadata...'
      call debugmsg

      call mpi_ssend(m,metasize,mpi_character,wt,tag_metadata,intercomm,ierr)
      write (msg,'(2a)') trim(debughdr),' Failed to send metadata.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a,i0,a)') trim(debughdr),' Sent    wt ',wt,' metadata.'
      call debugmsg

      deallocate(m)

! Send variable data.

      write (msg,'(2a,i0,a)') trim(debughdr),' Sending wt ',wt,&
        ' variable data...'
      call debugmsg

      if (segments.gt.1) then
        call mpi_ssend(d,datasize,mpi_real,wt,tag_data,intercomm,ierr)
        deallocate(d)
      else
        call mpi_ssend(head%segment,datasize,mpi_real,wt,tag_data,intercomm,&
          ierr)
      end if

      write (msg,'(2a)') trim(debughdr),' Failed to send variable data.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a,i0,a)') trim(debughdr),' Sent    wt ',wt,' variable data.'
      call debugmsg

    end if associated_head

    write (msg,'(2a,i0,a)') trim(debughdr),' wt=',wt,' exit.'
    call debugmsg

  end subroutine flush_one

#endif /* SERIAL */

! icosio_in_* README
!
! The required format for input files is the same as for 3D model-array output
! files: two records per file, where the first record contains a header, and the
! second the actual data. The header is simply discarded during input, so its
! length and contents are irrelevant.
!
! Step 1: Stage the input requests:
!
!   call icosio_in_request(its,filename_1,size_1)
!   call icosio_in_request(its,filename_2,size_2)
!   ...
!   call icosio_in_request(its,filename_n,size_n)
!
! Integer 'its' is used to construct debug and error messages. Its value is
! arbitrary, but a timestep value is useful. The 'filename_*' strings specify
! the files to be read, and the 'size_*' integers specify the number of REAL
! values in the data record of each file.
!
! Step 2: Start asynchronous input by IO tasks:
!
!   call icosio_in_start(its)
!
! After this call, the IO tasks proceed, asynchronously, to read the specified
! files. The compute tasks can go about other business, if they have any.
!
! Step 3: Complete the inputs to set variables from read-in data:
!
!   call icosio_in_receive(its,var_1)
!   call icosio_in_receive(its,var_2)
!   ...
!   call icosio_in_receive(its,var_n)
!
! Here, 'var_1' will receive the data read from 'filename_1', 'var_2' from
! 'filename_2', etc. The receive calls must be made in the same order as the
! request calls were made. On compute tasks, these calls block until the IO
! tasks can provide the requested data.
!
! Step 4: Complete input processing:
!
!   call icosio_in_stop(its)
!
! This call closes the input frame. IO tasks then wait for new commands.

!-------------------------------------------------------------------------------
  subroutine icosio_in_receive_2dr4(its,var)
!-------------------------------------------------------------------------------

! Receive a segment of a globally-sized variable from the write tas! that was
! previously asked to read the file containing that variable. This is a
! public routine -- part of the icosio API. It is called only by compute tasks.

    integer,intent(in)::its
    real(kind=4),intent(inout)::var(:)

    character(len=22)::this='icosio_in_receive_2dr4'

#include "icosio_in_receive.decl"

    levels=1

#include "icosio_in_receive.exec"

! After execution of the main (included) body code, single-task runs will have
! an allocated "glbvar" array containing the input data, which can simply by
! transferred to the model-variable arugment. Parallel runs will each have an
! allocted "segment" array, containing the input data they own, which can
! likewise be transferred to "var".

    if (single) then
      var=glbvar
      deallocate(glbvar)
#ifndef SERIAL
    else if (allocated(segment)) then
      var=segment
      deallocate(segment)
#endif /* SERIAL */
    end if

  end subroutine icosio_in_receive_2dr4

!-------------------------------------------------------------------------------
  subroutine icosio_in_receive_3dr4(its,var)
!-------------------------------------------------------------------------------

! See documentation comment for icosio_in_receive_2dr4.

    integer,intent(in)::its
    real(kind=4),intent(inout)::var(:,:)

    character(len=22)::this='icosio_in_receive_3dr4'

#include "icosio_in_receive.decl"

    levels=size(var,1)

#include "icosio_in_receive.exec"

! For 3D model input data, reshape the 2D input-data array to conform to the
! shape of the model array variable before transferring it.

    if (single) then
      var=reshape(glbvar,(/levels,nip/))
      deallocate(glbvar)
#ifndef SERIAL
    else if (allocated(segment)) then
      var=reshape(segment,(/levels,interior_size()/))
      deallocate(segment)
#endif /* SERIAL */
    end if

  end subroutine icosio_in_receive_3dr4

!-------------------------------------------------------------------------------
  subroutine icosio_in_request(its,filename,data_size)
!-------------------------------------------------------------------------------

! Issue a request to a write task to read in a file. Requests are issued to
! write tasks in round-robin order. This is a public routine -- part of the
! icosio API. It is called only by compute tasks.

    character(len=*),intent(in)::filename
    integer,intent(in)::data_size,its

    character(len=100)::debughdr
    character(len=17)::this='icosio_in_request'

#ifndef SERIAL
    integer::cmd(cmdsize),filename_len,wt
#endif /* SERIAL */

    call setdebughdr(debughdr,this,its)

    call check_setup_called(this,its)

#ifdef SERIAL
    call ignore_unused(kind(data_size))
#endif /* SERIAL */

    if (using_write_tasks) then

#ifndef SERIAL

    filename_len=len(filename)

! Select a write task to read this file.

      wtindex=mod(wtindex,nwt)
      wt=wtindex
      wtindex=wtindex+1

! Make a note of the fact that this write task was selected for input duty, as
! well as which write task was assigned to the first requested file. This allows
! us to (a) send a 'start reading' command to all readers when icosio_in_start
! is called; and (b) match up file readers to input variables when, later, calls
! to icosio_in_receive_* are made.

      if (first_reader.lt.0) then
        first_reader=wt
        allocate(readers(0:nwt-1),stat=ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to allocate readers.'
        call die_if(ierr.ne.0,ierr,msg)
        readers=.false.
      end if
      readers(wt)=.true.

! Set up the command to send to the selected reader.

      cmd(1)=cmd_in_request
      cmd(2)=filename_len
      cmd(3)=data_size
      cmd(4)=-1 ! not used

      if (i_am_compute_root) then

! Send frame command with filename length to selected write task.

        write (msg,'(2a,i0,a)') trim(debughdr),' Sending wt ',wt,&
          ' input-request command.'
        call debugmsg

        call mpi_ssend(cmd,cmdsize,mpi_integer,wt,tag_cmd,intercomm,ierr)
        write (msg,'(2a)') trim(debughdr),&
          ' Failed to send input-request command.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        write (msg,'(2a,i0,a)') trim(debughdr),'    Sent wt ',wt,&
          ' input-request command.'
        call debugmsg

! Send filename to selected write task.

        write (msg,'(2a,i0,3a)') trim(debughdr),' Sending wt ',wt,&
          ' filename "',filename,'".'
        call debugmsg

        call mpi_ssend(str_to_arr(filename),filename_len,mpi_char,wt,&
          tag_in_filename,intercomm,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to send filename_array.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        write (msg,'(2a,i0,3a)') trim(debughdr),'    Sent wt ',wt,&
          ' filename "',filename,'".'
        call debugmsg

      end if

#endif /* SERIAL */

    else ! not using write tasks...

      if (i_am_compute_root) then

        call stringlist_add(input_filenames,trim(filename))

      end if

    end if

    in_request_called=.true.

  end subroutine icosio_in_request

!-------------------------------------------------------------------------------
  subroutine icosio_in_start(its)
!-------------------------------------------------------------------------------

! Inform write tasks that were previously issued requests to read files that
! they should now read those files. This is a public routine -- part of the
! icosio API. It is called only by compute tasks.

    integer,intent(in)::its

    character(len=100)::debughdr
    character(len=15)::this='icosio_in_start'
#ifndef SERIAL
    integer::cmd(cmdsize),wt
#endif /* SERIAL */

    call setdebughdr(debughdr,this,its)

    call check_setup_called(this,its)

    write (msg,'(4a)') trim(debughdr),&
      ' icosio_in_request must be called before ',this,'.'
    call die_if(.not.in_request_called,ierr,msg)

    if (using_write_tasks) then

#ifndef SERIAL

      cmd(1)=cmd_in_start
      cmd(2)=-1 ! not used
      cmd(3)=-1 ! not used
      cmd(4)=-1 ! not used

      if (i_am_compute_root) then

        do wt=0,nwt-1

          if (readers(wt)) then

            write (msg,'(2a,i0,a)') trim(debughdr),' Sending wt ',wt,&
              ' input-start command.'
            call debugmsg

            call mpi_ssend(cmd,cmdsize,mpi_integer,wt,tag_cmd,intercomm,ierr)
            write (msg,'(2a)') trim(debughdr),&
              ' Failed to send input-start command.'
            call die_if(ierr.ne.mpi_success,ierr,msg)

            write (msg,'(2a,i0,a)') trim(debughdr),'    Sent wt ',wt,&
              ' input-start command.'
            call debugmsg

          end if

        end do

      end if

#endif /* SERIAL */

    else ! not using write tasks...

      if (i_am_compute_root) then

        write (msg,'(2a)') trim(debughdr),' No input filenames found.'
        call die_if(.not.associated(input_filenames),ierr,msg)

      end if

    end if

    in_start_called=.true.

  end subroutine icosio_in_start

!-------------------------------------------------------------------------------
  subroutine icosio_in_stop(its)
!-------------------------------------------------------------------------------

! Inform write tasks that previously read in files that they should now dispose
! of read-in data, clean up, and prepare to respond to a new command. This is a
! public routine -- part of the icosio API. It is called only by compute tasks.

    integer,intent(in)::its

    character(len=100)::debughdr
    character(len=14)::this='icosio_in_stop'
#ifndef SERIAL
    integer::cmd(cmdsize),wt
#endif /* SERIAL */

    call setdebughdr(debughdr,this,its)

    call check_setup_called(this,its)

    write (msg,'(4a)') trim(debughdr),&
      ' icosio_in_receive must be called before ',this,'.'
    call die_if(.not.in_receive_called,ierr,msg)

    if (using_write_tasks) then

#ifndef SERIAL

      cmd(1)=cmd_in_stop
      cmd(2)=-1 ! not used
      cmd(3)=-1 ! not used
      cmd(4)=-1 ! not used

      if (i_am_compute_root) then

        do wt=0,nwt-1

          if (readers(wt)) then

            write (msg,'(2a,i0,a)') trim(debughdr),' Sending wt ',wt,&
              ' input-stop command.'
            call debugmsg

            call mpi_ssend(cmd,cmdsize,mpi_integer,wt,tag_cmd,intercomm,ierr)
            write (msg,'(2a)') trim(debughdr),&
              ' Failed to send input-stop command.'
            call die_if(ierr.ne.mpi_success,ierr,msg)

            write (msg,'(2a,i0,a)') trim(debughdr),'    Sent wt ',wt,&
              ' input-stop command.'
            call debugmsg

          end if

        end do

      end if

      deallocate(readers)
      in_receive_called=.false.
      in_request_called=.false.
      in_start_called=.false.
      first_reader=-1
      reader=-1

#endif /* SERIAL */

    else ! not using write tasks...

      if (i_am_compute_root) then

        if (associated(input_filenames)) call stringlist_destroy(input_filenames)
        nullify(current_filename)

      end if

    end if

  end subroutine icosio_in_stop

!-------------------------------------------------------------------------------
  subroutine icosio_out_2dr4(its,time,varname,var,filename,header,scalefactor,&
    accum_start)
!-------------------------------------------------------------------------------

! A wrapper around icosio_out_3dr4 for shape conformance.

    character(len=*),intent(in)::filename,header,varname
    integer,intent(in),optional::accum_start
    integer,intent(in)::its,time
    real(kind=4),intent(in),optional::scalefactor
    real(kind=4),intent(in)::var(:)

    call icosio_out_3dr4(its,time,varname,reshape(var,(/1,size(var,1)/)),&
      filename,header,scalefactor,accum_start)

  end subroutine icosio_out_2dr4

!-------------------------------------------------------------------------------
  subroutine icosio_out_2dr8(its,time,varname,var,filename,header,scalefactor,&
    accum_start)
!-------------------------------------------------------------------------------

! Convert REAL*8 var to REAL*4 and pass on to icosio_out_2dr4.

    character(len=*),intent(in)::filename,header,varname
    integer,intent(in),optional::accum_start
    integer,intent(in)::its,time
    real(kind=8),intent(in),optional::scalefactor
    real(kind=8),intent(in)::var(:)

    real(kind=4)::scalefactor_r4

    scalefactor_r4=1.
    if (present(scalefactor)) scalefactor_r4=scalefactor

    call icosio_out_2dr4(its,time,varname,real(var,4),filename,header,&
      scalefactor_r4,accum_start)

  end subroutine icosio_out_2dr8

!-------------------------------------------------------------------------------
  subroutine icosio_out_3dr4(its,time,varname,var,filename,header,scalefactor,&
    accum_start)
!-------------------------------------------------------------------------------

! Receive a field-data array and either write it to disk directly, or buffer it
! for a later transmisstion to write task(s). This is a public routine -- part
! of the icosio API. It is called only by compute tasks.

    character(len=*),intent(in)::filename,header,varname
    integer,intent(in),optional::accum_start
    integer,intent(in)::its,time
    real(kind=4),intent(in),optional::scalefactor
    real(kind=4),intent(in),target::var(:,:)

    character(len=100)::debughdr
    character(len=15)::this='icosio_out_3dr4'
    integer::accum_start_local
    real(kind=4),pointer::glbvar(:,:)=>null()
    real(kind=4)::scalefactor_local

    call setdebughdr(debughdr,this,its)

    call check_setup_called(this,its)

! Set default values for scalefactor (factor to scale input by for grib output)
! and accum_start (accumulation start) if not provided by caller.

    scalefactor_local=1.
    if (present(scalefactor)) scalefactor_local=scalefactor
    accum_start_local=-1
    if (present(accum_start)) accum_start_local=accum_start

! If write tasks are in use, buffer data for eventual overlapped write;
! otherwise, write data directly.

    if (using_write_tasks) then

      call buffer_var(its,varname,var,filename,header,time,scalefactor_local,&
        accum_start_local)

    else ! when not using write tasks...

! Collect inv_perm on the compute root if fixed grid-order reordering is needed.

      if (i_am_compute_task) then
        if (permute) then
          if (.not.associated(inv_perm_global)) then
            call collect_inv_perm(its)
          end if
        end if
      end if

! If no reordering is required in single mode, save memory by writing from the
! passed-in array directly. If reordering *is* required in single mode, make a
! copy of the passed-in array first. In parallel mode, collect the distributed
! array onto the compute root.

      if (single) then
        if (permute) then
          allocate(glbvar(size(var,1),size(var,2)),stat=ierr)
          write (msg,'(2a)') trim(debughdr),' Failed to allocate glbvar.'
          call die_if(ierr.ne.0,ierr,msg)
          glbvar=var ! copy
        else
          glbvar=>var ! point
        end if
#ifndef SERIAL
      else
        call collect_var(its,var,glbvar)
#endif /* SERIAL */
      end if

      if (single.or.i_am_compute_root) then
#ifndef NOGRIB
        call var_to_disk(    &
          filename,header,   &
          its,glbvar,        &
          varname,           &
          accum_start_local, &
          scalefactor_local  &
          )
#else
        call var_to_disk(    &
          filename,          &
          header,            &
          its,               &
          glbvar,            &
          varname            &
          )
#endif /* NOGRIB */
      end if

      if ((single.and.permute).or.(.not.single.and.i_am_compute_root)) then
        deallocate(glbvar)
      end if
      if (associated(glbvar)) then
        nullify(glbvar)
      end if

    end if

  end subroutine icosio_out_3dr4

!-------------------------------------------------------------------------------
  subroutine icosio_out_3dr8(its,time,varname,var,filename,header,scalefactor,&
    accum_start)
!-------------------------------------------------------------------------------

! Convert REAL*8 var to REAL*4 and pass on to icosio_out_3dr4.

    character(len=*),intent(in)::filename,header,varname
    integer,intent(in),optional::accum_start
    integer,intent(in)::its,time
    real(kind=8),intent(in),optional::scalefactor
    real(kind=8),intent(in)::var(:,:)

    real(kind=4)::scalefactor_r4

    scalefactor_r4=1.
    if (present(scalefactor)) scalefactor_r4=scalefactor

    call icosio_out_3dr4(its,time,varname,real(var,4),filename,header,&
      scalefactor_r4,accum_start)

  end subroutine icosio_out_3dr8

!-------------------------------------------------------------------------------
  subroutine icosio_out_stop(its)
!-------------------------------------------------------------------------------

! Reset the list of filenames written to during this output interval. If write
! tasks are not in use, simply return; otherwise, call flush_all(). This is a
! public routine -- part of the icosio API. It is called only by compute tasks.

    integer,intent(in)::its

    character(len=100)::debughdr
    character(len=15)::this='icosio_out_stop'

    call setdebughdr(debughdr,this,its)

    call check_setup_called(this,its)

    if (associated(output_filenames)) call stringlist_destroy(output_filenames)

    if (nwt.eq.0) return

#ifndef SERIAL
    call flush_all(its)
#endif /* SERIAL */

  end subroutine icosio_out_stop

#ifndef SERIAL

!-------------------------------------------------------------------------------
  subroutine icosio_run
!-------------------------------------------------------------------------------

! Main routine for write tasks. Allocate dynamic buffers, receive some one-time
! data, then enter a loop waiting for and responding to frame commands, which
! trigger input, output or shutdown actions. This is a public routine -- part of
! the icosio API. It is called only by write tasks.

    character(len=100)::debughdr
    character(len=10)::this='icosio_run'
    integer::cmd(cmdsize),ct,data_size,expected,filename_len,&
      i_params(num_i_params),ipn,is(1),lb,levels,mpitype,&
      mpistatus(mpi_status_size),offset,received,segment_size,tag,ub
    logical::l_params(num_l_params),receive_active=.false.,running=.true.
    type(in_buffer_node),pointer::next,node

    call setdebughdr(debughdr,this,0)

    write (msg,'(2a)') trim(debughdr),' Entry.'
    call debugmsg

    call check_setup_called(this,0)

! Get sizes of intrinsic types for communication.

    call comm_type_sizes

! Receive integer configuration parameters from compute root, broadcast to all
! write tasks, then unpack.

    if (i_am_write_root) then

      write (msg,'(2a)') trim(debughdr),&
        ' Receiving integer configuration parameters from compute root...'
      call debugmsg

      call mpi_recv(i_params,num_i_params,mpi_integer,0,tag_integer_parameters,&
        intercomm,mpistatus,ierr)
      write (msg,'(2a)') trim(debughdr),&
        ' Write root failed to receive integer configuration parameters.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a)') trim(debughdr),&
        ' Received  integer configuration parameters from compute root.'
      call debugmsg

    end if

    write (msg,'(2a)') trim(debughdr),&
      ' Integer configuration parameters broadcast: start.'
    call debugmsg

    call mpi_bcast(i_params,num_i_params,mpi_integer,0,intracomm,ierr)
    write (msg,'(2a)') trim(debughdr),&
      ' Failed to broadcast integer configuration parameters.'
    call die_if(ierr.ne.mpi_success,ierr,msg)

    write (msg,'(2a)') trim(debughdr),&
      ' Integer configuration parameters broadcast: done.'
    call debugmsg

    lun = i_params(1)
    nip = i_params(2)

! Receive logical configuration parameters from compute root, broadcast to all
! write tasks, then unpack.

    if (i_am_write_root) then

      write (msg,'(2a)') trim(debughdr),&
        ' Receiving logical configuration parameters from compute root...'
      call debugmsg

      call mpi_recv(l_params,num_l_params,mpi_logical,0,tag_logical_parameters,&
        intercomm,mpistatus,ierr)
      write (msg,'(2a)') trim(debughdr),&
        ' Write root failed to receive logical configuration parameters.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a)') trim(debughdr),&
        ' Received  logical configuration parameters from compute root.'
      call debugmsg

    end if

    write (msg,'(2a)') trim(debughdr),&
      ' Logical configuration parameters broadcast: start.'
    call debugmsg

    call mpi_bcast(l_params,num_l_params,mpi_logical,0,intracomm,ierr)
    write (msg,'(2a)') trim(debughdr),&
      ' Failed to broadcast logical configuration parameters.'
    call die_if(ierr.ne.mpi_success,ierr,msg)

    write (msg,'(2a)') trim(debughdr),&
      ' Logical configuration parameters broadcast: done.'
    call debugmsg

    binout      = l_params(1)
    gribout     = l_params(2)
    permute     = l_params(3)
    print_diags = l_params(4)

! Report configuration parameters.

    if (debugmsg_on) call show_config

! If the model is running in an "only print diagnostics" mode, either stop or
! return, depending on the client/server io settting.

    if (print_diags) then
      if (client_server_io) then
        call mpi_finalize(ierr)
        stop
      else
        return
      end if
    end if

! Allocate input buffers.

    allocate(in_buffers)
    write (msg,'(2a)') trim(debughdr),' Failed to allocate input buffers.'
    call die_if(ierr.ne.0,ierr,msg)

! Allocate output buffers.

    allocate(out_buffers(0:nct-1),stat=ierr)
    write (msg,'(2a)') trim(debughdr),' Failed to allocate output buffers.'
    call die_if(ierr.ne.0,ierr,msg)

! Allocate and receive interior sizes from compute task(s). The root write
! task collects the interior sizes and populates its buffer, then broadcasts to
! the other write tasks, if any.

    allocate(interior_sizes(0:nct-1),stat=ierr)
    write (msg,'(2a)') trim(debughdr),' Failed to allocate interior_sizes.'
    call die_if(ierr.ne.0,ierr,msg)

    if (i_am_write_root) then

      do ct=0,nct-1

        write (msg,'(2a,i0,a)') trim(debughdr),&
          ' Receiving interior size from ct ',ct,'...'
        call debugmsg

        expected=1
        mpitype=mpi_integer
        tag=tag_interior_sizes

        call mpi_recv(is,expected,mpitype,ct,tag,intercomm,mpistatus,ierr)
        write (msg,'(2a,i0,a)') trim(debughdr),&
          ' Failed to receive interior size from ',ct,'.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_get_count(mpistatus,mpitype,received,ierr)
        write (msg,'(2a,2(i0,a,i0,a))') trim(debughdr),' is expected ',&
          expected,' elements, received ',received,'.'
        call die_if(received.ne.expected,ierr,msg)

        interior_sizes(ct)=is(1)

        write (msg,'(2a,2(i0,a))') trim(debughdr),&
          ' Received  interior size from ct ',ct,': ',interior_sizes(ct),'.'
        call debugmsg

      end do

    end if

    if (nwt.gt.1) then

      write (msg,'(2a)') trim(debughdr),' Interior sizes broadcast entry.'
      call debugmsg

      call mpi_bcast(interior_sizes,nct,mpi_integer,0,intracomm,ierr)
      write (msg,'(2a)') trim(debughdr),' Failed to broadcast interior sizes.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a)') trim(debughdr),' Interior sizes broadcast exit.'
      call debugmsg

    end if

! Allocate and receive from compute task(s) inv_perm and perm, if needed. The
! compute root sends a control flag to the root write task (which broadcasts to
! all other write tasks, if any) to say whether or not inv_perm and perm will be
! sent.

    if (i_am_write_root) then

      write (msg,'(2a)') trim(debughdr),&
        ' Receiving perm control flag from compute root...'
      call debugmsg

      call mpi_recv(permute,1,mpi_logical,0,tag_permute_control,intercomm,&
        mpistatus,ierr)
      write (msg,'(2a)') trim(debughdr),&
        ' Write root failed to receive permute flag.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a)') trim(debughdr),&
        ' Received  perm control flag from compute root.'
      call debugmsg

    end if

    if (nwt.gt.1) then

      write (msg,'(2a)') trim(debughdr),' perm control flag broadcast entry.'
      call debugmsg

      call mpi_bcast(permute,1,mpi_logical,0,intracomm,ierr)
      write (msg,'(2a)') trim(debughdr),' Failed to broadcast permute flag.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a)') trim(debughdr),' perm control flag broadcast exit.'
      call debugmsg

    end if

    if_permute: if (permute) then ! we ARE receiving [inv_]perm...

      allocate(inv_perm_global(nip),stat=ierr)
      write (msg,'(2a)') trim(debughdr),' Failed to allocate inv_perm_global.'
      call die_if(ierr.ne.0,ierr,msg)

      allocate(perm_global(nip),stat=ierr)
      write (msg,'(2a)') trim(debughdr),' Failed to allocate perm_global.'
      call die_if(ierr.ne.0,ierr,msg)

      write_root: if (i_am_write_root) then

        offset=1

        do ct=0,nct-1

          expected=interior_sizes(ct)
          mpitype=mpi_integer

! inv_perm

          tag=tag_inv_perm_data

          write (msg,'(2a,i0,a)') trim(debughdr),&
            ' Receiving inv_perm segment from ct ',ct,'...'
          call debugmsg

          call mpi_recv(inv_perm_global(offset:offset-1+interior_sizes(ct)),&
            expected,mpitype,ct,tag,intercomm,mpistatus,ierr)
          write (msg,'(2a,i0,a)') trim(debughdr),&
            ' Failed to receive inv_perm_global segment from ',ct,'.'
          call die_if(ierr.ne.mpi_success,ierr,msg)

          call mpi_get_count(mpistatus,mpitype,received,ierr)
          write (msg,'(2a,2(i0,a))') trim(debughdr),&
            ' inv_perm_global expected ',expected,' elements, received ',&
            received,'.'
          call die_if(received.ne.expected,ierr,msg)

          write (msg,'(2a,i0,a)') trim(debughdr),&
            ' Received  inv_perm segment from ct ',ct,'.'
          call debugmsg

! perm

          tag=tag_perm_data

          write (msg,'(2a,i0,a)') trim(debughdr),&
            ' Receiving perm segment from ct ',ct,'...'
          call debugmsg

          call mpi_recv(perm_global(offset:offset-1+interior_sizes(ct)),&
            expected,mpitype,ct,tag,intercomm,mpistatus,ierr)
          write (msg,'(2a,i0,a)') trim(debughdr),&
            ' Failed to receive perm_global segment from ',ct,'.'
          call die_if(ierr.ne.mpi_success,ierr,msg)

          call mpi_get_count(mpistatus,mpitype,received,ierr)
          write (msg,'(2a,2(i0,a))') trim(debughdr),' perm_global expected ',&
            expected,' elements, received ',received,'.'
          call die_if(received.ne.expected,ierr,msg)

          write (msg,'(2a,i0,a)') trim(debughdr),&
            ' Received  perm segment from ct ',ct,'.'
          call debugmsg

          offset=offset+interior_sizes(ct)

        end do

      end if write_root

      if (nwt.gt.1) then

! Broadcast inv_perm_global from write root to other write tasks.

        write (msg,'(2a)') trim(debughdr),' inv_perm_global broadcast entry.'
        call debugmsg

        call mpi_bcast(inv_perm_global,nip,mpi_integer,0,intracomm,ierr)
        write (msg,'(2a)') trim(debughdr),&
          ' Failed to broadcast inv_perm_global.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        write (msg,'(2a)') trim(debughdr),' inv_perm_global broadcast exit.'
        call debugmsg

! Broadcast perm_global from write root to other write tasks.

        write (msg,'(2a)') trim(debughdr),' perm_global broadcast entry.'
        call debugmsg

        call mpi_bcast(perm_global,nip,mpi_integer,0,intracomm,ierr)
        write (msg,'(2a)') trim(debughdr),&
          ' Failed to broadcast perm_global.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        write (msg,'(2a)') trim(debughdr),' perm_global broadcast exit.'
        call debugmsg

      end if

    end if if_permute

! Main loop

    do while (running)

! Receve a frame command from the compute root.

      write (msg,'(2a)') trim(debughdr),&
        ' Waiting for frame command from compute root.'
      call debugmsg

      expected=cmdsize
      mpitype=mpi_integer
      tag=tag_cmd

      call mpi_recv(cmd,expected,mpitype,0,tag,intercomm,mpistatus,ierr)
      write (msg,'(2a)') trim(debughdr),&
        ' Failed to receive frame command from compute root.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      call mpi_get_count(mpistatus,mpitype,received,ierr)
      write (msg,'(2a,2(i0,a))') trim(debughdr),' cmd expected ',expected,&
        ' elements, received ',received,'.'
      call die_if(received.ne.expected,ierr,msg)

! Handle frame commands.

      if (cmd(1).eq.cmd_in_receive) then

! Receive input.

        write (msg,'(2a)') trim(debughdr),' Received input-receive command.'
        call debugmsg

        if (receive_active) then
          node=>node%next
        else
          node=>in_buffers%head
          receive_active=.true.
        end if

        lb=1

        do ct=0,nct-1

          levels=node%data_size/nip

          segment_size=levels*interior_sizes(ct)

          ub=lb+segment_size-1

          write (msg,'(4a,3(i0,a))') trim(debughdr),' Sending "',&
            arr_to_str(node%filename,node%filename_len),'" elements ',lb,'-',&
            ub,' to ct ',ct,'.'
          call debugmsg

          call mpi_ssend(node%data(lb:ub),segment_size,mpi_real,ct,&
            tag_in_receive,intercomm,ierr)
          write (msg,'(2a,i0,a)') trim(debughdr),&
            ' Failed to send variable data to ct ',ct,'.'
          call die_if(ierr.ne.mpi_success,ierr,msg)

          write (msg,'(4a,3(i0,a))') trim(debughdr),'    Sent "',&
            arr_to_str(node%filename,node%filename_len),'" elements ',lb,'-',&
            ub,' to ct ',ct,'.'
          call debugmsg

          lb=ub+1

        end do

      else if (cmd(1).eq.cmd_in_request) then

! Request input.

        write (msg,'(2a)') trim(debughdr),' Received input-request command.'
        call debugmsg

        filename_len=cmd(2)
        data_size=cmd(3)

        if (.not.associated(in_buffers%head)) then
          allocate(in_buffers%head,stat=ierr)
          write (msg,'(2a)') trim(debughdr),&
            ' Failed to allocate in_buffers%head.'
          call die_if(ierr.ne.0,ierr,msg)
          in_buffers%tail=>in_buffers%head
        else
          allocate(in_buffers%tail%next,stat=ierr)
          write (msg,'(2a)') trim(debughdr),&
            ' Failed to allocate in_buffers%tail%next.'
          call die_if(ierr.ne.0,ierr,msg)
          in_buffers%tail=>in_buffers%tail%next
        end if

        allocate(in_buffers%tail%filename(filename_len),stat=ierr)
        write (msg,'(2a)') trim(debughdr),&
          ' Failed to allocate in_buffers%tail%filename.'
        call die_if(ierr.ne.0,ierr,msg)

        in_buffers%tail%filename_len=filename_len
        in_buffers%tail%data_size=data_size

        expected=filename_len
        mpitype=mpi_character
        tag=tag_in_filename

        write (msg,'(2a)') trim(debughdr),&
          ' Receiving input filename from compute root.'
        call debugmsg

        call mpi_recv(in_buffers%tail%filename,expected,mpitype,0,tag,&
          intercomm,mpistatus,ierr)
        write (msg,'(2a)') trim(debughdr),&
          ' Failed to receive input filename from compute root.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_get_count(mpistatus,mpitype,received,ierr)
        write (msg,'(2a,2(i0,a))') trim(debughdr),' filename expected ',&
          expected,' elements, received ',received,'.'
        call die_if(received.ne.expected,ierr,msg)

        write (msg,'(4a)') trim(debughdr),&
          '  Received input filename from compute root: ',&
          arr_to_str(in_buffers%tail%filename,in_buffers%tail%filename_len),'.'
        call debugmsg

      else if (cmd(1).eq.cmd_in_start) then

! Start input.

        write (msg,'(2a)') trim(debughdr),' Received input-start command.'
        call debugmsg

        node=>in_buffers%head

        do while (associated(node))

          write (msg,'(4a)') trim(debughdr),' Reading ',&
            arr_to_str(node%filename,node%filename_len),'.'
          call debugmsg

          open (lun,file=arr_to_str(node%filename,node%filename_len),&
            form='unformatted',status='old',iostat=ierr)
          write (msg,'(4a)') trim(debughdr),' Failed to open file ',&
            arr_to_str(node%filename,node%filename_len),'.'
          call die_if(ierr.ne.0,ierr,msg)

          read (lun,iostat=ierr)
          write (msg,'(2a)') trim(debughdr),' Failed to read header.'
          call die_if(ierr.ne.0,ierr,msg)

          allocate(node%data(node%data_size),stat=ierr)
          write (msg,'(2a)') trim(debughdr),' Failed to allocate node%data.'
          call die_if(ierr.ne.0,ierr,msg)

          read (lun,iostat=ierr) node%data
          write (msg,'(2a)') trim(debughdr),' Failed read to data.'
          call die_if(ierr.ne.0,ierr,msg)

          levels=node%data_size/nip

          call permute_input(node%data,levels)

          close (lun)

          write (msg,'(4a)') trim(debughdr),'    Read ',&
            arr_to_str(node%filename,node%filename_len),'.'
          call debugmsg

          node=>node%next

        end do

        write (msg,'(2a)') trim(debughdr),' Finished reading files.'
        call debugmsg

      else if (cmd(1).eq.cmd_in_stop) then

! Shut down.

        write (msg,'(2a)') trim(debughdr),' Received stop command.'
        call debugmsg

        node=>in_buffers%head
        do while (associated(node))
          if (associated(node%filename)) deallocate(node%filename)
          if (associated(node%data)) deallocate(node%data)
          next=>node%next
          deallocate(node)
          node=>next
        end do
        nullify(in_buffers%head)
        nullify(in_buffers%tail)
        nullify(node)

        receive_active=.false.

      else if (cmd(1).eq.cmd_out) then

! Output.

        write (msg,'(2a)') trim(debughdr),' Received output command.'
        call debugmsg

        call die_if(cmd(2).lt.0,ierr,'Frame command segment count < 0.')

        call setdebughdr(debughdr,this,cmd(3))

        write (msg,'(2a,i0,a)') trim(debughdr),&
          ' Received frame command from compute root for ',cmd(2),' segments.'
        call debugmsg

        call process_output_frame(cmd(2),cmd(3),cmd(4))

      else if (cmd(1).eq.cmd_stop) then

        write (msg,'(2a)') trim(debughdr),' Stop command received.'
        call debugmsg

        running=.false.

        if (client_server_io) then

          write (msg,'(2a)') trim(debughdr),&
            ' Client-server mode enabled, stopping...'
          call debugmsg

          call mpi_finalize(ierr)
          stop

        else

          return

        end if

      else

        write (msg,'(2a,i0,a)') trim(debughdr),' Unrecognized command value: ',&
          cmd(1),'.'
        call die(ierr,msg)

      end if

    end do

  end subroutine icosio_run

#endif /* SERIAL */

!-------------------------------------------------------------------------------
  subroutine icosio_setup(binout_in,client_server_io_in,comm_in,debugmsg_on_in,&
    gribout_in,i_am_write_task_in,ips_in,ipe_in,inv_perm_in,lun_in,perm_in,&
    permute_in,print_diags_in,using_write_tasks_in)
!-------------------------------------------------------------------------------

! Sets values known by the calling model and required by icosio. This is a
! public routine -- part of the icosio API. It is called by both compute and
! write tasks.

! Required arguments

    integer,intent(in)::comm_in
    logical,intent(in)::client_server_io_in,debugmsg_on_in,i_am_write_task_in

! Optional arguments (see below for default values)

    integer,intent(in),optional::ips_in,ipe_in,lun_in
    integer,intent(in),allocatable,optional,target::inv_perm_in(:)
    integer,intent(in),allocatable,optional,target::perm_in(:)
    logical,intent(in),optional::binout_in,gribout_in,&
      permute_in,print_diags_in,using_write_tasks_in

! Locals

    character(len=100)::debughdr
    character(len=12)::this='icosio_setup'
    integer::key,nipsum

#ifndef SERIAL
    integer::allcomm
    logical::have_intercomm,high
#endif /* SERIAL */

    call setdebughdr(debughdr,this,0)

#ifdef SERIAL
    call ignore_unused(kind(comm_in))
    call ignore_unused(kind(key))
    call ignore_unused(kind(nipsum))
#endif /* SERIAL */

! Set module variables.
!
! A number of necessary parameters are deduced from three settings: comm,
! i_am_write_task, and using_write_tasks. If write tasks are in use, comm must
! be an intercommunicator between the compute and write tasks; otherwise, it
! must be an intracommunicator for the compute tasks. The two logical values
! have the obvious meanings.
!
! If client/server IO is disabled, write tasks return from the icosio_setup()
! call when they receive a shutdown command from the compute root (i.e. after
! the last output frame has been processed. Otherwise, they call mpi_finalize().

    client_server_io=client_server_io_in
    debugmsg_on=debugmsg_on_in
    i_am_write_task=i_am_write_task_in

#ifndef SERIAL
    comm=comm_in
#endif /* SERIAL */

! Override optional-argument defaults, when possible.

    if (present(binout_in)) binout=binout_in
    if (present(gribout_in)) gribout=gribout_in
    if (present(lun_in)) lun=lun_in
    if (present(print_diags_in)) print_diags=print_diags_in

    if (present(ipe_in)) then
      ipe=ipe_in
    else
      write (msg,'(2a)') trim(debughdr),&
        ' Argument "ipe_in" is mandatory for compute tasks.'
      call die_if(.not.i_am_write_task,ierr,msg)
    end if

    if (present(ips_in)) then
      ips=ips_in
    else
      write (msg,'(2a)') trim(debughdr),&
        ' Argument "ips_in" is mandatory for compute tasks.'
      call die_if(.not.i_am_write_task,ierr,msg)
    end if

    if (present(permute_in)) then
      permute=permute_in
    else
      write (msg,'(2a)') trim(debughdr),&
        ' Argument "permute_in" is mandatory for compute tasks.'
      call die_if(.not.i_am_write_task,ierr,msg)
    end if

    if (present(using_write_tasks_in)) then
      using_write_tasks=using_write_tasks_in
    else
      write (msg,'(2a)') trim(debughdr),&
        ' Argument "using_write_tasks_in" is mandatory for compute tasks.'
      call die_if(.not.i_am_write_task,ierr,msg)
    end if

    if (present(inv_perm_in)) inv_perm=>inv_perm_in
    write (msg,'(2a)') trim(debughdr),&
      ' "inv_perm_in" required when "permute_in" is true.'
    call die_if(permute.and..not.associated(inv_perm),ierr,msg)

    if (present(perm_in)) perm=>perm_in
    write (msg,'(2a)') trim(debughdr),&
      ' "perm_in" required when "permute_in" is true.'
    call die_if(permute.and..not.associated(perm),ierr,msg)

#ifdef SERIAL

    perm_global => perm_in

#else

    intercomm=mpi_comm_null ! default value
    intracomm=mpi_comm_null ! default value

! Override defaults.

    if (i_am_write_task) then
      tasktype='wt'
      using_write_tasks=.true.
    end if

    if (using_write_tasks) then
      call mpi_comm_test_inter(comm,have_intercomm,ierr)
      write (msg,'(2a)') trim(debughdr),&
        ' Pass icosio an intercommunicator when write task(s) enabled.'
      call die_if(.not.have_intercomm,ierr,msg)
      call mpi_comm_rank(comm,me,ierr)
      if (i_am_write_task) then
        high=.true.
        key=1
      else
        high=.false.
        key=0
      end if
      call mpi_intercomm_merge(comm,high,allcomm,ierr)
      call mpi_comm_split(allcomm,key,me,intracomm,ierr)
      intercomm=comm
    else if (comm.ne.mpi_comm_null) then
      call mpi_comm_test_inter(comm,have_intercomm,ierr)
      write (msg,'(2a)') trim(debughdr),&
        ' Pass icosio an intracommunicator when write task(s) disabled.'
      call die_if(have_intercomm,ierr,msg)
      call mpi_comm_rank(comm,me,ierr)
      intracomm=comm
    end if

    if (intracomm.ne.mpi_comm_null) then
      call mpi_comm_set_errhandler(intracomm,mpi_errors_are_fatal,ierr)
      write (msg,'(2a,i0,a)') trim(debughdr),&
        ' MPI_Comm_set_errhandler for intracomm returned ',ierr,'.'
      call die_if(ierr.ne.mpi_success,ierr,msg)
    end if

    if (intercomm.ne.mpi_comm_null) then
      call mpi_comm_set_errhandler(intercomm,mpi_errors_are_fatal,ierr)
      write (msg,'(2a,i0,a)') trim(debughdr),&
        ' MPI_Comm_set_errhandler for intercomm returned ',ierr,'.'
      call die_if(ierr.ne.mpi_success,ierr,msg)
    end if

    if (i_am_write_task) then
      call mpi_comm_remote_size(intercomm,nct,ierr)
      call mpi_comm_size(intracomm,nwt,ierr)
      i_am_compute_root=.false.
      i_am_compute_task=.false.
      serial=.false.
      if (nwt.gt.1) then
        single=.false.
      end if
      if (me.eq.0) then
        i_am_write_root=.true.
      end if
    else
      if (intracomm.ne.mpi_comm_null) then
        call mpi_comm_size(intracomm,nct,ierr)
        serial=.false.
        if (nct.ne.1) then
          single=.false.
        end if
        if (me.ne.0) then
          i_am_compute_root=.false.
        end if
      end if
      if (intercomm.ne.mpi_comm_null) then
        call mpi_comm_remote_size(intercomm,nwt,ierr)
      end if
    end if

#endif /* SERIAL */

! Deduce nip from interior sizes of compute tasks.

    if (i_am_compute_task) then
      nip=interior_size()
#ifndef SERIAL
      if (.not.single) then
        call mpi_allreduce(nip,nipsum,1,mpi_integer,mpi_sum,intracomm,ierr)
        write (msg,'(2a,i0,a)') trim(debughdr),&
          ' MPI_Allreduce on nip returned ',ierr,'.'
        call die_if(ierr.ne.mpi_success,ierr,msg)
        nip=nipsum
      end if
#endif /* SERIAL */
    end if

! Report configuration parameters.

    if (i_am_compute_task.and.debugmsg_on) call show_config

! Record that icosio_setup has run.

    icosio_setup_called=.true.

#ifndef SERIAL

! Compute tasks should prepare write tasks.

    if (using_write_tasks.and.i_am_compute_task) call prep_write_tasks

! If I am a write task, enter client-server mode if selected.

    if (i_am_write_task.and.client_server_io) call icosio_run

#endif /* SERIAL */

  end subroutine icosio_setup

!-------------------------------------------------------------------------------
  subroutine icosio_stop(its)
!-------------------------------------------------------------------------------

! If write tasks are not in use, simply return; otherwise, tell write tasks to
! shut down. This is a public routine -- part of the icosio API. It is called
! only by compute tasks.

    integer,intent(in)::its

#ifndef SERIAL
    character(len=100)::debughdr
    character(len=11)::this='icosio_stop'
    integer::cmd(cmdsize),wt
#endif /* SERIAL */

#ifdef SERIAL

    call ignore_unused(kind(its))

#else

    call setdebughdr(debughdr,this,its)

    call check_setup_called(this,its)

    if (nwt.eq.0) return

    if (i_am_compute_root) then
      write (msg,'(2a)') trim(debughdr),' Sending stop commands...'
      call debugmsg
      cmd(1)=cmd_stop
      cmd(2)=0
      cmd(3)=its
      cmd(4)=0
      do wt=0,nwt-1
        write (msg,'(2a,i0,a)') trim(debughdr),&
          ' Sending stop command to wt ',wt,'...'
        call debugmsg
        call mpi_ssend(cmd,cmdsize,mpi_integer,wt,tag_cmd,intercomm,ierr)
        write (msg,'(2a,i0,a)') trim(debughdr),' Sent    stop command to wt ',&
          wt,'.'
        call debugmsg
      end do
      write (msg,'(2a)') trim(debughdr),' Sent stop commands.'
      call debugmsg
    end if

#endif /* SERIAL */

  end subroutine icosio_stop

!-------------------------------------------------------------------------------
  subroutine permute_input(var,levels)
!-------------------------------------------------------------------------------

! Perform in situ permutation of input data to match model curve. This is a
! private routine -- not part of the icosio API. It is called by both compute
! and write tasks.

    integer,intent(in)::levels
    real,intent(inout)::var(:)

    integer::fetch_from,ipn
    real,allocatable::column(:)

    allocate(column(levels))
    do ipn=1,nip-1
      fetch_from=perm_global(ipn)
      do while (fetch_from.lt.ipn)
        fetch_from=perm_global(fetch_from)
      end do
      if (fetch_from.ne.ipn) then
        column=var((ipn-1)*levels+1:ipn*levels)
        var((ipn-1)*levels+1:ipn*levels)=&
          var((fetch_from-1)*levels+1:fetch_from*levels)
        var((fetch_from-1)*levels+1:fetch_from*levels)=column
      end if
    end do
    deallocate(column)

  end subroutine permute_input

#ifndef SERIAL

!-------------------------------------------------------------------------------
  subroutine prep_write_tasks
!-------------------------------------------------------------------------------

! Send some one-time initialization data to write tasks. Allocate some arrays
! needed for communication with write tasks. This is a private routine -- not
! part of the icosio API. It is called only by compute tasks.

    character(len=100)::debughdr
    character(len=16)::this='prep_write_tasks'
    integer::i_params(num_i_params)
    logical::l_params(num_l_params)

    call setdebughdr(debughdr,this,0)

    write (msg,'(2a)') trim(debughdr),' Entry.'
    call debugmsg

! Get sizes of intrinsic types for communication.

    call comm_type_sizes

! Allocate input buffers.

    allocate(in_buffers)
    write (msg,'(2a)') trim(debughdr),' Failed to allocate input buffers.'
    call die_if(ierr.ne.0,ierr,msg)

! Allocate output buffers.

    allocate(out_buffers(0:nwt-1),stat=ierr)
    write (msg,'(2a)') trim(debughdr),' Failed to allocate output buffers.'
    call die_if(ierr.ne.0,ierr,msg)

! Pack and send integer configuration parameters to write root.

    i_params=-1 ! default value

    if (i_am_compute_root) then

      i_params(1)=lun
      i_params(2)=nip

      write (msg,'(2a)') trim(debughdr),&
        ' Sending integer configuration parameters to write root...'
      call debugmsg

      call mpi_ssend(i_params,num_i_params,mpi_integer,0,&
        tag_integer_parameters,intercomm,ierr)
      write (msg,'(2a)') trim(debughdr),&
        ' Failed to send integer configuration parameters.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a)') trim(debughdr),&
        ' Sent    integer configuration parameters to write root.'
      call debugmsg

    end if

! Pack and send logical configuration parameters to write root.

    l_params=.false. ! default value

    if (i_am_compute_root) then

      l_params(1)=binout
      l_params(2)=gribout
      l_params(3)=permute
      l_params(4)=print_diags

      write (msg,'(2a)') trim(debughdr),&
        ' Sending logical configuration parameters to write root...'
      call debugmsg

      call mpi_ssend(l_params,num_l_params,mpi_logical,0,&
        tag_logical_parameters,intercomm,ierr)
      write (msg,'(2a)') trim(debughdr),&
        ' Failed to send logical configuration parameters.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a)') trim(debughdr),&
        ' Sent    logical configuration parameters to write root.'
      call debugmsg

    end if

! Calculate interior points owned by this compute task and send to the root
! write task, which will broadcast to all other write tasks, if any.

    write (msg,'(2a,i0,a)') trim(debughdr),' Sending interior size ',&
      interior_size(),' to write root...'
    call debugmsg

    call mpi_ssend(interior_size(),1,mpi_integer,0,tag_interior_sizes,&
      intercomm,ierr)
    write (msg,'(2a)') trim(debughdr),' Failed to send interior sizes.'
    call die_if(ierr.ne.mpi_success,ierr,msg)

    write (msg,'(2a,i0,a)') trim(debughdr),' Sent    interior size ',&
      interior_size(),' to write root.'
    call debugmsg

! If permute is true, send inv_perm to the write task(s). The compute root
! informs the write task(s) whether or not inv_perm will be sent.

    if (i_am_compute_root) then

      write (msg,'(2a)') trim(debughdr),&
        ' Sending inv_perm control flag to write root...'
      call debugmsg

      call mpi_ssend(permute,1,mpi_logical,0,tag_permute_control,intercomm,&
        ierr)
      write (msg,'(2a)') trim(debughdr),' Failed to send inv_perm control flag.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a)') trim(debughdr),&
        ' Sent    inv_perm control flag to write root.'
      call debugmsg

    end if

! Send local inv_perm segment, if needed, to the root write task, which then
! broadcasts the complete inv_perm array to all other write tasks, if any.

    if (permute) then

      write (msg,'(2a)') trim(debughdr),&
        ' Sending inv_perm segment to write root...'
      call debugmsg

      call mpi_ssend(inv_perm(ips:ipe),interior_size(),mpi_integer,0,&
        tag_inv_perm_data,intercomm,ierr)
      write (msg,'(2a)') trim(debughdr),' Failed to send inv_perm segment.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a)') trim(debughdr),&
        ' Sent    inv_perm segment to write root.'
      call debugmsg

      write (msg,'(2a)') trim(debughdr),&
        ' Sending perm segment to write root...'
      call debugmsg

      call mpi_ssend(perm(ips:ipe),interior_size(),mpi_integer,0,&
        tag_perm_data,intercomm,ierr)
      write (msg,'(2a)') trim(debughdr),' Failed to send perm segment.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      write (msg,'(2a)') trim(debughdr),&
        ' Sent    perm segment to write root.'
      call debugmsg

    end if

    write (msg,'(2a)') trim(debughdr),' Exit.'
    call debugmsg

  end subroutine prep_write_tasks

#endif /* SERIAL */

!-------------------------------------------------------------------------------
  integer function interior_size()
!-------------------------------------------------------------------------------

    interior_size=ipe-ips+1

  end function interior_size

#ifndef SERIAL

!-------------------------------------------------------------------------------
  subroutine process_output_frame(segments,its,metasize)
!-------------------------------------------------------------------------------

! Receive metadata and field data from compute tasks, then write to disk. This
! is a private routine -- not part of the icosio API. It is called only by write
! tasks.

    integer,intent(in)::its,segments,metasize

    character(len=100)::debughdr
    character(len=20)::this='process_output_frame'
    character,allocatable::m(:)
    integer::ct,ct_ipe(0:nct-1),ct_ipns,ct_ips(0:nct-1),datasize(0:nct-1),&
      expected,ierr,ivar,mpitype,mpistatus(mpi_status_size),pos,received,tag
    real,allocatable::d(:)
    type(out_buffer_node),pointer::node
    type(glbvar_ptrs)::glbvars(segments)

    call setdebughdr(debughdr,this,its)

    write (msg,'(2a,i0,a)') trim(debughdr),' Segments=',segments,' entry.'
    call debugmsg

! Receive metadata.

! The size of the received metadata will be the size of one metadata set times
! the number of segments in this frame. All metadata sets are the same size, so
! allocate the receive buffer once and re-use for receives from all compute
! tasks.

    allocate(m(metasize),stat=ierr)
    write (msg,'(2a)') trim(debughdr),' Failed to allocate metadata buffer.'
    call die_if(ierr.ne.0,ierr,msg)

    write (msg,'(2a)') trim(debughdr),' Receiving metadata...'
    call debugmsg

! Receive one metadata set from each compute task.

    metadata_cts_loop: do ct=0,nct-1

      datasize(ct)=0

      write (msg,'(2a,i0,a)') trim(debughdr),' Receiving metadata from ct ',&
        ct,'...'
      call debugmsg

      expected=metasize
      mpitype=mpi_character
      tag=tag_metadata

      call mpi_recv(m,expected,mpitype,ct,tag,intercomm,mpistatus,ierr)
      write (msg,'(2a,i0,a)') trim(debughdr),&
        ' Failed to receive metadata from ',ct,'.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      call mpi_get_count(mpistatus,mpitype,received,ierr)
      write (msg,'(2a,2(i0,a))') trim(debughdr),' m expected ',expected,&
        ' elements, received ',received,'.'
      call die_if(received.ne.expected,ierr,msg)

      write (msg,'(2a,i0,a)') trim(debughdr),' Received  metadata from ct ',&
        ct,'.'
      call debugmsg

! Unpack metadata set for each segment and buffer in list.

      pos=0

      metadata_vars_loop: do ivar=1,segments

        call append_to_list(its,ct) ! append an empty buffer node

        node=>out_buffers(ct)%tail

        write (msg,'(2a,2(i0,a))') trim(debughdr),' Unpacking metadata ',ivar,&
          ' from ct ',ct,'...'
        call debugmsg

! Note that m is of type character, so that m is a byte count, which is what is
! required for the second argument to mpi_unpack().

        call mpi_unpack(m,metasize,pos,node%scalefactor,1,mpi_real,intercomm,&
          ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to unpack scalefactor.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_unpack(m,metasize,pos,node%accum_start,1,mpi_integer,&
          intercomm,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to unpack accum_start.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_unpack(m,metasize,pos,node%filename_len,1,mpi_integer,&
          intercomm,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to unpack filename_len.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_unpack(m,metasize,pos,node%header_len,1,mpi_integer,&
          intercomm,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to unpack header_len.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_unpack(m,metasize,pos,node%levels,1,mpi_integer,intercomm,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to unpack levels.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_unpack(m,metasize,pos,node%segment_size,1,mpi_integer,&
          intercomm,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to unpack segment_size.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_unpack(m,metasize,pos,node%time,1,mpi_integer,intercomm,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to unpack time.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        call mpi_unpack(m,metasize,pos,node%varname_len,1,mpi_integer,&
          intercomm,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to unpack varname_len.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        allocate(node%filename(node%filename_len),stat=ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to allocate filename member.'
        call die_if(ierr.ne.0,ierr,msg)
        call mpi_unpack(m,metasize,pos,node%filename,node%filename_len,&
          mpi_character,intercomm,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to unpack filename.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        allocate (node%header(node%header_len),stat=ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to allocate header member.'
        call die_if(ierr.ne.0,ierr,msg)
        call mpi_unpack(m,metasize,pos,node%header,node%header_len,&
          mpi_character,intercomm,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to unpack header.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        allocate(node%varname(node%varname_len),stat=ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to allocate varname member.'
        call die_if(ierr.ne.0,ierr,msg)
        call mpi_unpack(m,metasize,pos,node%varname,node%varname_len,&
          mpi_character,intercomm,ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to unpack varname.'
        call die_if(ierr.ne.mpi_success,ierr,msg)

        write (msg,'(2a,2(i0,a),2a)') trim(debughdr),' Unpacked  metadata ',&
          ivar,' from ct ',ct,' for ',&
          arr_to_str(node%varname,node%varname_len),'.'
        call debugmsg

        write (msg,'(2a,3(i0,a),2a)') trim(debughdr),' Received  metadata ',&
          ivar,' from ct ',ct,': allocating ',node%segment_size,' bytes for ',&
          arr_to_str(node%varname,node%varname_len),'.'
        call debugmsg

! Update the storage requirement for the eventual variable-data receive.

        datasize(ct)=datasize(ct)+node%segment_size

! Allocate a contiguous buffer for the global-sized version of the variable now
! under consideration. Sufficient information is available from the metadata
! provided by the compute root (ct 0): The levels for this variable will be the
! same across compute tasks, as will its second (nip) dimension.

        if (ct.eq.0) then
          allocate(glbvars(ivar)%ptr(node%levels,nip),stat=ierr)
          write (msg,'(2a,i0,a)') trim(debughdr),&
            ' Failed to allocate glbvars(',ivar,').'
          call die_if(ierr.ne.0,ierr,msg)
          write (msg,'(2a,3(i0,a))') trim(debughdr),' Allocated glbvars(',&
            ivar,') at ',size(glbvars(ivar)%ptr,1),' x ',&
            size(glbvars(ivar)%ptr,2),'.'
          call debugmsg
        end if

      end do metadata_vars_loop

! Calculate and record ips and ipe for this compute task. Here, after the loop
! over varaibles, 'node' still points to last output buffer node for this
! compute task, which contains sufficient information to calculte ips and ipe.

      ct_ipns=node%segment_size/node%levels
      if (ct.eq.0) then
        ct_ips(ct)=1
        ct_ipe(ct)=ct_ipns
      else
        ct_ips(ct)=ct_ipe(ct-1)+1
        ct_ipe(ct)=ct_ipe(ct-1)+ct_ipns
      end if

    end do metadata_cts_loop

    write (msg,'(2a)') trim(debughdr),' Received  metadata.'
    call debugmsg

    deallocate(m)

! Receive variable data.

    write (msg,'(2a)') trim(debughdr),' Receiving variable data...'
    call debugmsg

! Receive one variable data set from each compute task.

    vardata_cts_loop: do ct=0,nct-1

      write (msg,'(2a,i0,a)') trim(debughdr),&
        ' Receiving variable data from ct ',ct,'...'
      call debugmsg

      expected=datasize(ct)
      mpitype=mpi_real
      tag=tag_data

! If we are receiving segments of multiple variables from the current compute
! task, allocate an intermediate receive buffer. Otherwise, receive directly
! into the final buffer.

      if (segments.gt.1) then
        allocate(d(datasize(ct)),stat=ierr)
        write (msg,'(2a)') trim(debughdr),' Failed to allocate data buffer.'
        call die_if(ierr.ne.0,ierr,msg)
        call mpi_recv(d,expected,mpitype,ct,tag,intercomm,mpistatus,ierr)
      else
        call mpi_recv(glbvars(1)%ptr(:,ct_ips(ct):ct_ipe(ct)),expected,mpitype,&
          ct,tag_data,intercomm,mpistatus,ierr)
      end if

      write (msg,'(2a,i0,a)') trim(debughdr),&
        ' Failed to receive variable data from ',ct,'.'
      call die_if(ierr.ne.mpi_success,ierr,msg)

      call mpi_get_count(mpistatus,mpitype,received,ierr)
      write (msg,'(2a,2(i0,a))') trim(debughdr),&
        ' d expected ',expected,' elements, received ',received,'.'
      call die_if(received.ne.expected,ierr,msg)

      write (msg,'(2a,i0,a)') trim(debughdr),&
        ' Received  variable data from ct ',ct,'.'
      call debugmsg

! If only one segment was received, it was received directly into glbvars.
! Otherwise, unpack variable data into linked list nodes. Note that datasize
! represents number of reals, but MPI requires number of bytes in the second
! argument to mpi_unpack(), thus the scaling by size_r.

      if (segments.gt.1) then

        pos=0
        node=>out_buffers(ct)%head

        write (msg,'(2a)') trim(debughdr),' Failed to unpack data.'

        vardata_vars_loop: do ivar=1,segments

          call mpi_unpack(                              &
            d,                                          &
            datasize(ct)*size_r,                        &
            pos,                                        &
            glbvars(ivar)%ptr(:,ct_ips(ct):ct_ipe(ct)), &
            node%segment_size,                          &
            mpi_real,                                   &
            intercomm,                                  &
            ierr                                        &
            )

          write (msg,'(4a)') trim(debughdr),&
            ' Failed to unpack data segment for ',&
            arr_to_str(node%varname,node%varname_len),'.'
          call die_if(ierr.ne.mpi_success,ierr,msg)

          write (msg,'(2a,3(i0,a))') trim(debughdr),&
            ' Unpacked data segment into glbvars(',ivar,') for ipns ',&
            ct_ips(ct),'-',ct_ipe(ct),'.'
          call debugmsg

          node=>node%next

        end do vardata_vars_loop

        deallocate(d)

      end if

    end do vardata_cts_loop

    write (msg,'(2a)') trim(debughdr),' Received  variable data.'
    call debugmsg

! Flush the accumulated data to disk.

    call write_vars_to_disk(its,glbvars)

! Clean up.

    deallocate_glbvars: do ivar=1,segments
      if (associated(glbvars(ivar)%ptr)) then
        deallocate(glbvars(ivar)%ptr)
        nullify(glbvars(ivar)%ptr)
      end if
      write (msg,'(2a,i0,a)') trim(debughdr),' Deallocated glbvars(',ivar,').'
      call debugmsg
    end do deallocate_glbvars

    call clear_output_list(its)

! Sign off.

    write (msg,'(2a,i0,a)') trim(debughdr),' Segments=',segments,' exit.'
    call debugmsg

  end subroutine process_output_frame

#endif /* SERIAL */

!-------------------------------------------------------------------------------
  subroutine setdebughdr(debughdr,caller,its)
!-------------------------------------------------------------------------------

    character(len=*),intent(out)::debughdr
    character(len=*),intent(in)::caller
    integer,intent(in)::its

    write (debughdr,'(3a,2(a,i0))') trim(caller),' (',tasktype,' ',me,&
      '): its=',its

  end subroutine setdebughdr

!-------------------------------------------------------------------------------
  subroutine show_config
!-------------------------------------------------------------------------------

! Report configuration parameters. This is a private routine -- not part of the
! icosio API. It is called by both compute and write tasks.

    character(len=100)::debughdr
    character(len=11)::this='show_config'
    character*80::fmti='(a,a25,i0)'
    character*80::fmtl='(a,a25,l1)'

    call setdebughdr(debughdr,this,0)

    write (msg,fmtl) trim(debughdr),'binout = ',binout
    call debugmsg
    write (msg,fmtl) trim(debughdr),'client_server_io = ',client_server_io
    call debugmsg
    write (msg,fmtl) trim(debughdr),'debugmsg_on = ',debugmsg_on
    call debugmsg
    write (msg,fmtl) trim(debughdr),'gribout = ',gribout
    call debugmsg
    write (msg,fmtl) trim(debughdr),'i_am_compute_root = ',i_am_compute_root
    call debugmsg
    write (msg,fmtl) trim(debughdr),'i_am_compute_task = ',i_am_compute_task
    call debugmsg
    write (msg,fmtl) trim(debughdr),'i_am_write_root = ',i_am_write_root
    call debugmsg
    write (msg,fmtl) trim(debughdr),'i_am_write_task = ',i_am_write_task
    call debugmsg
    write (msg,fmti) trim(debughdr),'ipe = ',ipe
    call debugmsg
    write (msg,fmti) trim(debughdr),'ips = ',ips
    call debugmsg
    write (msg,fmti) trim(debughdr),'lun = ',lun
    call debugmsg
    write (msg,fmti) trim(debughdr),'nct = ',nct
    call debugmsg
    write (msg,fmti) trim(debughdr),'nip = ',nip
    call debugmsg
    write (msg,fmti) trim(debughdr),'nwt = ',nwt
    call debugmsg
    write (msg,fmtl) trim(debughdr),'permute = ',permute
    call debugmsg
    write (msg,fmtl) trim(debughdr),'print_diags = ',print_diags
    call debugmsg
    write (msg,fmtl) trim(debughdr),'serial = ',serial
    call debugmsg
    write (msg,fmtl) trim(debughdr),'single = ',single
    call debugmsg
    write (msg,fmtl) trim(debughdr),'using_write_tasks = ',using_write_tasks
    call debugmsg

  end subroutine show_config

!-------------------------------------------------------------------------------
  function str_to_arr(str)
!-------------------------------------------------------------------------------

! Return a character array corresponding to the supplied character string. This
! is a private routine -- not part of the icosio API. It may called by either
! compute or write tasks.

    character(len=*),intent(in)::str

    character::str_to_arr(len(str))
    integer::i

    do i=1,len(str)
      str_to_arr(i)=str(i:i)
    end do

  end function str_to_arr

!-------------------------------------------------------------------------------
  subroutine var_to_disk(filename,header,its,glbvar,varname,accum_start,&
    scalefactor)
!-------------------------------------------------------------------------------

! Write a global-size model field array to disk. This is a private routine --
! not part of the icosio API. It is called by both compute and write tasks.

    character(len=*),intent(in)::filename,header,varname
    integer,intent(in),optional::accum_start
    integer,intent(in)::its
    real,intent(in),optional::scalefactor
    real,intent(inout)::glbvar(:,:)

    character(len=100)::debughdr
    character(len=11)::this='var_to_disk'
    integer::ipn,level,levels,moved_to
    logical::append
    real::tmp

    call setdebughdr(debughdr,this,its)

    write (msg,'(2a)') trim(debughdr),' Entry.'
    call debugmsg

#ifdef SERIAL
    call ignore_unused(kind(accum_start))
    call ignore_unused(kind(scalefactor))
#endif /* SERIAL */

    levels=size(glbvar,1)

! Decide whether to append to the output file: If the current filename has been
! seen in this frame, append; otherwise create a new file, perhaps overwriting
! an existing file (as may be the case for restart runs).

    append=stringlist_has(output_filenames,trim(filename))

! Record the current filename if it is new this frame.

    if (.not.append) call stringlist_add(output_filenames,trim(filename))

! If binary history output is enabled:

    if (binout) then
      if (append) then
        open (lun,file=trim(filename),form='unformatted',iostat=ierr,&
          position='append')
      else
        open (lun,file=trim(filename),form='unformatted',iostat=ierr)
      end if
      write (msg,'(4a)') trim(debughdr),' Failed to open file ',&
        trim(filename),'.'
      call die_if(ierr.ne.0,ierr,msg)
      write (lun) str_to_arr(header)

      if (permute) then

        write (msg,'(2a)') trim(debughdr),&
          ' "inv_perm_global" must be associated when "permute" is true'
        call die_if(.not.associated(inv_perm_global),ierr,msg)
        write (msg,'(4a)') trim(debughdr),' Reordering ',trim(varname),'...'
        call debugmsg

        do ipn=1,nip-1
          moved_to=inv_perm_global(ipn)
          do while (moved_to.lt.ipn)
            moved_to=inv_perm_global(moved_to)
          end do
          if (moved_to.ne.ipn) then
            do level=1,levels
              tmp=glbvar(level,ipn)
              glbvar(level,ipn)=glbvar(level,moved_to)
              glbvar(level,moved_to)=tmp
            end do
          end if
        end do

        write (msg,'(4a)') trim(debughdr),'  Reordered ',trim(varname),'.'
        call debugmsg

      end if

      write (lun) glbvar
      write (msg,'(4a)') trim(debughdr),' Wrote binary field ',trim(varname),'.'
      call debugmsg
      close (lun)
    end if

#ifndef NOGRIB

! If grib output is enabled:

    if (gribout) then
      write (msg,'(3a)') trim(debughdr),' When gribout is true, optional ',&
        'arguments accum_start and scalefactor must be present.'
      call die_if(.not.present(accum_start).or..not.present(scalefactor),ierr,&
        msg)
      call post_write_field(glbvar,trim(varname),scalefactor,accum_start,&
        ierr)
      write (msg,'(2a)') trim(debughdr),' Bad return from post_write_field.'
      call die_if(ierr.ne.0,ierr,msg)
    end if

#endif /* NOGRIB */

    write (msg,'(2a)') trim(debughdr),' Exit.'
    call debugmsg

  end subroutine var_to_disk

!-------------------------------------------------------------------------------
  subroutine write_vars_to_disk(its,glbvars)
!-------------------------------------------------------------------------------

! Write to disk one output frame of all the variables for which this write task
! is responsible. If grib output is enabled, only one write task is allowed: The
! single write task opens the grib file, writes all necessary data to it, then
! closes it. This is a private routine -- not part of the icosio API. It is
! called only by write tasks.

    integer,intent(in)::its
    type(glbvar_ptrs),intent(inout)::glbvars(:)

    character(len=100)::debughdr
    character(len=18)::this='write_vars_to_disk'
    integer::ivar
    type(out_buffer_node),pointer::ptr

#ifndef NOGRIB
    integer::ierr,time
#endif /* NOGRIB */

    call setdebughdr(debughdr,this,its)

    write (msg,'(2a)') trim(debughdr),' Entry.'
    call debugmsg

    ptr=>out_buffers(0)%head

#ifndef NOGRIB

    time=ptr%time

! If gribout enabled, open the grib file.

    if (gribout) then
      call post_init_file(time,ierr)
      write (msg,'(2a)') trim(debughdr),' post_init_file failed.'
      call die_if (ierr.ne.0,ierr,msg)
      write (msg,'(2a,i0,a)') trim(debughdr),&
        ' post_init_file opened grib file time=',time,'.'
      call debugmsg
    end if

#endif /* NOGRIB */

! Loop over variables and have them written to disk.

    vars_loop: do ivar=1,size(glbvars)

#ifndef NOGRIB
      call var_to_disk(                            &
        arr_to_str(ptr%filename,ptr%filename_len), &
        arr_to_str(ptr%header,ptr%header_len),     &
        its,                                       &
        glbvars(ivar)%ptr,                         &
        arr_to_str(ptr%varname,ptr%varname_len),   &
        ptr%accum_start,                           &
        ptr%scalefactor                            &
        )
#else
      call var_to_disk(                            &
        arr_to_str(ptr%filename,ptr%filename_len), &
        arr_to_str(ptr%header,ptr%header_len),     &
        its,                                       &
        glbvars(ivar)%ptr,                         &
        arr_to_str(ptr%varname,ptr%varname_len)    &
        )
#endif /* NOGRIB */

      ptr=>ptr%next

    end do vars_loop

! Reset the list of filenames written to during this output interval.

    if (associated(output_filenames)) call stringlist_destroy(output_filenames)

#ifndef NOGRIB

! Close the output grib file.

    if (gribout) then
      call post_finalize_file(ierr)
      write (msg,'(2a,i0,a)') trim(debughdr),&
        ' post_finalize_file closed grib file time=',time,'.'
      call debugmsg
    end if

#endif

    write (msg,'(2a)') trim(debughdr),' Exit.'
    call debugmsg

  end subroutine write_vars_to_disk

end module icosio
