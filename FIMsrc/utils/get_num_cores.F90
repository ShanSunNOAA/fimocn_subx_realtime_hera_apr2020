program get_num_cores

  use module_wtinfo,only:wtinfo
  use fimnamelist, only : nprocs, cpn, mpipn, nthreads, root_own_node, &
                          num_write_tasks, max_write_tasks_per_node

  implicit none

  integer :: leftover                  ! left over tasks for a mod() calculation
  integer :: mwtpn                     ! max write tasks per node
  integer :: nct                       ! number of compute tasks
  integer :: num_nodes_wt              ! number of nodes filled w/ write tasks
  integer :: numtasks_batch            ! numtasks_mpirun modified to fill all nodes
  integer :: numcores_batch            ! numtasks_batch * mpipn
  integer :: numtasks_donothing        ! number of tasks which will be MPI do-nothing tasks
  integer :: numtasks_mpirun           ! tasks needed by mpirun cmd
  integer :: nwt                       ! number of write tasks
  integer :: tot_nodes                 ! number of nodes needed by batch environment

  logical :: compute_tasks_after_write_tasks ! whether there are compute tasks after write tasks

  call wtinfo ()
  nwt   = num_write_tasks
  mwtpn = max_write_tasks_per_node
  nct = nprocs

  if (mwtpn > mpipn) then
    write(6,*) 'get_num_cores: Max write tasks per node=', mwtpn, &
               ' exceeds MPI tasks per node=', mpipn
    stop 999
  end if

  numtasks_donothing = 0

! Initialize numtasks_mpirun to number of compute tasks, then add for root being
! on his own node, and write tasks

  compute_tasks_after_write_tasks = .false.
  numtasks_mpirun = nct
  if (root_own_node) then
    if (nwt > 0 .or. nct > 1) then ! fill rest of 1st node with idle tasks
      numtasks_mpirun = numtasks_mpirun + mpipn - 1
      numtasks_donothing = numtasks_donothing + mpipn - 1
    end if
    if (nct > 1) then
      compute_tasks_after_write_tasks = .true.
    end if
  else
    if (nct > mpipn) then    ! more than root node needed for compute tasks
      compute_tasks_after_write_tasks = .true.
    else if (nwt > 0) then   ! fill rest of 1st node with idle tasks
      numtasks_mpirun = numtasks_mpirun + (mpipn - nct)
      numtasks_donothing = numtasks_donothing + (mpipn - nct)
    end if
  end if

! Write tasks: first handle nodes filled with write tasks

  num_nodes_wt = nwt / mwtpn
  numtasks_mpirun = numtasks_mpirun + (num_nodes_wt * mpipn)
  numtasks_donothing = numtasks_donothing + num_nodes_wt * (mpipn - mwtpn)

! Last node with write tasks may have less than mwtpn

  leftover = mod (nwt, mwtpn)
  if (leftover > 0) then
    if (compute_tasks_after_write_tasks) then   
      ! Allocate the full node since compute tasks start at next empty node
      numtasks_mpirun = numtasks_mpirun + mpipn
      numtasks_donothing = numtasks_donothing + (mpipn - leftover)
    else
      ! No compute tasks to follow: Just account for the remaining write tasks
      numtasks_mpirun = numtasks_mpirun + leftover
      num_nodes_wt = num_nodes_wt + 1
    end if
  end if

! Check to ensure numbers were calculated correctly

  if (nct + nwt + numtasks_donothing /= numtasks_mpirun) then
    write(*,*) 'get_num_cores failure:', nct, nwt, numtasks_donothing, numtasks_mpirun
    stop 999
  end if

  write (*,'(a,i0)') 'num_tasks_mpirun:', numtasks_mpirun
  write (*,'(a,i0)') 'num_tasks_per_node:', mpipn

! numtasks_batch considers all nodes to be full.
! This is critical for jaguar* because the core count
! requested on the #PBS line must be a multiple of cpn

  numtasks_batch = numtasks_mpirun
  leftover = mod (numtasks_batch, mpipn)
  if (leftover /= 0) then
    numtasks_batch = numtasks_batch + (mpipn - leftover)
  end if

  write (*,'(a,i0)') 'num_tasks_batch:', numtasks_batch

  numcores_batch = numtasks_batch * nthreads
  write (*,'(a,i0)') 'num_cores_batch:', numcores_batch

  tot_nodes = numtasks_batch / mpipn
  write (*,'(a,i0)') 'tot_nodes:', tot_nodes
  write (*,'(a,i0)') 'num_nodes_wt:', num_nodes_wt
  write (*,'(a,i0)') 'num_tasks_donothing:', numtasks_donothing
  write (*,'(a,i0)') 'num_threads:', nthreads
  write (*,'(a,i0)') 'cores_per_node:', cpn

! root_own_node is needed in some of the scripts: print in a way that grep 
! will be sure to find it

  if (root_own_node) then
    write(*,'(a)')'root_own_node:TRUE'
  else
    write(*,'(a)')'root_own_node:FALSE'
  end if
end program get_num_cores
