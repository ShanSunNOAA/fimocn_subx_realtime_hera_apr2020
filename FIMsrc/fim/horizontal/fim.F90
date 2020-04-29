!*********************************************************************
!  icosahedral flow-following model
!  Authors: Alexander E. MacDonald & Jin-Luen Lee  11/12/05
!  Lead Developer:  J. L. Lee
!  Design Team:  J.L. Lee, R. Bleck, A. E. MacDonald, S. Benjamin
!  Computational design: A. E. MacDonald, J. Middlecoff, D. Schaffer
!*********************************************************************

program fim

  use fimnamelist,       only: mpipn, cpn, root_own_node, nthreads
  use global_bounds,     only: myrank
  use module_core_setup, only: base_comm, iam_fim_task, iam_write_task
  use module_fim_init,   only: init

#ifndef SERIAL
  use icosio,            only: icosio_run
#endif

  implicit none

!TODO:  Strictly speaking, the MPI timers called from StartTimer should not
!TODO:  be called prior to MPI_INIT.  For a serial build, MPI_INIT is never
!TODO:  called.  In practice this has not been a problem.  If it becomes
!TODO:  a problem fix it by using a different timer inside StartTimer and
!TODO:  IncrementTimer for a serial build.

#include <gptl.inc>

  integer :: ret
  logical :: setaffinity = .false.     ! Belongs in a namelist
  integer :: pin2core = 1              ! 1=>pin to core, 0=>pin to range of cores
  integer, external :: set_affinity

!sms$start

  call gptlprocess_namelist ('GPTLnamelist', 77, ret) ! process gptl namelist, if present
  ret = gptlinitialize ()                            ! 1-time gptl init function
  ret = gptlstart ('main')                           ! start timing the main program

! NOTE: Executable SMS directives must not be placed before the call to init!
! NOTE: This includes writes or prints without !SMS$ignore because they
!       generate SMS code.

  call init

  if (iam_fim_task) then
! Only set affinity when OpenMP is enabled for the entire model
#ifdef WHOLE_MODEL_OMP
! Set setaffinity=.true. only with caution!
    if (setaffinity) then
      if (root_own_node) then
        if (myrank == 0) then
! Affinity = entire node for rank 0 when root_own_node=.true.
          write(6,*)'Rank 0 setting affinity with cpn, pin2core=',cpn,pin2core
          ret = set_affinity (myrank, cpn, cpn, pin2core)
        else
! Since rank 1 starts a node when root_own_node=.true., set affinity as
! though ranks start at 1 not 0
          if (myrank < 2*mpipn) then
            write(6,*)'Rank ',myrank,' setting affinity with cpn, nthreads, pin2core=', &
            cpn,nthreads,pin2core
          end if
          ret = set_affinity (myrank-1, cpn, nthreads, pin2core)
        end if
      else
! root_own_node=.false.: the simplest case
        if (myrank == 0) then
          write(6,*)'Rank 0 setting affinity with cpn, nthreads, pin2core=', &
                    cpn,nthreads,pin2core
        end if
        ret = set_affinity (myrank, cpn, nthreads, pin2core)
      end if
      if (ret /= 0) then
        write(6,*)'FIM error: set_affinity gave non-zero return value=', ret
        stop
      end if
    end if
#endif
! Print affinity for 1st 2 nodes
    if (myrank < 2*mpipn) then
      call print_affinity (myrank)
    end if

    call run ()
  else if (iam_write_task) then
#ifndef SERIAL
    call icosio_run ()
#endif
  end if

  ret = gptlstop ('main')

  if (iam_fim_task) then
!SMS$INSERT call sms__unstructured_print_timers
    call finalize
    print *,'Program exited normally'
  end if

!sms$stop (base_comm)

end program fim
