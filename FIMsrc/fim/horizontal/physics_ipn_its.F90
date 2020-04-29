! Module contains public and private routines to enable retrieving ipn and its
! value from within GFS physics. It is more complicated now with thread-enabled
! GFS physics because ipn has a different value for each thread.
!
! Jim Rosinski, April, 2015

module physics_ipn_its
  use fimnamelist, only: nthreads  ! number of OMP threads

  implicit none

  private
  public :: physics_init_ipn, physics_set_ipn_its, physics_get_ipn_its, physics_reset_ipn

  integer, parameter :: nofalse = 16      ! size to prevent false sharing
  integer, allocatable :: ipn_thread(:,:) ! private array holds current ipn for each thread
  integer :: its = -1                     ! private copy of its
  logical :: init_done = .false.          ! verify init routine called before any worker routines

contains
  
! physics_init_ipn: Allocate private array of size "nthreads" to hold per-thread ipn value
  subroutine physics_init_ipn ()
    if (init_done) then
!SMS$IGNORE BEGIN
      write(6,*)'physics_init_ipn: init_done is true: Doing nothing'
!SMS$IGNORE END
      return
    end if

    allocate (ipn_thread(nofalse,0:nthreads-1))
    ipn_thread(:,:) = -1
    init_done = .true.
  end subroutine physics_init_ipn

! physics_set_ipn_its: Store the current value of ipn for this thread
  subroutine physics_set_ipn_its (ipn, itsin)
    integer, intent(in) :: ipn
    integer, intent(in) :: itsin
    integer :: mythread
#ifdef WHOLE_MODEL_OMP
    integer, external :: omp_get_thread_num
#endif    
    if (.not. init_done) then
!SMS$IGNORE BEGIN
      write(6,*)'physics_set_ipn_its: physics_init_ipn has not been called: Doing nothing!'
!SMS$IGNORE END
      return
    end if

    mythread = 0
#ifdef WHOLE_MODEL_OMP
    mythread = omp_get_thread_num ()
    if (mythread >= nthreads) then
!SMS$IGNORE BEGIN
      write(6,*)'physics_set_ipn_its: mythread=',mythread,' nthreads=',nthreads
!SMS$IGNORE END
      return
    end if
#endif
    ipn_thread(1,mythread) = ipn
    its = itsin
  end subroutine physics_set_ipn_its

! physics_get_ipn_its: Retrieve its and thread-local value of ipn
  subroutine physics_get_ipn_its (ipnout, itsout)
    integer, intent(out) :: ipnout
    integer, intent(out) :: itsout
    integer :: mythread
#ifdef WHOLE_MODEL_OMP
    integer, external :: omp_get_thread_num
#endif    

    if (.not. init_done) then
!SMS$IGNORE BEGIN
      write(6,*)'physics_get_ipn_its: physics_init_ipn has not been called: Doing nothing!'
!SMS$IGNORE END
      return
    end if

    mythread = 0
#ifdef WHOLE_MODEL_OMP
    mythread = omp_get_thread_num ()
    if (mythread >= nthreads) then
!SMS$IGNORE BEGIN
      write(6,*)'physics_get_ipn_its: mythread=',mythread,' nthreads=',nthreads
!SMS$IGNORE END
      return
    end if
#endif
    ipnout = ipn_thread(1,mythread)
    itsout = its
  end subroutine physics_get_ipn_its

! physics_reset_ipn: Reset ipn to bad (out-of-scope) value 
  subroutine physics_reset_ipn (ipn)
    integer, intent(in) :: ipn
    integer :: mythread
#ifdef WHOLE_MODEL_OMP
    integer, external :: omp_get_thread_num
#endif    

    if (.not. init_done) then
!SMS$IGNORE BEGIN
      write(6,*)'physics_reset_ipn: physics_init_ipn has not been called: Doing nothing!'
!SMS$IGNORE END
      return
    end if

    mythread = 0
#ifdef WHOLE_MODEL_OMP
    mythread = omp_get_thread_num ()
    if (mythread >= nthreads) then
!SMS$IGNORE BEGIN
      write(6,*)'physics_reset_ipn: mythread=',mythread,' nthreads=',nthreads
!SMS$IGNORE END
      return
    end if
#endif
    ipn_thread(1,mythread) = -1  ! Reset to out-of-scope
  end subroutine physics_reset_ipn
end module physics_ipn_its
