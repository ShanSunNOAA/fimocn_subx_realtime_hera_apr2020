module global_bounds
  use module_control, only: nip
  implicit none

  public

  integer :: ihe     = -1 ! upper halo index
  integer :: ihs     = -1 ! lower halo index
  integer :: ime     = -1 ! upper memory bound
  integer :: ims     = -1 ! lower memory bound
  integer :: ipe     = -1 ! upper owned index
  integer :: ipe_set = -1 ! upper owned index (default)
  integer :: ips     = -1 ! lower owned index
  integer :: ips_set = -1 ! lower owned index (default)
  integer :: npes    = 1  ! number of compute tasks (1 in serial mode)
  integer :: myrank  = 0  ! my rank (0 in serial mode)

  logical          :: parallel_build = .false. ! set to true by SMS below
  logical, private :: set_global_bounds_called = .false.

contains

  subroutine set_global_bounds ()
!sms$distribute(dh,1) begin
    integer, allocatable :: dummy(:)
!sms$distribute end

    allocate (dummy(nip))
    ips_set = 1
    ipe_set = nip
!sms$to_local (dh:<1,ips_set:lbound>,<1,ipe_set:ubound>) begin
    ips = ips_set
    ipe = ipe_set
!sms$to_local end
!sms$ignore begin
    ims = lbound(dummy,1)
    ime = ubound(dummy,1)
!sms$ignore end
    deallocate (dummy)
!sms$insert parallel_build=.true.
!sms$comm_size(npes)
!sms$comm_rank(myrank)
    set_global_bounds_called = .true.
  end subroutine set_global_bounds

  subroutine set_global_bounds_halo ()
    integer :: ipn

    if (.not.set_global_bounds_called) then
      write (6,*) 'set_global_bounds_halo: must call set_global_bounds first!'
      stop
    end if
    ihs = ipe + 1
    ihe = ipe
!sms$parallel (dh,ipn) begin
!sms$halo_comp(<1,1>) begin
    do ipn=1,nip
      continue
    enddo
!sms$halo_comp end
!sms$parallel end
    ihe = ipn - 1
!sms$comm_size(npes)
!sms$comm_rank(myrank)
  end subroutine set_global_bounds_halo
end module global_bounds
