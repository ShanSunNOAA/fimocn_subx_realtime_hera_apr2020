! Band-aid routine callable by GFS physics routines which is not
! contained in a module. It is needed because it needs to live
! in dynamics (so it can use the module), and thus GFS physics
! routines which call it cannot use it!

subroutine physics_get_ipn_its_mype_nomodule (ipn, its, mype)
  use global_bounds,   only: myrank
  use physics_ipn_its, only: physics_get_ipn_its

  integer, intent(out) :: ipn
  integer, intent(out) :: its
  integer, intent(out) :: mype

  call physics_get_ipn_its (ipn, its)
  mype = myrank
end subroutine physics_get_ipn_its_mype_nomodule
