module module_wrf_phy_run

  implicit none

#include <gptl.inc>

contains

!*********************************************************************
!       "Run" method for WRF physics
!*********************************************************************

subroutine wrf_phy_run(its)
use fimnamelist            ,only: nts, itsStart
use module_wrfphysics      ,only: wrf_physics
use module_initial_chem_namelists, only: cu_physics, mp_physics

!  Declare dummy arguments
integer, intent(in) :: its

!  Declare local variables:
integer :: ret

ret = gptlstart ('wrf_phy_run')

  !...........................................................
  ! Advance the physics component by one time step unless this 
  ! is the last (nts+1) iteration.  
  ! This complexity is required for the NCEP ESMF approach 
  ! in which single-phase DYN and PHY components alternate 
  ! execution during each time step.  
  !
if (its < itsStart+nts) then

  !...........................................................
  ! call WRF phyics
  !...........................................................
  if ((.not.mp_physics==0).or.(.not.cu_physics==0)) then
!   write(6,*)'call wrfphysics ',mp_physics,cu_physics
    call wrf_physics(its)
  endif

endif

ret = gptlstop ('wrf_phy_run')

return
end subroutine wrf_phy_run
end module module_wrf_phy_run
