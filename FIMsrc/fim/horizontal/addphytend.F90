!*********************************************************************
        subroutine addphytend(u_tdcy_phy,v_tdcy_phy,trc_tdcy_phy, &
                              us3d,vs3d,tr3d)
!       T. Smirnova            - April 2012
!*********************************************************************

use module_constants
use module_control  ,only: dt,nvlp1,nip,ntra,ntrb
use fimnamelist     ,only: nvl,addtend
! <begin rant>
! pr3d should be an intent(in) argument but it is contrary to ESMF architecture 
! to make it available to be read from the DYN import state (because import 
! state should be write-only).  So we must get it from the DYN module directly. 
! This illustrates a fundamental flaw in the ESMF component design.  Doing 
! this in the "ESMF way" would require separating the computation of 
! tr3d(:,:,1) into to parts with the part requring pr3d being done in the 
! DYN component.  Or a copy of pr3d could be stored as state in the coupler 
! (requring care to keep the copy in sync every time a changes was made).  
! Or a multi-phase component interaction could be used to make pr3d available 
! which negates the whole value of a component architecture (you try to 
! "couple" a 3-phase component with another component that uses 5 phases).  
! None of these approaches are practical because they would take a simple 
! equation now implemented in two lines of code and make its implementation 
! very complex.  Changes would be difficult to implement and prone to error.  
! <end rant>
use module_variables,only: pr3d
USE MACHINE         ,only: kind_evod

implicit none

!  Declare dummy arguments
!SMS$DISTRIBUTE(dh,2) BEGIN
real*8, intent(in) :: u_tdcy_phy  (nvl,nip)           ! physics forcing of u
real*8, intent(in) :: v_tdcy_phy  (nvl,nip)           ! physics forcing of v
real*8, intent(in) :: trc_tdcy_phy(nvl,nip,ntra+ntrb) ! physics forcing of tracers
real,  intent(out) :: us3d(nvl,nip)
real,  intent(out) :: vs3d(nvl,nip)
real,  intent(out) :: tr3d(nvl,nip,ntra+ntrb)
!SMS$DISTRIBUTE END

! Local Variables
integer ipn,ivl
  real    :: press3d
  real    :: rocp1,rocpr

      rocp1=rd/cp+1.
      rocpr=cp/rd

!SMS$PARALLEL (dh,ipn) BEGIN
!    add physics tendencies to dynamic arrays

    if(addtend) then
      do ipn=1,nip
        do ivl=1,nvl
          us3d(ivl,ipn) = us3d(ivl,ipn) + u_tdcy_phy(ivl,ipn)*dt
          vs3d(ivl,ipn) = vs3d(ivl,ipn) + v_tdcy_phy(ivl,ipn)*dt
! mid-layer pressure
          press3d=((pr3d(ivl,ipn)**rocp1-pr3d(ivl+1,ipn)**rocp1)/     &
                 ((pr3d(ivl,ipn)       -pr3d(ivl+1,ipn)       )*rocp1))**rocpr
!tgs - test_p          press3d=0.5*(pr3d(ivl,ipn)+pr3d(ivl+1,ipn))
          tr3d(ivl,ipn,1) = tr3d(ivl,ipn,1) + trc_tdcy_phy(ivl,ipn,1)*dt            &
            *(p1000/press3d)**(rd/cp)*(1.+0.6078*max(REAL(qvmin,kind_evod), &
             (tr3d(ivl,ipn,2)+ trc_tdcy_phy(ivl,ipn,2)*dt)))
!tgs - if temperature tendency computed in terms of potential temp then use the following 2 lines
!          tr3d(ivl,ipn,1) = tr3d(ivl,ipn,1) + trc_tdcy_phy(ivl,ipn,1)*dt &
!             *(1.+0.6078*max(REAL(qvmin,kind_evod),(tr3d(ivl,ipn,2)+ trc_tdcy_phy(ivl,ipn,2)*dt)))
          tr3d(ivl,ipn,2) = max(REAL(qvmin,kind_evod),(tr3d(ivl,ipn,2)+ trc_tdcy_phy(ivl,ipn,2)*dt)) ! qv
          tr3d(ivl,ipn,3) = max(REAL(qwmin,kind_evod),(tr3d(ivl,ipn,3) + trc_tdcy_phy(ivl,ipn,3)*dt)) ! qc
          tr3d(ivl,ipn,4) = max(REAL(qvmin,kind_evod),(tr3d(ivl,ipn,4) + trc_tdcy_phy(ivl,ipn,4)*dt)) ! oz
          ! Save DYN fields to avoid INOUT INTENT for these arrays in 
          ! cpl_phy_to_dyn() which breaks ESMF data model.  See rant above.  
        enddo
      end do
    endif ! addtend
!SMS$PARALLEL END

return
end subroutine addphytend
