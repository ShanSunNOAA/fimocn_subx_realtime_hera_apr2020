Module module_set_wrfphys
USE module_initial_chem_namelists
CONTAINS
SUBROUTINE set_wrfphys (mp_physics)
USE module_wrfphysvars
USE module_control, only: ntra,ntrb
USE module_wrf_control, only: num_moist,num_chem
implicit none
integer, intent(in) :: mp_physics
integer :: itest
!
! microphysics
!
if(mp_physics.eq.4)then
   if(num_moist.ne.5)then
      write(6,*) ' num_moist is not equal 5 '
      stop
   endif
   itest=ntra+ntrb-num_moist-num_chem+2
   if(itest.ne.ntra)then
      write(6,*) ' ntra ist falsch'
      write(6,*) ' ntra,num_moist,num_chem = ',ntra,num_moist,num_chem
      stop
   endif
   p_qv=1
   f_qv=.true.
   p_qc=2
   f_qc=.true.
   p_qr=3
   f_qr=.true.
   p_qi=4
   f_qi=.true.
   p_qs=5
   f_qs=.true.
else if(mp_physics.eq.2)then
   if(num_moist.ne.6)then
      write(6,*) ' num_moist is not equal 6 '
      stop
   endif
   itest=ntra+ntrb-num_moist-num_chem+2
   if(itest.ne.ntra)then
      write(6,*) ' ntra ist falsch'
      stop
   endif
   p_qv=1
   f_qv=.true.
   p_qc=2
   f_qc=.true.
   p_qr=3
   f_qr=.true.
   p_qi=4
   f_qi=.true.
   p_qs=5
   f_qs=.true.
   p_qg=6
   f_qg=.true.
else if(mp_physics.eq.0)then
   if(num_moist.ne.3)then
      write(6,*) ' num_moist is not equal 3 '
!     stop
   endif
   itest=ntra+ntrb-num_moist-num_chem+3
   if(itest.ne.ntra)then
      write(6,*) ' ntra ist falsch',ntra,num_moist,num_chem
!     stop
   endif
! qc and qi are only dependent on T, no separate transport
   p_qv=1
   f_qv=.true.
   p_qc=2
   f_qc=.true.
   p_qi=3
   f_qi=.true.
endif
END SUBROUTINE set_wrfphys

END MODULE module_set_wrfphys
