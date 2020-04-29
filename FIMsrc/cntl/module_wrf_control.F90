module module_wrf_control
  USE module_initial_chem_namelists,only:chem_opt,mp_physics
implicit none

!********************************************************************
!	This module specifies control variables for wrf physics 
!       and chemistry.  
!********************************************************************

! standard WRF index bounds
integer            ::  ims    ! =1
integer            ::  ime    ! =1
integer            ::  ids    ! =1
integer            ::  ide    ! =1
integer            ::  its    ! =1
integer            ::  ite    ! =1
integer            ::  jms    ! =1
integer            ::  jme    ! =nip
integer            ::  jds    ! =1
integer            ::  jde    ! =nip
integer            ::  jts    ! =1
integer            ::  jte    ! =nip
integer            ::  kms    ! =1
integer            ::  kme    ! =nvl+1
integer            ::  kds    ! =1
integer            ::  kde    ! =nvl+1
integer            ::  kts    ! =1
integer            ::  kte    ! =nvl

! WRF physics and chem parameters
integer, parameter :: numgas=1        ! # of tracers for gas phase chemistry
integer :: DoForce                    ! Timestep interval for averagin forcing for GF
integer :: DoForce_ave_length=820    ! Timestep interval n seconds for averagin forcing for GF
logical :: AveCuforcingNow

!
! GFS physics _ gocart very light for fim
!
!integer, parameter :: num_moist=2+1
!integer, parameter :: num_chem=13
!integer, parameter :: num_emis_ant = 6
!
! following for Lin et al. + regular GOCART
!
!integer, parameter :: num_moist=6+1
!integer, parameter :: num_chem=18
!integer, parameter :: num_emis_ant = 6
!integer, parameter :: num_emis_vol = 0
!
! volcanic ash only
!integer, parameter :: num_moist=2+1
!integer, parameter :: num_chem=23
!integer, parameter :: num_emis_ant = 6
!integer, parameter :: num_emis_vol = 10
!
! light gocart + reduced volcanic ash only (4 size bins)
!!!!!!!!!! REG TEST SETUP !!!!
!integer, parameter :: num_moist=2+1
!integer, parameter :: num_chem=17
!integer, parameter :: num_emis_ant = 6
!integer, parameter :: num_emis_vol = 4
 integer            :: num_chem,num_moist,num_emis_ant,num_emis_vol
!integer            :: num_moist =0
!integer            :: num_emis_ant =0
!integer            :: num_emis_vol =0
! volcanic ash only (4 size bins)
!integer, parameter :: num_moist=2+1
!integer, parameter :: num_chem=4
!integer, parameter :: num_emis_ant = 0
!integer, parameter :: num_emis_vol = 4
! Pure GOCART (volcanic ash included in p25 and p10)
!!!!! CURRENT REAL-TIME
!integer, parameter :: num_moist=2+1
!integer, parameter :: num_chem=19
!integer, parameter :: num_emis_ant = 6
!integer, parameter :: num_emis_vol = 4
! Pure GOCART (volcanic ash included in p25 and p10) + mp_phys=4 (qv,qc,qr,qi,qs)
!integer, parameter :: num_moist=4+1
!integer, parameter :: num_chem=19
!integer, parameter :: num_emis_ant = 6
!integer, parameter :: num_emis_vol = 4
!
!
integer            :: num_emis_seas_bb=0
integer		   :: num_emis_seas_ant=0
integer, parameter :: nbands=14
integer, parameter :: nbandlw=16
integer, parameter :: num_soil_layers=4
integer, parameter :: num_scalar=1
integer, parameter :: nvl_gocart=55  ! number of input levels from gocart file
integer, parameter :: num_ext_coef = 5
integer, parameter :: num_bscat_coef = 3
integer, parameter :: num_asym_par = 3

! namelist variables
! not yet used
integer            :: ChemistryInterval = 0  ! Interval in seconds to call chemistry, 0 => every time step
!Control variables calculated in init.F90 from namelist variables
integer            :: CallChemistry          ! Timestep interval to call chemistry
integer            :: CallBiom               ! Timestep interval to call biomass burning plumerise

contains

!subroutine wrf_control(nvarp,ntrb,num_chem,num_moist,num_emis_ant,num_emis_vol)
!  integer, intent(inout) :: nvarp,ntrb,num_chem,num_moist,num_emis_ant,num_emis_vol
subroutine wrf_control(nvarp,ntrb)
  integer, intent(inout) :: nvarp,ntrb
  ! add WRF variables to ntrb-dimensioned arrays
  if(mp_physics.eq.0 .and. chem_opt .eq. 317)then
     num_chem=17
     num_moist=3
     num_emis_ant=6
     num_emis_vol=4
     ntrb=ntrb+num_moist+num_chem-3           ! # of tracers + num_moist-3, num_chem - no ice variable transported
!    nvarp=nvarp+num_chem+num_moist - 3
  elseif (chem_opt .eq. 300 .and. mp_physics.eq.0 )then
     num_chem=19
     num_moist=3
     num_emis_ant=6
     num_emis_vol=0
     num_emis_seas_bb=6
     num_emis_seas_ant=6
     ntrb=ntrb+num_moist+num_chem-3           ! # of tracers + num_moist-3, num_chem - no ice variable transported
!    nvarp=nvarp+num_chem+num_moist - 2       ! nvarp includes all variables for pressure level output
  elseif (chem_opt .eq. 317 .and. mp_physics.eq.3 )then
     num_chem=17
     num_moist=5
     num_emis_ant=6
     num_emis_vol=4
     ntrb=ntrb+num_moist+num_chem-2           ! # of tracers + num_moist-2, num_chem - transport of qi
!    nvarp=nvarp+num_chem+num_moist - 2       ! nvarp includes all variables for pressure level output
  endif
end subroutine wrf_control

end module module_wrf_control
