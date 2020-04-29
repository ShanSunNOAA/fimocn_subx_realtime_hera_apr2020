module module_wrf_phy_init

  implicit none

#include <gptl.inc>

contains

  subroutine wrf_phy_init
!*********************************************************************
!       Loads the initial variables and constants for the WRF physics 
!       component.  
!*********************************************************************

    use module_control    ,only: nip
    use fimnamelist       ,only: nvl
    use module_wrf_control,only: ids,ide,jds,jde,kds,kde, &
                                 ims,ime,jms,jme,kms,kme, &
                                 its,ite,jts,jte,kts,kte
    use module_wrfphysvars
    use module_wrfphys_alloc,only: wrfphys_alloc
    use module_wrf_share    ,only: wrf_set_array_bounds
    use module_wrf_variables,only: phys3dwrf,phys3dwrf_ave,phys2dwrf,exch,pb2d

!TODO:  clean up duplication in these modules and add "only"
! for chemistry and WRF physics namelists:
    use module_chem_namelist_defaults
! contains config_flags
    use module_initial_chem_namelists  !, only: cu_physics, mp_physics
! TBH:  Ignore these so PPP doesn't have to translate them
    use module_cu_gd        ,only: gdinit
    use module_species_decs
    use module_set_wrfphys
    use module_ra_rrtmg_lw, only: rrtmg_lwinit
    use module_ra_rrtmg_sw, only: rrtmg_swinit
    use fimnamelist,        only: cu_physics, mp_physics, chem_opt, &
                                  ra_sw_physics, ra_lw_physics

! Local variables
    integer :: ret
    real*8 :: twrf_phy_init

!SMS$insert integer :: mype

!sms$comm_rank(mype)

    ret = gptlstart ('wrf_phy_init')

    call wrf_set_array_bounds(nvl,nip,                 &
                              ids,ide,jds,jde,kds,kde, &
                              ims,ime,jms,jme,kms,kme, &
                              its,ite,jts,jte,kts,kte)

    print *,'wrfphysics, cu_physics, mp_physics = ',cu_physics,mp_physics
    print *,'chemwrf, chem_opt = ',chem_opt

    allocate( phys3dwrf(nvl,nip,11)) ! WRF Physics diagnostic variable to store tendencies for microphys and cu
    allocate( phys3dwrf_ave(nvl,nip,6)) ! some timeavraged WRF Physics diagnostic variable to store tendencies for cu
                                     ! 1 = rqvcu
                                     ! 2 = rqvbl
                                     ! 3 = rqvf
                                     ! 4 = rthcu
                                     ! 5 = rthbl
                                     ! 6 = rthra
                                     ! 7 = rthf
                                     ! 8 = rqccu
                                     ! 9 = rqrcu
                                     ! 10 = rqscu
                                     ! 11 = rqicu
    allocate( phys2dwrf(nip,8))      ! Physics diagnostic variable

    phys3dwrf = 0.
    phys3dwrf_ave = 0.
    phys2dwrf = 0.

    if (chem_opt == 0 .and. mp_physics == 0 .and. cu_physics == 0 .and. ra_sw_physics == 0) then
      ret = gptlstop ('wrf_phy_init')
      return
    end if

    config_flags%mp_physics = mp_physics
    config_flags%cu_physics = cu_physics
    config_flags%ra_sw_physics = ra_sw_physics
    config_flags%ra_lw_physics = ra_lw_physics
    call set_wrfphys (mp_physics)
!
! we need the wrfphys variables for wrfchem
!
    if ((.not.mp_physics==0).or.(.not.cu_physics==0).or.(.not.chem_opt==0)) then
      write(0,*)'allocatewrfphys variables'
      call wrfphys_alloc
    end if   ! if ((.not.mp_physics==0).or.(.not.cu_physics==0)) then
!TODO:  move to wrfphys_alloc() ??  
!TODO:  avoid allocation when these variables are not used
    allocate( exch(nvl,nip)) ! 
    exch = 0.
    allocate( pb2d(nip)) ! 
    pb2d = 0.
    if (ra_sw_physics > 0 )then
     write(0,*)'initialize wrf rrtmg'
     call rrtmg_swinit
     call rrtmg_lwinit
    endif

    if (.not.cu_physics==0) then
!TODO:  this really only applies to GD cumulus scheme -- improve "if" 
!TODO:  logic and generalize call
      CALL gdinit(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,        &
                  MASS_FLUX,1004.6855,.false.,                &
                  0,0,0,                                      &
                  RTHFTEN, RQVFTEN,                           &
                  APR_GR,APR_W,APR_MC,APR_ST,APR_AS,          &
                  APR_CAPMA,APR_CAPME,APR_CAPMI,              &
                  .false.,                                    &
                  ids, ide, jds, jde, kds, kde,               &
                  ims, ime, jms, jme, kms, kme,               &
                  its, ite, jts, jte, kts, kte               )
    endif

    ret = gptlstop ('wrf_phy_init')
    ret = gptlget_wallclock ('wrf_phy_init', 0, twrf_phy_init)  ! The "0" is thread number
    print"(' WRF PHYSICS INIT time:',F10.0)", twrf_phy_init

    return

  end subroutine wrf_phy_init
end module module_wrf_phy_init
