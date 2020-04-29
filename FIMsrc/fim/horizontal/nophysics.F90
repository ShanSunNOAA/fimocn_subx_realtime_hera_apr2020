!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! "nophysics()" sets to zero fields normally output from GFS physics. Dynamics fields for wind, 
! temperature, and constituents are not mentioned here because they go through the coupler and will 
! be unmodified from previous dynamics values if nophysics is called instead of phy_run.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nophysics ()

  use gfs_physics_internal_state_mod, only: gis_phy
  use global_bounds, only: ims, ime, ips, ipe

  implicit none

  integer :: ipn  ! horizontal index
  integer :: k    ! vertical index

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (k) schedule (static)
  do ipn=ips,ipe
    gis_phy%flx_fld%geshem(ipn,1) = 0.
    gis_phy%flx_fld%rainc(ipn,1) = 0.
    gis_phy%sfc_fld%tsea(ipn,1) = 0.
    gis_phy%sfc_fld%uustar(ipn,1) = 0.
    gis_phy%flx_fld%hflx(ipn,1) = 0.
    gis_phy%flx_fld%evap(ipn,1) = 0.
    gis_phy%sfc_fld%sheleg(ipn,1) = 0.
    gis_phy%sfc_fld%canopy(ipn,1) = 0.
    gis_phy%sfc_fld%hice(ipn,1) = 0.
    gis_phy%sfc_fld%fice(ipn,1) = 0.
    do k=1,gis_phy%lsoil
      gis_phy%sfc_fld%stc(k,ipn,1) = 0.
      gis_phy%sfc_fld%smc(k,ipn,1) = 0.
    end do
    gis_phy%flx_fld%sfcdsw(ipn,1) = 0.
    gis_phy%flx_fld%sfcdlw(ipn,1) = 0.
    gis_phy%flx_fld%sfcusw(ipn,1) = 0.  !hli 09/2014
    gis_phy%flx_fld%sfculw(ipn,1) = 0.  !hli 09/2014
    gis_phy%flx_fld%topdsw(ipn,1) = 0.  !hli 09/2014
    gis_phy%flx_fld%topusw(ipn,1) = 0.  !hli 09/2014
    gis_phy%flx_fld%topulw(ipn,1) = 0.  !hli 09/2014
    gis_phy%sfc_fld%t2m(ipn,1) = 0.
    gis_phy%sfc_fld%q2m(ipn,1) = 0.
    gis_phy%flx_fld%u10m(ipn,1) = 0.
    gis_phy%flx_fld%v10m(ipn,1) = 0.
    do k=1,gis_phy%nmtvr
      gis_phy%hprime(k,ipn,1) = 0.
    end do
    gis_phy%fluxr(1,ipn,1) = 0.        ! 1st element is rlut. Rest not copied to dynamics
  end do
!$OMP END PARALLEL DO
!sms$ignore end

  return
end subroutine nophysics
