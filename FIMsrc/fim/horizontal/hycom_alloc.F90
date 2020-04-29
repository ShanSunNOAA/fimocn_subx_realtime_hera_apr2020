module hycom_alloc
!	This module allocates variables used in ocean model
contains

subroutine hyc_alloc
use module_control, only: npp,nip
use hycom_variables
use hycom_constants,only: wet,odepth,dpth_edg,nuphill,uphill,dnhill,	&
	meshsz,thkmx,tfrez,theta,salmin,numflds,dpmin,huuge,land_spval
use hycom_control, only: numtr,thruflsecs,bdcurrsecs,nsecs,iocnmx,	&
	subcycle_ocn,trcfrq
use fimnamelist ,only: kdm,itest
use global_bounds, only: ims, ime

implicit none

integer :: i

   print *,'entering hyc_init ...'

print *,'allocating ocean arrays ... nip,kdm,numtr =',nip,kdm,numtr

! Allocate arrays defined in hycom_variables and hycom_constants

allocate (theta(kdm),salmin(kdm),dpmin(kdm))	! target densities, min saln, min layer thickness

allocate (utrop      (3,nip))		! barotropic u velocity
allocate (vtrop      (3,nip))		! barotropic v velocity
allocate (ptrop      (3,nip))		! barotropic pressure
allocate (ubforc       (nip))		! forcing of barotropic u velocity
allocate (vbforc       (nip))		! forcing of barotropic v velocity
allocate (dpth_edg (npp,nip))		! water depth on cell edges
allocate (btropfx  (npp,nip))		! forcing of barotropic mass flux 
allocate (ustarb       (nip))		! bottom friction velocity
allocate (hmixl        (nip))		! mixed layer depth (m)

allocate (pbot(nip))			! initial bottom pressure
allocate (gzbot(nip))			! bottom geopotential
allocate (odepth(nip),wet(nip))		! water depth, land(=0)/sea(>0) mask
allocate (meshsz(nip))			! scale for horiz.mesh size
allocate (thkmx(nip))			! max. ice thickness (m)
allocate (tfrez(nip))			! freezing temp
allocate (psikk(nip),spvkk(nip))	! bottom Montg.pot. & specif.vol.
allocate (dnhill(nip))			! lookup table of downhill neighbors
allocate (nuphill(nip))			! no.of cells draining into ipn

allocate (uclin    (kdm,nip,2))		! baroclinic u field
allocate (vclin    (kdm,nip,2))		! baroclinic v field
allocate (pres   (kdm+1,nip))		! interface pressure
allocate (geop   (kdm+1,nip))		! geopotential
allocate (dp       (kdm,nip,2))		! layer thickness
allocate (dpinit   (kdm,nip))		! initial lyr.thknss for trcr transport
allocate (dpfinl   (kdm,nip))		! final lyr.thknss for trcr transport
allocate (montg    (kdm,nip))		! Montgomery potential
allocate (dens     (kdm,nip))		! water density anomaly
allocate (spcifv   (kdm,nip))		! specific volume anomaly
allocate (temp     (kdm,nip))		! temperature
allocate (saln     (kdm,nip))		! salinity
allocate (passv_tr (kdm,nip,numtr))	! passive tracer
allocate (srforc   (4,nip,numflds))	! atmospheric forcings
allocate (temp_ave (kdm,nip))		! archvintvl-avg'd temperature
allocate (saln_ave (kdm,nip))		! archvintvl-avg'd salinity
allocate (dens_ave (kdm,nip))		! archvintvl-avg'd density
allocate (uvel_ave (kdm,nip))		! archvintvl-avg'd u-velocity
allocate (vvel_ave (kdm,nip))		! archvintvl-avg'd v-velocity
allocate (  dp_ave (kdm,nip))		! archvintvl-avg'd layer thickness
allocate (sstndcy  (kdm,nip))		! SST tdcy, up to kdm indiv.processes

allocate (um_edg   (kdm,npp,nip))	! total u vel.on cell edges, mid-time
allocate (vm_edg   (kdm,npp,nip))	! total v vel.on cell edges, mid-time
allocate (un_edg   (kdm,npp,nip))	! total u vel.on cell edges, old
allocate (vn_edg   (kdm,npp,nip))	! total v vel.on cell edges, old
allocate (dp_edg  ( kdm,npp,nip))	! lyr thknss on cell edges, mid-time
allocate (pr_edg (kdm+1,npp,nip))	! intfc pres on cell edges, mid-time
allocate (uvsq_edg (kdm,npp,nip))	! velocity-squared on cell edges
allocate (geop_edg (kdm,npp,nip))	! geopotential on cell edges, mid-lyr
allocate (mont_edg (kdm,npp,nip))	! Montgomery potential on cell edges
allocate (spcv_edg (kdm,npp,nip))	! specific volume on cell edges
allocate (massflx  (kdm,npp,nip))	! lateral mass flux on cell edges
allocate (cumuflx  (kdm,npp,nip))	! time integral of massflx
allocate (mssfx_ave(kdm,npp,nip))	! archvintvl-avg'd mass flux
allocate (uphill       (npp,nip))	! cells draining river water into ipn

allocate (srfx         (nip))		! surface heat flux
allocate (ssht         (nip))		! sea surface height
allocate (pmne         (nip))		! precip minus evap, positive down
allocate (ustar        (nip))		! friction velocity
allocate (taux         (nip))		! wind stress, u component
allocate (tauy         (nip))		! wind stress, v component
allocate (curl         (nip))		! wind stress curl
allocate (prcp         (nip))		! precipitation, pos.down; m/s
allocate (uwnd         (nip))		! u surface wind speed; m/s
allocate (vwnd         (nip))		! v surface wind speed; m/s
allocate (wspd         (nip))		! wind surface wind speed; m/s
allocate (airt         (nip))		! 2m air temperature
allocate (vpmx         (nip))		! Specific Humidity; kg/kg

allocate (curl_ave     (nip))		! archvintvl-avg'd wind stress curl
allocate (ssht_ave     (nip))		! archvintvl-avg'd sea surface height
allocate (srfx_ave     (nip))		! archvintvl-avg'd heat flux, pos down
allocate (pmne_ave     (nip))		! archvintvl-avg'd p-minus-e, pos down
allocate (qf2d_ave     (nip))		! archvintvl-avg'd latent heat
allocate (hf2d_ave     (nip))		! archvintvl-avg'd sensible heat
allocate (airt_ave     (nip))		! archvintvl-avg'd air temperature
allocate (vpmx_ave     (nip))		! archvintvl-avg'd vapor mixing ratio
allocate (prcp_ave     (nip))		! archvintvl-avg'd precip
allocate (swdn_ave     (nip))		! archvintvl-avg'd shortwave down; W/m^2
allocate (lwdn_ave     (nip))		! archvintvl-avg'd longwave down; W/m^2
allocate ( rad_ave     (nip))		! archvintvl-avg'd net radiation; W/m^2
allocate (wspd_ave     (nip))		! archvintvl-avg'd wind speed; m/s
allocate (taux_ave     (nip))		! archvintvl-avg'd eastwd wind stress
allocate (tauy_ave     (nip))		! archvintvl-avg'd northwd wind stress
allocate (hmixl_ave    (nip))		! archvintvl-avg'd mixed layer depth
allocate (temice_ave   (nip))		! archvintvl-avg'd ice temperature
allocate (covice_ave   (nip))		! archvintvl-avg'd ice coverage
allocate (thkice_ave   (nip))		! archvintvl-avg'd ice thickness

allocate (srfx_bcl     (nip))		! bclin_frq-avg'd heat flux, pos down
allocate (sw2d_bcl     (nip))		! bclin_frq-avg'd shortwave, pos down
allocate (pmne_bcl     (nip))		! bclin_frq-avg'd p-minus-e, pos down
allocate (ustar_bcl    (nip))		! bclin_frq-avg'd friction velocity
allocate (ustarb_bcl   (nip))		! bclin_frq-avg'd bottom friction vel.
allocate (rivflo       (nip))		! river water load
allocate (covice       (nip))		! sea ice coverage; rel.units
allocate (thkice       (nip))		! sea ice thickness; m
allocate (thksno       (nip))		! snow depth over ice; m
allocate (temice       (nip))		! ice surface temperature; C
allocate (ticeol       (nip))		! previous ice surface temperature; C
allocate (salnow       (nip))		! observed SSS for the present time step
allocate (airtmn       (nip))		! annual min. air tempeature
allocate (srfxcum      (nip))		! integral of surface heatflux; pos down (W/m2)
allocate (pmnecum      (nip))		! integral of P minus E; pos down (m/s)
allocate (tdpold       (nip))		! old SST*dp for tendency calculation

! kpp related:
allocate(jerlov(nip),akpar(nip))
!allocate(zgrid(kdm+1,nip),vcty(kdm+1,nip),difs(kdm+1,nip),dift(kdm+1,nip),  &
!               ghats(kdm+1,nip))
!allocate(hmonob(nip),dpbl(nip),dpbbl(nip),dpmold(nip),tmix(nip),smix(nip), &
!         thmix(nip),umix(nip),vmix(nip))
!allocate(buoflx(nip),bhtflx(nip),mixflx(nip),sswflx(nip))

!$OMP PARALLEL DO
    do i=ims,ime
      srfx(i)   = land_spval
      ssht(i)   = land_spval
      pmne(i)   = land_spval
      ustar(i)  = land_spval
      taux(i)   = land_spval
      tauy(i)   = land_spval
      curl(i)   = land_spval
      prcp(i)   = land_spval
      uwnd(i)   = land_spval
      vwnd(i)   = land_spval
      wspd(i)   = land_spval
      airt(i)   = land_spval
      vpmx(i)   = land_spval
      salnow(i) = land_spval
      airtmn(i) = huuge	! cannot use land_spval=-1.e33, as mininum is used to reset the initial ice temp
    end do
!$OMP END PARALLEL DO

return
end subroutine hyc_alloc

end module hycom_alloc
