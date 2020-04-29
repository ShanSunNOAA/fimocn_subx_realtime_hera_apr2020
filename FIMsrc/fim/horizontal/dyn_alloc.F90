!*********************************************************************
! module_dyn_alloc: allocates variables used in dynamics
!   Middlecoff         October       2008
!*********************************************************************
module module_dyn_alloc

contains

  subroutine dyn_alloc
    use module_control,only: nvlp1,nvlp,npp,nip,nabl,ntra,ntrb,nd,nvlsig
    use fimnamelist   ,only: nvl,readrestart
    use module_wrf_control,only: nbands
    use module_constants
    use module_variables
    use module_sfc_variables
    use global_bounds, only: ims, ime
!use module_chem_variables,only: sscal,ext_cof,asymp,extlw_cof
    use infnan, only: inf, negint

    implicit none

    integer :: ipn

! Allocate variables from module_constants

!..................................................................
!	Sec. 1  Math and Physics Constants
!..................................................................

!..................................................................
!	Sec. 2.  Grid Descriptive Variables
!..................................................................

    allocate(dpsig (nvl))           ! list of minimum layer thickness (Pa)
    allocate(thetac(nvl))           ! target theta for hybgen
    allocate(dfflen(nvl))           ! diffusion length scale (m)
    allocate(biharw(nvl))           ! Laplacian/biharm blending coefficient
    allocate(smoo_coeff(nvlp1))     ! interface smoothing

    if (readrestart) then
      allocate(sigak(nvlsig+1))           ! sigma def for GFS coord
      allocate(sigbk(nvlsig+1))           ! "
    end if

! velocity transform constants for projection from cell edges
    allocate(cs(4,npp,nip),sn(4,npp,nip))
    dpsig(:) = inf
    thetac(:) = inf
    dfflen(:) = inf
    smoo_coeff(:) = inf
!$OMP PARALLEL DO
    do ipn=ims,ime
      cs(:,:,ipn) = inf
      sn(:,:,ipn) = inf
    end do
!$OMP END PARALLEL DO

! Variables to describe the icos grid in xy (local stereographic)
    allocate(sidevec_c(nd,npp,nip)) ! side vectors projected from center
    allocate(sidevec_e(nd,npp,nip)) ! side vectors projected from edge
    allocate(sideln   (   npp,nip)) ! the length of side vectors (m)
    allocate(rsideln  (   npp,nip)) ! reciprocal of "sideln" (m**-1)
    allocate(rprox_ln (   npp,nip)) ! reciprocal of distance cell cent to prox pts
    allocate(area     (       nip)) ! the area of cell polygon (m**2)
    allocate(rarea    (       nip)) ! reciprocal of the "area"
    allocate(grndspd  (       nip)) ! ground speed due to earth rotation
    allocate(corner_xy(npp,nd,nip)) ! side vectors projected from edge
    allocate(hfcoef(6, npp,nip))  ! interpolation weights 

!$OMP PARALLEL DO
    do ipn=ims,ime
      sidevec_c(:,:,ipn) = inf
      sidevec_e(:,:,ipn) = inf
      sideln(:,ipn) = inf
      rsideln(:,ipn) = inf
      rprox_ln(:,ipn) = inf
      area(ipn) = inf
      rarea(ipn) = inf
      corner_xy(:,:,ipn) = inf
      hfcoef(:,:,ipn) = inf
    end do
!$OMP END PARALLEL DO

!.................................................................
!	Sec. 3.  Neighbor Lookup Tables etc.
!.................................................................

    allocate(prox       (npp,nip))  ! Holds index of proximity points
    allocate(proxs      (npp,nip))  ! Holds index of proximity sides
    allocate(nprox      (    nip))  ! Holds number of proximity points
    allocate(perm       (    nip))  ! translation table local -> global IPN
    allocate(inv_perm   (    nip))  ! inverse of 'perm'

    allocate (idxldg(npp,3,nip),idxtrl(npp,3,nip),			&
              wgtldg(npp,3,nip),wgttrl(npp,3,nip))

!$OMP PARALLEL DO
    do ipn=ims,ime
      prox(:,ipn) = negint
      proxs(:,ipn) = negint
      nprox(ipn) = negint
      inv_perm(ipn) = negint
    end do
!$OMP END PARALLEL DO

! nedge holds the number of edges valid at each grid cell on this task.  
! For a serial ! case, nedge == nprox.  
! For a parallel case, nedge == nprox on "interior" cells 
! and, nedge < nprox on "halo" cells.  
    allocate(nedge      (    nip))

! permedge stores a look-up table for edge indexes.  
! For a serial case, permedge does nothing:  
!   permedge(:,ipn) = 1, 2, 3, ... nprox(ipn)
! For a parallel case, permedge does nothing on "interior" cells.  
! For a parallel case, permedge skips "missing" edges on "halo" cells.  
    allocate(permedge   (npp,nip))

!$OMP PARALLEL DO
    do ipn=ims,ime
      nedge(ipn) = negint
      permedge(:,ipn) = negint
      perm(ipn) = negint
    end do
!$OMP END PARALLEL DO

!.....................................................................
!	Sec. 4.  Geo Variables:
!.....................................................................

allocate(corio(nip))			! Coriolis acceleration 
allocate(lat(nip),lon(nip))		! lat and lon in radians
allocate(deg_lat(nip),deg_lon(nip))	! lat and lon in degrees

!$OMP PARALLEL DO
    do ipn=ims,ime
      corio(ipn) = inf
      lat(ipn) = inf
      lon(ipn) = inf
      deg_lat(ipn) = inf
      deg_lon(ipn) = inf
    end do
!$OMP END PARALLEL DO

! Allocate variables from module_variables

!.....................................................................
!	Sec. 1.  3D Primary Variables
!.....................................................................
!  State variables at center point of cell for 3D grid:

!  Layer variables are defined in the middle of the layer

    allocate( us3d(nvl,nip))      ! zonal wind (m/s)
    allocate( vs3d(nvl,nip))      ! meridional wind (m/s)
    allocate( ws3d(nvl,nip))      ! vertical wind (Pa/s)
    allocate( dp3d(nvl,nip))      ! press.diff. between coord levels (Pa)
    allocate( dpinit(nvl,nip))    ! lyr thknss for class B tracer transport
    allocate( mp3d(nvl,nip))      ! Montgomery Potential (m^2/s^2)
    allocate( tk3d(nvl,nip))      ! temperature, kelvin
    allocate( dcudt(nvl,nip))     ! convective tendency from GF for GWD
    allocate( dcudq(nvl,nip))     ! convective tendency from GF for GWD
    allocate( dcudu(nvl,nip))     ! convective tendency from GF for GWD
    allocate( dcudv(nvl,nip))     ! convective tendency from GF for GWD
    allocate( relvor(nvl,nip))    ! relative vorticity (s^-1)
    allocate( potvor(nvl,nip))    ! potential vorticity
    allocate( pvsnap(  5,nip))    ! PV for SNAP project
    allocate( tr3d(nvl,nip,ntra+ntrb)) ! 1=pot.temp, 2=water vapor, 3=cloud water, 4=ozone
    allocate( trdp(nvl,nip,ntra+ntrb)) ! (tracer times dp3d ) for tracer transport eq.
    allocate( rh3d(nvl,nip))      ! relative humidity from 0 to 1
    allocate( qs3d(nvl,nip))      ! saturation specific humidity
    allocate( pw2d(nip))          ! precipitable water
    allocate( pq2d(nip))          ! vertically integrated hydrometeor condensate

    ALLOCATE(wt_int_rev(nvl,nip)) ! inverse weights
    ALLOCATE(k_int_rev(nvl,nip))  ! k levs for inverse weights 

!$OMP PARALLEL DO
    do ipn=ims,ime
      us3d(:,ipn) = inf
      vs3d(:,ipn) = inf
      ws3d(:,ipn) = inf
      dp3d(:,ipn) = inf
      dpinit(:,ipn) = inf
      mp3d(:,ipn) = inf
      tk3d(:,ipn) = inf
      dcudt(:,ipn) = 0.
      dcudq(:,ipn) = 0.
      dcudu(:,ipn) = 0.
      dcudv(:,ipn) = 0.
      relvor(:,ipn) = inf
      potvor(:,ipn) = inf
      tr3d(:,ipn,:) = inf
      trdp(:,ipn,:) = inf
      rh3d(:,ipn) = inf
      qs3d(:,ipn) = inf
      pw2d(ipn) = inf
      wt_int_rev(:,ipn) = inf
      k_int_rev(:,ipn) = -1
    end do
!$OMP END PARALLEL DO

!  Level variables defined at layer interfaces
    allocate( pr3d(nvlp1,nip))    ! pressure (pascal)
    allocate( ex3d(nvlp1,nip))    ! exner function
    allocate( ph3d(nvlp1,nip))    ! geopotential (=gz), m^2/s^2
    allocate( sdot(nvlp1,nip))    ! mass flux across interfaces, sdot*(dp/ds)

!$OMP PARALLEL DO
    do ipn=ims,ime
      pr3d(:,ipn) = inf
      ex3d(:,ipn) = inf
      ph3d(:,ipn) = inf
      sdot(:,ipn) = inf
    end do
!$OMP END PARALLEL DO

!..................................................................
!	Sec. 2. Edge Variables
!..................................................................
!  Variables carried at the midpoints of the 6(5) sides of each cell
    allocate( u_edg   (nvl,npp,nip))    ! u on edge
    allocate( v_edg   (nvl,npp,nip))    ! v on edge
    allocate( dp_edg  (nvl,npp,nip))    ! dp on edge
    allocate( trc_edg (nvl,npp,nip,ntra+ntrb))! tracers on edge
    allocate( lp_edg  (nvl,npp,nip))    ! mid-layer pressure on edge
    allocate( geop_edg(nvl,npp,nip))    ! geopotential on edge
    allocate( mont_edg(nvl,npp,nip))    ! montg.potential on edge
    allocate( uvsq_edg(nvl,npp,nip))    ! velocity-squared on edge
    allocate( massfx  (nvl,npp,nip,3))  ! mass fluxes on edge
    allocate( cumufx  (nvl,npp,nip))    ! time-integrated mass flx on edge
    allocate( flxavg  (nvl,npp,nip))    ! time-integrated flx across transects

!$OMP PARALLEL DO
    do ipn=ims,ime
      u_edg(:,:,ipn) = inf
      v_edg(:,:,ipn) = inf
      dp_edg(:,:,ipn) = inf
      trc_edg(:,:,ipn,:) = inf
      lp_edg(:,:,ipn) = inf
      geop_edg(:,:,ipn) = inf
      mont_edg(:,:,ipn) = inf
      uvsq_edg(:,:,ipn) = inf
      massfx(:,:,ipn,:) = inf
      cumufx(:,:,ipn) = inf
    end do
!$OMP END PARALLEL DO

!.....................................................................
! 	Sec. 3. Forcing (tendency) Variables
!.....................................................................
    allocate( u_tdcy  (nvl,nip,nabl))     ! forcing of u
    allocate( v_tdcy  (nvl,nip,nabl))     ! forcing of v
    allocate( dp_tdcy (nvl,nip,nabl))     ! forcing of dp
    allocate( dpl_tdcy(nvl,nip,nabl))     ! forcing dp, low order
    allocate( trc_tdcy(nvl,nip,nabl,ntra+ntrb))! forcing of tracers
    allocate( trl_tdcy(nvl,nip,nabl,ntra+ntrb))! forcing of tracers, low order

!$OMP PARALLEL DO
    do ipn=ims,ime
      u_tdcy(:,ipn,:) = inf
      v_tdcy(:,ipn,:) = inf
      dp_tdcy(:,ipn,:) = inf
      dpl_tdcy(:,ipn,:) = inf
      trc_tdcy(:,ipn,:,:) = inf
      trl_tdcy(:,ipn,:,:) = inf
    end do
!$OMP END PARALLEL DO

!....................................................................
!       Sec. 4. Misc. arrays
!....................................................................
    allocate(  work2d(nip  ))
    allocate( iwork2d(nip  ))
    allocate( psrf   (nip  ))   ! surface pressure
    allocate(th_pvsrf(nip  ))   ! pot.temperature on pot.vorticity surface
    allocate(pr_pvsrf(nip  ))   ! pressure on pot.vorticty surface
    allocate(us_pvsrf(nip  ))   ! u component on pot.vorticty surface
    allocate(vs_pvsrf(nip  ))   ! v component on pot.vorticty surface
    allocate( ptdcy  (nip,2))   ! sfc.pres.tdcy at 2 consec. time levels
    allocate( conv_act(nip))    ! counter for conv activity over last  time steps
    allocate( worka  (nvl,nip)) ! 3d work array
    allocate( workb  (nvl,nip)) ! 3d work array
    allocate( q_min_out (nvl,nip)) ! 3d work array
    allocate( q_max_out (nvl,nip)) ! 3d work array
    allocate( grad(nvlp1,2,nip))! work array to store gradient of variables
    allocate( diaga(nvl,nip))   ! diagnostic, for output, fill with anything
    diaga=0.
    allocate( diagb(nvl,nip))   ! diagnostic, for output, fill with anything
    diagb=0.

!TODO Initialize diag* arrays to inf instead of zero. Init to zero currently required or
!TODO da and db fields will be wrong on history files.
!$OMP PARALLEL DO
    do ipn=ims,ime
      work2d(ipn) = inf
      iwork2d(ipn) = negint
      psrf(ipn) = inf
      ptdcy(ipn,:) = inf
      worka(:,ipn) = inf
      workb(:,ipn) = inf
      q_min_out(:,ipn) = 0.
      q_max_out(:,ipn) = 0.
      grad(:,:,ipn) = 0.
    end do
!$OMP END PARALLEL DO

! Allocate variants of physics variables from module module_sfc_variables
! These are single-precision copies of physics variables passed from physics 
! to dynamics via the coupler.  
!JR These 5 things were moved from output.F90 so they can be written to the restart file.
    allocate(rsds(nip),rsds_ave(nip))   ! radiation short-wave downward at surface, snapshot/averaged over ArchvIntvl (W/m2)
    allocate(rlds(nip),rlds_ave(nip))   ! radiation  long-wave downward at surface, snapshot/averaged over ArchvIntvl (W/m2)
    allocate(rsus(nip),rsus_ave(nip))   ! radiation short-wave   upward at surface, snapshot/averaged over ArchvIntvl (W/m2)
    allocate(rlus(nip),rlus_ave(nip))   ! radiation  long-wave   upward at surface, snapshot/averaged over ArchvIntvl (W/m2)
    allocate(rsdt(nip),rsdt_ave(nip))   ! radiation short-wave downward        TOA, snapshot/averaged over ArchvIntvl (W/m2)
    allocate(rsut(nip),rsut_ave(nip))   ! radiation short-wave   upward        TOA, snapshot/averaged over ArchvIntvl (W/m2)
    allocate(rlut(nip),rlut_ave(nip))   ! radiation  long-wave   upward        TOA, snapshot/averaged over ArchvIntvl (W/m2)
    allocate(qf2d(nip),  qf_ave(nip))	! latent heatflux, snapshot/averaged over ArchvIntvl (W/m2)
    allocate(hf2d(nip),  hf_ave(nip))	! sensible heatflux, snapshot/averaged over ArchvIntvl (W/m2)
    allocate (cld_hi(nip),cld_md(nip),cld_lo(nip),cld_tt(nip),cld_bl(nip))   ! high/mid/low/total/boundary cloud cover
    allocate (cld_hi_ave(nip),cld_md_ave(nip),cld_lo_ave(nip),cld_tt_ave(nip),cld_bl_ave(nip)) ! averaged over ArchvIntvl
    allocate(rn2d0(nip))
    allocate(rc2d0(nip))
    allocate(rg2d0(nip))
    allocate(rn2d6_0(nip))
    allocate(rc2d6_0(nip))
    allocate(rg2d6_0(nip))
    allocate(rn2d(nip))		! accumulated total precipitation/rainfall
    allocate(sn2d(nip))		! accumulated total snowfall
    allocate(rc2d(nip))		! accumulated convective precipitation/rainfall 
    allocate(rain(nip))		! m/step
    allocate(runoff(nip))	! m/step
    allocate(ts2d(nip))		! skin temperature
    allocate(sst_prev(nip))	! skin temperature
    allocate(sst_next(nip))	! skin temperature
    allocate(us2d(nip))		! friction velocity/equivalent momentum flux
    allocate(sheleg2d(nip))
    allocate(tg32d(nip))
    allocate(zorl2d(nip))
    allocate(vfrac2d(nip))
    allocate(vtype2d(nip))
    allocate(stype2d(nip))
    allocate(cv2d(nip))
    allocate(cvb2d(nip))
    allocate(cvt2d(nip))
    allocate(alvsf2d(nip))
    allocate(alvwf2d(nip))  
    allocate(alnsf2d(nip))
    allocate(alnwf2d(nip))
    allocate(f10m2d(nip))
    allocate(ts_ave(nip),t2_ave(nip),q2_ave(nip),td_ave(nip),sm_ave(nip),runoff_ave(nip),slp_ave(nip),wsp80_ave(nip))
    allocate(ustr_ave(nip),vstr_ave(nip),t2max(nip),t2min(nip))
    allocate(u10m(nip),v10m(nip),u10_ave(nip),v10_ave(nip),g3p_ave(nvlp,nip,7))
    allocate(facsf2d(nip))
    allocate(facwf2d(nip))
    allocate(uustar2d(nip))
    allocate(alb2d(nip))
    allocate(strs2d(nip))
    allocate(ffmm2d(nip))
    allocate(ffhh2d(nip))
    allocate(srflag2d(nip))
    allocate(snwdph2d(nip))
    allocate(shdmin2d(nip))
    allocate(shdmax2d(nip))
    allocate(slope2d(nip))
    allocate(snoalb2d(nip))
    allocate(canopy2d(nip))
    allocate(hice2d(nip),hice_ave(nip))
    allocate(fice2d(nip),fice_ave(nip))
    allocate(fice2d_prev(nip))
    allocate(fice2d_next(nip))
    allocate(slc2d(nip))
    allocate(t2m2d(nip))
    allocate(q2m2d(nip))
    allocate(slmsk2d(nip))
    allocate(tprcp2d(nip))  ! precip rate (1000*kg/m**2)
    allocate(st3d(4,nip))   ! soil temperature
    allocate(sm3d(4,nip))   ! total soil moisture
    allocate(slc3d(4,nip))  ! soil liquid content
    allocate(hprm2d(14,nip))

!$OMP PARALLEL DO
    do ipn=ims,ime
! replace inf by 0 to prevent inf in the initial output when physics is NOT called yet
      rsds(ipn) = 0.; rlds(ipn) = 0.; rsus(ipn) = 0.; rlus(ipn) = 0. 
      rsdt(ipn) = 0.; rsut(ipn) = 0.; rlut(ipn) = 0.; runoff(ipn) = 0.
      qf2d(ipn) = 0.; hf2d(ipn) = 0.
      cld_hi(ipn) = 0.; cld_md(ipn) = 0.; cld_lo(ipn) = 0.; cld_tt(ipn) = 0.; cld_bl(ipn) = 0.

      rsds_ave(ipn) = 0.; rlds_ave(ipn) = 0.; rsus_ave(ipn) = 0.; rlus_ave(ipn) = 0. 
      rsdt_ave(ipn) = 0.; rsut_ave(ipn) = 0.; rlut_ave(ipn) = 0.
        qf_ave(ipn) = 0.; hf_ave(ipn) = 0.
      cld_hi_ave(ipn) = 0.; cld_md_ave(ipn) = 0.; cld_lo_ave(ipn) = 0.; cld_tt_ave(ipn) = 0.; cld_bl_ave(ipn) = 0.
      conv_act(ipn)=0

      rn2d0(ipn) = inf
      rc2d0(ipn) = inf
      rg2d0(ipn) = inf
      rn2d6_0(ipn) = inf
      rc2d6_0(ipn) = inf
      rg2d6_0(ipn) = inf
      rn2d(ipn) = inf
      sn2d(ipn) = inf
      rc2d(ipn) = inf
      rain(ipn) = inf
      ts2d(ipn) = inf
      sst_prev(ipn) = inf
      sst_next(ipn) = inf
      us2d(ipn) = inf
      sheleg2d(ipn) = inf
      tg32d(ipn) = inf
      zorl2d(ipn) = inf
      vfrac2d(ipn) = inf
      vtype2d(ipn) = inf
      stype2d(ipn) = inf
      cv2d(ipn) = inf
      cvb2d(ipn) = inf
      cvt2d(ipn) = inf
      alvsf2d(ipn) = inf
      alvwf2d(ipn) = inf
      alnsf2d(ipn) = inf
      alnwf2d(ipn) = inf
      f10m2d(ipn) = inf
      u10m(ipn) = inf
      v10m(ipn) = inf
      facsf2d(ipn) = inf
      facwf2d(ipn) = inf
      uustar2d(ipn) = inf
      alb2d(ipn) = inf
      strs2d(ipn) = inf
      ffmm2d(ipn) = inf
      ffhh2d(ipn) = inf
      srflag2d(ipn) = inf
      snwdph2d(ipn) = inf
      shdmin2d(ipn) = inf
      shdmax2d(ipn) = inf
      slope2d(ipn) = inf
      snoalb2d(ipn) = inf
      canopy2d(ipn) = inf
      hice2d(ipn) = inf
      hice_ave(ipn) = inf
      fice2d(ipn) = inf
      fice_ave(ipn) = inf
      fice2d_prev(ipn) = inf
      fice2d_next(ipn) = inf
      slc2d(ipn) = inf
      t2m2d(ipn) = inf
      q2m2d(ipn) = inf
      slmsk2d(ipn) = inf
      tprcp2d(ipn) = inf
      st3d(:,ipn) = inf
      sm3d(:,ipn) = inf
      slc3d(:,ipn) = inf
      hprm2d(:,ipn) = inf
      cld_hi    (ipn) = inf
      cld_md    (ipn) = inf
      cld_lo    (ipn) = inf
      cld_tt    (ipn) = inf
      cld_bl    (ipn) = inf
      cld_hi_ave(ipn) = inf
      cld_md_ave(ipn) = inf
      cld_lo_ave(ipn) = inf
      cld_tt_ave(ipn) = inf
      cld_bl_ave(ipn) = inf
    end do
!$OMP END PARALLEL DO

! chem arrays

! TODO We'd rather not allocate these arrays if chem is off. But Lahey does not
! TODO allow passing unallocated arrays as subroutine parameters. So, work needs
! TODO to be done to (maybe) have these arrays as optional parameters and wrap
! TODO accesses in "if (present())" conditionals -- or some other approach. For
! TODO now, we have to allocate these.

!if(aer_ra_feedback == 1 ) then
    allocate(sscal(nvl,nip,nbands))
    allocate(ext_cof(nvl,nip,nbands))
    allocate(asymp(nvl,nip,nbands))
    allocate(extlw_cof(nvl,nip,16))
    allocate(aod2d(nip))      ! aerosol optical depth
!$OMP PARALLEL DO
    do ipn=ims,ime
      sscal(:,ipn,:)     = 1.
      ext_cof(:,ipn,:)   = 0.
      asymp(:,ipn,:)     = 0.
      extlw_cof(:,ipn,:) = 0.
      aod2d(ipn)         = 0.
    end do
!$OMP END PARALLEL DO
!endif

    return
  end subroutine dyn_alloc
end module module_dyn_alloc
