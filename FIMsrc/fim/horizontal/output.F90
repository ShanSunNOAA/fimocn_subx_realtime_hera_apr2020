module module_output

  use module_constants
  use module_control,           only: dt,nip,ntra,ntrb,nvar2d,nvarsig,nvarp,    &
                                      nvlsig,nvlp,nvlp1,ArchvStep,hrs_in_month, &
                                      itsDFI,ArchvStep6h,ArchvStep0
  use fimnamelist,              only: FixedGridOrder,nvl,ArchvTimeUnit,         &
                                      PrintDiags,restart_freq,yyyymmddhhmm,     &
                                      enkfio_out,readrestart,itsstart,digifilt, &
                                      cycle_freq,enkf_diag,ArchvIntvl,atmonly
  use module_variables,         only: potvor,th_pvsrf,pr_pvsrf,us_pvsrf,	&
				      vs_pvsrf,pvsnap,flxavg
  use module_core_setup,        only: use_write_tasks
  use module_op_diag,           only: op_diag
  use module_outFMTed,          only: outFMTed
  use module_outqv,             only: outqv
  use module_outqv_mn_lat,      only: outqv_mn_lat
  use module_outqv_mn_lat_land, only: outqv_mn_lat_land
  use module_outvar_enkf_spectral,       only: outvar_enkf_spectral
  use module_outvar_enkf_sfcio,          only: outvar_enkf_sfcio
  use module_enkf_io, only : write_enkfio
  use module_printMAXMIN,       only: printMAXMIN
  use findmaxmin1
  use mdul_pvsurf
  use icosio_wrapper,           only: maybe_write
  use module_wrf_control,       only: num_moist,num_chem
  use global_bounds,            only: ims, ime, ips, ipe
  use module_sfc_variables,     only: rsds_ave,rlds_ave,rsus_ave,rlus_ave,rsdt_ave,rsut_ave,rlut_ave,	&
				qf_ave,hf_ave, rn2d6_0,rc2d6_0,rg2d6_0,ts_ave,t2_ave,q2_ave,		&
				td_ave,sm_ave,runoff_ave,slp_ave,wsp80_ave,g3p_ave,ustr_ave,		&
				vstr_ave,u10_ave,v10_ave,fice2d,fice_ave,hice_ave,runoff,t2max,t2min,	&
				cld_hi,cld_md,cld_lo,cld_tt,cld_bl,					&
				cld_hi_ave,cld_md_ave,cld_lo_ave,cld_tt_ave,cld_bl_ave,strs2d
  use stencilprint
  use mdul_transects,           only: atm_flxsum

  implicit none

#include <gptl.inc>

contains

!*********************************************************************
!	output
!	Output program for fim global model
!	Alexander E. MacDonald  12/27/2004
!	J. Lee                  September, 2005
!*********************************************************************

subroutine output (its, nts,				& ! index time step, final timestep
                   us3d, vs3d, dp3d, sdot,		& ! west wind, south wind, delta pres, vert.velocity
                   pr3d, ex3d, mp3d,			& ! pressure, Exner, mont pot,
                   tr, rh3d, relvor, ws3d,		& ! tracers, humidity, vorticity, omega
                   chem_opt,diaga, diagb,		& ! diagnostic arrays
                   ph3d, tk3d, rn2d, sn2d, rc2d, pw2d, pq2d,  &
                   ts2d, us2d, rsds, rlds, rsus, rlus,	&
                   rsdt, rsut, rlut, qf2d, hf2d,	&
                   st3d, sm3d, t2m2d, q2m2d, canopy2d,	& ! geopotential, accumulated precip/rainfall
                   fice2d, hice2d, sheleg2d, slmsk2d,	&
		   u10m, v10m, rn2d0, rc2d0, rg2d0,	&
                   TimingBarriers, curr_write_time,	&
                   wt_int_rev,k_int_rev)

! Arguments
    integer, intent(in) :: its      ! current time step
    integer, intent(in) :: nts      ! final time step
    integer, intent(in) :: chem_opt ! for chem pressure level files we need to know chem_opt
    real, intent(in)    :: us3d(nvl,  ims:ime)
    real, intent(in)    :: vs3d(nvl,  ims:ime)
    real, intent(in)    :: dp3d(nvl,  ims:ime)
    real, intent(in)    :: sdot(nvlp1,ims:ime)
    real, intent(in)    :: pr3d(nvlp1,ims:ime)
    real, intent(in)    :: mp3d(nvl,  ims:ime)
    real, intent(in)    :: diaga(nvl, ims:ime)
    real, intent(in)    :: diagb(nvl, ims:ime)
    real, intent(in)    :: relvor(nvl,ims:ime)
    real, intent(in)    :: ws3d(nvl  ,ims:ime)
    real, intent(in)    :: ph3d(nvlp1,ims:ime)
    real, intent(in)    :: ex3d(nvlp1,ims:ime)
    real, intent(in)    :: tr(nvl,    ims:ime,ntra+ntrb)
    real, intent(inout) :: rh3d(nvl,  ims:ime)
    real, intent(inout) :: tk3d(nvl,  ims:ime)
    real, intent(in)    :: rn2d(ims:ime)
    real, intent(inout) :: sn2d(ims:ime)
    real, intent(in)    :: rc2d(ims:ime)
    real, intent(in)    :: u10m(ims:ime)
    real, intent(in)    :: v10m(ims:ime)
!JR Moved these 5 things to arg list so they can be written to the restart file.
    real, intent(inout) :: rn2d0(ims:ime)
    real, intent(inout) :: rc2d0(ims:ime)
    real, intent(inout) :: rg2d0(ims:ime)
    real, intent(inout) :: pw2d(ims:ime)
    real, intent(inout) :: pq2d(ims:ime)
    real, intent(in)    :: ts2d(ims:ime)
    real, intent(in)    :: us2d(ims:ime)
    real, intent(in)    :: hf2d(ims:ime)
    real, intent(in)    :: qf2d(ims:ime)
    real, intent(in)    :: rsds(ims:ime)
    real, intent(in)    :: rlds(ims:ime)
    real, intent(in)    :: rsus(ims:ime)
    real, intent(in)    :: rlus(ims:ime)
    real, intent(in)    :: rsdt(ims:ime)
    real, intent(in)    :: rsut(ims:ime)
    real, intent(in)    :: rlut(ims:ime)
    real, intent(in)    :: st3d(4,ims:ime)
    real, intent(in)    :: sm3d(4,ims:ime)
    real, intent(in)    :: t2m2d(ims:ime)
    real, intent(in)    :: q2m2d(ims:ime)
!JR ADDED ITEMS FROM MODULE_SFC_VARIABLES SO EVERYTHING COMES FROM INPUT ARG LIST
!JR RATHER THAN SOME USED FROM MODULE
    real, intent(in) :: canopy2d(ims:ime)
    real, intent(in) :: fice2d(ims:ime)
    real, intent(in) :: hice2d(ims:ime)
    real, intent(in) :: sheleg2d(ims:ime)
    real, intent(in) :: slmsk2d(ims:ime)

    logical, intent(in)    :: TIMINGBARRIERS
    integer, intent(inout) :: CURR_WRITE_TIME ! MOST RECENT TIME VARS. WERE WRITTEN

    real, intent(inout) :: wt_int_rev(nvl,ims:ime)! weights for reverse interpolation from
! gfs sig levels back to hybrid levels
! for fim ensemble data assimila\tion
    integer, intent(inout) :: k_int_rev(nvl,ims:ime)! k levs for reverse interpolation from 
!GFS sig levels back to hybrid levels

! LOCAL VARIABLES
    real :: td3d(nvl,ims:ime)
    real :: mslp(ims:ime)
    real :: rg2d(ims:ime)
    real :: rn_xh(ims:ime), rc_xh(ims:ime), rg_xh(ims:ime), sn_xh(ims:ime)
    real :: rn_6h(ims:ime), rc_6h(ims:ime), rg_6h(ims:ime), sn_6h(ims:ime)
    real :: g3p(nvlp,ims:ime,nvarp)
    !      (nvarp=6 + NUM_CHEM)
    !      1=HEIGHT,2=TEMP,3=RH (W.R.T. WATER),4=U WIND,5=V WIND
    real :: g3p_chem(nvlp,ims:ime,num_chem+2)
!   real :: g3p_chem(nvlp,ims:ime,num_chem)
    !      (nvarp=6 + NUM_CHEM)
    !      1=HEIGHT,2=TEMP,3=RH (W.R.T. WATER),4=U WIND,5=V WIND
    real, allocatable :: g3sig (:,:,:)
    real :: g2d(ims:ime,nvar2d)
    !       ADDITIONAL DIAGNOSTIC 2D VARIABLES FROM OP_DIAG.F90
    real :: spd10m_dif(ims:ime)
    real :: t2m_dif(ims:ime)
    real :: dum(ims:ime)

    integer :: time, ipn, laststep,ichemstart
    integer :: accum_start     ! value of "time" from previous output call
    integer :: accum_start_6h  ! value of "time" 6 hours ago 
    real*8  :: toutput, tmaxmin
    real :: es,esln,pkap,wspd,iextent
    real :: icex_nh, icev_nh, icex_sh, icev_sh  ! NH_ice_ext/NH_ice_vol/SH_ice_ext/SH_ice_vol

    integer :: ret,k
    logical, parameter:: ave_intvl=.true.	! T: write fluxes averaged over "ArchvIntvl"
    logical           :: ave6hr=.false.         ! T: write averages over 6hr
    logical	      :: output2d=.false.	! T: write 2d & g3p fields only to save diskspace
    logical	      :: dayave9z		! T: snapshot daily ave is done by 0z,6z,12z and 18z
    character*80 stamp

    real, external :: its2time

    allocate (g3sig(nvlsig+1,ims:ime,nvarsig))
    output2d=.false.
    dayave9z=.false.
!   if (ArchvTimeUnit=='hr' .or. ArchvTimeUnit=='dy') output2d=.true.
    if (ArchvTimeUnit=='hr' .and. ArchvIntvl == 6) ave6hr=.false.	! no need to repeat 
    if (ArchvTimeUnit=='hr' .and. ArchvIntvl == 24) dayave9z=.true.
    time = nint(its2time(its))

! there is extra varibles in nvarp, so chem starts one later than ntra!!
! need to recover original nvarp, assuming there is some sort of chemistry in this run
! so FROM FIMsrc/cntl/module_wrf_control
!   if(num_moist.lt.4)then
!    nvarp=nvarp+num_chem+num_moist - 3
! else
!    nvarp=nvarp+num_chem+num_moist - 2       ! nvarp includes all variables for pressure level output
! endif
! THEREFORE:

!    ichemstart=(nvarp -num_chem-num_moist +3) +num_moist - 3 ! as of sep 2012, this should be 6
!   if (num_moist > 3) then
!     ichemstart=(nvarp-num_chem-num_moist +2) + num_moist - 2
!   end if

!SB - 1/26/2014 - chem vars on isobaric levels are now in separate g3p_chem array
!      to avoid awkward indexing
    ichemstart=0

!   write(stamp,'(a,i4)')  ' step=',its
!     call stencl(cld_bl_ave,1,1000.,'clbl output_write'//trim(stamp))

    if (its-ArchvStep0 == itsstart-1 .and. its-ArchvStep0 >1) then

!sms$ignore begin
!$OMP PARALLEL DO schedule (static)
      do ipn=ips,ipe
          rsds_ave(ipn) = 0.; rlds_ave(ipn) = 0.; rsus_ave(ipn) = 0.; rlus_ave(ipn) = 0.
          rsdt_ave(ipn) = 0.; rsut_ave(ipn) = 0.; rlut_ave(ipn) = 0.
            qf_ave(ipn) = 0.;   hf_ave(ipn) = 0.
            ts_ave(ipn) = 0.;   t2_ave(ipn) = 0.;   q2_ave(ipn) = 0.; fice_ave(ipn) = 0.; hice_ave(ipn) = 0.
            td_ave(ipn) = 0.;   sm_ave(ipn) = 0.
       runoff_ave (ipn) = 0.
         wsp80_ave(ipn) = 0.
          ustr_ave(ipn) = 0.; vstr_ave(ipn) = 0.
           u10_ave(ipn) = 0.;  v10_ave(ipn) = 0.
        cld_hi_ave(ipn) = 0.; cld_md_ave(ipn) = 0.; cld_lo_ave(ipn) = 0.; cld_tt_ave(ipn) = 0.; cld_bl_ave(ipn) = 0.
             t2min(ipn) =  9999.; t2max(ipn) = -9999.
         g3p_ave(:,ipn,:) = 0.
           slp_ave(ipn) = 0.
	   rn_xh(ipn) = 0.; rc_xh(ipn) = 0.; rg_xh(ipn) = 0.; sn_xh(ipn) = 0.; !sn2d(ipn) = 0.
	   rn_6h(ipn) = 0.; rc_6h(ipn) = 0.; rg_6h(ipn) = 0.; sn_6h(ipn) = 0.
             rn2d0(ipn) = rn2d(ipn); rc2d0(ipn) = rc2d(ipn);	!no need to save rn2d0.rc2d0 in the restart file; note rg2d is local
             rg2d0(ipn) = rn2d(ipn) - rc2d(ipn)
           rn2d6_0(ipn) = rn2d(ipn); rc2d6_0(ipn)=rc2d(ipn);
           rg2d6_0(ipn) = rg2d0(ipn)
      end do
!$OMP END PARALLEL DO
!sms$ignore end
      return	! cannot (and no need to) rewrite binary at itsstart due to the lack of avg field in restart
    end if


    if (its-ArchvStep0 == 0) then

      print *,'chk ini its,ArchvStep0,ArchvStep,ArchvStep6h=',its,ArchvStep0,ArchvStep,ArchvStep6h
      if (ArchvStep < ArchvStep6h) stop 'wrong ArchvStep: ArchvStep cannot be shorter than ArchvStep6h'

!sms$ignore begin
!$OMP PARALLEL DO schedule (static)
      do ipn=ips,ipe
! load 0hr snapshot as values at 0hr
	rsds_ave(ipn) = 0.; rlds_ave(ipn) = 0.; rsus_ave(ipn) = 0.; rlus_ave(ipn) = 0.
	rsdt_ave(ipn) = 0.; rsut_ave(ipn) = 0.; rlut_ave(ipn) = 0.
	  qf_ave(ipn) = 0.; hf_ave(ipn) = 0.
          ts_ave(ipn) = 0.; t2_ave(ipn) = 0.; q2_ave(ipn) = 0.; fice_ave(ipn) = 0.; hice_ave(ipn) = 0.
          td_ave(ipn) = 0.; sm_ave(ipn) = 0.
     runoff_ave (ipn) = 0.
       wsp80_ave(ipn) = 0.
	ustr_ave(ipn) = 0.; vstr_ave(ipn) = 0.
         u10_ave(ipn) = 0.; v10_ave(ipn) = 0.
	   rn_xh(ipn) = 0.; rc_xh(ipn) = 0.; rg_xh(ipn) = 0.; sn_xh(ipn) = 0.; sn2d(ipn) = 0.
	   rn_6h(ipn) = 0.; rc_6h(ipn) = 0.; rg_6h(ipn) = 0.; sn_6h(ipn) = 0.
	   rn2d0(ipn) = 0.; rc2d0(ipn) = 0.; rg2d0(ipn) = 0.
         rn2d6_0(ipn) = 0.; rc2d6_0(ipn)=0.; rg2d6_0(ipn)=0.
	   t2max(ipn) = -9999.; t2min(ipn) =  9999.
       g3p_ave(:,ipn,:) = 0.
	 slp_ave(ipn) = 0.
      cld_hi_ave(ipn) = 0.; cld_md_ave(ipn) = 0.; cld_lo_ave(ipn) = 0.; cld_tt_ave(ipn) = 0.; cld_bl_ave(ipn) = 0.
	     dum(ipn) = 0.
      end do
!$OMP END PARALLEL DO
!sms$ignore end
      open(429,file='ice_hemi.txt',form='formatted',status='unknown')
      icex_nh=0.; icev_nh=0.; icex_sh=0.; icev_sh=0.; 
    end if	! its-ArchvStep0 == 0
!
!
      if (ave_intvl) then	
! need to skip the do loop below on the 1st time step after restart
!sms$ignore begin
!$OMP PARALLEL DO schedule (static)
      do ipn=ips,ipe
! --- accumulate fields averaged over "ArchvIntvl"
        rsds_ave(ipn) =   rsds_ave(ipn) + rsds(ipn)
        rlds_ave(ipn) =   rlds_ave(ipn) + rlds(ipn)
        rsus_ave(ipn) =   rsus_ave(ipn) + rsus(ipn)
        rlus_ave(ipn) =   rlus_ave(ipn) + rlus(ipn)
        rsdt_ave(ipn) =   rsdt_ave(ipn) + rsdt(ipn)
        rsut_ave(ipn) =   rsut_ave(ipn) + rsut(ipn)
        rlut_ave(ipn) =   rlut_ave(ipn) + rlut(ipn)
        fice_ave(ipn) =   fice_ave(ipn) + fice2d(ipn)	! ice fraction     over ArchvTimeUnit
        hice_ave(ipn) =   hice_ave(ipn) + hice2d(ipn)	! ice thickness    over ArchvTimeUnit
          qf_ave(ipn) =     qf_ave(ipn) + qf2d(ipn)
          hf_ave(ipn) =     hf_ave(ipn) + hf2d(ipn)
          ts_ave(ipn) =     ts_ave(ipn) + ts2d(ipn)	! t_skin  averaged over ArchvTimeUnit
          t2_ave(ipn) =     t2_ave(ipn) + t2m2d(ipn)	! 2m temp averaged over ArchvTimeUnit
          q2_ave(ipn) =     q2_ave(ipn) + q2m2d(ipn)	! 2m humidity      over ArchvTimeUnit
            pkap  = 0.5*(ex3d(1,ipn) + ex3d(2,ipn))/cp
            es    = (p1000*pkap**(cp/rd)-20.)*(q2m2d(ipn) + 1.e-8)/(ratio_h20_dry+(q2m2d(ipn) + 1.e-8))	! pres(1)-20. => 2m
            esln  = log(abs(es) + 1.e-6)
          td_ave(ipn) = td_ave(ipn)+(35.86*esln - 4947.2325)/(esln - 23.6837)   ! Teten's equation
          sm_ave(ipn) =     sm_ave(ipn) + (sm3d(1,ipn)*0.1+sm3d(2,ipn)*0.3+sm3d(3,ipn)*0.6+sm3d(4,ipn)*1.0)/2.0
      runoff_ave(ipn) = runoff_ave(ipn) + runoff(ipn)*area(ipn)/dt	! m3/sec
                 wspd = sqrt(us3d(1,ipn)**2+vs3d(1,ipn)**2)		! m/s
       if (wspd.gt.1.e-5) then                 ! scaled by ice coverage
        ustr_ave(ipn) =   ustr_ave(ipn) + strs2d(ipn)*us3d(1,ipn)*(1.-.9*fice2d(ipn))/wspd
        vstr_ave(ipn) =   vstr_ave(ipn) + strs2d(ipn)*vs3d(1,ipn)*(1.-.9*fice2d(ipn))/wspd
       end if
         u10_ave(ipn) =    u10_ave(ipn) +   u10m(ipn)
         v10_ave(ipn) =    v10_ave(ipn) +   v10m(ipn)
      cld_hi_ave(ipn) = cld_hi_ave(ipn) + cld_hi(ipn)
      cld_md_ave(ipn) = cld_md_ave(ipn) + cld_md(ipn)
      cld_lo_ave(ipn) = cld_lo_ave(ipn) + cld_lo(ipn)
      cld_tt_ave(ipn) = cld_tt_ave(ipn) + cld_tt(ipn)
      cld_bl_ave(ipn) = cld_bl_ave(ipn) + cld_bl(ipn)
           t2max(ipn) =  max(t2max(ipn),   t2m2d(ipn))
           t2min(ipn) =  min(t2min(ipn),   t2m2d(ipn))

      end do
!$OMP END PARALLEL DO
!sms$ignore end
      end if	! ave_intvl


!    if (mod(its,ArchvStep6h) == 0 ) then
!      accum_start_6h = max(0, its-6)
    if (mod(its-ArchvStep0,min(ArchvStep,ArchvStep6h)) == 0 ) then
      accum_start_6h = max(0, its-1)
!sms$ignore begin
!$OMP PARALLEL DO schedule (static)
        DO ipn=ips,ipe
           rg2d(ipn) = rn2d(ipn) - rc2d(ipn)
          rn_6h(ipn) = rn2d(ipn) - rn2d6_0(ipn)
          rc_6h(ipn) = rc2d(ipn) - rc2d6_0(ipn)
          rg_6h(ipn) = rg2d(ipn) - rg2d6_0(ipn)
          rn2d6_0(ipn) = rn2d(ipn)
          rc2d6_0(ipn) = rc2d(ipn)
          rg2d6_0(ipn) = rg2d(ipn)
        ENDDO
!$OMP END PARALLEL DO
!sms$ignore end
    end if	! mod(its,min(ArchvStep,ArchvStep6h)) == 0

!---------------------------------------------------------------------------------
!         Output diagnostics of additional variables - op_diag
!---------------------------------------------------------------------------------
      ! Calculate various 3-d and 2-d diagnostic variables for outputting below.
      ! These are generally multivariate diagnostics, thus not do-able by the
      ! scalar FIMpost, which can only do horizontal interpolation (icos to
      ! lat/lon) one variable at a time.

! op_diag is called every 6hrs or ArchvStep, whatever is smaller for day & hr, or ArchvStep for "month"
    if ((mod(its-ArchvStep0,min(ArchvStep,ArchvStep6h)) == 0 .and. (ArchvTimeUnit=='dy' .or. ArchvTimeUnit=='hr')) &
      .or.  mod(its,ArchvStep) == 0) then	! op_diag is called every 6hr for day & hr archive, and monthly for mo archive
      call op_diag (its, time,				&
                    us3d, vs3d, dp3d, sdot, pr3d,	&
                    ex3d, mp3d, tr, ws3d,		&
                    ph3d, rn2d, rc2d, ts2d, us2d,	&
                    hf2d, qf2d, rsds, rlds, st3d,	&
                    tr(:,:,2), t2m2d, dum, rn_6h,	&
                    ! Below are output variables from op_diag
                    tk3d, rh3d, td3d, pw2d, pq2d, mslp,	&
                    g3p, g3p_chem, g3sig, g2d, t2m_dif,	&
                    dum, sn_6h, wt_int_rev, k_int_rev	&	! snow is accumulated at 6hr interval
                    )
      if (dayave9z .and. mod(its-ArchvStep0,ArchvStep) == 0) then
        if (its-ArchvStep0 > 0) then
          g3p_ave(:,:,:) = g3p_ave(:,:,:)/ArchvStep*min(ArchvStep6h,ArchvStep)
          slp_ave(:)     = slp_ave(:)    /ArchvStep*min(ArchvStep6h,ArchvStep)
        else 
          g3p_ave(:,:,:) = g3p(:,:,:)
          slp_ave(:)     = mslp(:)
        end if
        call maybe_write (its, 'hgtP', g3p_ave(:,:,1), nvlp)
        call maybe_write (its, 'tmpP', g3p_ave(:,:,2), nvlp)
        call maybe_write (its, 'qv3P', g3p_ave(:,:,3), nvlp)
        call maybe_write (its, 'up3P', g3p_ave(:,:,4), nvlp)
        call maybe_write (its, 'vp3P', g3p_ave(:,:,5), nvlp)
        call maybe_write (its, 'vv3P', g3p_ave(:,:,6), nvlp)
        call maybe_write (its, 'qc3P', g3p_ave(:,:,7), nvlp)
        call maybe_write (its, 'ms2D', slp_ave,   1, twodfile=.true.)
        g3p_ave(:,:,:)=0.
        slp_ave(:)=0.
      end if

!sms$ignore begin
!$OMP PARALLEL DO schedule (static)
      do ipn=ips,ipe
! sn_6h, g3p & mslp are output of op_diag, need this at its=0
        sn_xh   (ipn) =    sn_6h(ipn)	
        sn2d    (ipn) =     sn2d(ipn)  + sn_xh(ipn)	
       wsp80_ave(ipn) = wsp80_ave(ipn) +   g2d(ipn,7)
      g3p_ave(:,ipn,:)= g3p_ave(:,ipn,:)+g3p(:,ipn,:)
         slp_ave(ipn) =   slp_ave(ipn) +  mslp(ipn)
      end do
!$OMP END PARALLEL DO
!sms$ignore end

!     call stencl(cld_bl_ave,1,1.,'cld_bl_ave aft op_diag'//trim(stamp))
    end if	!mod(its,...) == 0

!---------------------------------------------------------------------------------
!JR Write history info every so often (mod(its,archvstep) == 0)
    if (mod(its-ArchvStep0,ArchvStep) == 0) then
      accum_start = curr_write_time
      curr_write_time = time

!sms$ignore begin
!$OMP PARALLEL DO schedule (static)
      do ipn=ips,ipe
        rg2d(ipn) = rn2d(ipn) - rc2d(ipn)
! Calculate x-hour interval precip (difference from total precip at the last output time)
        rn_xh(ipn) = rn2d(ipn) - rn2d0(ipn)
        rc_xh(ipn) = rc2d(ipn) - rc2d0(ipn)
        rg_xh(ipn) = rg2d(ipn) - rg2d0(ipn)
        rn2d0(ipn) = rn2d(ipn)
        rc2d0(ipn) = rc2d(ipn)
        rg2d0(ipn) = rg2d(ipn)
        spd10m_dif(ipn) = sqrt(u10m(ipn)**2 + v10m(ipn)**2) - &
                          sqrt(us3d(1,ipn)**2 + vs3d(1,ipn)**2)
      end do
!$OMP END PARALLEL DO
!sms$ignore end
     
      if (ips.eq.1) print *,'(subr.output) archiving data, time =',time

      if (ave_intvl .and. its-ArchvStep0 > 0) then
        icex_nh=0.; icev_nh=0.; icex_sh=0.; icev_sh=0.; 
!sms$ignore begin
!$OMP PARALLEL DO SCHEDULE (static) REDUCTION(+:icex_nh,icev_nh,icex_sh,icev_sh)
      do ipn=ips,ipe
        rsds_ave(ipn) =   rsds_ave(ipn) / ArchvStep
        rlds_ave(ipn) =   rlds_ave(ipn) / ArchvStep
        rsus_ave(ipn) =   rsus_ave(ipn) / ArchvStep
        rlus_ave(ipn) =   rlus_ave(ipn) / ArchvStep
        rsdt_ave(ipn) =   rsdt_ave(ipn) / ArchvStep
        rsut_ave(ipn) =   rsut_ave(ipn) / ArchvStep
        rlut_ave(ipn) =   rlut_ave(ipn) / ArchvStep
        fice_ave(ipn) =   fice_ave(ipn) * 100. / ArchvStep	! fraction => %
        hice_ave(ipn) =   hice_ave(ipn) / ArchvStep
          qf_ave(ipn) =     qf_ave(ipn) / ArchvStep
          hf_ave(ipn) =     hf_ave(ipn) / ArchvStep
          ts_ave(ipn) =     ts_ave(ipn) / ArchvStep
          t2_ave(ipn) =     t2_ave(ipn) / ArchvStep
          q2_ave(ipn) =     q2_ave(ipn) / ArchvStep
          td_ave(ipn) =     td_ave(ipn) / ArchvStep
          sm_ave(ipn) =     sm_ave(ipn) / ArchvStep
      runoff_ave(ipn) = runoff_ave(ipn) / ArchvStep
        ustr_ave(ipn) =   ustr_ave(ipn) / ArchvStep
        vstr_ave(ipn) =   vstr_ave(ipn) / ArchvStep
         u10_ave(ipn) =    u10_ave(ipn) / ArchvStep
         v10_ave(ipn) =    v10_ave(ipn) / ArchvStep
        if (ArchvTimeUnit=='dy' .or. ArchvTimeUnit=='hr') then
          wsp80_ave(ipn) =  wsp80_ave(ipn) / ArchvStep &
                          * min(ArchvStep6h,ArchvStep)
          if (.not. dayave9z) then
           g3p_ave(:,ipn,:)= g3p_ave(:,ipn,:)/ArchvStep &
                         * min(ArchvStep6h,ArchvStep)
           slp_ave(ipn) =    slp_ave(ipn) / ArchvStep &	! depend on how often op_diag is called; no ave is needed for mo archive
                         * min(ArchvStep6h,ArchvStep)
          end if
        end if
        cld_hi_ave(ipn)= cld_hi_ave(ipn) / ArchvStep
        cld_md_ave(ipn)= cld_md_ave(ipn) / ArchvStep 
        cld_lo_ave(ipn)= cld_lo_ave(ipn) / ArchvStep
        cld_tt_ave(ipn)= cld_tt_ave(ipn) / ArchvStep
        cld_bl_ave(ipn)= cld_bl_ave(ipn) / ArchvStep

        if (fice2d(ipn) .lt. .15) then
          iextent=0.
        else
          iextent=1.
        end if

        if (deg_lat(ipn).gt.0.) then
          icex_nh=icex_nh+iextent    *area(ipn)		    ! ice extent in NH
          icev_nh=icev_nh+fice2d(ipn)*area(ipn)*hice2d(ipn) ! ice volume in NH
        else
          icex_sh=icex_sh+iextent    *area(ipn)		    ! ice extent in SH
          icev_sh=icev_sh+fice2d(ipn)*area(ipn)*hice2d(ipn) ! ice volume in SH
        end if
      end do
!$OMP END PARALLEL DO
!sms$ignore end
      end if	!ave_intvl
!!SMS$REDUCE(icex_nh,icev_nh,icex_sh,icev_sh,sum)
      write(429,*) its,icex_nh,icev_nh,icex_sh,icev_sh

      call pvsurf(its, relvor, dp3d, tr, potvor, th_pvsrf, pr_pvsrf,	&
                  us_pvsrf, vs_pvsrf,pvsnap)

    if (atmonly) then	! to reduce output in the coupled run

      call outqv        (us3d,       1,   deg_lat, deg_lon,     1., 'U wind comp at k=1 - max/min')
      call outqv_mn_lat (us3d,       1,   deg_lat, deg_lon, 50.,1., 'U wind comp at k=1 - max/min')
      call outqv        (u10m,       1,   deg_lat, deg_lon,     1., 'U10M - max/min')
      call outqv_mn_lat (u10m,       1,   deg_lat, deg_lon, 50.,1., 'U10M - max/min')
      call outqv        (vs3d,       1,   deg_lat, deg_lon,     1., 'V wind comp at k=1 - max/min')
      call outqv_mn_lat (vs3d,       1,   deg_lat, deg_lon, 50.,1., 'V wind comp at k=1 - max/min')
      call outqv        (v10m,       1,   deg_lat, deg_lon,     1., 'V10M - max/min')
      call outqv_mn_lat (v10m,       1,   deg_lat, deg_lon, 50.,1., 'V10M - max/min')
      call outqv        (spd10m_dif, 1,   deg_lat, deg_lon,     1., 'spd_dif(wind-10M- wind(k=1)) - max/min')
      call outqv_mn_lat (spd10m_dif, 1,   deg_lat, deg_lon, 50.,1., 'spd_dif(wind-10M- wind(k=1)) - max/min')
      call outqv        (t2m2d,      1,   deg_lat, deg_lon,     1., 'T2M - max/min')
      call outqv_mn_lat (t2m2d,      1,   deg_lat, deg_lon, 50.,1., 'T2M - max/min')
      call outqv        (t2m_dif,    1,   deg_lat, deg_lon,     1., 'T2M - T(k=1)) - max/min')
      call outqv_mn_lat (t2m_dif,    1,   deg_lat, deg_lon, 50.,1., 'T2M - T(k=1)) - max/min')
      call outqv        (q2m2d,      1,   deg_lat, deg_lon,     1., 'Q2M - max/min')
      call outqv_mn_lat (q2m2d,      1,   deg_lat, deg_lon, 50.,1., 'Q2M - max/min')
      call outqv        (g2d(:,7),   1,   deg_lat, deg_lon,     1., 'Wind at 80 m - max/min')
      call outqv_mn_lat (g2d(:,7),   1,   deg_lat, deg_lon, 50.,1., 'Wind at 80 m - max/min')

      call outqv        (tr(:,:,4), nvl, deg_lat, deg_lon,     1.E5, 'O3 - max/min')
      call outqv_mn_lat (tr(:,:,4), nvl, deg_lat, deg_lon, 50.,1.E5, 'O3 - max/min')

      call outqv        (sheleg2d,   1,   deg_lat, deg_lon,     1., 'Snow water equivalent - max/min')
      call outqv_mn_lat (sheleg2d,   1,   deg_lat, deg_lon, 50.,1., 'Snow water equivalent - max/min')
      call outqv        (canopy2d,   1,   deg_lat, deg_lon,     1., 'Canopy water - max/min')
      call outqv_mn_lat (canopy2d,   1,   deg_lat, deg_lon, 30.,1., 'Canopy water - max/min')
      call outqv_mn_lat_land (canopy2d, 1, deg_lat, deg_lon, 30., slmsk2d, 1, 1., ' Canopy water - mean - land only')

      call outqv_mn_lat_land (g2d(:,5),1,deg_lat, deg_lon, 30., slmsk2d, 1, 1., ' RH w.r.t. PW - mean - land only')
      call outqv_mn_lat_land (pw2d,    1,deg_lat, deg_lon, 30., slmsk2d, 1, 1., ' PW - mean - land only')
      call outqv_mn_lat_land (pq2d,    1,deg_lat, deg_lon, 30., slmsk2d, 1, 1., ' Integ condensate - mean - land only')

      call outqv        (rn_xh,      1,   deg_lat, deg_lon,     1., 'Rain - since last output - max/min')
      call outqv        (rc_xh,      1,   deg_lat, deg_lon,     1., 'Rain-conv - since last output - max/min')
      call outqv        (rg_xh,      1,   deg_lat, deg_lon,     1., 'Rain-grid-scale - since last output - max/min')
      if (mod(its,ArchvStep6h) == 0) then
        call outqv       (rn_6h, 1, deg_lat, deg_lon, 1., 'Rain -6h-total  - max/min')
        call outqv       (sn_6h, 1, deg_lat, deg_lon, 1., 'Snow -6h-total  - max/min')
      endif
      call outqv       (rn2d, 1, deg_lat, deg_lon, 1., 'Rain -run-total  - max/min')
!     call outqv       (sn2d, 1, deg_lat, deg_lon, 1., 'Snow -run-total  - max/min')
      ! call outqv       (rc2d0, 1, deg_lat, deg_lon, 1., 'Rain0-conv  - max/min')
      if (its == itsStart-1 .or. (its==itsDFI-1.and.digifilt)) then
        call outqv        (hice2d,   1,   deg_lat, deg_lon,     1., 'Sea ice-hice - max/min')
        call outqv_mn_lat (hice2d,   1,   deg_lat, deg_lon, 50.,1., 'Sea ice-hice - max/min')
        call outqv        (fice2d,   1,   deg_lat, deg_lon,     1., 'Sea ice-fice - max/min')
        call outqv_mn_lat (fice2d,   1,   deg_lat, deg_lon, 50.,1., 'Sea ice-fice - max/min')
      end if
    end if	! atmonly=T

      IF (PrintDiags) THEN
        call outFMTed (its,pr3d,ph3d,us3d,vs3d,dp3d,mp3d,relvor,tr,rh3d,tk3d,ws3d, &
                       st3d,sm3d,rn2d,pw2d,pq2d,ts2d,us2d,hf2d,qf2d,rsds,rlds,time)
      else

        IF (enkfio_out .AND. .NOT. use_write_tasks .AND. .NOT. FixedGridOrder) THEN
! write only fields needed for EnKF DA cycle to two files.
! 3d dynamical and surface fields. Only supported if no write task, and FixedGridOrder
! = .false.
          CALL outvar_enkf_spectral(time,g3sig,ph3d(1,:))
          CALL outvar_enkf_sfcio(time,ph3d(1,:))

          IF ( time == cycle_freq ) CALL write_enkfio ()

        ENDIF			! enkfio_out

        IF (.NOT. enkfio_out .OR. enkf_diag) THEN

!JR Redo order for testing convenience to match specification in FIMnamelist ("wgrib" will give 
!JR same order). Optional argument "scalefactor" is scaling factor to be applied when GRIB file 
!JR is written. Optional argument "accum_start" specifies an accumulation start time different 
!JR from default when GRIB file is written.

          call maybe_write (its, 'rn2D', rn2d, 1, accum_start=0, twodfile=.true.)
!         call maybe_write (its, 'sn2D', sn2d, 1, accum_start=0, twodfile=.true.)
          call maybe_write (its, 'rc2D', rc2d, 1, accum_start=0, twodfile=.true.)
          call maybe_write (its, 'rg2D', rg2d, 1, accum_start=0, twodfile=.true.)
          call maybe_write (its, 'r12D', rn_xh,1, accum_start=accum_start, twodfile=.true.)
          call maybe_write (its, 'r22D', rc_xh,1, accum_start=accum_start, twodfile=.true.)
          call maybe_write (its, 'r32D', rg_xh,1, accum_start=accum_start, twodfile=.true.)
          call maybe_write (its, 's12D', sn_xh,1, accum_start=accum_start, twodfile=.true.)
 
          if (mod(its,ArchvStep6h) == 0 .and. ave6hr .and. ArchvTimeUnit.ne.'mo') then
            call maybe_write (its, 'r42D', rn_6h,1,accum_start=accum_start_6h,twodfile=.true.)
            call maybe_write (its, 's62D', sn_6h,1,accum_start=accum_start_6h,twodfile=.true.)
            call maybe_write (its, 'r52D', rc_6h,1,accum_start=accum_start_6h,twodfile=.true.)
            call maybe_write (its, 'r62D', rg_6h,1,accum_start=accum_start_6h,twodfile=.true.)
          endif

          call maybe_write (its, 'pw2D', pw2d, 1, twodfile=.true.)
          call maybe_write (its, 'pq2D', pq2d, 1, twodfile=.true.)
          call maybe_write (its, 'us2D', us2d, 1, twodfile=.true.)
          call maybe_write (its, 'sa2D', sheleg2d ,1, twodfile=.true.)
          call maybe_write (its, 'cb2D', g2d(:,1),1, twodfile=.true.)
          call maybe_write (its, 'ct2D', g2d(:,3),1, twodfile=.true.)
          call maybe_write (its, 'rp2D', g2d(:,5),1, twodfile=.true.)
          call maybe_write (its, 'cn2D', canopy2d,  1, twodfile=.true.)
          call maybe_write (its, 'lmsk', slmsk2d,   1, twodfile=.true.)
          call maybe_write (its, 'st3D', st3d,      4, twodfile=.true.)

          if (ave_intvl) then
            call maybe_write (its, 'rsds', rsds_ave,  1, twodfile=.true.)
            call maybe_write (its, 'rlds', rlds_ave,  1, twodfile=.true.)
            call maybe_write (its, 'rsus', rsus_ave,  1, twodfile=.true.)
            call maybe_write (its, 'rlus', rlus_ave,  1, twodfile=.true.)
            call maybe_write (its, 'rsdt', rsdt_ave,  1, twodfile=.true.)
            call maybe_write (its, 'rsut', rsut_ave,  1, twodfile=.true.)
            call maybe_write (its, 'rlut', rlut_ave,  1, twodfile=.true.)
            call maybe_write (its, 'fice', fice_ave , 1, twodfile=.true.)
            call maybe_write (its, 'hice', hice_ave , 1, twodfile=.true.)
            call maybe_write (its, 'hfls', qf_ave,    1, twodfile=.true.)
            call maybe_write (its, 'hfss', hf_ave,    1, twodfile=.true.)
            call maybe_write (its, 'ts2D', ts_ave,    1, twodfile=.true.)
            call maybe_write (its, 't22D', t2_ave,    1, twodfile=.true.)
            call maybe_write (its, 'q22D', q2_ave,    1, twodfile=.true.)
            call maybe_write (its, 'td2D', td_ave,    1, twodfile=.true.)
            call maybe_write (its, 'sm2D', sm_ave,    1, twodfile=.true.)
            call maybe_write (its, 'runo', runoff_ave,1, twodfile=.true.)
            call maybe_write (its, 'w080', wsp80_ave, 1, twodfile=.true.)
            call maybe_write (its, 'ustr', ustr_ave,  1, twodfile=.true.)
            call maybe_write (its, 'vstr', vstr_ave,  1, twodfile=.true.)
            call maybe_write (its, 'u12D', u10_ave,   1, twodfile=.true.)
            call maybe_write (its, 'v12D', v10_ave,   1, twodfile=.true.)
            call maybe_write (its, 'tmax', t2max,     1, twodfile=.true.)
            call maybe_write (its, 'tmin', t2min,     1, twodfile=.true.)

           if (.not. dayave9z) then
            call maybe_write (its, 'hgtP', g3p_ave(:,:,1), nvlp)
            call maybe_write (its, 'tmpP', g3p_ave(:,:,2), nvlp)
            call maybe_write (its, 'qv3P', g3p_ave(:,:,3), nvlp)
            call maybe_write (its, 'up3P', g3p_ave(:,:,4), nvlp)
            call maybe_write (its, 'vp3P', g3p_ave(:,:,5), nvlp)
            call maybe_write (its, 'vv3P', g3p_ave(:,:,6), nvlp)
            call maybe_write (its, 'qc3P', g3p_ave(:,:,7), nvlp)
            call maybe_write (its, 'ms2D', slp_ave,   1, twodfile=.true.)
           end if
          else
            call maybe_write (its, 'rsds', rsds,      1, twodfile=.true.)
            call maybe_write (its, 'rlds', rlds,      1, twodfile=.true.)
            call maybe_write (its, 'rsus', rsus,      1, twodfile=.true.)
            call maybe_write (its, 'rlus', rlus,      1, twodfile=.true.)
            call maybe_write (its, 'rsdt', rsdt,      1, twodfile=.true.)
            call maybe_write (its, 'rsut', rsut,      1, twodfile=.true.)
            call maybe_write (its, 'rlut', rlut,      1, twodfile=.true.)
            call maybe_write (its, 'hfls', qf2d,      1, twodfile=.true.)
            call maybe_write (its, 'hfss', hf2d,      1, twodfile=.true.)
            call maybe_write (its, 'w080', g2d(:,7),  1, twodfile=.true.)
            call maybe_write (its, 'ts2D', ts2d,      1, twodfile=.true.)
            call maybe_write (its, 't22D', t2m2d,     1, twodfile=.true.)
            call maybe_write (its, 'q22D', q2m2d,     1, twodfile=.true.)
            call maybe_write (its, 'u12D', u10m,      1, twodfile=.true.)
            call maybe_write (its, 'v12D', v10m,      1, twodfile=.true.)
            call maybe_write (its, 'ms2D', mslp,      1, twodfile=.true.)
            call maybe_write (its, 'fice', fice2d,    1, twodfile=.true.)
            call maybe_write (its, 'hice', hice2d,    1, twodfile=.true.)
            call maybe_write (its, 'sm3D', sm3d,      4, twodfile=.true.)

            call maybe_write (its, 'hgtP', g3p(:,:,1), nvlp)
            call maybe_write (its, 'tmpP', g3p(:,:,2), nvlp)
            call maybe_write (its, 'qv3P', g3p(:,:,3), nvlp)
            call maybe_write (its, 'up3P', g3p(:,:,4), nvlp)
            call maybe_write (its, 'vp3P', g3p(:,:,5), nvlp)
            call maybe_write (its, 'vv3P', g3p(:,:,6), nvlp)
            call maybe_write (its, 'qc3P', g3p(:,:,7), nvlp)
          end if	! ave_intvl

          call maybe_write (its, 'thpv', th_pvsrf , 1, twodfile=.true.)
          call maybe_write (its, 'prpv', pr_pvsrf , 1, twodfile=.true.)
          call maybe_write (its, 'uspv', us_pvsrf , 1, twodfile=.true.)
          call maybe_write (its, 'vspv', vs_pvsrf , 1, twodfile=.true.)
          call maybe_write (its, 'clhi', cld_hi_ave,1, twodfile=.true.)
          call maybe_write (its, 'clmd', cld_md_ave,1, twodfile=.true.)
          call maybe_write (its, 'cllo', cld_lo_ave,1, twodfile=.true.)
          call maybe_write (its, 'cltt', cld_tt_ave,1, twodfile=.true.)
! cloud cover is not bit-wise identical after restart
!         call maybe_write (its, 'clbl', cld_bl_ave,1, twodfile=.true.)
          call maybe_write (its, 'snp1', pvsnap(1,:),1,twodfile=.true.)
          call maybe_write (its, 'snp2', pvsnap(2,:),1,twodfile=.true.)
          call maybe_write (its, 'snp3', pvsnap(3,:),1,twodfile=.true.)
          call maybe_write (its, 'snp4', pvsnap(4,:),1,twodfile=.true.)
          call maybe_write (its, 'snp5', pvsnap(5,:),1,twodfile=.true.)

        IF (.NOT. output2d) THEN         ! write out 3d fields
        
            call maybe_write (its, 'pr3D', pr3d, nvlp1)
            call maybe_write (its, 'ph3D', ph3d, nvlp1, scalefactor=1./9.8)
            call maybe_write (its, 'tk3D', tk3d, nvl)
            call maybe_write (its, 'td3D', td3d, nvl)
            call maybe_write (its, 'ws3D', ws3d, nvl)
            call maybe_write (its, 'rh3D', rh3d, nvl)
            call maybe_write (its, 'us3D', us3d, nvl)
            call maybe_write (its, 'vs3D', vs3d, nvl)

!JR   These arent specified in default FIMnamelist

            call maybe_write (its, 'dp3D', dp3d,      nvl)
            call maybe_write (its, 'mp3D', mp3d,      nvl)
            call maybe_write (its, 'th3D', tr(:,:,1), nvl)
            call maybe_write (its, 'qv3D', tr(:,:,2), nvl, scalefactor=1000.)
            call maybe_write (its, 'qw3D', tr(:,:,3), nvl, scalefactor=1000.)
            call maybe_write (its, 'oz3D', tr(:,:,4), nvl, scalefactor=1000.)
            call maybe_write (its, 'vo3D', relvor,    nvl)
            call maybe_write (its, 'pv3D', potvor,    nvl)
            call maybe_write (its, 'da3D', diaga,     nvl)
!           call maybe_write (its, 'da3D', ext_cof(:,:,10),nvl)	! a test
            call maybe_write (its, 'db3D', diagb,     nvl)
        end if   ! (.NOT. output2d)

!sms$ignore begin
!$OMP PARALLEL DO schedule (static)
      do ipn=ips,ipe
	rsds_ave(ipn) = 0.; rlds_ave(ipn) = 0.; rsus_ave(ipn) = 0.; rlus_ave(ipn) = 0.
	rsdt_ave(ipn) = 0.; rsut_ave(ipn) = 0.; rlut_ave(ipn) = 0.
	  qf_ave(ipn) = 0.; hf_ave(ipn) = 0.
	  ts_ave(ipn) = 0.; t2_ave(ipn) = 0.; q2_ave(ipn) = 0.; fice_ave(ipn)=0.; hice_ave(ipn)=0.
          td_ave(ipn) = 0.; sm_ave(ipn) = 0.
      runoff_ave(ipn) = 0.
       wsp80_ave(ipn) = 0.
	ustr_ave(ipn) = 0.; vstr_ave(ipn) = 0.
         u10_ave(ipn) = 0.;  v10_ave(ipn) = 0.
	   rn_xh(ipn) = 0.; rc_xh(ipn) = 0.; rg_xh(ipn) = 0.; sn_xh(ipn) = 0.
	   rn_6h(ipn) = 0.; rc_6h(ipn) = 0.; rg_6h(ipn) = 0.; sn_6h(ipn) = 0.
	   t2max(ipn) = -9999.; t2min(ipn) =  9999.
      if (.not. dayave9z) then
	 g3p_ave(:,ipn,:) = 0.
	 slp_ave(ipn) = 0.
      end if
      cld_hi_ave(ipn) = 0.; cld_md_ave(ipn) = 0.; cld_lo_ave(ipn) = 0.; cld_tt_ave(ipn) = 0.; cld_bl_ave(ipn) = 0.
      end do
!$OMP END PARALLEL DO
!sms$ignore end

! so far only outputting pressure level chem data for GOCART option
!
          IF (num_chem > 0 .AND. chem_opt == 300)THEN
            CALL maybe_write (its, 'so2P', g3p_chem(:,:,1), nvlp)
            CALL maybe_write (its, 'slfP', g3p_chem(:,:,2), nvlp)
            CALL maybe_write (its, 'dmsP', g3p_chem(:,:,3), nvlp)
            CALL maybe_write (its, 'msaP', g3p_chem(:,:,4), nvlp)
            CALL maybe_write (its, 'p25P', g3p_chem(:,:,5), nvlp)
            CALL maybe_write (its, 'bc1P', g3p_chem(:,:,6), nvlp)
            CALL maybe_write (its, 'bc2P', g3p_chem(:,:,7), nvlp)
            CALL maybe_write (its, 'oc1P', g3p_chem(:,:,8), nvlp)
            CALL maybe_write (its, 'oc2P', g3p_chem(:,:,9), nvlp)
            CALL maybe_write (its, 'd1sP', g3p_chem(:,:,10), nvlp)
            CALL maybe_write (its, 'd2sP', g3p_chem(:,:,11), nvlp)
            CALL maybe_write (its, 'd3sP', g3p_chem(:,:,12), nvlp)
            CALL maybe_write (its, 'd4sP', g3p_chem(:,:,13), nvlp)
            CALL maybe_write (its, 'd5sP', g3p_chem(:,:,14), nvlp)
            CALL maybe_write (its, 's1sP', g3p_chem(:,:,15), nvlp)
            CALL maybe_write (its, 's2sP', g3p_chem(:,:,16), nvlp)
            CALL maybe_write (its, 's3sP', g3p_chem(:,:,17), nvlp)
            CALL maybe_write (its, 's4sP', g3p_chem(:,:,18), nvlp)
            CALL maybe_write (its, 'p10P', g3p_chem(:,:,19), nvlp)
!           CALL maybe_write (its, 'slfP', g3p(:,:,ichemstart+2), nvlp)
!           CALL maybe_write (its, 'dmsP', g3p(:,:,ichemstart+3), nvlp)
!           CALL maybe_write (its, 'msaP', g3p(:,:,ichemstart+4), nvlp)
!           CALL maybe_write (its, 'p25P', g3p(:,:,ichemstart+5), nvlp)
!           CALL maybe_write (its, 'bc1P', g3p(:,:,ichemstart+6), nvlp)
!           CALL maybe_write (its, 'bc2P', g3p(:,:,ichemstart+7), nvlp)
!           CALL maybe_write (its, 'oc1P', g3p(:,:,ichemstart+8), nvlp)
!           CALL maybe_write (its, 'oc2P', g3p(:,:,ichemstart+9), nvlp)
!           CALL maybe_write (its, 'd1sP', g3p(:,:,ichemstart+10), nvlp)
!           CALL maybe_write (its, 'd2sP', g3p(:,:,ichemstart+11), nvlp)
!           CALL maybe_write (its, 'd3sP', g3p(:,:,ichemstart+12), nvlp)
!           CALL maybe_write (its, 'd4sP', g3p(:,:,ichemstart+13), nvlp)
!           CALL maybe_write (its, 'd5sP', g3p(:,:,ichemstart+14), nvlp)
!           CALL maybe_write (its, 's1sP', g3p(:,:,ichemstart+15), nvlp)
!           CALL maybe_write (its, 's2sP', g3p(:,:,ichemstart+16), nvlp)
!           CALL maybe_write (its, 's3sP', g3p(:,:,ichemstart+17), nvlp)
!           CALL maybe_write (its, 's4sP', g3p(:,:,ichemstart+18), nvlp)
!           CALL maybe_write (its, 'p10P', g3p(:,:,ichemstart+19), nvlp)
          END IF		! (num_chem > 0 .AND. chem_opt == 300)
           
        END IF			! (.NOT. enkfio_out .OR. enkf_diag)
         
      END IF			! (PrintDiags)

      ret = gptlstart ('maxmin')
      call printMAXMIN(its,nvl,nip,ntra+ntrb,tr,dp3d)
      ret = gptlstop ('maxmin')
      if (ips.eq.1) print *,'(subr.output) archiving completed'

!     call stencl(cld_bl_ave,1,1000.,'clbl output_write'//trim(stamp))
    end if			! (mod(its,ArchvStep) == 0)

    if (digifilt) then
      laststep = itsDFI + nts
    else
      laststep = itsStart + nts
    end if

!TODO: Fix this to print max, min times and keep separate from OMP
!    if (its == laststep-1) then
!      ret = gptlget_wallclock ('Output', 0, toutput)  ! The "0" is thread number
!      ret = gptlget_wallclock ('maxmin', 0, tmaxmin)  ! The "0" is thread number
!      print "(' OUTPUT time, maxmin time:',2F10.1)", toutput, tmaxmin
!    endif

    print *,'(output) starting atm_flxsum calculation'
    if (its.ge.itsStart-1+nts) call atm_flxsum(itsStart-1,nts,flxavg)
    print *,'(output) atm_flxsum calculation completed'

    deallocate (g3sig)
    return

  end subroutine output

end module module_output
