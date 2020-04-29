module module_do_physics_one_step

!*********************************************************************
!     do_physics_one_step
!	Calculates column forcing for global fim
!	12/21/2005 - Alexander E. MacDonald     - original version
!	05/01/2006 - Jian-Wen Bao               - modified for GFS physics
!       04/14/2008 - Stan Benjamin, John Brown  - modifications
!                      for introduction of virtual pot temp for temp prog 
!                      variable instead of previous non-virtual pot temp
!       02/26/2009 - Tom Henderson       - moved here from physics.F90 to more 
!                                          closely match new GFS r3038
!       07/21/2009 - Jian-wen Bao        - change to random number generator
!                      for xkt2 (cloud-top height) instead of previous 0.6 constant
!                      This follows NCEP's use of this random number generator
!                      (Mersenne twister) and appears to qualitatively improve
!                      forecasts in tropics.
!       06/08/2015 - Jim Rosinski - Rewrote to enable OpenMP of loop over calls
!                      to gloopr and gbphys. Extracted code within that giant
!                      loop and subsumed into a single routine named "column_physics".
!                      That new routine is itself "contained" within subroutine
!                      do_physics_one_step, enabling "column_physics" to see all
!                      input arguments and local variables in do_physics_one_step.
!                      Also: Extracted ozone computation code into separate ozone 
!                      initialization and computation routines "ozone_init" and 
!                      "ozone_calc" (contained within the module). Previously ozone
!                      data were read every physics time step, now this is only
!                      done once per run.
!                      Also: moved the call to rad_initialize here from gloopr. This
!                      call should actually be moved even higher in the call tree to
!                      avoid "if (first) then" logic.
!*********************************************************************

  use global_bounds,  only: ips, ipe, ihe, myrank
  use module_control, only: nip, itsdfi, dospptdt
  use fimnamelist,    only: nvl, digifilt, ishal_cnv, ideep_cnv, sppt, coupled
  use machine,        only: kind_phys, kind_rad, kind_evod
  use infnan,         only: inf

  implicit none

#include <gptl.inc>

  integer, parameter :: IM = 1
  integer, parameter :: IX = 1
! Ozone stuff
  integer, parameter :: ko3 = 80
  integer, parameter :: pcoeff = 4    ! ozone levels in climatology
  integer, parameter :: latsoz = 35
  integer, parameter :: oztime = 13
  integer :: IDAT(8)
  integer :: timeoz
  integer :: latsozp
  real*8 :: poz(KO3)
  real ::  pl_time(oztime)
  real, save :: ozplin(latsoz,ko3,pcoeff,oztime)

!nwang added for sppt
  integer, save :: pat_lun

!sms$distribute (dh,1) begin
  real, allocatable :: pattern(:) 
!sms$distribute end

!sms$distribute (dh,3) begin
  real*8, allocatable :: ozplout(:,:,:) ! Dims will be (ko3,pcoeff,nip)
!sms$distribute end

contains

  subroutine do_physics_one_step (its, dtp, kdt, phour,deltim,          &
                 gfs_ps, gfs_dp, gfs_dpdt, gfs_p, gfs_u, gfs_v,         &
                 gfs_t, gfs_q, gfs_oz, gfs_cld,                         &
                 sfc_fld, flx_fld,                                      &
                 XLON,XLAT,COSZDG,                                      &
                 HPRIME, SWH, HLW, SWHC, HLWC, FLUXR, SFALB, SLAG, SDEC,&
                 CDEC,                                                  &
                 phy_f3d, phy_f2d, NBLCK,                               &
                 CLDCOV,                                                &
                 LATS_NODE_R, NMTVR, num_p3d, num_p2d, NFXR,            &
                 lsoil, CallRadiation,                                  &
                 yyyymmddhhmm,                                          &
                 skip_cu_physics, skip_mp_physics, skip_chem,sscal,     &
                 ext_cof,asymp,extlw_cof)

    USE gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data, Flx_Var_Data
    use layout1, only: me
! jbao orig USE resol_def,		ONLY : levs,levp1
    USE resol_def,		ONLY : levs,levp1,levr, ntcw, ntoz, ncld, ntrac_resol_def=>ntrac
    USE module_variables,	ONLY : diaga,diagb,conv_act
    USE module_constants,	ONLY : prox,nprox
    USE findmaxmin1
    use module_wrf_control, only: nbands
    USE module_wrf_variables, ONLY:phys3dwrf,phys2dwrf
    USE module_wrfphysvars, ONLY: RTHRATENSW,RTHRATENLW
    USE mersenne_twister
!JR Moved namelist_def here from gloopr to enable calling rad_initialize
    use namelist_def, only: crick_proof, ccnorm, norad_precip, ras, &
                            ldiag3d, lssav
    use fimnamelist,  only: ictm, isol, ico2, iaer, iaer_mdl, ialb, iems, &
                            iovr_sw, iovr_lw, isubc_sw, isubc_lw, n3dflxtvd
    use date_def, only: idate
    use module_do_physics_one_step_chem,only:do_physics_one_step_chem
    use physics_ipn_its, only: physics_init_ipn

! Dimension and type external variables:

    integer              , intent(in   ) :: its
    real (kind=kind_phys), intent(in   ) :: dtp
    integer              , intent(in   ) :: kdt
    real (kind=kind_rad) , intent(in   ) :: phour
! for may 2014 jbao
    real (kind=kind_rad) , intent(in   ) :: deltim
! for may 2014 jbao
    TYPE(Sfc_Var_Data)   , intent(inout) :: sfc_fld
    TYPE(Flx_Var_Data)   , intent(inout) :: flx_fld
    integer              , intent(in   ) :: NBLCK
    integer              , intent(in   ) :: LATS_NODE_R
    integer              , intent(in   ) :: NMTVR
    integer              , intent(in   ) :: num_p3d
    integer              , intent(in   ) :: num_p2d
    integer              , intent(in   ) :: NFXR
    integer              , intent(in   ) :: lsoil
    integer              , intent(in   ) :: CallRadiation
    CHARACTER(len=12)    , intent(in   ) :: yyyymmddhhmm
    logical              , intent(in   ) :: skip_cu_physics,skip_mp_physics,skip_chem

!sms$distribute (dh,1) begin
    real*8, intent(in   ) :: gfs_ps(nip)
    real*8, intent(in   ) :: gfs_dp(nip,nvl)
    real*8, intent(in   ) :: gfs_dpdt(nip,nvl)
    real*8, intent(in   ) :: gfs_p(nip,nvl)
    real*8, intent(inout) :: gfs_u(nip,nvl)
    real*8, intent(inout) :: gfs_v(nip,nvl)
    real*8, intent(inout) :: gfs_t(nip,nvl)
    real*8, intent(inout) :: gfs_q(nip,nvl)
    real*8, intent(inout) :: gfs_oz(nip,nvl)
    real*8, intent(inout) :: gfs_cld(nip,nvl)
    real :: field(nip)
    real (kind=kind_rad) , intent(in   ) :: XLON(nip,LATS_NODE_R)
    real (kind=kind_rad) , intent(in   ) :: XLAT(nip,LATS_NODE_R)
    real (kind=kind_rad) , intent(inout) :: COSZDG(nip,LATS_NODE_R)
    real (kind=kind_rad) , intent(inout) :: SWH(nip,nvl,NBLCK,LATS_NODE_R)
    real (kind=kind_rad) , intent(inout) :: HLW(nip,nvl,NBLCK,LATS_NODE_R)
! jbao gfs 2014 update for clear sky
    real (kind=kind_rad) , intent(inout) :: SWHC(nip,nvl,NBLCK,LATS_NODE_R)
    real (kind=kind_rad) , intent(inout) :: HLWC(nip,nvl,NBLCK,LATS_NODE_R)
! jbao gfs 2014 update for clear sky
    real (kind=kind_rad) , intent(inout) :: SFALB(nip,LATS_NODE_R)
    real (kind=kind_evod), intent(inout) :: SLAG(nip,LATS_NODE_R)
    real (kind=kind_evod), intent(inout) :: SDEC(nip,LATS_NODE_R)
    real (kind=kind_evod), intent(inout) :: CDEC(nip,LATS_NODE_R)
    real (kind=kind_rad) , intent(inout) :: phy_f3d(nip,nvl,NBLCK,LATS_NODE_R,num_p3d)
    real (kind=kind_rad) , intent(inout) :: phy_f2d(nip,LATS_NODE_R,num_p2d)
!sms$distribute end

!sms$distribute (dh,2) begin
    real (kind=kind_rad) , intent(inout) :: HPRIME(NMTVR,nip,LATS_NODE_R)
    real (kind=kind_rad) , intent(inout) :: FLUXR(NFXR,nip,LATS_NODE_R)
    real (kind=kind_rad) , intent(inout) :: CLDCOV(nvl,nip,LATS_NODE_R)
    real, intent(inout) :: sscal(nvl,nip,nbands)
    real, intent(inout) :: ext_cof(nvl,nip,nbands)
    real, intent(inout) :: asymp(nvl,nip,nbands)
    real, intent(inout) :: extlw_cof(nvl,nip,16)
!sms$distribute end

!Parameters and arrays used in the GFS physics
!----------------------------------------------------------------------
! TODO:  move nrcm to module resol_def
! Local variables

    integer, parameter :: ntrac = 3   ! Number of tracers
    integer, parameter :: nrcm = 1
    integer :: k, ipn
    integer, save :: lonf,latg
    
    real(kind=kind_phys) :: dtf,FHOUR,solhr

! TBH:  These arrays are used to avoid recomputation of values 
! TBH:  between calls to GLOOPR and GBPHYS.  
    real(kind=kind_phys) :: PRSI_S(IX,NVL+1)

    integer :: hour
    logical :: CallRadiationNow
    logical :: DoSPPT
    integer :: year,month,day

    integer iter,nvars,edg,smpass
    real sum
    
    integer :: ret  ! return code from GPTL routines
    real    :: dx

    logical, save :: sas_shal 
    logical, save :: first = .true.   !JR TODO: Get rid of this!
!  ---  iflip indicates model vertical index direction:
!     integer, parameter :: iflip = 0    ! vertical profile index from top to bottom
!!    integer, parameter :: iflip = 1    ! vertical profile index from bottom to top
! jbao gfs 2014 from 2012 gloopr
    integer, parameter :: iflip = 1      ! vertical profile index from bottom to top
    real (kind=kind_phys) :: si_loc(nvl+1)
#ifdef NAG
! INF is not permitted in an initialisation expression, NAG does this with -nan
    real (kind=kind_phys), save :: cdmbgwd(2)
#else
    real (kind=kind_phys), save :: cdmbgwd(2) = (/inf,inf/)  ! Init to bad value
#endif

    ret = gptlstart ('do_physics_one_step')

    ntrac_resol_def = ntrac
    levs  = nvl    !jfm
    levr  = levs
    levp1 = levs+1 !JFM
    smpass = 0
    dtf = dtp

    idat(:) = 0
! get date info from the date string
    READ(UNIT=yyyymmddhhmm(1:4), FMT='(I4)') year
    READ(UNIT=yyyymmddhhmm(5:6), FMT='(I2)') month
    READ(UNIT=yyyymmddhhmm(7:8), FMT='(I2)') day
    READ(UNIT=yyyymmddhhmm(9:10), FMT='(I2)') hour

    fhour = phour
    solhr = REAL(mod(REAL(fhour+hour),24.0),kind_phys)

    idat(1) = year
    idat(2) = month
    idat(3) = day
    idat(5) = hour
!jbao, jr: Bugfix: Set idate elements--previously default initial value of 0
!jbao, jr: was used
    idate(1) = idat(5)
    idate(2) = idat(3)
    idate(3) = idat(2)
    idate(4) = idat(1)

! Calculate grid length, used later to calculate resolution-related dimension
! parameters 
!     used for gravity wave drag
    dx = sqrt (5.e14 / float(nip))               ! 27,500m (27.5km) at g8

    CallRadiationNow = (mod(kdt,CallRadiation) == 0 .or. kdt == 1 .or.	&
		        (kdt==itsDFI.and.digifilt))
    if(CallRadiationNow)write(0,*)'calling radiation, kdt = ',kdt

    DoSPPT = (mod(kdt,dospptdt) == 0 .or. kdt == 1) 
    if(sppt .and. DOSPPT)write(0,*)'applying SPPT, kdt = ',kdt

!----------------------------------------------------------------------
!    Loop begins over all horizontal grid points
!----------------------------------------------------------------------

!SMS$PARALLEL(dh,ipn) BEGIN 
!sms$compare_var(gfs_ps, "do_physics.F90 - gfs_ps0 ")
!sms$compare_var(gfs_dp, "do_physics.F90 - gfs_dp0 ")
!sms$compare_var(gfs_dpdt, "do_physics.F90 - gfs_dpdt0 ")
!sms$compare_var(gfs_p, "do_physics.F90 - gfs_p0 ")
!sms$compare_var(gfs_oz, "do_physics.F90 - gfs_oz0 ")
!sms$compare_var(gfs_cld, "do_physics.F90 - gfs_cld0 ")
!sms$compare_var(gfs_u, "do_physics.F90 - gfs_u0 ")
!sms$compare_var(gfs_v, "do_physics.F90 - gfs_v0 ")
!sms$compare_var(gfs_t, "do_physics.F90 - gfs_t0 ")
!sms$compare_var(gfs_q, "do_physics.F90 - gfs_q0 ")
!SMS$PARALLEL END

!JR Moved rad_initialize call to here because it is called
!JR only once per model run. For now use the horrible
!JR if(first)then; first=.false. approach
!JR At least that allows the ipn loop below to be threaded.
!TODO: Move this call up the tree, probably to run.F90
!JR Preserve for now what appears to be a bug: Each MPI task
!JR uses si_loc for its starting ipn value. Why this does
!JR not cause a test suite failure I am not sure

    if (first) then
      first = .false.
      ret = gptlstart ('dpos_first')

      if (n3dflxtvd /= ntrac) then
        write(6,*)'do_physics_one_step: namelist var n3dflxtvd, ntrac=',n3dflxtvd, ntrac
        stop
      end if

      call physics_init_ipn ()         ! For retrieving ipn in physics
      call ozone_init ()
      call ozone_calc (lats_node_r, xlat, phour)

! nwang addeed to sppt
      if (sppt) then
        call sppt_init ()
      
        if (DoSPPT) then
!SMS$SERIAL(<pattern,out> : default=ignore) BEGIN
          read(pat_lun) pattern
!SMS$SERIAL END
        endif
      endif
 
      lonf = int(40075000. / dx)   !  40,075km = circumference of earth at equator
      latg =  lonf / 2
! 20150414 - SGB - Per Fanglin Yang, use /2.0,0.25/ for T1534
!    but                                 /0.25,2.0/ for T574
!    with 2015 version of GFS.
! 20151119 - SGB - value recommended >25km give poorer results
!   per 2015 10-month FIM experiment - FIMRETRO_FIM_r5137_cdmbgwd
!sb   if (dx.ge.25000.) then
! change for gfs input 2014
!sb     cdmbgwd(1) = 0.25
!sb     cdmbgwd(2) = 2.0 
!sb    else

        cdmbgwd(1) = 2.0
        cdmbgwd(2) = 0.25

!sb   end if
      me   = myrank

      PRSI_S(1,1) = gfs_ps(ips)
      do k=1,nvl
        PRSI_S(1,k+1) = PRSI_S(1,k) - gfs_dp(ips,k)
      end do
      do k=1,nvl+1
        si_loc(k) = prsi_s(1,k)/prsi_s(1,1)
      end do

!jbao 2014 added sas_shal as in 2012 gloopr for radiation initialization
!     sas_shal = sashal .and. (.not. ras)	! this is replaced by ishal_cnv

      ntoz = 2   !JR Is this used?
      ntcw = 3   !JR This comes from resol_def. Used to be parameter in dpos, but unused below gbphys
      call rad_initialize                                              &
!  ---  inputs:
          ( si_loc,levr,ictm,isol,ico2,iaer,iaer_mdl,ialb,iems,ntcw,   &
            num_p3d,ntoz,iovr_sw,iovr_lw,isubc_sw,isubc_lw,            &
            ishal_cnv,ideep_cnv,crick_proof,ccnorm,norad_precip,idate, &
            iflip,me )

!JR Do initial iteration outside of threaded loop to accommodate 
!JR numerous "if (first) then" constructs inside the GFS physics
      call column_physics (ips)
      ret = gptlstop ('dpos_first')

!$OMP PARALLEL DO SCHEDULE (guided)
      do ipn=ips+1,ipe
        call column_physics (ipn)
      end do
!$OMP END PARALLEL DO

    else   ! Not first time through

      if (sppt) then
        if (DoSPPT) then
!SMS$SERIAL(<pattern,out> : default=ignore) BEGIN
          read(pat_lun) pattern
!SMS$SERIAL END
        endif
      endif
 
      call ozone_calc (lats_node_r, xlat, phour)
!JR Do initial iteration outside of threaded loop due to dataset reading
!JR e.g. aerosols which may happen in the middle of the run.
!JR Another solution might be to instead install "if(mythread==0)" inside
!JR all the places inside GFS physics where dataset reading may happen.
      call column_physics (ips)
!$OMP PARALLEL DO SCHEDULE (guided)
      do ipn=ips+1,ipe
        call column_physics (ipn)
      end do
!$OMP END PARALLEL DO
    end if

!SMS$PARALLEL (dh,ipn) BEGIN
!sms$compare_var(gfs_u, "do_physics.F90 - gfs_u1 ")
!sms$compare_var(gfs_v, "do_physics.F90 - gfs_v1 ")
!sms$compare_var(gfs_t, "do_physics.F90 - gfs_t1 ")
!sms$compare_var(gfs_q, "do_physics.F90 - gfs_q1 ")
!sms$compare_var(gfs_ps, "do_physics.F90 - gfs_ps1 ")
!sms$compare_var(gfs_dp, "do_physics.F90 - gfs_dp1 ")
!sms$compare_var(gfs_dpdt, "do_physics.F90 - gfs_dpdt1 ")
!sms$compare_var(gfs_p, "do_physics.F90 - gfs_p1 ")
!sms$compare_var(gfs_oz, "do_physics.F90 - gfs_oz1 ")
!sms$compare_var(gfs_cld, "do_physics.F90 - gfs_cld1 ")
!SMS$PARALLEL END

!
! do this for many levels
!
    if ( skip_cu_physics .and. smpass > 0 ) then
!SMS$PARALLEL (dh,ipn) BEGIN
      do nvars=2,7
        if (nvars == 4) cycle
        do k=1,nvl-5
          field(:) = phys3dwrf(k,:,nvars)
          do iter=1,smpass
!     if(nvars.eq.5 .and. k.eq.5 .and. iter.eq.1)call findmxmn1(field,nip,'smoothr')
!SMS$EXCHANGE(field)
            do ipn=1,nip
              sum = field(ipn)*.5*nprox(ipn)
              do edg=1,nprox(ipn)
                sum = sum + field(prox(edg,ipn))
              end do
              field(ipn) = sum/(1.5*nprox(ipn))
            end do
          end do                ! iter
          phys3dwrf(k,:,nvars) = field(:)
        end do
      end do
!SMS$PARALLEL END
    end if
    ret = gptlstop ('do_physics_one_step')
    return

  CONTAINS   ! Note this "contains" is inside the subroutine, not the module

!JR Make new threaded driver routine "contain"ed inside do_physics_one_step so that all input arguments
!JR to do_physics_one_step are also visible here
    subroutine column_physics (ipn)
      use fimnamelist, only: PrintIpnDiag
      use module_initial_chem_namelists, only: chem_opt,ra_lw_physics
      use physics_ipn_its, only: physics_set_ipn_its, physics_reset_ipn

      integer, intent(in) :: ipn
  
      logical :: vrbos
      integer :: k
      integer :: ret

      real(kind=kind_phys) :: prsi_s(ix,nvl+1)
      real(kind=kind_phys) :: prsl_s(ix,nvl)
      real(kind=kind_phys) :: prsl(ix,nvl)
      real(kind=kind_phys) :: prsi(ix,nvl+1)
      real(kind=kind_phys) :: tg(ix,nvl)
      real(kind=kind_phys) :: qg(ix,nvl,ntrac)  !JR TODO: change from parameter to resol_def value
      real(kind=kind_phys) :: qg1(ix,nvl,2)
      real(kind=kind_phys) :: ug(ix,nvl)
      real(kind=kind_phys) :: vg(ix,nvl)
      real(kind=kind_phys) :: tisfc(im)
      real(kind=kind_phys) :: forcet(ix,nvl)
      real(kind=kind_phys) :: forceq(ix,nvl)

      real(kind=kind_rad) :: sscal_rad(nvl,nbands)
      real(kind=kind_rad) :: ext_cof_rad(nvl,nbands)
      real(kind=kind_rad) :: asymp_rad(nvl,nbands)
      real(kind=kind_rad) :: extlw_cof_rad(nvl,16)

      real(kind=kind_phys) :: phy_f3d_loc(ix,nvl,num_p3d)
      real(kind=kind_phys) :: phy_f2d_loc(ix,num_p2d)
      real(kind=kind_phys) :: sfcemis(im)
      real(kind=kind_phys) :: hlwc_loc(ix,nvl)
      real(kind=kind_phys) :: acv(im)
      real(kind=kind_phys) :: acvb(im)
      real(kind=kind_phys) :: acvt(im)

!JR Define arrays that prevent copyin/copyout
      real(kind=kind_phys) :: vvel(ix,nvl)
      real(kind=kind_phys) :: swh_loc(ix,nvl)
      real(kind=kind_phys) :: hlw_loc(ix,nvl)
      real(kind=kind_phys) :: swhc_loc(ix,nvl)

      real*8 :: prdout(ix,ko3,pcoeff)      !JR I think this should be kind_phys not real*8
      real(kind=kind_phys) :: prslk(ix,nvl)
      real(kind=kind_phys) :: pgr(im)
      real(kind=kind_phys) :: suntim(im)
      real(kind=kind_phys) :: phil (ix,nvl)
      real(kind=kind_phys) :: phii (ix,nvl+1)
      real(kind=kind_phys) :: prsik(ix,nvl+1)
      real(kind=kind_phys) :: dpshc(im)
      real(kind=kind_phys) :: crtrh(3)
      real(kind=kind_phys) :: sinlat(im)
      real(kind=kind_phys) :: coslat(im)
      real(kind=kind_phys) :: prsshc
      real(kind=kind_phys) :: cubot
      real(kind=kind_phys) :: cutop
      real(kind=kind_phys) :: rqtk(im)
      real(kind=kind_phys) :: evbsa(im)
      real(kind=kind_phys) :: evcwa(im)
      real(kind=kind_phys) :: transa(im)
      real(kind=kind_phys) :: sbsnoa(im)
      real(kind=kind_phys) :: snowca(im)
      real(kind=kind_phys) :: soilm(im)
      real(kind=kind_phys) :: snohfa(im)
      real(kind=kind_phys) :: spfhmin(im)
      real(kind=kind_phys) :: spfhmax(im)
      real(kind=kind_phys) :: oro(im)
      real(kind=kind_phys) :: upd_mf(ix,nvl)
      real(kind=kind_phys) :: dwn_mf(ix,nvl)
      real(kind=kind_phys) :: det_mf(ix,nvl)
      real(kind=kind_phys) :: srunoff(im)
      real(kind=kind_phys) :: weasd(im)
!TODO:  Move these next 4 vars to new module d3d_def and (maybe) use d3d_zero to set
      real(kind=kind_rad) :: dt3dt(ix,nvl,6)
      real(kind=kind_rad) :: dq3dt(ix,nvl,9)
      real(kind=kind_rad) :: du3dt(ix,nvl,4)
      real(kind=kind_rad) :: dv3dt(ix,nvl,4)
      real(kind=kind_phys) :: dt_cool(im) !JR I/O arg to gbphys
      real(kind=kind_phys) :: flgmin(2)   !JR Input arg to gbphys
      real(kind=kind_phys) :: ccwf(2)     !JR Input arg to gbphys
      real(kind=kind_phys) :: I_Qrain(im) !JR I/O arg to gbphys
      real(kind=kind_phys) :: ifd(im)     !JR Input arg to gbphys
    
! Outputs from gbphys
      real(kind=kind_phys) :: gt0(ix,nvl)
      real(kind=kind_phys) :: gu0(ix,nvl)
      real(kind=kind_phys) :: gv0(ix,nvl)
      real(kind=kind_phys) :: gq0(ix,nvl,ntrac)
      real(kind=kind_phys) :: dkt(im,nvl-1)
      real(kind=kind_phys) :: sdiaga(im,nvl,11)
      real(kind=kind_phys) :: sdiagb(im,nvl,11)
      real(kind=kind_phys) :: zlvl(im)
      real(kind=kind_phys) :: t1(im)
      real(kind=kind_phys) :: q1(im)
      real(kind=kind_phys) :: u1(im)
      real(kind=kind_phys) :: v1(im)
      real(kind=kind_phys) :: CHH(IM)
      real(kind=kind_phys) :: CMM(IM)
      real(kind=kind_phys) :: DLWSFCI(IM)
      real(kind=kind_phys) :: ULWSFCI(IM)
      real(kind=kind_phys) :: DSWSFCI(IM)
      real(kind=kind_phys) :: USWSFCI(IM)
      real(kind=kind_phys) :: DTSFCI(IM)
      real(kind=kind_phys) :: DQSFCI(IM)
      real(kind=kind_phys) :: GFLUXI(IM)
      real*8 :: epi(im)
      real(kind=kind_phys) :: SMCWLT2(IM)
      real(kind=kind_phys) :: SMCREF2(IM)
      real(kind=kind_phys) :: sr(im)           ! in gsmddrive, initialize to 0??
      real(kind=kind_phys) :: DLWSFC_cc_dummy(im)
      real(kind=kind_phys) :: ULWSFC_cc_dummy(im)
      real(kind=kind_phys) :: SWSFC_cc_dummy(im)
      real(kind=kind_phys) :: XMU_cc_dummy(im)
      real(kind=kind_phys) :: DLW_cc_dummy(im)
      real(kind=kind_phys) :: DSW_cc_dummy(im)
      real(kind=kind_phys) :: SNW_cc_dummy(im)
      real(kind=kind_phys) :: LPREC_cc_dummy(im)
      real(kind=kind_phys) :: DUSFC_cc_dummy(im)
      real(kind=kind_phys) :: DVSFC_cc_dummy(im)
      real(kind=kind_phys) :: DTSFC_cc_dummy(im)
      real(kind=kind_phys) :: DQSFC_cc_dummy(im)
      real(kind=kind_phys) :: PRECR_cc_dummy(im)
      real(kind=kind_phys) :: Tref(im)
      real(kind=kind_phys) :: z_c(im)
      real(kind=kind_phys) :: c_0(im)
      real(kind=kind_phys) :: c_d(im)
      real(kind=kind_phys) :: w_0(im)
      real(kind=kind_phys) :: w_d(im)
      real(kind=kind_phys) :: dtdtr(ix,nvl)  ! jbao 2014 gfs phys add clear sky
! End outputs from gbphys

!JR Static variables that do not change
      integer, save :: jcap = 126 
      integer, save :: lsm  = 1
      integer, save :: lat = 1
      integer, save :: nlons(im) = (/200/)
      integer, save :: nst_fcst = 0
      integer, save :: thermodyn_id = 3
      integer, save :: sfcpress_id = 1
!     integer, save :: ncw(2) = (/20,120/)   ! from gfs input
      integer       :: ncw(2)                                
      integer, save :: global_lats_r(1) = 1
      integer, save :: lonsperlar(1) = 1 !JFM

      logical, save :: cnvgwd = .true.   ! change for gfs input 2014
      logical	    :: mom4ice = .false.
      logical, save :: trans_trac = .true. !  change for gfs input 2014
! jbao cal_pre: test gfsnamelist 2014 has true and there are
!                    are differences between latest version of 2014
!                    and original 2014 we got .false.
      logical, save :: cal_pre = .true.
      logical, save :: gen_coord_hybrid = .false.
      logical, save :: pre_rad = .false.
      logical, save :: lssav_cc_dummy = .false.
      logical, save :: lsfwd = .true.
      logical, save :: flipv = .true.     !  change gfs input 2014
      logical, save :: pdfcld = .false. ! from compns.f
      logical, save :: shcnvcw = .false. ! from compns.f
      logical, save :: redrag = .true. ! from paraconfigs
      logical, save :: hybedmf = .true.  ! .true. from paraconfigs
      logical, save :: dspheat = .true.  ! from paraconfigs

      real(kind=kind_phys), save :: rcs2(im) = (/1./)
      real(kind=kind_phys), save :: clstp = 1110.0
      real(kind=kind_phys), save :: fscav(ntrac-2) = 0.
! From Evelyns gfs runs change for gfs input 2014
! change for gfs input 2014
      real(kind=kind_phys), save :: dlqf(2) = (/0.0,0.0/)
      real(kind=kind_phys), save :: xkzm_m = 3.0
      real(kind=kind_phys), save :: xkzm_h = 1.0
      real(kind=kind_phys), save :: xkzm_s = 0.2
!     real(kind=kind_phys), save :: psautco(2) = (/6.0E-4,3.0E-4/) !  update from gfs 2014 input
      real(kind=kind_phys)       :: psautco(2) 
      real(kind=kind_phys), save :: prautco(2) = (/1.0E-4,1.0E-4/) !  update from gfs 2014 input
      real(kind=kind_phys), save :: evpco = 2.0E-5
      real(kind=kind_phys), save :: wminco(2) = (/1.0E-5,1.0E-5/)
      real(kind=kind_phys), save :: xt(im) = 0.
      real(kind=kind_phys), save :: xs(im) = 0.
      real(kind=kind_phys), save :: xu(im) = 0.
      real(kind=kind_phys), save :: xv(im)  = 0.
      real(kind=kind_phys), save :: xz(im) = 0.
      real(kind=kind_phys), save :: zm(im) = 0.
      real(kind=kind_phys), save :: xtts(im) = 0.
      real(kind=kind_phys), save :: xzts(im) = 0.
      real(kind=kind_phys), save :: d_conv(im) = 0.
      real(kind=kind_phys), save :: cgwf(2) = (/0.5,0.05/)   ! jbao 2014
      real(kind=kind_phys), save :: sup = 1.1                ! from compns.f
      real(kind=kind_phys), save :: adjtrc(ntrac) = 1.0    ! initialized as 1 in gloopb.f in gfs???
      real(kind=kind_phys), save :: oro_uf(im) = 0.
      real(kind=kind_phys), save :: rann(ix,nrcm) = 0.
      integer			 :: cactiv(ix)

      integer :: ipnout, itsout
      real :: pat_ipn
      real :: gfsq_ipn, taperoff

      ret = gptlstart ('column_physics')
      call physics_set_ipn_its (ipn, its)      ! Enable physics routines to call physics_get_ipn_its_mype
      vrbos = ipn == PrintIpnDiag
      prsi_s(1,1) = gfs_ps(ipn)
      if (ishal_cnv+ideep_cnv > 2) then
        cactiv(1)=conv_act(ipn)
      else
        cactiv(1)=0.
      end if

      do k=1,nvl
        prsi_s(1,k+1) = prsi_s(1,k) - gfs_dp(ipn,k)
        !tgs - BUG!!!    PRSL_S  (1,k) = 0.5*(PRSI_S(1,k)+PRSI_S(1,k+1))
        prsl_s(1,k) = gfs_p(ipn,k) * 0.001 ! units [cbar] 
        !TODO:  can this copy be avoided, eliminating tg?  
        tg(1,k)   = gfs_t(ipn,k)
        !TODO:  can these copies be avoided, eliminating qg?  
        qg(1,k,1) = gfs_q(ipn,k)
        qg(1,k,2) = gfs_oz(ipn,k)
        qg(1,k,3) = gfs_cld(ipn,k)
        !
        ! Copy data from FIM arrays to set up for GBPHYS call (below)
        !
        !TODO:  can this copy be avoided, eliminating ug and vg?  
        ug(1,k) = gfs_u(ipn,k)
        vg(1,k) = gfs_v(ipn,k)
      end do
    
      !----------------------------------------------------------------------
      if (CallRadiationNow .and. ra_lw_physics == 0) then ! Call radiation
      !----------------------------------------------------------------------
        fluxr(:,ipn,:) = 0.0
        do k=1,nvl
          !TODO:  can these copies be avoided, eliminating qg1?  
          qg1(1,k,1) = gfs_oz(ipn,k)
          qg1(1,k,2) = gfs_cld(ipn,k)
! activate the following line for radiation coupling
          qg1(1,k,2) = gfs_cld(ipn,k) + diaga(k,ipn) + diagb(k,ipn) 
          prsl(1,k) = prsl_s(1,k)
        end do

        do k=1,nvl+1
          prsi(1,k) = prsi_s(1,k)
        end do

        tisfc(:) = sfc_fld%tsea(ipn,1)
      
        if (chem_opt > 0) then
          sscal_rad(:,:) = sscal(:,ipn,:)
          asymp_rad(:,:) = asymp(:,ipn,:)
          ext_cof_rad(:,:) = ext_cof(:,ipn,:)
          extlw_cof_rad(:,:) = extlw_cof(:,ipn,:)
        else
          sscal_rad(:,:) = 999.
        end if

!----------------------------------------------------------------------
!   Call GFS radiation (longwave and shortwave)
!----------------------------------------------------------------------
!JR Copy things manually to avoid copyin/copyout
!JR If any of these things are output only, the copy can be eliminated
        phy_f3d_loc(1,:,:) = 0.0 ! phy_f3d(ipn,:,1,1,:)
        phy_f2d_loc(1,:) = 0.0 ! phy_f2d(ipn,1,:)
! jbao 2012 update
        sfcemis(:) = 1.
! jbao 2014 update: for passing the Lahey compiler
        slag(ipn,1)    = 0.
        sdec(ipn,1)    = 0.
        cdec(ipn,1)    = 0.
        swh(ipn,:,1,1) = 0.
        swh_loc(1,:)   = 0.
        hlw(ipn,:,1,1) = 0.
        hlw_loc(1,:)   = 0.
! jbao 2014

        ret = gptlstart ('gloopr')
        CALL GLOOPR                                           &
           (1,global_lats_r, lonsperlar,                  & ! jbao add ncld
           phour,deltim,                                       & ! jbao fcst hour
           XLON(ipn,1),XLAT(ipn,1),COSZDG(ipn,1),              &
           flx_fld%COSZEN(ipn,1),                              & ! COSZEN is output and used in gbphys
           sfc_fld%SLMSK(ipn,1),                               &
           sfc_fld%SHELEG(ipn,1),                              & ! jbao new gfs 2014 sheleg
           sfc_fld%SNCOVR(ipn,1),                              & ! jbao new gfs needs sncovr
           sfc_fld%SNOALB(ipn,1),                              & ! jbao new gfs nees snoalb
           sfc_fld%ZORL(ipn,1),                                &
           sfc_fld%TSEA(ipn,1),                                &
           HPRIME(1,ipn,1),                                    &
           SFALB(ipn,1),                                       & ! SFALB is set in gloopr but then set to 0 before call gbphys
           sfc_fld%ALVSF(ipn,1),                               &
           sfc_fld%ALNSF(ipn,1),                               &
           sfc_fld%ALVWF(ipn,1),                               &
           sfc_fld%ALNWF(ipn,1),                               &
           sfc_fld%FACSF(ipn,1),                               &
           sfc_fld%FACWF(ipn,1),sfc_fld%CV(ipn,1),             &
           sfc_fld%CVT(ipn,1),                                 &
           sfc_fld%CVB(ipn,1),swh_loc,swhc_loc,                & !jbao 2014 update swhc
           hlw_loc,hlwc_loc,                                   & !jbao 2014 update 
           flx_fld%SFCNSW(ipn,1),                              &
           flx_fld%SFCDLW(ipn,1),                              & ! SWH,HLW,SFCNSW,SFCDLW output and used in gbphys
           sfc_fld%FICE(ipn,1) ,                               &
           TISFC,                                              & ! jbao new gfs needs tisfc
           flx_fld%SFCDSW(ipn,1),                              & ! FOR SEA-ICE - XW Nov04, SFCDSW output and used in gbphys
           sfcemis,                                            & !  gloopr 2012
           flx_fld%SFCUSW(ipn,1),                              & ! hli 09/2014
           flx_fld%SFCULW(ipn,1),                              & ! hli 09/2014
           flx_fld%TOPDSW(ipn,1),                              & ! hli 09/2014
           flx_fld%TOPUSW(ipn,1),                              & ! hli 09/2014
           flx_fld%TOPULW(ipn,1),                              & ! hli 09/2014
           flx_fld%TSFLW(ipn,1),FLUXR(:,ipn,1),                & ! jbao new gfs does not need cldcov
           phy_f3d_loc,phy_f2d_loc,SLAG(ipn,1) ,SDEC(ipn,1),   & ! jbao 2014 new gfs needs phy_f2d 
           CDEC(ipn,1),KDT,                                    & 
           0.0D0, prsl(1,1),prsi(1,1),prslk(1,1),tg(1,1),      &
           qg(1,1,1),qg1(1,1,1),                               &
           yyyymmddhhmm)
        ret = gptlstop ('gloopr')
        
!JR Copy things manually to avoid copyin/copyout
!JR If any of these things are input only, the copy can be eliminated
        swh(ipn,:,1,1) = swh_loc(1,:)
        hlw(ipn,:,1,1) = hlw_loc(1,:)
        phy_f3d(ipn,:,1,1,:) = phy_f3d_loc(1,:,:)
! jbao gfs 2014 update
        swhc(ipn,:,1,1) = swhc_loc(1,:)
        hlwc(ipn,:,1,1) = hlwc_loc(1,:)
        phy_f2d(ipn,1,:) = phy_f2d_loc(1,:)
        !if (chem_opt.eq.0) then
        !  sscal(:,ipn,:)=sscal_rad(:,:)
        !  asymp(:,ipn,:)=asymp_rad(:,:)
        !  ext_cof(:,ipn,:)=ext_cof_rad(:,:)
        !  extlw_cof(:,ipn,:)=extlw_cof_rad(:,:)
        !endif
      else if (ra_lw_physics > 0 ) then
!
! need to get the following from wrf physics
!  swh and hlw are tendencies from radiation-wrf routines. swh includes
!  everything
!
        swh(ipn,:,1,1) = rthratensw(1,:,ipn)
        hlw(ipn,:,1,1) = rthratenlw(1,:,ipn)
      else
! initialize for restart
        swhc(ipn,:,1,1) = 0.
        hlwc(ipn,:,1,1) = 0.
      end if

      !
      ! Local set up for GBPHYS
      !
      do k=1,nvl
! jbao may 2012
        prsl(1,k) = 1000.0*prsl_s(1,k)
! jbao may 2012
      end do
      do k=1,nvl+1
! jbao may 2012
        prsi(1,k) = 1000.0*prsi_s(1,k)
      end do
      pgr(:) = 1000.0*prsi_s(1,1)

!  change for gfs input 2014
      suntim(:) = 0.0
      phil(:,:) = 0.0
      phii(:,:) = 0.0
      prsik(:,:) = 0.0
      dpshc(1) = 0.3 * prsi(1,1) ! jbao new GFS physics as of Feb 2010
      if (dx.lt.25000.) then
        crtrh(:) = 0.9  ! 2014 from gfs input 0.85
        ncw(1) = 20   ! from gfs input
        ncw(2) = 120
        psautco(1) = 6.0E-4
        psautco(2) = 3.0E-4
      else
        crtrh(:) = 0.85 ! 2014 from gfs input 0.85
!       crtrh(3) = 0.85 ! testing
!       crtrh(2) = 0.85 ! testing
!       crtrh(1) = 0.80 ! testing
        ncw(1) = 50   ! from gfs input
        ncw(2) = 150
        psautco(1) = 4.0E-4
        psautco(2) = 4.0E-4
      end if

! change for gfs input 2014
      tisfc(:)  = sfc_fld%tsea(ipn,1)
      sinlat(:) = sin(xlat(ipn,1))
      coslat(:) = cos(xlat(ipn,1))
      prsshc    = prsi_s(1,1)
      prdout(1,:,:) = ozplout(:,:,ipn)
      flx_fld%PSMEAN(ipn,1) = prsshc
      flx_fld%PSURF(ipn,1)  = prsshc
!TODO:  call flx_init() instead as in GFS r3038 do_physics_one_step.f and 
!TODO:  remove some of these statements...  
      flx_fld%GESHEM(ipn,1) = 0.0
      flx_fld%RAINC(ipn,1)  = 0.0
      cubot                 = 0.0
      cutop                 = 0.0
      rqtk(:)               = 0.0  !JR output from gbphys: maybe no init needed?
!
! if running wrfphysics, these are data needed (?)
!
      if (skip_cu_physics) then
        flx_fld%RAINC(ipn,1)  = phys2dwrf (ipn,6)
        cubot                 = phys2dwrf (ipn,7)
        cutop                 = phys2dwrf (ipn,8)
!       if(flx_fld%RAINC(ipn,1).gt.0)write(6,*)'do_phys',flx_fld%RAINC(ipn,1), &
!                                    cubot,cutop
      end if

      flx_fld%DUSFC(ipn,1)  = 10.0
      flx_fld%DVSFC(ipn,1)  = 10.0
      flx_fld%DTSFC(ipn,1)  = 1.0
      flx_fld%DQSFC(ipn,1)  = 0.0
      flx_fld%GFLUX(ipn,1)  = 0.0
      flx_fld%RUNOFF(ipn,1) = 0.0
      flx_fld%EP(ipn,1)     = 0.0
      flx_fld%CLDWRK(ipn,1) = 0.0
      flx_fld%DUGWD(ipn,1)  = 0.0
      flx_fld%DVGWD(ipn,1)  = 0.0
      DT3DT(:,:,:)     = 0.0
      DQ3DT(:,:,:)     = 0.0
      DU3DT(:,:,:)     = 0.0
      DV3DT(:,:,:)     = 0.0
      EVBSA(:)  = 0.0
      EVCWA(:)  = 0.0
      TRANSA(:) = 0.0
      SBSNOA(:) = 0.0
      SNOWCA(:) = 0.0
      SNOHFA(:) = 0.0
      SPFHMAX(:) = 0.0
      SPFHMIN(:) = 1.e10
      oro(:) = 0.0
      upd_mf(:,:) = 0.0
      dwn_mf(:,:) = 0.0
      det_mf(:,:) = 0.0
      sdiaga(:,:,:) = 0.0
      sdiagb(:,:,:) = 0.0
      srunoff(:) = 0.0
! jbao may 2012
      soilm(:) = 0.0
      weasd(:) = sfc_fld%SHELEG(ipn,:)
      sfcemis(:) = 1.
   
!----------------------------------------------------------------------
!   Call all other (non-radiation) GFS physics parameterizations
!----------------------------------------------------------------------
!JR Copy things manually to avoid copyin/copyout
!JR If any of these things are output only, the copy can be eliminated
!sam test orig  vvel(1,:)          = 1000.0*gfs_dpdt(ipn,:)
      vvel(1,:)          = 1000.0*gfs_dpdt(ipn,:)
      swh_loc(1,:)       = swh(ipn,:,1,1)
      hlw_loc(1,:)       = hlw(ipn,:,1,1)
! jbao new for 2014
      swhc_loc(1,:)       = swhc(ipn,:,1,1)
      hlwc_loc(1,:)       = hlwc(ipn,:,1,1)
! gg : advective forcing for GF
      forcet(1,:)         = phys3dwrf(:,ipn,7)
      forceq(1,:)         = phys3dwrf(:,ipn,3)
      if (ishal_cnv+ideep_cnv >2) cactiv(1)=conv_act(ipn)

! jbao new for 2014
      phy_f3d_loc(1,:,:) = phy_f3d(ipn,:,1,1,:)
      phy_f2d_loc(1,:)   = phy_f2d(ipn,1,:)

!SMS$IGNORE BEGIN
      if (vrbos) then
        print *,ipn, &    
           'SMC=',(sfc_fld%SMC(k,ipn,1),k=1,LSOIL),'STC=',(sfc_fld%STC(k,ipn,1),k=1,LSOIL),  &
           'SOILTYP=',sfc_fld%STYPE(ipn,1),  &
           'VEGTYP=',sfc_fld%vtype(ipn,1),'SLMSK=',sfc_fld%slmsk(ipn,1), &
           'VEGFRAC=',sfc_fld%vfrac(ipn,1),'CANOPY=',sfc_fld%CANOPY(ipn,1), &
           'ZORL=',sfc_fld%ZORL(ipn,1),'SNOW COVER=',sfc_fld%sncovr(ipn,1)
      end if
!SMS$IGNORE END

!JR Provide initial values for these unused things
      acv(:)     = inf
      acvb(:)    = inf
      acvt(:)    = inf
      flgmin(:)  = inf
      ccwf(:)    = inf
      I_Qrain(:) = inf
      dt_cool(:) = inf
      ifd(:)     = inf
      if (coupled) then
        mom4ice = .true.
      else
        mom4ice = .false.
      end if


      ret = gptlstart ('gbphys')
      call GBPHYS(cactiv,ipn,im,ix,nvl,lsoil,ntrac,ncld,ntoz,ntcw,           &
           nmtvr,nrcm,ko3,lonf,latg,jcap,num_p3d,num_p2d,                    &
           kdt,lat,me,pcoeff,nlons,ncw,flgmin,crtrh,cdmbgwd,                 &
           ccwf,dlqf,clstp,cgwf,dtp,dtf,fhour,solhr,                         &
           slag(ipn,1),sdec(ipn,1),cdec(ipn,1),sinlat,coslat,pgr,ug,vg,      &
           tg,qg,vvel,prsi,prsl,prslk,prsik,phii,phil,                       &
           rann,prdout,poz,dpshc,HPRIME(:,ipn,1),XLON(ipn,1),XLAT(ipn,1),    &
           sfc_fld%SLOPE(ipn,1),sfc_fld%SHDMIN(ipn,1),sfc_fld%SHDMAX(ipn,1), &
           sfc_fld%SNOALB(ipn,1),sfc_fld%TG3(ipn,1),sfc_fld%SLMSK(ipn,1),    &
           sfc_fld%VFRAC(ipn,1),sfc_fld%VTYPE(ipn,1),sfc_fld%STYPE(ipn,1),   &
           sfc_fld%UUSTAR(ipn,1),sfc_fld%STRESS(ipn,1),                      &
           oro,oro_uf,flx_fld%COSZEN(ipn,1),                                 &
           flx_fld%SFCDSW(ipn,1),flx_fld%SFCNSW(ipn,1),                      &
           flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1),sfcemis,SFALB(ipn,1),  &
           SWH_loc,swhc_loc,HLW_loc,hlwc_loc,ras,ideep_cnv,pre_rad,          &
           ldiag3d,lssav,lssav_cc_dummy,                                     &
           xkzm_m,xkzm_h,xkzm_s,psautco,prautco,evpco,wminco,sup,            &
!          pdfcld,shcnvcw,flipv,cnvgwd,shal_cnv,                             &	! shal_cnv is replaced by ishal_cnv
           pdfcld,shcnvcw,flipv,cnvgwd,ishal_cnv,                            &
           redrag,hybedmf,dspheat,cal_pre,                                   &
           mom4ice,trans_trac,nst_fcst,fscav,               		     &
           thermodyn_id, sfcpress_id, gen_coord_hybrid, adjtrc,              &
!  --  input/outputs:
           sfc_fld%HICE(ipn,1),sfc_fld%FICE(ipn,1),TISFC,sfc_fld%TSEA(ipn,1),&
           sfc_fld%TPRCP(ipn,1),sfc_fld%CV(ipn,1),sfc_fld%CVB(ipn,1),        &
           sfc_fld%CVT(ipn,1), sfc_fld%SRFLAG(ipn,1),                        &
           sfc_fld%SNWDPH(ipn,1),WEASD,sfc_fld%SNCOVR(ipn,1),		     &
           sfc_fld%ZORL(ipn,1),sfc_fld%CANOPY(ipn,1),                        &
           sfc_fld%FFMM(ipn,1),sfc_fld%FFHH(ipn,1),sfc_fld%F10M(ipn,1),      &
           srunoff,evbsa,evcwa,snohfa,transa,sbsnoa,                         &
           snowca,SOILM,flx_fld%TMPMIN(ipn,1),flx_fld%TMPMAX(ipn,1),         &
           flx_fld%DUSFC(ipn,1),flx_fld%DVSFC(ipn,1),flx_fld%DTSFC(ipn,1),   &
           flx_fld%DQSFC(ipn,1),flx_fld%GESHEM(ipn,1),flx_fld%GFLUX(ipn,1),  &
           flx_fld%DLWSFC(ipn,1),flx_fld%ULWSFC(ipn,1),suntim,               &
           flx_fld%RUNOFF(ipn,1),flx_fld%EP(ipn,1),flx_fld%CLDWRK(ipn,1),    &
           flx_fld%DUGWD(ipn,1),flx_fld%DVGWD(ipn,1),                        &
           flx_fld%PSMEAN(ipn,1),flx_fld%RAINC(ipn,1),spfhmin,spfhmax,       &
           dt3dt,dq3dt,du3dt,dv3dt,acv(1),acvb(1),acvt(1),                   &
           sfc_fld%SLC(:,ipn,1),sfc_fld%SMC(:,ipn,1),sfc_fld%STC(:,ipn,1),   &
           upd_mf,dwn_mf,det_mf,phy_f3d_loc,phy_f2d_loc,             	     &
           dlwsfc_cc_dummy,ulwsfc_cc_dummy,dtsfc_cc_dummy,swsfc_cc_dummy,    &
           dusfc_cc_dummy, dvsfc_cc_dummy, dqsfc_cc_dummy,precr_cc_dummy,    &
           xt,xs,xu,xv,xz,zm,xtts,xzts,d_conv,ifd,dt_cool,I_Qrain,           &
!  ---  outputs:
           gt0,gq0,gu0,gv0,sfc_fld%T2M(ipn,1),sfc_fld%Q2M(ipn,1),            &
           flx_fld%U10M(ipn,1),flx_fld%V10M(ipn,1),                          &
           zlvl,flx_fld%PSURF(ipn,1),flx_fld%HPBL(ipn,1),                    &
           flx_fld%PWAT(ipn,1),t1,q1,u1,v1,                                  &
           chh,cmm,dlwsfci,ulwsfci,dswsfci,uswsfci,                          &
           dtsfci,dqsfci,gfluxi,epi,smcwlt2,smcref2,                         &
           sr,dkt,                                                           &
           xmu_cc_dummy,dlw_cc_dummy,dsw_cc_dummy,snw_cc_dummy,              &
           lprec_cc_dummy,tref, z_c, c_0, c_d, w_0, w_d,                     &
           rqtk   , dtdtr, flx_fld%EVAP(ipn,1), flx_fld%HFLX(ipn,1),         &
           skip_cu_physics,skip_mp_physics,cubot,cutop,                      &
           forcet,forceq,sdiaga,sdiagb)
      ret = gptlstop ('gbphys')
!JR Copy things manually to avoid copyin/copyout.
!JR If any of these things were input only, the copy can be eliminated
      swh(ipn,:,1,1)       = swh_loc(1,:)
      hlw(ipn,:,1,1)       = hlw_loc(1,:)
      phy_f3d(ipn,:,1,1,:) = phy_f3d_loc(1,:,:)
      phy_f2d(ipn,1,:)     = phy_f2d_loc(1,:)
      if (ishal_cnv + ideep_cnv >2 .and. cactiv(1) > 0)then
        conv_act(ipn) = conv_act(ipn)+1
      else
        conv_act(ipn) = 0
      endif

! fill diagnostic output array diagb with whatever you would like to see....
! you could also fill diaga, unless you would want to do this somewhere else
!
     do k=1,nvl-1
       diaga(k,ipn) = sdiaga(1,k,1)
       diagb(k,ipn) = sdiagb(1,k,1)
     enddo

! Set tendency arrays from physics

      do k=1,nvl
!SMS$IGNORE BEGIN
        if (gq0(1,k,1) < 0.) print*,'qv negative',gfs_q(ipn,k),gq0(1,k,1),ipn,k,prsl(1,k)
!SMS$IGNORE END
        if (.not. sppt) then
          gfs_u(ipn,k)   = gu0(1,k)
          gfs_v(ipn,k)   = gv0(1,k)
          gfs_t(ipn,k)   = gt0(1,k)
          gfs_q(ipn,k)   = gq0(1,k,1)
          gfs_oz(ipn,k)  = gq0(1,k,2)
          gfs_cld(ipn,k) = gq0(1,k,3)
        else
          if (DoSPPT) then
            if (k >= 3 .and. k <= 43) then
              taperoff = 1.0 
            else if (k > 45) then
              taperoff = 0.0
            else if (k == 1 .or. k == 45) then
              taperoff = 0.125 
            else if (k == 2 .or. k == 44) then
              taperoff = 0.5 
            endif
            pat_ipn = pattern(ipn)*taperoff+1.0
        
            gfs_u(ipn,k)   = gfs_u(ipn,k)+(gu0(1,k)-gfs_u(ipn,k))*pat_ipn
            gfs_v(ipn,k)   = gfs_v(ipn,k)+(gv0(1,k)-gfs_v(ipn,k))*pat_ipn
            gfs_oz(ipn,k)  = gfs_oz(ipn,k)+(gq0(1,k,2)-gfs_oz(ipn,k))*pat_ipn
            gfs_cld(ipn,k) = gfs_cld(ipn,k)+(gq0(1,k,3)-gfs_cld(ipn,k))*pat_ipn
            gfsq_ipn = gfs_q(ipn,k)+(gq0(1,k,1)-gfs_q(ipn,k))*pat_ipn
            if (gfsq_ipn > 0.0) then
              gfs_q(ipn,k)   = gfsq_ipn
              gfs_t(ipn,k)   = gfs_t(ipn,k)+(gt0(1,k)-gfs_t(ipn,k)-dtdtr(1,k))*pat_ipn+dtdtr(1,k)
            else
              gfs_q(ipn,k)   = gq0(1,k,1)
              gfs_t(ipn,k)   = gt0(1,k)
            endif
          else
            gfs_u(ipn,k)   = gu0(1,k)
            gfs_v(ipn,k)   = gv0(1,k)
            gfs_t(ipn,k)   = gt0(1,k)
            gfs_q(ipn,k)   = gq0(1,k,1)
            gfs_oz(ipn,k)  = gq0(1,k,2)
            gfs_cld(ipn,k) = gq0(1,k,3)
          endif
        endif
      end do

! jbao 2014 
      sfc_fld%SHELEG(ipn,:) = weasd(:)
! jbao 2014

      if (chem_opt > 0 .or. skip_cu_physics  .or. skip_mp_physics) then
        call do_physics_one_step_chem(kdt,ipn,skip_chem,skip_cu_physics,skip_mp_physics,dkt,dq3dt,dt3dt)
      endif
      call physics_reset_ipn (ipn)      ! This thread no longer working on this ipn
      ret = gptlstop ('column_physics')
    end subroutine column_physics
  end subroutine do_physics_one_step

  subroutine ozone_init ()
    use units, only: getunit, returnunit

    integer :: unitno, ios
    integer :: ioz, koz, noz
    integer :: pl_coeff_oz
    integer :: levozp
    integer :: ipn
    integer :: ret
    real :: pl_lat(latsoz)
    real :: pl_pres(ko3)
    real :: tempin(latsoz)

    ret = gptlstart ('ozone_init')
    if (allocated (ozplout)) then
      write(6,*)'ozone_init: ozplout is already allocated!'
      stop
    end if

    allocate (ozplout(ko3,pcoeff,nip))
!JWB: need to refne...   if(CallRadiationNow) then ! ozone sources and sinks update
    ozplin(:,:,:,:) = inf
!JR Use OMP here for first-touch. Initialize the unused halo as well.
!$OMP PARALLEL DO
    do ipn=ips,ihe
      ozplout(:,:,ipn) = inf
    end do
!$OMP END PARALLEL DO
    poz(:) = inf

!SMS$SERIAL(<ozplin,poz,timeoz,pl_time,latsozp,out> : default=ignore) BEGIN
    unitno = getunit (76)
    if (unitno < 0) then
      print*,'ozone_init: getunit failed for ozone phys input (big endian). Stopping'
      stop
    end if

    open (unitno, file='global_o3prdlos.f77', form='unformatted', action='read', iostat=ios)
    if (ios.ne.0) then
      print*,'ozone_init: open failed for ozone (big endian). Stopping'
      stop
    end if

    rewind (unitno)
    read (unitno,iostat=ios) pl_coeff_oz, latsozp, levozp, timeoz, pl_lat, pl_pres, pl_time
    if (ios /= 0) then
      write(6,*)'ozone_init: error reading ozone header. Stopping'
      stop
    end if

    if (levozp /= ko3) then
      write(6,*)'ozone_init: levozp=',levozp,' not = ko3=',ko3
      stop
    end if
      
    if (pl_coeff_oz /= pcoeff) then
      write(6,*)'ozone_init: pl_coeff_oz=',pl_coeff_oz,' not = pcoeff=',pcoeff
      stop
    end if

    if (timeoz /= oztime-1) then
      write(6,*)'ozone_init: timeoz=',timeoz,' not = oztime=',oztime
      stop
    end if

    do ioz=1,timeoz
      do noz=1,pcoeff
        do koz=1,ko3
          read(unitno,iostat=ios) tempin
          if (ios /= 0) then
            write(6,*)'do_physics_one_step: error reading ozone info.'
            write(6,*)'ioz,noz,koz=',ioz,noz,koz,' Stopping.'
            stop
          end if
          ozplin(:,koz,noz,ioz) = tempin(:)
        end do
      end do
    end do

    do koz=1,ko3
      poz(koz) = log(100.0*pl_pres(koz))
    end do

    close (unitno)
    call returnunit (unitno)
!SMS$SERIAL END
!sms$compare_var(ozplin, "do_physics.F90 - ozplin ")
    ret = gptlstop ('ozone_init')
  end subroutine ozone_init

  subroutine ozone_calc (lats_node_r, xlat, phour)
    integer, intent(in) :: lats_node_r
    real (kind=kind_rad), intent(in) :: phour
!sms$distribute (dh,1) begin
    real (kind=kind_rad), intent(in) :: xlat(nip,lats_node_r)
!sms$distribute end

    real, parameter:: con_pi = 3.14159265
    real, parameter :: raddeg=180.0/con_pi

    integer :: ilatoz
    integer :: ioz,koz,noz
    integer :: ipn
    integer :: j,n1,n2
    integer :: jdat(8)
    integer :: jdow, jdoy, jday
    real*8 :: rinc(5), rjday
    real*8 :: tx1, tx2
    integer :: ret

    ret = gptlstart ('ozone_calc')
    RINC=0.
    RINC(2)=PHOUR
    CALL W3MOVDAT(RINC,IDAT,JDAT)
      
    jdow = 0
    jdoy = 0
    jday = 0
    call w3doxdat(jdat,jdow,jdoy,jday)
!  print *,'rjday, jdoy ,jdat(5)',jday, jdoy ,jdat(5)
    rjday = jdoy + jdat(5) / 24.
    IF (RJDAY .LT. PL_time(1)) RJDAY = RJDAY+365.
 
    n2 = timeoz + 1
    do j=1,timeoz
      if (rjday .lt. pl_time(j)) then
        n2 = j
        exit
      endif
    end do
    n1 = n2 - 1
    if (n1 <= 0)     n1 = n1 + timeoz
    if (n2 > timeoz) n2 = n2 - timeoz

    tx1 = (pl_time(n2) - rjday) / (pl_time(n2) - pl_time(n1))
    tx2 = 1.0 - tx1

!$OMP PARALLEL DO PRIVATE (ilatoz, noz, koz)
    do ipn=ips,ipe
      ilatoz = min(latsozp, int((90.0+xlat(ipn,1)*raddeg)/5.0) + 1 )
      do noz=1,pcoeff
        do koz=1,ko3
          ozplout(koz,noz,ipn) = tx1*ozplin(ilatoz,koz,noz,n1) + tx2*ozplin(ilatoz,koz,noz,n2)
        end do
      end do
    end do
!$OMP END PARALLEL DO

!sms$compare_var(ozplout, "do_physics.F90 - ozplout ")
    ret = gptlstop ('ozone_calc')
  end subroutine ozone_calc

!nwang added subroutine sppt_init
  subroutine sppt_init ()
    use units, only: getunit, returnunit

    integer :: ios
    allocate (pattern(nip))

!sms$serial begin
    pat_lun = getunit ()
    if (pat_lun < 0) then
      print*,'sppt_init: getunit failed for pattern data input. Stopping'
      stop
    end if

    open (pat_lun, file='patterns_icos.dat', form='unformatted', action='read', iostat=ios)
    if (ios.ne.0) then
      print*,'sppt_init: open failed for pattern data input. Stopping'
      stop
    end if
!sms$serial end

  end subroutine sppt_init

end module module_do_physics_one_step
