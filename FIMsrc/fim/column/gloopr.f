       subroutine gloopr (lats_nodes_r,global_lats_r,
     &     lonsperlar,
     &     phour,deltim,
     &     xlon,xlat,coszdg,coszen,
     &     slmsk,snwdph,sncovr,snoalb,zorl,tsea,hprime,
     &     sfalb,
     &     alvsf,alnsf ,alvwf ,alnwf,facsf ,facwf,cv,cvt ,
     &     cvb  ,swh,swhc,hlw,hlwc,sfcnsw,sfcdlw,
     &     fice ,tisfc, sfcdsw, sfcemis,                ! for sea-ice - xw nov04
     &     sfcusw, sfculw, topdsw, topusw, topulw,      ! hli 09/2014
     &     tsflw,fluxr, phy_f3d,phy_f2d,slag,sdec,cdec,kdt,
     &     global_times_r,prsl,prsi,prslk,gt,gr,gr1,
     &     yyyymmddhhmm)
!!
!
!   program modification log:
!      apr 2012  - y-t hou     -  modified radiation calling interval
!             parameters fhswr/fhlwr.  changed units from hour to second
!             (in multiple of model time step length), thus radiation can
!             be invoked more frequently than the minimum of 1 hr limit.
!                              -  moved radiation initialization part to
!             the beginning of model initialization process. in the time
!             loop, added a subroutine "radupdate" to updae radiation related
!             external data sets (such as aerosols, gases, astronomy, etc.)
!             moved computation of coszen to inside subrouine "grrad"
!      nov 2012  - y-t hou     - modified many physics/radiation control
!             parameters through module 'physpara' which contains the use of
!             module 'machine'.
!      may 2013 s. mooorthi - removed fpkapx
!      jul 2013 r sun  - to use pdf cloud and convective cld cover and water 

      use physpara
      use physcons, fv => con_fvirt, rerth => con_rerth, rocp =>con_rocp

      use module_radiation_driver,   only : radupdate, grrad

      use module_radsw_parameters,  only : topfsw_type, sfcfsw_type
      use module_radlw_parameters,  only : topflw_type, sfcflw_type
      USE FUNCPHYS             ,     ONLY : fpkap
!
!! ---  for optional spectral band heating outputs
!!    use module_radsw_parameters,   only : nbdsw
!!    use module_radlw_parameters,   only : nbdlw
!
!JR Keep resol_def for now to ensure correct settings. BUT, num_p3d exists
!JR in 3 places! It is in resol_def, gis_phy, and module_control
!TODO FIX THIS!
      use resol_def, only: num_p3d, nfxr, levr, ntcw, ncld, ntoz, latr,
     &                     lonr, levs, ntrac, jcap1, lonrx, lonr, lots,
     &                     lotd, levp1
      use layout1
! jbao not used in fim      use layout_grid_tracers
!JR Comment out gg_def and use locally-defined sinlat_r and coslat_r.
!JR Necessary since otherwise in threaded mode the allocate fails across
!JR multiple threads
!      use gg_def
!end jbao for fim 
      use vert_def
      use date_def
      use namelist_def
      use fimnamelist, only : fhswr, fhlwr, isubc_lw, isubc_sw,
     &                        n3dzhaocld, n3dcldpdf
      use coordinate_def					! hmhj
      use tracer_const						! hmhj
! jbao old gloopr cldcov passed in
!      use d3d_def , only : cldcov
!
      use mersenne_twister, only : random_setseed, random_index,
     &                             random_stat
!
      implicit none
!     include 'mpif.h'
!

      real (kind=kind_phys), parameter :: qmin =1.0e-10
      real (kind=kind_evod), parameter :: typical_pgr = 95.0
      real (kind=kind_evod), parameter :: cons0 = 0.0

!  --- ...  inputs:
      integer, intent(in) ::  lats_nodes_r, 
     &                       global_lats_r(latr), lonsperlar(latr)

      real (kind=kind_phys), dimension(lonr,lats_node_r), intent(in) :: &
     &                       xlon, xlat, slmsk, snwdph, zorl, tsea,     &
     &                       alvsf, alnsf, alvwf, alnwf, facsf, facwf,  &
     &                       cv, cvt, cvb, fice, tisfc, sncovr, snoalb

      real (kind=kind_phys), intent(in) ::                              &
     &                    HPRIME(    1,lonr,lats_node_r), phour,deltim
!
!  --- ...  input and output:
      real (kind=kind_phys), intent(inout) ::                           &
     &                    phy_f3d(lonr,levs,4,lats_node_r),             &
     &                    phy_f2d(lonr,3,lats_node_r)
      integer              ipt_ls                                       ! hmhj
      real(kind=kind_evod) reall                                        ! hmhj
      real(kind=kind_evod) rlcs2(jcap1)                                 ! hmhj

      real (kind=kind_phys), intent(inout) ::                           &
!jbao change 33 to match nfxr 2014  
     &                    fluxr (nfxr,LONR,LATS_NODE_R)

! orig     &                    fluxr (lonr,nfxr,lats_node_r)
! jbao orig ,nfxr isn't defined til later but=27  and added cldcov not
! 27, but levs  &                    fluxr (NFXR,LONR,LATS_NODE_R)
      real (kind=kind_phys) CLDCOV(LEVS,lonr,lats_node_r)


      integer, intent(in) :: kdt

!  --- ...  outputs:
! jbao old gloopr had global_times_r defined as the following:
             real(kind=kind_evod) global_times_r

!jbao orig      real(kind=kind_evod), intent(inout) ::                            &
!jbao orig     &                    global_times_r(latr,nodes)

      real (kind=kind_phys), intent(inout) ::                           &
     &                    swh(lonr,levs,lats_node_r),                   &
     &                    swhc(lonr,levs,lats_node_r),                  &
     &                    hlwc(lonr,levs,lats_node_r),                  &
     &                    hlw(lonr,levs,lats_node_r)

      real (kind=kind_phys),dimension(lonr,lats_node_r), intent(inout)::&
     &                    coszdg, coszen, sfcnsw, sfcdlw, tsflw,        &
     &                    sfcdsw, sfalb, sfcemis

!hli radiation variables 09/2014
      real (kind=kind_phys),dimension(lonr,lats_node_r), intent(inout)::&
     &                    sfcusw, sfculw, topdsw, topusw, topulw

      real (kind=kind_phys), intent(inout) :: slag, sdec, cdec

!! --- ...  optional spectral band heating rates
!!    real (kind=kind_phys), optional, intent(out) ::                   &
!!   &                 htrswb(ngptc,levs,nbdsw,nblck,lats_node_r),      &
!!   &                 htrlwb(ngptc,levs,nbdlw,nblck,lats_node_r)

!  --- ...  locals:
      real(kind=kind_evod) ::                                           &
     &                    for_gr_r_1(LONRX,LOTS,LATS_DIM_R),            &
     &                    dyn_gr_r_1(lonrx,lotd,lats_dim_r),            &
     &                    for_gr_r_2(LONR ,LOTS,LATS_DIM_R),            &
     &                    dyn_gr_r_2(lonr ,lotd,lats_dim_r)

      real(kind=kind_phys) :: prsl(ngptc,levs),  prslk(ngptc,levs),     &
     &                        prsi(ngptc,levp1), prsik(ngptc,levp1),    &
     &                        hlwc_v(ngptc,levs), swhc_v(ngptc,levs),   &
     &                        hlw_v(ngptc,levs), swh_v(ngptc,levs)

      real (kind=kind_phys) ::                                          &
     &                      gt(ngptc,levs),   gq(ngptc),                &
! jbao orig old has ngptc,levs,2 ntrac -1     &
! gr(NGPTC,levs), gr1(NGPTC,levs,NTRAC-1),  &
     &                      gr(ngptc,levs),   gr1(ngptc,levs,2),        &
     &                      gtv(ngptc,levs)

      real (kind=kind_phys), allocatable ::  sumq(:,:), xcp(:,:),       &
     &                                       gtvx(:,:), gtvy(:,:),      &
     &                                                  gd(:,:),        &
!    &                                       vvel(:,:), gd(:,:),        &
     &                                       gu(:,:),   gv1(:,:),       &
     &                                       gphi(:),   glam(:)

      real (kind=kind_phys) :: f_ice(ngptc,levs), f_rain(ngptc,levs),   &
     &                         r_rime(ngptc,levs)

      real (kind=kind_phys) :: deltaq(ngptc,levs),cnvw(ngptc,levs)
      real (kind=kind_phys) :: cnvc(ngptc,levs ), vvel(ngptc,levs)

      real (kind=kind_phys) :: cldcov_v(ngptc,levs), fluxr_v(ngptc,33)
      real (kind=kind_phys) :: flgmin_v(ngptc), work1, work2

      real (kind=kind_phys) :: coszen_v(ngptc), coszdg_v(ngptc),        &
     &                         sinlat_v(ngptc), coslat_v(ngptc)

      real (kind=kind_phys) :: rinc(5), dtsw, dtlw, solcon, raddt, solhr

      integer :: njeff, lon, lan, lat, lons_lat, lmax
      integer :: idat(8), jdat(8)

!  ---  variables of instantaneous calculated toa/sfc radiation fluxes
      type (topfsw_type), dimension(ngptc) :: topfsw
      type (sfcfsw_type), dimension(ngptc) :: sfcfsw

      type (topflw_type), dimension(ngptc) :: topflw
      type (sfcflw_type), dimension(ngptc) :: sfcflw

!  ---  variables used for random number generator (thread safe mode)
      type (random_stat) :: stat
      integer :: numrdm(lonr*latr*2), ixseed(lonr,lats_node_r,2)
      integer :: ipseed, icsdlw(ngptc), icsdsw(ngptc)
      integer, parameter :: ipsdlim = 1.0e8      ! upper limit for random seeds

!  --- ...  control parameters:
!           (some of the them may be moved into model namelist)

!  ---  ictm=yyyy#, controls time sensitive external data (e.g. co2, solcon, aerosols, etc)
!     integer, parameter :: ictm =   -2 ! same as ictm=0, but add seasonal cycle from
!                                       ! climatology. no extrapolation.
!     integer, parameter :: ictm =   -1 ! use user provided external data set for the
!                                       ! forecast time. no extrapolation.
!     integer, parameter :: ictm =    0 ! use data at initial cond time, if not
!                                       ! available, use latest, no extrapolation.
!!    integer, parameter :: ictm =    1 ! use data at the forecast time, if not
!                                       ! available, use latest and extrapolation.
!     integer, parameter :: ictm =yyyy0 ! use yyyy data for the forecast time,
!                                       ! no further data extrapolation.
!     integer, parameter :: ictm =yyyy1 ! use yyyy data for the fcst. if needed, do
!                                       ! extrapolation to match the fcst time.

!  ---  isol controls solar constant data source
!!    integer, parameter :: isol = 0   ! use prescribed solar constant
!     integer, parameter :: isol = 1   ! use varying solar const with 11-yr cycle

!  ---  ico2 controls co2 data source for radiation
!     integer, parameter :: ico2 = 0   ! prescribed global mean value (old opernl)
!!    integer, parameter :: ico2 = 1   ! use obs co2 annual mean value only
!     integer, parameter :: ico2 = 2   ! use obs co2 monthly data with 2-d variation

!  ---  ialb controls surface albedo for sw radiation
!!    integer, parameter :: ialb = 0   ! use climatology alb, based on sfc type
!     integer, parameter :: ialb = 1   ! use modis derived alb (to be developed)

!  ---  iems controls surface emissivity and sfc air/ground temp for lw radiation
!        ab: 2-digit control flags. a-for sfc temperature;  b-for emissivity
!!    integer, parameter :: iems = 00  ! same air/ground temp; fixed emis = 1.0
!!    integer, parameter :: iems = 01  ! same air/ground temp; varying veg typ based emis
!!    integer, parameter :: iems = 10  ! diff air/ground temp; fixed emis = 1.0
!!    integer, parameter :: iems = 11  ! diff air/ground temp; varying veg typ based emis

!  ---  iaer  controls aerosols scheme selections
!     integer, parameter :: iaer = abc --> abc are three digits integer numbers
!                                          to control aerosol effect
!     a: stratospheric-volcanic forcing;  b: lw radiation;  c: sw radiation.
!      a=0: no stratospheric-volcanic aerosol effect;   =1: include the effect.
!      b=0: no lw tropospheric aerosols; =1: lw compute 1 bnd; =2: lw compute multi bnd.
!      c=0: no sw tropospheric aerosols; =1: sw compute multi band.

!  ---  iovr controls cloud overlapping method in radiation:
!     integer, parameter :: iovr_sw = 0  ! sw: random overlap clouds
!!    integer, parameter :: iovr_sw = 1  ! sw: max-random overlap clouds

!     integer, parameter :: iovr_lw = 0  ! lw: random overlap clouds
!!    integer, parameter :: iovr_lw = 1  ! lw: max-random overlap clouds

!  ---  isubc controls sub-column cloud approximation in radiation:
!!    integer, parameter :: isubc_sw = 0 ! sw: without sub-col clds approx
!     integer, parameter :: isubc_sw = 1 ! sw: sub-col clds with prescribed seeds
!     integer, parameter :: isubc_sw = 2 ! sw: sub-col clds with random seeds

!!    integer, parameter :: isubc_lw = 0 ! lw: without sub-col clds approx
!     integer, parameter :: isubc_lw = 1 ! lw: sub-col clds with prescribed seeds
!     integer, parameter :: isubc_lw = 2 ! lw: sub-col clds with random seeds

!
!    the following parameters are from gbphys

      real (kind=kind_phys), parameter :: dxmax=-16.118095651,          &
     &                dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)

      integer :: i, j, k, n, nt
      integer :: kdtphi,kdtlam,ks                                ! hmhj
      integer :: dbgu

      logical :: change
!JR Removed "first" because it is not thread-safe

!  ---  for debug test use
!     real (kind=kind_phys) :: temlon, temlat, alon, alat
      integer :: ipt
      logical :: lprnt

!  jbao 2012 update
      CHARACTER(len=12) :: yyyymmddhhmm
      INTEGER year, month, day, hour, minute
!  end jbao
!  ---  timers:
!     real*8 :: rtc, timer1, timer2
!
!===> *** ...  begin here
!
!$$$      integer                lots,lotd,lota
!$$$      parameter            ( lots = 5*levs+1*levh+3 )
!$$$      parameter            ( lotd = 6*levs+2*levh+0 )
!$$$      parameter            ( lota = 3*levs+1*levh+1 )

      integer              ksd,ksplam,kspphi,ksq,ksr,kst
      integer              ksu,ksv,ksz,item,jtem
! jbao add for may 2014
      real(kind=kind_phys) :: sinlat_r(1), coslat_r(1)
!JR Variables that need thread-local storage
      integer :: nstp                         ! Hoisted from radiation_astronomy.f
      real(kind=kind_phys) :: anginc          ! Hoisted from radiation_astronomy.f
      real(kind=kind_phys) :: sindec          ! Hoisted from radiation_astronomy.f
      real(kind=kind_phys) :: cosdec          ! Hoisted from radiation_astronomy.f
      real(kind=kind_phys) :: sollag          ! Hoisted from radiation_astronomy.f
      

!     real(kind=kind_evod) spdlat(levs,lats_dim_r)
!moor real(kind=kind_phys) slk(levs)
!     real(kind=kind_evod) spdmax_node (levs)
!     real(kind=kind_evod) spdmax_nodes(levs,nodes)
!!
!!................................................................
!!  syn(1, 0*levs+0*levh+1, lan)  ze
!!  syn(1, 1*levs+0*levh+1, lan)  di
!!  syn(1, 2*levs+0*levh+1, lan)  te
!!  syn(1, 3*levs+0*levh+1, lan)  rq
!!  syn(1, 3*levs+1*levh+1, lan)  q
!!  syn(1, 3*levs+1*levh+2, lan)  dpdlam
!!  syn(1, 3*levs+1*levh+3, lan)  dpdphi
!!  syn(1, 3*levs+1*levh+4, lan)  uln
!!  syn(1, 4*levs+1*levh+4, lan)  vln
!!................................................................
!!  dyn(1, 0*levs+0*levh+1, lan)  d(t)/d(phi)
!!  dyn(1, 1*levs+0*levh+1, lan)  d(rq)/d(phi)
!!  dyn(1, 1*levs+1*levh+1, lan)  d(t)/d(lam)
!!  dyn(1, 2*levs+1*levh+1, lan)  d(rq)/d(lam)
!!  dyn(1, 2*levs+2*levh+1, lan)  d(u)/d(lam)
!!  dyn(1, 3*levs+2*levh+1, lan)  d(v)/d(lam)
!!  dyn(1, 4*levs+2*levh+1, lan)  d(u)/d(phi)
!!  dyn(1, 5*levs+2*levh+1, lan)  d(v)/d(phi)
!!................................................................
!!  anl(1, 0*levs+0*levh+1, lan)  w     dudt
!!  anl(1, 1*levs+0*levh+1, lan)  x     dvdt
!!  anl(1, 2*levs+0*levh+1, lan)  y     dtdt
!!  anl(1, 3*levs+0*levh+1, lan)  rt    drdt
!!  anl(1, 3*levs+1*levh+1, lan)  z     dqdt
!!................................................................
!!
!$$$      parameter(ksz     =0*levs+0*levh+1,
!$$$     x          ksd     =1*levs+0*levh+1,
!$$$     x          kst     =2*levs+0*levh+1,
!$$$     x          ksr     =3*levs+0*levh+1,
!$$$     x          ksq     =3*levs+1*levh+1,
!$$$     x          ksplam  =3*levs+1*levh+2,
!$$$     x          kspphi  =3*levs+1*levh+3,
!$$$     x          ksu     =3*levs+1*levh+4,
!$$$     x          ksv     =4*levs+1*levh+4)
!!
!$$$      parameter(kdtphi  =0*levs+0*levh+1,
!$$$     x          kdrphi  =1*levs+0*levh+1,
!$$$     x          kdtlam  =1*levs+1*levh+1,
!$$$     x          kdrlam  =2*levs+1*levh+1,
!$$$     x          kdulam  =2*levs+2*levh+1,
!$$$     x          kdvlam  =3*levs+2*levh+1,
!$$$     x          kduphi  =4*levs+2*levh+1,
!$$$     x          kdvphi  =5*levs+2*levh+1)
!!
!$$$      parameter(kau     =0*levs+0*levh+1,
!$$$     x          kav     =1*levs+0*levh+1,
!$$$     x          kat     =2*levs+0*levh+1,
!$$$     x          kar     =3*levs+0*levh+1,
!$$$     x          kap     =3*levs+1*levh+1)
!!
!$$$      integer   p_gz,p_zem,p_dim,p_tem,p_rm,p_qm
!$$$      integer   p_ze,p_di,p_te,p_rq,p_q,p_dlam,p_dphi,p_uln,p_vln
!$$$      integer   p_w,p_x,p_y,p_rt,p_zq
!$$$      parameter(p_gz  = 0*levs+0*levh+1,  !      gze/o(lnte/od,2),
!$$$     x          p_zem = 0*levs+0*levh+2,  !     zeme/o(lnte/od,2,levs),
!$$$     x          p_dim = 1*levs+0*levh+2,  !     dime/o(lnte/od,2,levs),
!$$$     x          p_tem = 2*levs+0*levh+2,  !     teme/o(lnte/od,2,levs),
!$$$     x          p_rm  = 3*levs+0*levh+2,  !      rme/o(lnte/od,2,levh),
!$$$     x          p_qm  = 3*levs+1*levh+2,  !      qme/o(lnte/od,2),
!$$$     x          p_ze  = 3*levs+1*levh+3,  !      zee/o(lnte/od,2,levs),
!$$$     x          p_di  = 4*levs+1*levh+3,  !      die/o(lnte/od,2,levs),
!$$$     x          p_te  = 5*levs+1*levh+3,  !      tee/o(lnte/od,2,levs),
!$$$     x          p_rq  = 6*levs+1*levh+3,  !      rqe/o(lnte/od,2,levh),
!$$$     x          p_q   = 6*levs+2*levh+3,  !       qe/o(lnte/od,2),
!$$$     x          p_dlam= 6*levs+2*levh+4,  !  dpdlame/o(lnte/od,2),
!$$$     x          p_dphi= 6*levs+2*levh+5,  !  dpdphie/o(lnte/od,2),
!$$$     x          p_uln = 6*levs+2*levh+6,  !     ulne/o(lnte/od,2,levs),
!$$$     x          p_vln = 7*levs+2*levh+6,  !     vlne/o(lnte/od,2,levs),
!$$$     x          p_w   = 8*levs+2*levh+6,  !       we/o(lnte/od,2,levs),
!$$$     x          p_x   = 9*levs+2*levh+6,  !       xe/o(lnte/od,2,levs),
!$$$     x          p_y   =10*levs+2*levh+6,  !       ye/o(lnte/od,2,levs),
!$$$     x          p_rt  =11*levs+2*levh+6,  !      rte/o(lnte/od,2,levh),
!$$$     x          p_zq  =11*levs+3*levh+6)  !      zqe/o(lnte/od,2)
!!
! NAG compiler gives "Wrong data type of value for INF" if set to real

      integer            :: iinf
      real               :: inf
      equivalence (iinf,inf)
      data iinf/z'7f800000'/    ! Infinity

!JR Initialize new stack variables to infinity
      nstp = huge(nstp)
      anginc = inf
      sindec = inf
      cosdec = inf
      sollag = inf

! jbao initialize variables as in gloopr 2012
      f_ice = 0.0
      f_rain = 0.0
      r_rime =0.0
      cldcov =0.0

      vvel = 0.0 !JFM and Bao
! jbao 2012
        icsdlw=0
        icsdsw=0
        numrdm=0
        ixseed=0
        ipseed=0
!        stat=0.0
! jbao 2012

        fluxr = 0.0
        fluxr_v = 0.0
! jbao new gfs 2014 physics intialize radiation

      idat = 0
!  get date info from the date string
      READ(UNIT=yyyymmddhhmm(1:4), FMT='(I4)') year
      READ(UNIT=yyyymmddhhmm(5:6), FMT='(I2)') month
      READ(UNIT=yyyymmddhhmm(7:8), FMT='(I2)') day
      READ(UNIT=yyyymmddhhmm(9:10), FMT='(I2)') hour
      READ(UNIT=yyyymmddhhmm(11:12), FMT='(I2)') minute
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

      rinc = 0.
      rinc(2) = phour  ! hli 09/2014
      call w3movdat(rinc, idat, jdat)

!JR turn of entire if(first) stuff because it is not thread-safe
!JR      if (first) then
!JR        if( thermodyn_id.eq.3 ) then
!JR          if (.not. allocated(xcp)) allocate (xcp(ngptc,levr))
!JR          if (.not. allocated(sumq)) allocate (sumq(ngptc,levr))
!JR        endif
!JR
!JR        if( ntcw <= 0 ) then
!JR          if( gen_coord_hybrid .and. vertcoord_id == 3.) then
!JR            if (.not. allocated(gtvx)) allocate (gtvx(ngptc,levs))
!JR            if (.not. allocated(gtvy)) allocate (gtvy(ngptc,levs))
!JR          endif
!JR
!JR          if (.not. allocated(gu))   allocate (gu(ngptc,levs))
!JR          if (.not. allocated(gv1))  allocate (gv1(ngptc,levs))
!JR          if (.not. allocated(gd))   allocate (gd(ngptc,levs))
!         if (.not. allocated(vvel)) allocate (vvel(ngptc,levs))
!JR          if (.not. allocated(gphi)) allocate (gphi(ngptc))
!JR          if (.not. allocated(glam)) allocate (glam(ngptc))
!JR        endif

! jbao 2014 new gfs to add radiation initialization

!JR Moved rad_initialize call to do_physics_one_step because as written, it is not thread-safe
!JR        first = .false.
!JR
!JR      endif         ! end_if_first


!JR         allocate (sinlat_r(1),coslat_r(1))
         sinlat_r(1) = sin(xlat(1,1)) ! jbao
         coslat_r(1) = cos(xlat(1,1)) ! jbao

!  --- ...  update time varying radiation related quantities

      dtsw  = fhswr                  ! fhswr is in sec
      dtlw  = fhlwr                  ! fhlwr is in sec

      if ( phour > 0.0 ) then
        solhr = mod(float(jdat(5)),24.0) ! hour after 00z at current fcst time
      else
        solhr = idate(1)                 ! initial time
      endif

      call radupdate                                                    &
!  ---  inputs:
     &     ( idat, jdat, dtsw, deltim, lsswr, me,                       &
!  ---  outputs:
     &       slag, sdec, cdec, solcon,                                  &
     &       nstp, anginc, sindec, cosdec, sollag )

!  --- ...  generate 2-d random seeds array for sub-grid cloud-radiation

      if ( isubc_lw==2 .or. isubc_sw==2 ) then
        ipseed = mod(nint(100.0*sqrt(phour*3600)), ipsdlim) + 1 + ipsd0

        call random_setseed                                             &
!  ---  inputs:
     &     ( ipseed,                                                    &
!  ---  outputs:
     &       stat                                                       &
     &      )
        call random_index                                               &
!  ---  inputs:
     &     ( ipsdlim,                                                   &
!  ---  outputs:
     &       numrdm, stat                                               &
     &     )

        do k = 1, 2
          do j = 1, lats_node_r
            lat = global_lats_r(ipt_lats_node_r-1+j)

            do i = 1, lonr
              ixseed(i,j,k) = numrdm(i+(lat-1)*lonr+(k-1)*latr)
            enddo
          enddo
        enddo
      endif

!  --- ...  starting latitude loop

      do lan = 1, lats_node_r

        lat = global_lats_r(ipt_lats_node_r-1+lan)
        lons_lat = lonsperlar(lat)

!       print *,' lan=',lan

!$omp parallel do schedule(dynamic,1) private(lon,i,j,k)
!$omp+private(vvel,gu,gv1,gd,gt,gr,gr1,gq,gphi,glam)
!$omp+private(gtv,gtvx,gtvy,sumq,xcp)
!$omp+private(cldcov_v,fluxr_v,f_ice,f_rain,r_rime)
!$omp+private(deltaq,cnvw,cnvc)
!$omp+private(coszen_v,coszdg_v,sinlat_v,coslat_v)
!$omp+private(prslk,prsl,prsik,prsi,flgmin_v,hlw_v,swh_v)
!$omp+private(hlwc_v,swhc_v)
!$omp+private(njeff,n,item,jtem,ks,work1,work2)
!$omp+private(icsdsw,icsdlw,topfsw,sfcfsw,topflw,sfcflw)
!$omp+private(lprnt,ipt,dbgu)

        do lon = 1, lons_lat, ngptc

! jbao newgfs as previous gloopr for fim  set njeff and istrt=1
!          njeff   = min(ngptc,lons_lat-lon+1)
           njeff=1
! jbao newgfs as previous gloopr for fim  set njeff and istrt=1

          dbgu = 300 + lon
!         print *,' ngptc=',ngptc,' dbgu=',dbgu
!         write(dbgu,*)' dbgu=',dbgu,' lon=',lon,' lan=',lan
          lprnt = .false.

!  --- ...  for debug test
!         alon = 139.20
!         alat = 74.71
!         alon = 97.5
!         alat = -6.66
          ipt = 1
!         do i = 1, njeff
!           item = lon + i - 1
!           temlon = xlon(item,lan) * 57.29578
!           if (temlon < 0.0) temlon = temlon + 360.0
!           temlat = xlat(item,lan) * 57.29578
!           lprnt = abs(temlon-364.5) < 0.5 .and. abs(temlat-65.7) < 0.1&
!    &             .and. kdt > 13
!           lprnt = abs(temlon-alon) < 0.1 .and. abs(temlat-alat) < 0.1 &
!    &             .and. kdt > 0
!           if ( lprnt ) then
!             ipt = item
!             print *,' ipt=',ipt,' lon=',lon,' lan=',lan
!             exit
!           endif
!         enddo
!         lprnt = .false.


          do j = 1, njeff
            sinlat_v(j) = sinlat_r(lat)
            coslat_v(j) = coslat_r(lat)
          enddo
          do k=1,levs
            do j=1,njeff
              prslk(j,k)    = fpkap(prsl(j,k)*1000.0)
!              cldcov_v(j,k) = cldcov(k,istrt+j-1,lan)
            enddo
          enddo


          if (n3dfercld >0) then
            do k = 1, levr
              do j = 1, njeff
                jtem = lon+j-1
                f_ice (j,k) = phy_f3d(jtem,k,1,lan)
                f_rain(j,k) = phy_f3d(jtem,k,2,lan)
                r_rime(j,k) = phy_f3d(jtem,k,3,lan)
              enddo
            enddo

          else
            do j=1,njeff
              flgmin_v(j) = 0.0
            enddo
          endif

          do k = 1, levr
            do j = 1, njeff
              if(n3dzhaocld >0 .and. n3dcldpdf >0) then
                jtem = lon+j-1
                deltaq(j,k) = phy_f3d(jtem,k,n3dzhaocld+1,lan)
                cnvw(j,k)   = phy_f3d(jtem,k,n3dzhaocld+2,lan)
                cnvc(j,k)   = phy_f3d(jtem,k,n3dzhaocld+3,lan)
              else
                deltaq(j,k) = 0.0
                cnvw(j,k)   = 0.0
                cnvc(j,k)   = 0.0
              endif
            enddo
          enddo

!  --- ...  assign random seeds for sw and lw radiations

          if ( isubc_lw==2 .or. isubc_sw==2 ) then
            do j = 1, njeff
              icsdsw(j) = ixseed(lon+j-1,lan,1)
              icsdlw(j) = ixseed(lon+j-1,lan,2)
            enddo
          endif

!  --- ...  calling radiation driver

!         lprnt = me .eq. 0 .and. kdt .ge. 120
!         if (lprnt) then
!         if (kdt .gt. 85) then
!         print *,' calling grrad for me=',me,' lan=',lan,' lat=',lat   &
!    &,           ' num_p3d=',num_p3d
!         if (lan == 47) print *,' gt=',gt(1,:)
!         if (kdt > 3) call mpi_quit(5555)
!
!         if (lprnt) print *,' in grrad tsea=',tsea(ipt,lan)
          call grrad                                                    &
!  ---  inputs:
     &     ( prsi,prsl,prslk,gt,gr,gr1,vvel,slmsk(lon,lan),             &
     &       xlon(lon,lan),xlat(lon,lan),tsea(lon,lan),                 &
     &       snwdph(lon,lan),sncovr(lon,lan),snoalb(lon,lan),           &
     &       zorl(lon,lan),hprime(lon,1,lan),                           &
     &       alvsf(lon,lan),alnsf(lon,lan),alvwf(lon,lan),              &
     &       alnwf(lon,lan),facsf(lon,lan),facwf(lon,lan),              &
     &       fice(lon,lan),tisfc(lon,lan),                              &
     &       sinlat_v,coslat_v,solhr,jdat,solcon,                       &
     &       cv(lon,lan),cvt(lon,lan),cvb(lon,lan),                     &
     &       f_ice,f_rain,r_rime,flgmin_v,                              &
     &       icsdsw,icsdlw,ntcw-1,ncld,ntoz-1,ntrac-1,nfxr,             &
     &       dtlw,dtsw,lsswr,lslwr,lssav,                               &
     &       ngptc,njeff,levr, me, lprnt,ipt,kdt,deltaq,sup,cnvw,cnvc,  &
!  ---  outputs:
     &       swh_v,topfsw,sfcfsw,sfalb(lon,lan),coszen_v,coszdg_v,      &
     &       hlw_v,topflw,sfcflw,tsflw(lon,lan),                        &
     &       sfcemis(lon,lan),cldcov_v,                                 &
     &       nstp, anginc, sindec, cosdec, sollag,                      &
!  ---  input/output:
     &       fluxr_v,                                                   &
!    &       fluxr_v,dbgu                                               &
     &       htrlw0=hlwc_v,htrsw0=swhc_v                                &
     &       )
!         print *, 'after grrad lw',
!    &    minval(hlw_v),maxval(hlw_v),minval(hlwc_v),maxval(hlwc_v)
!         print *, 'after grrad sw',
!    &    minval(swh_v),maxval(swh_v),minval(swhc_v),maxval(swhc_v)

          do j = 1, njeff
            coszen(lon+j-1,lan) = coszen_v(j)
            coszdg(lon+j-1,lan) = coszdg_v(j)
          enddo

!  ---  check print
!         if (lprnt) print *,' returned from grrad for me=',me,' lan=', &
!    &      lan,' lat=',lat,' kdt=',kdt,' tsea=',tsea(ipt,lan)
!         print *,' end gloopr hlw=',hlw(lon,:,lan),' lan=',lan
!
!         if (lprnt) print *,' swh_vg=',swh_v(ipt,:)

          do k=1,nfxr
            do j=1,njeff
              fluxr(k,lon+j-1,lan) = fluxr_v(j,k)
            enddo
          enddo

!  --- ...  radiation fluxes and heating rates

          if (lsswr) then
            do i = 1, njeff
              jtem = lon + i - 1
              sfcdsw(jtem,lan) = sfcfsw(i)%dnfxc
              sfcusw(jtem,lan) = sfcfsw(i)%upfxc      ! hli
              topdsw(jtem,lan) = topfsw(i)%dnfxc      ! hli
              topusw(jtem,lan) = topfsw(i)%upfxc      ! hli
              sfcnsw(jtem,lan) = sfcfsw(i)%dnfxc - sfcfsw(i)%upfxc
            enddo
            do k=1,levr
              do j=1,njeff
                jtem = lon + j - 1
                swh(jtem,k,lan) = swh_v(j,k)
                swhc(jtem,k,lan) = swhc_v(j,k)
              enddo
            enddo

         endif    ! end if_lsswr_block
          if (lslwr) then
            do i = 1, njeff
              jtem = lon + i - 1
              sfcdlw(jtem,lan) = sfcflw(i)%dnfxc
              sfculw(jtem,lan) = sfcflw(i)%upfxc      ! hli
              topulw(jtem,lan) = topflw(i)%upfxc      ! hli
            enddo

            do k=1,levr
              do j=1,njeff
                jtem = lon + j - 1
                hlw(jtem,k,lan) = hlw_v(j,k)
                hlwc(jtem,k,lan) = hlwc_v(j,k)
              enddo
            enddo

          endif    ! end if_lslwr_block

!  ---  
!         if (lat == 45 .and. me == 0 .and. lon == 1) then
!           print *,' after grrad hlw_v=',hlw_v(1,:)
!           print *,' after grrad swh_v=',swh_v(1,:)
!         endif
!         if (lprnt) print *,' hlwg=',hlw(lon+ipt-1,:,lan)
!         if (lprnt) print *,' swhg=',swh(lon+ipt-1,:,lan)
!         if (lprnt) print *,' swh_vg=',swh_v(ipt,:)
 
!         print *,' completed grrad for lan=',lan,' istrt=',istrt
        enddo    ! do_lon_loop

      enddo    ! do_lan_loop

      return
      end subroutine gloopr
