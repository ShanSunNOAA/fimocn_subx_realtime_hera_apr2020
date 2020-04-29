module fimnamelist

  use module_initial_chem_namelists
  use units, only: getunit, returnunit

  implicit none

  character(13),private     :: nlfile = './FIMnamelist'       !
  integer,private           :: ierr,i                         !
  integer,private           :: u                              ! logical unit for namelist io
  integer,private,parameter :: number_of_namelists = 18       !

  public

  save

! Variable declarations and initialization

  integer :: nprocs ! Set depending on ComputeTasks

!!! chemwrf    declared in module_initial_chem_namelists

!!! CNTLnamelist

  integer :: glvl                  ! The grid level defined in the Makefile
  integer :: SubdivNum(20)     = 2 ! Subdivision numbers for each recursive refinement(2: bisection, 3: trisection, etc.)
  integer :: nvl                   ! Number of vertical native levels

!!! DIAGnamelist

  integer :: PrintDiagNoise    = 1       ! Hourly increment to print diagnostic gravity wave noise (-1 = none 0 = every step)
  integer :: PrintDiagProgVars = 6       ! Hourly increment to print diagnostic prognosis variables (-1 = none 0 = every step)
  logical :: PrintDiags        = .false. ! Print diagnostic messages?
  integer :: PrintIpnDiag      = -1      ! Global ipn value at which to print diagnostics (-1  = > no print)
  integer :: stencl_frst       = 0       ! lowest layer for stencl diagnostics
  integer :: stencl_last       = 0       ! uppermost layer for stencl diagnostics
  integer :: stencl_step       = 0       ! layer increment in stencl diagnostics

!!! gfsphys

! Set default values to bad values that can be tested later

  integer      :: ictm       = huge(ictm)
  integer      :: isol       = huge(isol)
  integer      :: ico2       = huge(ico2)
  integer      :: iaer       = huge(iaer)
  integer      :: iaer_mdl   = huge(iaer_mdl)
  integer      :: ialb       = huge(ialb)
  integer      :: iems       = huge(iems)
  integer      :: iovr_sw    = huge(iovr_sw)
  integer      :: iovr_lw    = huge(iovr_lw)
  integer      :: isubc_sw   = huge(isubc_sw)
  integer      :: isubc_lw   = huge(isubc_lw)
  integer      :: n2dzhaocld = huge(n2dzhaocld)
  integer      :: n3dzhaocld = huge(n3dzhaocld)
  integer      :: n2dcldpdf  = huge(n2dcldpdf)
  integer      :: n3dcldpdf  = huge(n3dcldpdf)
  integer      :: n3dflxtvd  = huge(n3dflxtvd)
  real(kind=8) :: fhswr      = huge(fhswr)
  real(kind=8) :: fhlwr      = huge(fhlwr)

!!! ISOBARICnamelist

  character(80) :: isobaric_levels_file = 'NO_SUCH_FILE'

!!! LANDnamelist

  real           :: landsmoothfact = 1.00            ! influence factor radius in grid lengths
  character(60)  :: landdatfile    = 'geo_em.d01.nc' ! file with land dat file
  character(120) :: landglvldir    = './'            !  dir of glvl.dat (must set in FIMnamelist)
  integer        :: niland         = 4001            ! number of grid points in the E-W (x or i) dir of the landdatfile
  integer        :: njland         = 2000            ! number of grid points in the N-S (y or j) dir of the landdatfile
  character(120) :: landdatdir     = '/scratch1/portfolios/BMC/wrfruc/smirnova/FIM_geo/static/'

!!! MODELnamelist

! For enkf:
! *.? stands for formats: b - binary, *.s - spectral *.g - NCEP gaussian
! For sfcio:
! ?*. stands for spatial: i - icos, g - gaussian, s - spectral

  logical        :: addtend                = .false.   ! Add physics tendencies using Tanya Smirnova''s method?
  logical        :: dophysics              = .true.    ! If false, physics is never called
  logical        :: digifilt               = .false.   ! Use digital filter?
  real           :: dt_reducer_numerator   = 9.        ! dt= dt * dt_reducer_numerator/dt_reducer_denominator
  real           :: dt_reducer_denominator = 9.        !
  integer        :: dtratio                = 6         ! tracer B / tracer A timestep ratio
  logical        :: enkfio_in              = .false.   ! Start by reading an EnKF analysis file? previouslu EnKFAnl
  logical        :: enkfio_out             = .false.   ! Do special IO for EnKF analysis?  previously EnKFIO
  character(120) :: incr_fname             = 'gincr.b' !  Name of the  increment file
  character(120) :: bckg_fname             = 'ibckg.b' !  Name of the the background file
  logical        :: enkf_diag              = .false.   ! to output standard fimout files for enkf diagnostics
  integer        :: cycle_freq             = 6         ! hours for cycle frequency
  integer        :: ens_member             = 1         ! member number
  logical, parameter :: eqwgt              = .true.    ! Use equal (1/3) hor.intpol.weights?
  integer        :: flux_conserv_schm      = 2         ! Use 1: FCT, 2: flux_limiter 
  integer, parameter :: janjic             = 0         ! 1: partial Janjic; 2: full Janjic
  logical        :: hiOrderFluxComp        = .true.    ! Use high order flux computation scheme
  logical        :: hiOrderUVComp          = .false.   ! Use high order u,v components computation scheme
  real           :: intfc_smooth           = 30.       ! diffusivity (m/s) for intfc smoothing
  integer        :: gz_smooth              = 2         ! # of gz smoothing passes in op_diag
  integer        :: nts                    = -999      ! # of time steps to run: init to bad value
  real           :: ptop                   = 30.       ! pres at top of model domain (Pa)
  real           :: thktop                 = 30.       ! min.thickness (Pa) of uppermost layer
  real           :: PVparam                = 2.        ! pot.vort. surface used for Rossby wave-breaking diagnostics
  logical        :: pure_sig               = .false.   ! Use pure sigma coord?
  real           :: rleigh                 = 2.        ! rayleigh damping time scale (days^-1)
  real           :: rlthresh               = 90.       ! min speed (m/s) for rayleigh damping
  real           :: rltaper                = 0.4       ! tapering coeff for rayleigh damping
  real           :: slak                   = 0.5       ! intfc movement retardation coeff.
  integer        :: miglim                 = 1         ! intfc migration limit (<0: no limit)
  integer        :: remap_optn             = 2         ! vert.advc optn: 1-PCM; 2-PLM; 3-PPM
  logical        :: smoothtend             = .false.   ! Smooth physics tendencies using Tanya Smirnova''s method?
  real           :: tfiltwin               = 5400.     ! length of digital filter window (secs)
  logical        :: UpdateSST              = .false.   ! Update SST field with interpolated monthly SST and Sea Ice Fraction?
  real           :: veldff_bkgnd           = 0.3       ! diffusion velocity (=diffusion/meshsz)
  real           :: veldff_boost           = 3.        ! veldff at model top (linear ramp-up over several layers)
  integer        :: biharm_frst            = 0         ! index range for using biharm.dissip
  integer        :: biharm_last            = 0         ! index range for using biharm.dissip
  integer        :: wts_type               = 3         ! type of digital filter window (1=Lanczos 2=Hamming 3=Dolph)
  character(80)  :: ocean_bcs_ltln         = 'ocean_bcs_ltln.360x180.dat'
  logical        :: sppt                   = .false.
  logical        :: fixedseed              = .false.
  character*256  :: atm_trnsecdir          = ' '       ! directory of atm.massflx transects

!!! OUTPUTnamelist

  integer      :: ArchvIntvl       = 6       ! archive interval in ArchvTimeUnit
  character(2) :: ArchvTimeUnit    = 'hr'    !  ts:timestep; mi:minute; hr:hour dy:day; mo:month
  logical      :: FixedGridOrder   = .true.  ! Always output in IJ order? (false => order determined by curve)
  integer      :: itsStart         = 1       ! index of starting time step (is and always has been 1)
  logical      :: PrintMAXMINtimes = .true.  ! Print MAX MIN routine times? (false => print for each PE)
  logical      :: readrestart      = .false. ! Start by reading a restart file?
  integer      :: restart_freq     = 999999  ! Interval (units=ArchvTimeUnit) to write restart file
  logical      :: TimingBarriers   = .false. ! Insert timed barriers to measure task skew? (slows down model)
  integer      :: TotalTime        = 5       ! total integration time in ArchvTimeUnit

!!! PHYSICSnamelist

  integer :: PhysicsInterval   = 0       ! Interval (seconds) to call physics (0 => every time step)
  integer :: RadiationInterval = 3600    ! Interval in seconds to call radiation (0 => every time step)
  integer :: SSTInterval       = -999    ! Interval in seconds to call sst (0 => every time step)
  integer :: num_p3d           = 4       ! 4 means call Zhao/Carr/Sundqvist Microphysics
  integer :: ishal_cnv         = 1	 ! shallow conv: 0=none; 1=SAS; 2=GF; 3=SAS&GF
  integer :: ideep_cnv         = 1	 ! deep    conv: 0=none; 1=SAS; 2=GF

!!! POSTnamelist

  integer,parameter         :: max_vars                    = 200     ! max number of output vars
  integer,parameter         :: max_pathlen                 = 512     ! max length of pathname
  integer,parameter         :: max_filelen                 = 32      ! max length of filename
  integer,parameter         :: max_varnamelen              = 16      ! max length of variable name
  character(max_pathlen)    :: post_datadir                = ''      ! path to input FIM data
  character(max_pathlen)    :: outputdir                   = ''      ! path to output GRIB or netCDF files
  character(1)              :: input                       = ''      ! unknown functionality not ported from pop
  character(1)              :: output                      = ''      ! unknown functionality not ported from pop
  character(max_varnamelen) :: mode                                  ! Unused namelist placebo (needed by pop)
  character(16)             :: output_fmt                            ! 'grib' or 'nc'
  logical                   :: multiple_output_files                 ! placebo: In pop a flag indicating multiple output files
  character(max_filelen)    :: gribtable                   = ''      ! name of grib table
  integer                   :: grid_id                               ! GRIB-specific var. related to horizontal resolution
  integer                   :: mx                          = -999    ! number of points in x on output lat-lon grid
  integer                   :: my                          = -999    ! number of points in y on output lat-lon grid
  logical                   :: latlonfld                   = .true.  ! flag indicates output grib file is lat-lon
  integer                   :: is                          = -999    ! Interpolation scheme. Only valid value is 1. Init to invalid
  integer                   :: vres                        = -999    ! number of vertical levels (init to bad value)
  character(max_varnamelen) :: var_list(max_vars)          = ' '     ! list of GRIB output variables
  integer                   :: nsmooth_var(max_vars)       = 0       ! number of smoothing invocations
  integer                   :: t1                          = -999    ! used only in pop.F90
  integer                   :: t2                          = -999    ! used only in pop.F90
  integer                   :: delta_t                     = -999    ! used only in pop.F90
  logical                   :: multiHrTu                   = .true.  ! using multi-hour time unit
  logical                   :: fimout                      = .true.  ! Write FIM-style output files directly?
  logical                   :: gribout                     = .false. ! whether to write GRIB files directly
  logical                   :: only_write_var_list_entries = .false. ! FIM binaries only written based on var_list

!!! PREPnamelist

  integer       :: curve               = 3              ! curve type for icosahedral grid points
  integer       :: gtype               = 0              ! grid type for icosahedral grid points
  integer       :: NumPostGrids        = 1              ! number of grids post-processing will create
  integer       :: PostGridIds(10)     = 228            ! post grid ids
  integer       :: NumCacheBlocksPerPE = 1              ! # of cache blocks per processor (only for ij-block curve)
  integer       :: topo_smoo           = 1              ! # of srf.hgt smoothing passes
  logical       :: alt_topo            = .false.        ! Use non-GFS surface height?
  logical       :: alt_land            = .false.        ! if true: use non-GFS surface height
  integer       :: atm_ic              = 2              ! 1 = GFS;         2 = CFSR
  integer       :: ocn_ic              = 2              ! 1 = climatology; 2 = CFSR
  character(80) :: aerosol_file        = 'NO_SUCH_FILE' !
  character(80) :: gfsltln_file        = 'NO_SUCH_FILE' !
  character(80) :: mtnvar_file         = 'NO_SUCH_FILE' !

!!! QUEUEnamelist

  character(8)   :: ComputeTasks = '10'            ! Number of compute tasks for FIM; 'S' means Serial
  character(8)   :: MaxQueueTime = '00:05:00'      ! Run time for the complete job (HH:MM:SS)
  character(160) :: PREPDIR      = 'nodir'         ! If exists, use for prep otherwise calculate prep
  character(160) :: FIMDIR       = 'nodir'         ! If exists, use for FIM otherwise calculate FIM
  character(160) :: DATADIR      = '/no_such_path' ! Location of gfsltln and global_mtnvar files
  character(160) :: DATADR2      = '/no_such_path' ! Location of the sanl file and the sfcanl file
  character(160) :: chem_datadir = '/no_such_path' ! Location of the chemistry data files

!!! SYSTEMnamelist

  character(80) :: mpiruncmd = 'false' ! MPI run command + num-cores switch: 'mpirun -n', 'mpiexec -np', etc.

!!! TIMEnamelist

  character(12) :: yyyymmddhhmm !  Forecast initial time

!!! TOPOnamelist

  character(120) :: topodatfile    = '/no_such_file' !  path to topo dat file
  character(120) :: topoglvldir    = './'            !  dir of glvl.dat
  integer        :: toponpass      =  0              ! # of passes of shuman smoother-desmoother of input 5' wrf topo grid
  real           :: toposmoothfact = 1.25            ! radius of influence factor in grid lengths

!!! wrfphysics declared in module_initial_chem_namelists

!!! WRITETASKnamelist

  integer :: cpn                       = -1      ! cores per node (user MUST specify)
  integer :: mpipn                     = -1      ! MPI tasks per node (default will be cpn)
  integer :: max_write_tasks_per_node  = 7       ! default
  integer :: num_write_tasks           = 0       ! default is no write tasks
  integer :: nthreads                  = 1       ! number of OpenMP threads (default 1)
  logical :: check_omp_consistency     = .false. ! default is do not check consistency of OMP settings
  logical :: abort_on_bad_task_distrib = .true.  ! default is to die if something fishy
  logical :: debugmsg_on               = .false. ! write-task debug message control
  logical :: root_own_node             = .true.  ! default is to put rank 0 on node by himself

!!! OCEANnamelist

  integer      :: kdm=20            ! default vertical ocean grid dimension
  integer      :: itest=-1          ! test pt, may differ from PrintIpnDiag
  integer      :: test_start=-1     ! start time for diagnostics
  integer      :: diag_intvl=730    ! time intvl for diagnostics
  integer      :: num_mthly_fields=720
  integer      :: num_daily_fields=21900
  integer      :: num_6hrly_fields=87600

  logical      :: ann_core=.false.  ! if true, use annual average instead of daily/monthly
  logical      :: do_rnoff=.true.   ! return precip to ocn via river runoff
  logical      :: do_radcor=.false. ! do radiation correction
  logical      :: do_pcpcor=.false. ! do precipitation correction
  logical      :: janjic_ocn=.false.! compute PGF using Janjic scheme
  logical      :: sss_rstor=.false. ! restore to observed monthly SSS
  logical      :: ocnonly=.false., atmonly=.false., coupled=.false.  ! only one can be true
  logical      :: ocntst=.false.    ! for test suite to test the case of ocnonly=T

  character*80 :: inicondir=' '     ! directory of initial conditions
  character*256:: ocn_trnsecdir=' ' ! directory of ocn.massflx transects
  real         :: smagcf = 0.3      ! nonlin.viscosity coeff. (nondim)
  real         :: veldff = 0.1      ! lat.momentum mix.coeff (m/s)
  real         :: temdff = 0.02     ! lat.t/s mixing coeff (m/s)
  real         :: thkdff = 0.05     ! thkns.diffu.coeff (m/s)
  real         :: biharm = 0.75     ! Lapl/biharm bolusflx blend (0=La,1=bi)
  real         :: bolfac = 1.0      ! factor for locally enhancing bolus flux
  real         :: diapyc=2.e-5      ! diapycnal diffusivity (m^2/s)
  real         :: diapyn=2.e-7      ! diapyc dfsivty x buoycy freq (m^2/s^2)
  real         :: ocnmx_factor_s=.5, ocnmx_factor_t=.5    !factor to reduce kpp difs/dift

! Namelist declarations

  namelist /chemwrf/           &
    emi_inname,                &
    fireemi_inname,            &
    emi_outname,               &
    fireemi_outname,           &
    input_chem_inname,         &
    input_chem_outname,        &
    frames_per_emissfile,      &
    frames_per_fireemissfile,  &
    io_style_emissions,        &
    io_form_emissions,         &
    bioemdt,                   &
    photdt,                    &
    chemdt,                    &
    ne_area,                   &
    kemit,                     &
    nmegan,                    &
    kfuture,                   &
    errosion_dim,              &
    chem_conv_tr,              &
    chem_opt,                  &
    gaschem_onoff,             &
    aerchem_onoff,             &
    wetscav_onoff,             &
    cldchem_onoff,             &
    vertmix_onoff,             &
    chem_in_opt,               &
    phot_opt,                  &
    drydep_opt,                &
    emiss_opt,                 &
    dust_opt,                  &
    dmsemis_opt,               &
    seas_opt,                  &
    bio_emiss_opt,             &
    biomass_burn_opt,          &
    plumerisefire_frq,         &
    emiss_inpt_opt,            &
    gas_bc_opt,                &
    gas_ic_opt,                &
    aer_bc_opt,                &
    aer_ic_opt,                &
    have_bcs_chem,             &
    aer_ra_feedback,           &
    aer_op_opt,                &
    ash_height,                &
    ash_mass,                  &
    tr_height,                 &
    tr_mass

  namelist /cntlnamelist/      &
    glvl,                      &
    SubdivNum,                 &
    nvl

  namelist /diagnamelist/      &
    PrintDiagNoise,            &
    PrintDiagProgVars,         &
    PrintDiags,                &
    PrintIpnDiag,              &
    stencl_frst,               &
    stencl_last,               &
    stencl_step

  namelist /gfsphys/           &
    ictm,                      &
    isol,                      &
    ico2,                      &
    iaer,                      &
    iaer_mdl,                  &
    ialb,                      &
    iems,                      &
    iovr_sw,                   &
    iovr_lw,                   &
    isubc_sw,                  &
    isubc_lw,                  &
    n2dzhaocld,                &
    n3dzhaocld,                &
    n2dcldpdf,                 &
    n3dcldpdf,                 &
    n3dflxtvd,                 &
    fhswr,                     &
    fhlwr

  namelist /isobaricnamelist/  &
    isobaric_levels_file

  namelist /landnamelist/      &
    landsmoothfact,            &
    landdatdir,                &
    landdatfile,               &
    landglvldir,               &
    niland,                    &
    njland

  namelist /modelnamelist/     &
    addtend,                   &
    dophysics,                 &
    digifilt,                  &
    dt_reducer_denominator,    &
    dt_reducer_numerator,      &
    dtratio,                   &
    enkfio_in,                 &
    enkfio_out,                &
    incr_fname,                &
    bckg_fname,                &
    enkf_diag,                 &
    cycle_freq,                &
    ens_member,                &
    flux_conserv_schm,         &
    hiOrderFluxComp,           &
    hiOrderUVComp,             &
    intfc_smooth,              &
    gz_smooth,                 &
    nts,                       &
    ocean_bcs_ltln,            &
    ptop,                      &
    PVparam,                   &
    pure_sig,                  &
    rleigh,                    &
    rlthresh,                  &
    rltaper,                   &
    slak,                      &
    miglim,                    &
    remap_optn,                &
    smoothtend,                &
    tfiltwin,                  &
    thktop,                    &
    UpdateSST,                 &
    veldff_bkgnd,              &
    veldff_boost,              &
    biharm_frst,               &
    biharm_last,               &
    sppt,                      &
    fixedseed,                 &
    wts_type,		       &
    atm_trnsecdir

  namelist /outputnamelist/    &
    ArchvIntvl,                &
    ArchvTimeUnit,             &
    FixedGridOrder,            &
    itsStart,                  &
    PrintMAXMINtimes,          &
    readrestart,               &
    restart_freq,              &
    TimingBarriers,            &
    TotalTime

  namelist /physicsnamelist/   &
    PhysicsInterval,           &
    RadiationInterval,         &
    SSTInterval,               &
    num_p3d,                   &
    ishal_cnv,                 &
    ideep_cnv

  namelist /postnamelist/      &
    post_datadir,              &
    outputdir,                 &
    input,                     &
    output,                    &
    output_fmt,                &
    multiple_output_files,     &
    gribtable,                 &
    grid_id,                   &
    mx,                        &
    my,                        &
    latlonfld,                 &
    is,                        &
    vres,                      &
    mode,                      &
    var_list,                  &
    nsmooth_var,               &
    t1,                        &
    t2,                        &
    delta_t,                   &
    multiHrTu,                 &
    gribout,                   &
    fimout,                    &
    only_write_var_list_entries

  namelist /prepnamelist/      &
    aerosol_file,              &
    topo_smoo,                 &
    alt_topo,                  &
    alt_land,                  &
    atm_ic, ocn_ic,            &
    curve,                     &
    gtype,                     &
    gfsltln_file,              &
    mtnvar_file,               &
    NumCacheBLocksPerPE,       &
    NumPostGrids,              &
    PostGridIds

  namelist /queuenamelist/     &
    ComputeTasks,              &
    MaxQueueTime,              &
    PREPDIR,                   &
    FIMDIR,                    &
    DATADIR,                   &
    DATADR2,                   &
    chem_datadir

  namelist /systemnamelist/    &
    mpiruncmd

  namelist /timenamelist/      &
    yyyymmddhhmm

  namelist /toponamelist/      &
    topodatfile,               &
    topoglvldir,               &
    toponpass,                 &
    toposmoothfact

  namelist /wrfphysics/        &
    mp_physics,                &
    gsfcgce_hail,              &
    gsfcgce_2ice,              &
    progn,                     &
    ra_lw_physics,             &
    ra_sw_physics,             &
    naer,                      &
    sf_sfclay_physics,         &
    sf_surface_physics,        &
    bl_pbl_physics,            &
    sf_urban_physics,          &
    cu_physics,                &
    num_urban_layers,          &
    cugd_avedx,                &
    imomentum,                 &
    radt,                      &
    clos_choice,               &
    num_land_cat,              &
    num_soil_cat,              &
    mp_zero_out,               &
    mp_zero_out_thresh,        &
    seaice_threshold,          &
    cu_rad_feedback,           &
    slope_rad,                 &
    topo_shading

  namelist /writetasknamelist/ &
    check_omp_consistency,     &
    abort_on_bad_task_distrib, &
    cpn,                       &
    mpipn,                     &
    nthreads,                  &
    debugmsg_on,               &
    max_write_tasks_per_node,  &
    num_write_tasks,           &
    root_own_node

  namelist /oceannamelist/ &
    ocn_trnsecdir,         &
    kdm,                   &
    smagcf,                &
    veldff,                &
    temdff,                &
    thkdff,                &
    biharm,                &
    bolfac,                &
    diapyc,                &
    diapyn,                &
    ocnmx_factor_s,	   &
    ocnmx_factor_t,	   &
    num_mthly_fields,      &
    num_daily_fields,      &
    num_6hrly_fields,      &
    diag_intvl,            &
    itest,                 &
    test_start,            &
    sss_rstor,             &
    inicondir,             &
    ann_core,              &
    ocnonly,               &
    atmonly,               &
    coupled,               &
    ocntst,                &
    do_rnoff,              &
    do_radcor,             &
    do_pcpcor,             &
    janjic_ocn

contains

! Routines for reading namelists in parallel and serial contexts

  subroutine readnl
    logical,save::first_call=.true.
#include "fimnamelist.exec"
  end subroutine readnl

  subroutine readnl_serial
    logical,save::first_call=.true.
!sms$ignore begin
#include "fimnamelist.exec"
!sms$ignore end
  end subroutine readnl_serial

end module fimnamelist
