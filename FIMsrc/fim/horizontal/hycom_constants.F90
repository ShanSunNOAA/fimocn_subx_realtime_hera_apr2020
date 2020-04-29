module hycom_constants
use fimnamelist,     only: glvl,kdm
use module_constants,only: grvity

  real,allocatable :: theta(:)		! target densities
  real,allocatable :: salmin(:)		! minimum salinity values
  real,allocatable :: dpmin(:)		! minimum layer thickness

!SMS$DISTRIBUTE (dh,1) BEGIN
  real   ,allocatable :: odepth(:)	! bathymetry
  real   ,allocatable :: meshsz(:)	! scale for horiz.mesh size
  real   ,allocatable :: thkmx(:)	! max. ice thickness (m)
  real   ,allocatable :: tfrez(:)	! freezing temp
  integer,allocatable :: wet(:)		! mask for ocean(>0)/land(=0) points
  integer,allocatable :: dnhill(:)	! lookup table of downhill neighbors
  integer,allocatable :: nuphill(:)	! no.of neighbors draining into ipn
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE (dh,2) BEGIN
  real   ,allocatable :: dpth_edg(:,:)	! water depth on cell edges
  integer,allocatable :: uphill(:,:)	! neighbors draining into ipn
!SMS$DISTRIBUTE END

  real           :: batrop	 	! barotropic time step
  real           :: dtleap		! leap frog time step
  integer        :: ArchvStepOcn	! archive steps for ocean
  integer        :: TotalStepOcn	! total time step

  real,parameter :: rho = 1000.		! nominal water density (kg/m^3)
  real,parameter :: rhoice = 917.	! ice density (kg/m^3)
  real,parameter :: thref = 1./rho	! nominal specific volume
  real,parameter :: airdns = 1.2	! nominal air density (kg/m^3)
  real,parameter :: stdsal = 34.7	! salinity used in virt.saltflx calc.
  real,parameter :: saldif = 22.	! salinity difference water--ice
  real,parameter :: onem  = rho*grvity	! 1 meter in prs units
  real,parameter :: tenm  = 10.*onem	! 10 m in prs units
  real,parameter :: tencm = 0.1*onem	! 10 cm in prs units
  real,parameter :: onecm = 0.01*onem	! 1 cm in prs units
  real,parameter :: onemm = 0.001*onem	! 1 mm in prs units
  real,parameter :: onemu = 1.e-6*onem	! 1 micron in prs units
  real,parameter :: r_onem = 1./onem	! reciprocal of onem
  real,parameter :: sigjmp = 0.001	! small density increment (kg/m^3)
  real,parameter :: dpthmn = 15.  	! minimum mixed layer depth (m)
  real,parameter :: drcoef = 3.e-3	! bottom drag coeff.
  real,parameter :: spcifh = 3990.	! specific heat of sea water (J/kg/deg)
  real,parameter :: evaplh = 2.47e6	! latent heat of evaporation (J/kg)
  real,parameter :: csubp  = 1005.7	! spcfc heat of air, const.pr.(J/kg/deg)
  real,parameter :: assluv = .25	! off-ctr time smoothing wgt for velo.
  real,parameter :: assldp = .25	! off-ctr time smoothing wgt for thknss

! --- state equation coefficients, Jackett et al. 2006, J. of atm & oceanic technology; T:[-2,30],S:[30,38]
  real, parameter ::		&
!
  c1= 5.019935E+00, c2=1.470290E-02, c3=7.904674E-01, c4=-6.861022E-03,	&
  c5=-3.051459E-03, c6=3.099320E-05, c7=3.050299E-05, c8=1.046036E-04,	&
  c9= 5.287567E-06, pref=1.e7	! for sigma_1

! c1= 9.929853E+00, c2=-1.990743E-02, c3=7.813083E-01, c4=-6.343472E-03,&
! c5=-2.883841E-03, c6=2.746757E-05,  c7=2.897260E-05, c8=1.160426E-04,	&
! c9= 5.007679E-06, pref=2.e7	! for sigma_2

! --- sub-coefficients for locally referenced sigma, a fit towards Jackett et al. 2006, J. of atm & oceanic technology
  real, parameter, dimension(9) ::  				& ! T:[-2,30],S:[30,38], P:[0:4000] in dbar
  alphap = (/-6.313707E-03, 5.066326E-02, 7.999489E-01,   &
             -7.403005E-03,-3.225958E-03, 3.471955E-05,   &
              3.210467E-05, 9.259891E-05, 5.580351E-06 /) &
 ,betap  = (/ 5.083339E-03,-3.664164E-05,-9.637964E-06,   &
              5.540090E-07, 1.778390E-07,-3.824842E-09,   &
             -1.636328E-09, 1.227641E-08,-2.989838E-10 /) &
 ,gammap = (/-5.762706E-08, 6.781354E-10, 1.588334E-10,   &
             -1.212103E-11,-3.390138E-12, 9.942362E-14,   &
              3.514305E-14,-2.772681E-13, 6.323460E-15 /)

! --- ocean namelist variables and their defaults:
  real           :: land_spval = -1.e33	! special value given to land points

  real,parameter :: epsil = 1.e-11
  real,parameter :: huuge = 1.e33	! uu avoids mixup with fortran intrinsic
  real,parameter :: slip = -1.		! +1 for free-slip, -1 for no-slip b.c.
  real,parameter :: temmin = -2.	! lowest allowed temperature (C)
  real,parameter :: KelvC = 273.15	! Celsius <=> Kelvin conversion
  real,parameter :: fusion = 334.e3	! latent heat of fusion (J/kg)
  real,parameter :: tmin0=-5., tmax0=33., smin0=1., smax0=42.   ! min/max T & S

! --- CORE forcing related, default is CORE2
      integer,parameter ::  numflds=9
!
      integer,parameter :: ixuwi=1	! index of u wind        in srforc
      integer,parameter :: ixvwi=2	! index of v wind        in srforc
      integer,parameter :: ixtem=3	! index of air temp.     in srforc
      integer,parameter :: ixvap=4	! index of vap.mix.ratio in srforc
      integer,parameter :: ixlng=5	! index of longwv.rad.   in srforc
      integer,parameter :: ixsho=6	! index of shortw.rad.   in srforc
      integer,parameter :: ixpcp=7	! index of precip        in srforc
      integer,parameter :: ixsss=8	! index of mxlyr.salin.  in srforc
      integer,parameter :: ixrff=9	! index of river runoff  in srforc
!!    integer,parameter :: ixsst=10	! index of mxlyr.temp.   in srforc
!!    integer,parameter :: ixice=11	! index of ice coverage  in srforc

  character(len=12) :: flnm_in(numflds)	! partial names of forcing files
  character(len=12) :: varname(numflds)	! variable names in netcdf files
  integer           :: ncid1(numflds)
  integer,parameter :: nsstpro=7	! max. # of sstbud calls in hycom_run
  character(len=16) :: sstproc(nsstpro) = ' (not assigned) '	! processes affecting SST

end module hycom_constants
