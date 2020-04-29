module hycom_variables
use hycom_control,only: nsecs

!SMS$DISTRIBUTE(dh,1) BEGIN

real,allocatable :: ubforc(:)		! forcing of barotropic u velocity
real,allocatable :: vbforc(:)		! forcing of barotropic v velocity

real,allocatable :: psikk (:)		! montg.pot. in bottom layer
real,allocatable :: spvkk (:)		! specific vol. in bottom layer
real,allocatable :: pbot  (:)		! initial bottom pressure
real,allocatable :: gzbot (:)		! bottom geopotential

real,allocatable :: srfx  (:)		! total heatflux, rad+snsbl+latnt; W/m^2
real,allocatable :: ssht  (:)		! sea surface height (m)
real,allocatable :: pmne  (:)		! fresh water flux, pos.down; m/s
real,allocatable :: ustar (:)		! u_star=sqrt(tau/rho); m/s
real,allocatable :: taux  (:)		! eastward wind stress; N/m^2
real,allocatable :: tauy  (:)		! northward wind stress; N/m^2
real,allocatable :: curl  (:)		! wind stress curl
real,allocatable :: prcp  (:)		! precipitation, pos.down; m/s
real,allocatable :: uwnd  (:)		! u surface wind speed; m/s
real,allocatable :: vwnd  (:)		! v surface wind speed; m/s
real,allocatable :: wspd  (:)		! total surface wind speed; m/s
real,allocatable :: airt  (:)		! 2m air temperature
real,allocatable :: vpmx  (:)		! Specific Humidity; kg/kg

real,allocatable :: curl_ave(:)		! archvintvl-avg'd wind stress curl
real,allocatable :: ssht_ave(:)		! archvintvl-avg'd sea surface height
real,allocatable :: srfx_ave(:)		! archvintvl-avg'd heatflux; W/m^2
real,allocatable :: pmne_ave(:)		! archvintvl-avg'd P minus E; m/s
real,allocatable :: qf2d_ave(:)		! archvintvl-avg'd latent heat; W/m^2
real,allocatable :: hf2d_ave(:)		! archvintvl-avg'd sensbl heat; W/m^2
real,allocatable :: airt_ave(:)		! archvintvl-avg'd air temperature
real,allocatable :: vpmx_ave(:)		! archvintvl-avg'd vapor mixing ratio
real,allocatable :: prcp_ave(:)		! archvintvl-avg'd precip
real,allocatable :: swdn_ave(:)		! archvintvl-avg'd shortwave down; W/m^2
real,allocatable :: lwdn_ave(:)		! archvintvl-avg'd longwave down; W/m^2
real,allocatable ::  rad_ave(:)		! archvintvl-avg'd net radiation; W/m^2
real,allocatable :: wspd_ave(:)		! archvintvl-avg'd wind speed; m/s
real,allocatable :: taux_ave(:)		! archvintvl-avg'd eastwd wind stress
real,allocatable :: tauy_ave(:)		! archvintvl-avg'd northwd wind stress
real,allocatable :: hmixl_ave(:)	! archvintvl-avg'd mixed layer depth (m)
real,allocatable :: temice_ave(:)	! archvintvl-avg'd ice temperature (degC)
real,allocatable :: covice_ave(:)	! archvintvl-avg'd ice coverage (%)
real,allocatable :: thkice_ave(:)	! archvintvl-avg'd ice thickness (m)

real,allocatable :: srfx_bcl  (:)	! bclin_frq-averaged heatflux; W/m^2
real,allocatable :: sw2d_bcl  (:)	! bclin_frq-averaged shortwave; W/m^2
real,allocatable :: pmne_bcl  (:)	! bclin_frq-averaged P minus E; m/s
real,allocatable :: ustar_bcl (:)	! bclin_frq-averaged ustar; m/s
real,allocatable :: ustarb_bcl(:)	! bclin_frq-averaged ustarb; m/s
real,allocatable :: rivflo  (:)		! river water load; kg
real,allocatable :: covice  (:)		! sea ice coverage; rel.units
real,allocatable :: thkice  (:)		! sea ice thickness; m
real,allocatable :: thksno  (:)		! snow depth on top of ice; m
real,allocatable :: temice  (:)		! ice surface temperature; C
real,allocatable :: ticeol  (:)		! previous ice surface temperature; C
real,allocatable :: ustarb  (:)		! bottom friction velocity; m/s
real,allocatable :: hmixl   (:)		! mixed layer depth (m)
real,allocatable :: akpar   (:)		! photosynthetically avail.radn.coeff.
real,allocatable :: salnow  (:)		! observed SSS for the present time step
real,allocatable :: airtmn  (:)		! annual min. air temperature
real,allocatable :: srfxcum (:)         ! intgl of srf heatflux; pos down (W/m2)
real,allocatable :: pmnecum (:)         ! integral of P minus E; pos down (m/s)
real,allocatable :: tdpold  (:)		! SST*dp used in tendency calculations
integer,allocatable :: jerlov(:)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE(dh,2) BEGIN

real,allocatable :: utrop    (:,:)	! barotropic u velocity
real,allocatable :: vtrop    (:,:)	! barotropic v velocity
real,allocatable :: ptrop    (:,:)	! barotropic pressure
real,allocatable :: uclin    (:,:,:)	! baroclinic u velocity
real,allocatable :: vclin    (:,:,:)	! baroclinic v velocity
real,allocatable :: dp       (:,:,:)	! layer thickness
real,allocatable :: dpinit   (:,:)	! initl.lyr.thknss for trcr transport
real,allocatable :: dpfinl   (:,:)	! final lyr.thknss for trcr transport
real,allocatable :: pres     (:,:)	! interface pressure
real,allocatable :: geop     (:,:)	! geopotential
real,allocatable :: montg    (:,:)	! Montgomery potential
real,allocatable :: dens     (:,:)	! water density anomaly
real,allocatable :: spcifv   (:,:)	! specific volume anomaly
real,allocatable :: temp     (:,:)	! temperature
real,allocatable :: saln     (:,:)	! salinity
real,allocatable :: passv_tr (:,:,:)	! passive tracers
real,allocatable :: btropfx  (:,:)	! forcing of barotropic mass flux
!real,allocatable :: zgrid(:,:),vcty(:,:),difs(:,:),dift(:,:),ghats(:,:)  ! kpp
real,allocatable :: srforc   (:,:,:)	! core forcings
real,allocatable :: sstndcy  (:,:)	! SST tendency due to various processes

real,allocatable :: temp_ave (:,:)	! archvintvl-avg'd temperature
real,allocatable :: saln_ave (:,:)	! archvintvl-avg'd salinity 
real,allocatable :: dens_ave (:,:)	! archvintvl-avg'd salinity 
real,allocatable :: uvel_ave (:,:)	! archvintvl-avg'd u velocity
real,allocatable :: vvel_ave (:,:)	! archvintvl-avg'd v velocity
real,allocatable ::   dp_ave (:,:)	! archvintvl-avg'd dp
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE(dh,3) BEGIN

real,allocatable :: um_edg   (:,:,:)	! total u vel.on cell edges, mid-time
real,allocatable :: vm_edg   (:,:,:)	! total v vel.on cell edges, mid-time
real,allocatable :: un_edg   (:,:,:)	! total u vel.on cell edges, old
real,allocatable :: vn_edg   (:,:,:)	! total v vel.on cell edges, old
real,allocatable :: dp_edg   (:,:,:)	! lyr thknss on cell edges, mid-time
real,allocatable :: pr_edg   (:,:,:)	! intfc pres on cell edges, mid-time
real,allocatable :: uvsq_edg (:,:,:)	! velocity-squared on cell edges
real,allocatable :: geop_edg (:,:,:)	! geopotential on cell edges (mid-lyr)
real,allocatable :: mont_edg (:,:,:)	! Montgomery potential on cell edges
real,allocatable :: spcv_edg (:,:,:)	! specific volume on cell edges
real,allocatable :: massflx  (:,:,:)	! lateral mass flux on cell edges
real,allocatable :: cumuflx  (:,:,:)	! mass flux time integral
real,allocatable :: mssfx_ave(:,:,:)	! archvintvl-averaged mass flux

!SMS$DISTRIBUTE END

real wm0,wm1,wm2,wm3		& ! interpol.weights, monthly frcg fields
    ,wd0,wd1,wd2,wd3		& ! interpol.weights, daily frcg fields
    ,ws0,ws1,ws2,ws3		  ! interpol.weights, 6-hrly frcg fields

integer  lm0,lm1,lm2,lm3	& ! pointers for monthly frcg time slots
	,ld0,ld1,ld2,ld3	& ! pointers for daily frcg time slots
	,ls0,ls1,ls2,ls3	& ! pointers for 6-hrly frcg time slots
        ,mth0,mth1,mth2,mth3,mcyc0,mcyc1,mcyc2,mcyc3	&
        ,day0,day1,day2,day3,six0,six1,six2,six3

integer :: nstep_ocn
real*8  :: ocnarea,totmass		! global ocean area (m^2) & mass (kg)
real*8  :: massglb0			! total ocean mass (kg) at t=0
real*8  :: saltini,heatini		! total oceanic salt/heat content at t=0
real*8  :: rhoglb0,temglb0,salglb0	! mean density, temp, saln at t=0
real*8  :: watcumu,pmecumu		! cumulative energy/freshwater fluxes
real*8  :: pcptot,radtot		! global precip/radiation integral
real*8  :: pcpcor=0.,radcor=0.		! anti-drift correction for prcp/rad'n
real*8  :: masscor=0.			! anti-drift correction for total mass

! --- chain of icos cells used for mass flux diagnostics:
  integer                 :: lgthset(nsecs,2)
  integer    ,allocatable :: cellset(:,:,:),edgeset(:,:,:)
  character*6,allocatable :: sense(:,:,:)
  real       ,allocatable :: ctrlat(:,:,:),ctrlon(:,:,:),transp(:,:,:,:)
end module hycom_variables
