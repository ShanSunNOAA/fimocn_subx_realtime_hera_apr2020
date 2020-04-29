module hycom_kpp_constants
! --- KPP related constants
  real,   parameter ::  dp00  =   3.0    ! deep    z-level spacing minimum thickness (m)
  logical,parameter ::  locsig= .false.  ! locally-referenced pot. density for stability
  real,   parameter ::  cb    =   3.e-3  ! coefficient of quadratic bottom friction
  real,   parameter ::  tmljmp=   0.2    ! equivalent temperature jump across mixed-layer (degC)
  real,   parameter ::  rigr  =   0.25   ! PWP:     critical gradient richardson number
  real,   parameter ::  ribc  =   0.65   ! PWP:     critical bulk     richardson number
  real,   parameter ::  rinfty=   0.7    ! KPP:     maximum  gradient richardson number (shear inst.)
  real,   parameter ::  ricr  =   0.45   ! KPP:     critical bulk     richardson number
  real,   parameter ::  bldmin=  20.0    ! KPP:     minimum surface boundary layer thickness
  real,   parameter ::  bldmax=1200.0    ! KPP:     maximum surface boundary layer thickness
  real,   parameter ::  cekman=   0.7    ! KPP/KT:  scale factor for Ekman depth
  real,   parameter ::  cmonob=   1.0    ! KPP:     scale factor for Monin-Obukov depth
  logical,parameter ::  bblkpp= .false.  ! KPP:     activate bottom boundary layer
  logical,parameter ::  shinst= .true.   ! KPP:     activate shear instability mixing
  logical,parameter ::  dbdiff= .true.   ! KPP:     activate double diffusion  mixing
  logical,parameter ::  nonloc= .true.   ! KPP:     activate nonlocal b. layer mixing
  logical,parameter ::  latdiw= .true.   ! K-PROF:  activate lat.dep. int.wave mixing
  logical,parameter ::  botdiw= .true.   ! GISS:    activate bot.enhan.int.wav mixing
  logical,parameter ::  difsmo= .false.  ! K-PROF:  activate horiz smooth diff coeffs
  real,   parameter ::  difm0 =  50.e-4  ! KPP:     max viscosity   due to shear instability (m**2/s)
  real,   parameter ::  difs0 =  50.e-4  ! KPP:     max diffusivity due to shear instability (m**2/s)
  real,   parameter ::  difmiw=   3.e-5  ! KPP:     background/inter :: nal wave viscosity       (m**2/s)
  real,   parameter ::  difsiw=   1.e-5  ! KPP:     background/inter :: nal wave diffusivity     (m**2/s)
  real,   parameter ::  dsfmax=  10.e-4  ! KPP:     salt fingering diffusivity factor        (m**2/s)
  real,   parameter ::  rrho0 =   1.9    ! KPP:     salt fingering rp=(alpha*delT)/(beta*delS)
  real,   parameter ::  cs0   =  98.96   ! KPP:     value for nonlocal flux term
  real,   parameter ::  cstar =  10.0    ! KPP:     value for nonlocal flux term
  real,   parameter ::  cv    =   0.0    ! KPP:     buoyancy frequency ratio (0.0 to use a fn. of N)
  real,   parameter ::  c11   =   5.0    ! KPP:     value for turb velocity scale
  real,   parameter ::  hblflg=   1      ! KPP:     b. layer inter :: polation flag (1=lin.,2=quad.)
  integer,parameter ::  jerlv0 =  1      ! KPP:     initial jerlov water type (1 to 5; 0 to use KPAR)

  integer, parameter :: nzehat=890
  integer, parameter :: nustar=192
  real, dimension (0:nzehat+1,0:nustar+1) ::   &
      wmt            & ! momentum velocity scale table
     ,wst              ! scalar   velocity scale table
  real :: vonk   =  0.4
  real :: zmin   = -0.4e-6
  real :: zmax   =  0.0
  real :: umin   =  0.0
  real :: umax   =  0.16
  real :: epsilon=  0.1
  real :: vtc, cg, dp0enh, deltaz, deltau, qdif0,qdifiw

! --- kpp variables
real :: betard(5),               & ! red  extinction coefficient
        betabl(5),               & ! blue extinction coefficient
        redfac(5)                  ! fract. of penetr. red light

integer, parameter :: lookup=762
integer :: irimax(-lookup:lookup),nb

real :: ribtbl(-lookup:lookup)                  &
       ,ridb(  -lookup:lookup)                  &
       ,slq2b( -lookup:lookup,-lookup:lookup)   &
       ,dri                                     &
       ,smb(   -lookup:lookup,-lookup:lookup)   &
       ,shb(   -lookup:lookup,-lookup:lookup)   &
       ,ssb(   -lookup:lookup,-lookup:lookup)   &
       ,back_ra_r(-39:117)                      &
       ,sisamax(  -39:117)                      &
       ,ra_rmax(  -39:117)                      &
       ,c_y_r0(   -39:117)                      &
       ,sm_r1(    -39:117)                      &
       ,sh_r1(    -39:117)                      &
       ,ss_r1(    -39:117)                      &
       ,slq2_r1(  -39:117)                      &
       ,b1,visc_cbu_limit,diff_cbt_limit        &
       ,theta_rcrp,theta_rcrn

end module hycom_kpp_constants

