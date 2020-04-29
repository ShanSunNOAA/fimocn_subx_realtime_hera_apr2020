program patterngenerator_icos
 use kinds, only: r_kind, r_single, r_double
 use shtns, only: shtns_init, spectogrd, nlm, areawts
 use patterngenerator, only: computevarspec, getnoise, rnorm,&
 patterngenerator_advance, computevargrid, set_random_seed, getvarspectrum,&
 patterngenerator_init,patterngenerator_destroy,random_pattern
 use slint, only: slint_init,bilinear_interp 
 implicit none
 integer, parameter :: nlons=384
 integer, parameter :: nlats=192
 integer, parameter :: ntrunc=126
 logical :: logit
 real(r_kind), parameter :: rerth=6.3712e6
 real(r_kind), parameter :: dt=3600.  ! pattern output interval (seconds)
 real(r_kind), parameter :: tau=6.*3600.  ! pattern decorrelation time scale (seconds)
! real(r_kind), parameter :: tau=3.*3600.  ! pattern decorrelation time scale (seconds)
 real(r_kind), parameter :: lengthscale=500.e3 ! pattern decorrelation length scale (meters)
 real(r_kind), parameter :: psistdev=0.8  ! amplitude of pattern
 integer ::ntime=3781         ! number of time steps to write out
 type(random_pattern) :: rpattern
 real(r_kind) pi
 integer n,iseed
 complex(r_kind), allocatable, dimension(:) :: psispec, noise
 real(r_kind), allocatable, dimension(:,:) :: psigrd
 real(r_kind) varspect(0:ntrunc),varspect1(0:ntrunc)

 real(r_kind), allocatable, dimension(:,:) :: gg_ltln, icos_ltln 
 real(r_kind), allocatable, dimension(:) :: icos_pat 
 integer :: nip, i
 character(16) :: header
 character(8) :: sntime

 integer :: ierr
 integer :: glvl                  ! The grid level defined in the Makefile
 integer :: SubdivNum(20)     = 2 ! Subdivision numbers for each recursive
 integer :: nvl                   ! Number of vertical native levels

 namelist /cntlnamelist/      &
   glvl,                      &
   SubdivNum,                 &
   nvl

 if (iargc() == 1) then
   call getarg(1,sntime)
   read(sntime,*) ntime
 endif

 logit=.true.  ! we use a logit transform (set to true) for SPPT
 pi = 4.*atan(1.0)

 open (20,file='FIMnamelist',status='old',action='read',iostat=ierr)
 if (ierr.ne.0) then
   write (*,'(a)') 'readnl: ERROR, could not open FIMnamelist file '
   stop
 endif
 read (20,nml=CNTLnamelist,iostat=ierr)
 if (ierr.eq.0) then   ! successful read
   write (*,CNTLnamelist)
 end if
 if (ierr.lt.0) then   
   write (*,'(a)') 'readnl: ERROR, CNTLnamelist not read, exiting...'
   stop
 endif

 nip = 1
 do i = 1, glvl
   nip = nip * SubdivNum(i)
 enddo
 nip = 10 * nip * nip + 2
print*, 'nip = ',nip

 call shtns_init(nlons,nlats,ntrunc)  ! initialize spectral transform subroutine
 allocate(psispec(nlm))   ! pattern in spectral space
 allocate(psigrd(nlons,nlats))  ! pattern in grid space

 if (ntime == 0) then
   iseed = 100     ! fix seed
   ntime = 500     ! enough patterns for FIM test suit
 else 
   iseed = 0  ! need to change to get different patterns
 endif

 call patterngenerator_init(lengthscale, dt, tau, psistdev, iseed, rpattern,&
                            nlons,nlats,ntrunc,2.*pi)
 print *,'phi = ',rpattern%phi
 print *,'lengthscale = ',rpattern%lengthscale
 print *,'tau = ',rpattern%tau

 call getnoise(rpattern,psispec)
 psispec = rpattern%varspectrum*psispec
 varspect = 0.
 open(7,file='patterns.dat',form='unformatted')
 do n=1,ntime
    call patterngenerator_advance(psispec,rpattern)
    call spectogrd(psispec,psigrd)
    ! logit transform
    if (logit) then
       psigrd = (2./(1.+exp(psigrd)))-1.
    endif
    write(7) psigrd
 enddo
 close(7)

! convert and write patterns to icos grid
  allocate(gg_ltln(nlons*nlats,2))
  open(7,file='gfsltln_384x192.dat',status='old',form='unformatted')
  read(7) gg_ltln(:,1), gg_ltln(:,2)
  close(7)

  open(7,file='icos_grid_info_level.dat',status='old',form='unformatted')
  read(7) header 
  read(7)
!  print*, header
!  read(header,"(6x,I2)") glvl
!  nip = 10*2**(2*glvl)+2
  allocate(icos_ltln(nip,2))
  read(7) icos_ltln(:,1), icos_ltln(:,2)
  close(7)

  allocate(icos_pat(nip))
  call slint_init(gg_ltln,nlons*nlats,icos_ltln,nip)

  open(7,file='patterns.dat',status='old',form='unformatted')
  open(8,file='patterns_icos.dat',form='unformatted')
  do n=1,ntime
    read(7) psigrd
    call bilinear_interp (psigrd, icos_pat)
    write(8) icos_pat
  enddo
  close(7)
  close(8)

end program patterngenerator_icos
