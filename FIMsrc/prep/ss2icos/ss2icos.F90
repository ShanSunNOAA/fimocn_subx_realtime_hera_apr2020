!------------------------------------------------------------------------------
! This file contains the subroutines that transform GFS analysis in spherical
! spectral coefficients to physical domain and interpolate them from Gaussian
! grid to an specified icosahedral grid.
! 
! History:
!   2006 - 2013: Jin Lee, Jacques Middlecoeff, and Ning Wang adapted spectral
!   transform code from NCEP subroutine ss2gg(), and developed a horizontal
!   interpolation package. Rainer Bleck developed vertical interpolation code 
!   to transform initial condition from sigma to hybrid-isentropic coordinate. 
!
!   2009: Rainer Bleck added solid rotation subroutine as a part of test code. 
!
!   2012: Jim Rosinski modified subroutines to added openmp directives to 
!   the code to significantly reduce the computation time.
!   
!   2014: Ning Wang rewrite and re-arrange the computation to reduce memory 
!   usage, in preparation for the use of high resolution initial condition.
!
!   2014: Jim Rosinski revised structure of OpenMP directives to retain ability
!   to parallelize the transform and interpolation procedure on a single node.
!   Also major cleanup of argument lists and code structure. Removed obsolete
!   features. Made ss2icos a module. All data and contained procedures are 
!   private, except for subroutine ss2icos.
!
!   Parallel algorithm: A single "chunk" of k-values is assigned to a set of
!   available worker threads (e.g. 12 k-values are processed by 12 threads on
!   zeus Westmere processors). When there are more k-values than threads,
!   multiple iterations are required to complete the transform procedure along
!   with interpolation to the icos grid.
!
!   2016: Rainer Bleck rewrote part of subr. ss2ggtopo to emphasize its
!   utility as a tool for reconciling a pressure/height column profile with
!   an independently prescribed surface height.
!------------------------------------------------------------------------------
module ss_gfs

  use headers,only: testcurveheader,testglvlheader
  use module_variables, only: us3d, vs3d, dp3d, mp3d, pr3d, ex3d, ph3d, tr3d
  use sigio_module
  use sptsnv, only: init_sp, end_sp, sptez_mod, sptezmv_mod
  use module_control,only: nvlp1,nip,ntra,ntrb
  use fimnamelist,   only: glvl,nvl,curve,alt_topo,topo_smoo,pure_sig,	&
                           yyyymmddhhmm, gfsltln_file, PrintIpnDiag,	&
                           nthreads, root_own_node, cpn, atm_ic
  use module_constants, only: rd, cp, qvmin, grvity, sigak, sigbk,	&
      raddeg, p1000, perm
  USE slint, ONLY: bilinear_init, bilinear_interp, bilinear_interp_uv, bilinear_interp_uv2
  use mdul_smoofld, only: smolight,smoheavy
  use stencilprint
  use findmaxmin2
  use infnan, only: negint
  use units, only: getunit, returnunit
  use global_bounds, only: myrank

  implicit none

  private
  public :: ss2icos

!sms$distribute(dh,1) begin
  real(4), allocatable :: hs_lev(:)               ! surface height (m)
  real(4), allocatable :: ps_lev(:)               ! surface pressure in pascals
  real(4), allocatable :: test(:)                 ! smoother test field
!sms$distribute end

! The value of nvp is not yet known, so these arrays must be allocatable

!sms$distribute(dh,2) begin
  real(4), allocatable :: t_lyr (:,:)         ! temperature in Kelvins
  real(4), allocatable :: qv_lyr(:,:)         ! specific humidity
  real(4), allocatable :: qc_lyr(:,:)         ! cloud condensate
  real(4), allocatable :: u_lyr (:,:)         ! zonal velocity
  real(4), allocatable :: v_lyr (:,:)         ! meridional velocity
  real(4), allocatable :: o3_lyr(:,:)         ! ozone mixing ratio
  real(4), allocatable :: p_lev (:,:)         ! interface pressure
  real(4), allocatable :: gz_lev(:,:)         ! geopotential height
!sms$distribute end

  real(4), allocatable :: ak(:)
  real(4), allocatable :: bk(:)
  
! Storage for GFS input on Gaussian grid -- may be used multiply as indicated by '/'
  real(4), allocatable :: f1(:,:)           ! surface pressure on gaussian grid
  real(4), allocatable :: f2(:,:)           ! surface height (m) on gaussian grid

! Chunking info
  integer, allocatable :: kbeg(:)   ! beginning k-index per chunk
  integer, allocatable :: kend(:)   ! ending k-index per chunk
  integer :: nchunks = negint       ! number of chunks (each of size at most nthreads)
  integer :: maxchunk = negint      ! max size of any chunk

! Gaussian grid
  integer :: imax = negint          ! size of lon dim on gaussian grid: init to bad value
  integer :: jmax = negint          ! size of lat dim on gaussian grid: init to bad value
  integer :: nvp = negint           ! size of vertical dim on gaussian grid: init to bad value
  integer :: maxwv = negint         ! max number of spectral waves for given truncation

  type(sigio_head), save :: head
  type(sigio_data), save :: data

  integer, parameter :: idrt = 4    ! needed by spectral transform code

#include <gptl.inc>

contains

  subroutine ss2icos ()
! read spherical data (GFS) and perform 2-step transform:
! --- (1) horizontal transform from spherical to icos grid
! --- (2) vertical transform from sigma to hybrid-isentropic coord.

! Local workspace

    character(len=14), parameter :: thisfunc = 'gsi2gauss2icos'
    character :: string*8
    integer :: ipn, k
    logical :: vrbos = .false.
    integer :: luglvl = negint           ! unit number for glvl.dat file

    integer(sigio_intkind) :: iret
    integer :: ret                    ! GPTL return code
#ifdef _OPENMP
    integer, external :: omp_set_num_threads
#endif

    ret = gptlstart (thisfunc)
! In OMP mode, if root_own_node is true then use all cores on the node. 
! No need to reset at the end of this routine either: might as well let 
! the root use all available cores.
! All this currently commented out because "omplace" on zeus has already
! assumed a different thread count and does the wrong thing.
!
! Note cpn wil be 0 in a serial run, thus the test on it below.
! TODO? Move this code into fimcore.F90? Problem is, myrank is not yet
! known there.
! TODO: Enable OpenMP in a serial run
#ifdef _OPENMP
!    if (root_own_node .and. myrank == 0 .and. cpn > 0) then
!      nthreads = cpn
!      ret = omp_set_num_threads (nthreads)
!    end if
#endif
! Slaves need to know nvp,imax,jmax for dimensioning
!SMS$SERIAL (<nvp,imax,jmax,maxwv,OUT> : default=ignore)  BEGIN
    ret = gptlprint_memusage ('start '//thisfunc)
    call initialize_ss ()
    nvp = head%levs			! # of layers in input grid
    maxwv = head%jcap
    ! set imax/jmax values based on jcap in sanl file since there were
    !   inconsistent values
    IF (maxwv == 254) THEN
      imax = 768 
      jmax = 384 
    ELSE IF (maxwv == 382) THEN
      imax = 1152 
      jmax = 576 
    ELSE IF (maxwv == 574) THEN
      imax = 1760 
      jmax = 880 
    ELSE IF (maxwv == 878) THEN
      imax = 2640 
      jmax = 1320 
    ELSE IF (maxwv == 1534) THEN
      imax = 3072 
      jmax = 1536 
    ENDIF
    write (*,'(a)'   ) 'Input spectral grid is as follows:'
    write (*,'(a,i0)') 'nvp (vertical levels)   = ',nvp
    write (*,'(a,i0)') 'imax (num longitudes)   = ',imax
    write (*,'(a,i0)') 'jmax (num latitudes)    = ',jmax
    write (*,'(a,i0)') 'maxwv (max wavenumbers) = ',maxwv
    ret = gptlstart ('init_sp')
    call init_sp (0, maxwv, 4, imax, jmax) ! 4 means gaussian grid
    ret = gptlstop ('init_sp')
!SMS$SERIAL END

    allocate (test(nip))
    allocate (hs_lev(nip))
    allocate (ps_lev(nip))
    allocate (t_lyr (nvp  ,nip))
    allocate (qv_lyr(nvp  ,nip))
    allocate (qc_lyr(nvp  ,nip))
    allocate (u_lyr (nvp  ,nip))
    allocate (v_lyr (nvp  ,nip))
    allocate (o3_lyr(nvp  ,nip))
    allocate (p_lev (nvp+1,nip))
    allocate (gz_lev(nvp+1,nip))
    allocate (ak(nvp+1))
    allocate (bk(nvp+1))

! Slaves need to know nchunks, maxchunk for dimensioning
!SMS$SERIAL (<nchunks,maxchunk,OUT> : default=ignore)  BEGIN
    allocate(f1(imax,jmax))
    allocate(f2(imax,jmax))

    ret = gptlstart ('bilinear_init')
    luglvl = getunit ()
    open (unit=luglvl, file="glvl.dat", status='old', form='unformatted')
    call TestGlvlHeader (luglvl, "glvl.dat", thisfunc, glvl)
    call TestCurveHeader(luglvl, "glvl.dat", thisfunc, curve)
    CALL bilinear_init(gfsltln_file, imax*jmax, luglvl, nip)
    close (luglvl)
    call returnunit (luglvl)
    ret = gptlstop ('bilinear_init')

!JR Recent PPP requires any variables with specific IN,OUT, or INGORE specs to be explicitly referenced
!JR in SERIAL regions, so pass them in to assign_work_to_threads even though they're module-local
    call assign_work_to_threads (nchunks, maxchunk)
!SMS$SERIAL END


!SMS$SERIAL (<hs_lev,ps_lev,OUT> : default=ignore)  BEGIN
    call transform_interp_1 (hs_lev, ps_lev)
!SMS$SERIAL end


!SMS$SERIAL (<qv_lyr,OUT> : default=ignore) BEGIN
    call transform_interp_2 (qv_lyr)
!SMS$SERIAL end

! avoid zero pressure at top, and preserve for now the ak, bk bug found by Ning
!SMS$SERIAL (default=ignore) BEGIN
    ak(:nvp) = head%ak(:nvp)
    ak(nvp+1) = max (head%ak(nvp+1), .2*head%ak(nvp))

    bk(:nvp) = head%bk(:nvp)
    bk(nvp+1) = 0.
!JR Need to explicitly restrict the range because size of head%ak is something like 101 (i.e. padded)
    head%ak(:nvp+1) = ak(:nvp+1)
    head%bk(:nvp+1) = bk(:nvp+1)
!SMS$SERIAL END

! These arrays are used elsewhere in FIM, so they will not be deallocated, and need an OUT clause
! on the SERIAL directive
    allocate (sigak(nvp+1))
    allocate (sigbk(nvp+1))

!SMS$SERIAL (<p_lev,sigak,sigbk,OUT> : default=ignore) BEGIN
    call transform_interp_3 (p_lev, sigak, sigbk)
!SMS$SERIAL END


!SMS$SERIAL (<gz_lev,t_lyr,OUT> : default=ignore)  BEGIN
    call transform_interp_4 (gz_lev, t_lyr)
!SMS$SERIAL end

!SMS$SERIAL (<u_lyr,v_lyr,OUT> : default=ignore)  BEGIN
    call transform_interp_5 (u_lyr, v_lyr)
!SMS$SERIAL END

    do k=1,nvp,7
      write (string,'(a,i3)') 'k=',k
      call findmxmn2 (u_lyr, nvp, nip, k, thisfunc//' u_lyr '//string)
      call findmxmn2 (v_lyr, nvp, nip, k, thisfunc//' v_lyr '//string)
    end do

    call stencl (u_lyr, nvp, 1., thisfunc//' -u- wind')
    call stencl (v_lyr, nvp, 1., thisfunc//' -v- wind')

!SMS$SERIAL (<o3_lyr,qc_lyr,OUT> : default=ignore)  BEGIN
    call transform_interp_6 (o3_lyr, qc_lyr)
!SMS$SERIAL END

!SMS$SERIAL (default=ignore)  BEGIN
    call sigio_axdata (data, iret)     ! deallocate array
    call end_sp ()                     ! deallocate arrays from init data

! also done with these
    deallocate (f1)
    deallocate (f2)
    deallocate (kbeg)
    deallocate (kend)
!SMS$SERIAL END
    ! horizontal interpolation done.

    if (alt_topo) then
     
!SMS$SERIAL (<hs_lev,INOUT> : default=ignore)  BEGIN
! get topo on icos grid
      ret = gptlstart ('mktopo')
      call mktopo (hs_lev, nip)
      ret = gptlstop ('mktopo')
!SMS$SERIAL END

    end if

! --- smooth topography
!   if (topo_smoo > 0) call smoheavy(hs_lev,1,topo_smoo,'srf.height')	! heavy
    if (topo_smoo > 0) call smolight(hs_lev,1,topo_smoo,'srf.height')	! light

! --- make surface pressure consistent with icos-grid surface height
    ret = gptlstart ('ss2ggtopo')
    call ss2ggtopo ()
    ret = gptlstop ('ss2ggtopo')

    call stencl (hs_lev, 1, 1., thisfunc//' surface height (m)')

!SMS$PARALLEL (dh,ipn) BEGIN
    do ipn=1,nip
      hs_lev(ipn) = hs_lev(ipn)*grvity	! surface height => surface geopot
      test(ipn)=0.
      if (ipn.eq.PrintIpnDiag) test(ipn)=1000.			! single-pt mtn
    end do
!SMS$PARALLEL END

!   call smoheavy(test,1,topo_smoo,'single-pt mtn test')		! heavy
    call smolight(test,1,topo_smoo,'single-pt mtn test')		! light
    call stencl(test,1,1.,'single-pt mtn test')

! Now do vertical interpolation.
    ret = gptlstart ('fimini')
    call fimini (nvp,hs_lev,ps_lev,gz_lev,p_lev,t_lyr,qv_lyr,	&
                 u_lyr,v_lyr,o3_lyr,qc_lyr,us3d,vs3d,		&
                 dp3d,mp3d,pr3d,ex3d,ph3d,tr3d)
    ret = gptlstop ('fimini')

    deallocate (hs_lev)
    deallocate (ps_lev)
    deallocate (t_lyr )
    deallocate (qv_lyr)
    deallocate (qc_lyr)
    deallocate (u_lyr )
    deallocate (v_lyr )
    deallocate (o3_lyr)
    deallocate (p_lev )
    deallocate (gz_lev)
    deallocate (ak)
    deallocate (bk)

    ret = gptlstop (thisfunc)
  end subroutine ss2icos

  subroutine initialize_ss ()
    character(len=14), parameter :: thisfunc = 'initialize_ss'
    CHARACTER(len=2 ) :: hh
    CHARACTER(len=9 ) :: jdate
    CHARACTER(len=80) :: sanlFile
    integer(sigio_intkind) :: lusig = negint            ! unit number for NCEP sigio interface
    integer :: ret
    
    if (atm_ic == 1) then                     ! GFS
      call GetJdate(yyyymmddhhmm,jdate)       ! Julian date conversion
      hh = yyyymmddhhmm(9:10)
      sanlFile = jdate // ".gfs.t" // hh // "z.sanl"
    elseif (atm_ic == 2) then                 ! CFSR
      sanlFile = "siganl.gdas." // yyyymmddhhmm(1:10)
    end if
    PRINT *,' read vertical levels from input data ',sanlFile

! sigio reads big-endian data so need to request a specific unit number known by surrounding 
! scriptery to be big-endian
    lusig = getunit (82)
    if (lusig < 0) then
      print*,thisfunc//': getunit failed for unit=82. stopping'
      stop
    end if

    ret = gptlstart ('sigio_srohdc')
! Read in all 3 dimensions of all spectral fields. At T1534 this is a HUGE amount of data
    call sigio_srohdc (lusig, sanlfile, head, data, ret)
    IF (ret /= 0) THEN
      PRINT '(a)',thisfunc//': error reading '//sanlFile
      STOP
!TBH:  this code hangs because errexit does not call MPI_ABORT
!    call errmsg('dyn_init: error reading '//sanlFile)
!    call errexit(2)
    END IF
    CALL returnunit (lusig)

    if (pure_sig .and. nvl /= head%levs) then
      call errmsg (thisfunc//': in "pure_sig" mode, nvl must match head%levs')
      print '(a,2i5)','nvl,head%levs =',nvl,head%levs
      STOP 
    end if
    if (head%idvc < 2) then
      call errmsg (thisfunc//': head%idvc < 2 is NOT SUPPORTED')
      call errmsg ('It is documented in sigio_module.F90 as OBSCOLESCENT')
      call errmsg ('The last FIM rev. which enabled idvc < 2 was r4356')
      call errmsg ('But probably idvc < 2 had bugs in it then')
      STOP
    end if
    ret = gptlstop ('sigio_srohdc')
    ret = gptlprint_memusage ('end '//thisfunc)
  end subroutine initialize_ss

  subroutine assign_work_to_threads (nchunks, maxchunk)
    integer, intent(out) :: nchunks
    integer, intent(out) :: maxchunk

    character(len=22) :: thisfunc = 'assign_work_to_threads'
    integer :: n

    if (nthreads < 1) then
      call errmsg (thisfunc//': do not yet know nthreads')
      stop
    end if

! Discover how to assign work to threads
    if (mod (nvp, nthreads) == 0) then
      nchunks = nvp / nthreads
    else
      nchunks = nvp / nthreads + 1
    end if
    maxchunk = min (nvp, nthreads)
    allocate (kbeg(nchunks))
    allocate (kend(nchunks))

    kbeg(1) = 1
    do n=1,nchunks-1
      kend(n  ) = kbeg(n) + nthreads - 1
      kbeg(n+1) = kend(n) + 1
    end do
    kend(nchunks) = nvp
   
! Print work distribution across threads
    write(6,*)thisfunc//': nchunks=',nchunks,' maxchunk=',maxchunk,' nthreads=',nthreads
    do n=1,nchunks
      write(6,*)thisfunc//': chunk=',n,': kbeg=',kbeg(n),' kend=',kend(n)
    end do
  end subroutine assign_work_to_threads

  subroutine transform_interp_1 (hs_lev, ps_lev)
    real(4), intent(inout) :: hs_lev(nip)
    real(4), intent(inout) :: ps_lev(nip)

    character(len=18), parameter :: thisfunc = 'transform_interp_1'
    integer :: i, j
    integer :: ret

    ret = gptlstart (thisfunc)

! Initialize bilinear interpolation
! perform spherical transform on surface height: hs->f2, and surface pressure ps->f1
! Then interpolate to icos grid.
! Cannot yet do these 2 in parallel with OpenMP because the 1st call to sptez does the initialization
    call sptez_mod (0, maxwv, idrt, imax, jmax, data%hs, f2, 1)
    CALL bilinear_interp (f2, hs_lev)

    call sptez_mod (0, maxwv, idrt, imax, jmax, data%ps, f1, 1)
!$OMP PARALLEL DO PRIVATE (i)
    do j=1,jmax
      do i=1,imax
        f1(i,j) = exp(f1(i,j))*1.e3		! convert ln(ps) in centibars to ps in Pa.
      end do
    end do
!$OMP END PARALLEL DO
    CALL bilinear_interp (f1, ps_lev)	! unit in pascal

    print 100,'min,max of srf.height on spherical grid:', minval(f2),maxval(f2)
    print 100,'min,max of srf.press  on spherical grid:', minval(f1),maxval(f1)
    print 100,'min,max of srf.height on icos grid:',      minval(hs_lev),maxval(hs_lev)
    print 100,'min,max of srf.press  on icos grid:',      minval(ps_lev),maxval(ps_lev)
100 format (a,2f13.2)
    ret = gptlstop (thisfunc)
  end subroutine transform_interp_1

  subroutine transform_interp_2 (qv_lyr)
    real(4), intent(inout) :: qv_lyr(nvp,nip)

    character(len=18), parameter :: thisfunc = 'transform_interp_2'
    integer :: k, ipn
    integer :: ret
    real(4) :: gauss2d(imax,jmax)
    real(4) :: icos2d(nip)                    ! work array on icos grid

    ret = gptlstart (thisfunc)
!$OMP PARALLEL DO PRIVATE(ipn,icos2d,gauss2d)
    do k=1,nvp
      call sptez_mod (0, maxwv, idrt, imax, jmax, data%q(:,k,1), gauss2d, 1)
      CALL bilinear_interp (gauss2d, icos2d)
      do ipn=1,nip
        qv_lyr(k,ipn) = max( icos2d(ipn), qvmin )
      end do
    end do
!$OMP END PARALLEL DO
    ret = gptlstop (thisfunc)
  end subroutine transform_interp_2

  subroutine transform_interp_3 (p_lev, sigak, sigbk)
    real(4), intent(out) :: p_lev(nvp+1,nip)
    real(4), intent(out) :: sigak(nvp+1)
    real(4), intent(out) :: sigbk(nvp+1)

    character(len=18), parameter :: thisfunc = 'transform_interp_3'
    integer :: k, j, i, ipn
    integer :: ret
    real(4) :: gauss2d(imax,jmax)
    real(4) :: icos2d(nip)                    ! work array on icos grid

    ret = gptlstart (thisfunc)
    write (*,'(/a)') 'GFS intfc.prs is defined as p(k) = ak(k) + bk(k) * surf.prs'
    sigak(:) = ak(:)
    sigbk(:) = bk(:)
!$OMP PARALLEL DO PRIVATE(ipn,j,i,icos2d,gauss2d)
    do k=1,nvp+1
      do j=1,jmax
        do i=1,imax
          gauss2d(i,j) = ak(k) + bk(k)*f1(i,j)
        end do
      end do
      CALL bilinear_interp (gauss2d, icos2d)
      do ipn=1,nip
        p_lev(k,ipn) = icos2d(ipn)
      end do
    end do
!$OMP END PARALLEL DO

    ret = gptlstop (thisfunc)
    ret = gptlprint_memusage ('end '//thisfunc)
  end subroutine transform_interp_3

  subroutine transform_interp_4 (gz_lev, t_lyr)
    real(4), intent(out) :: gz_lev(nvp+1,nip)
    real(4), intent(out) :: t_lyr(nvp,nip)

! Local workspace
    character(len=18), parameter :: thisfunc = 'transform_interp_4'
    real(4) :: gauss2d_lo
    real(4) :: gauss2d_up
    real(4) :: pl2d(imax,jmax,maxchunk)
    real(4) :: vtmp2d(imax,jmax,maxchunk)
    integer j,i,k,ipn,k3d,n
    integer :: ret
    real :: pid, piu
    real,parameter :: rocp  = 287.05/1004.6
    real,parameter :: rocp1 = rocp+1
    real,parameter :: rocpr = 1/rocp
    real(4) :: gauss2d(imax,jmax)
    real(4) :: icos2d(nip)                    ! work array on icos grid

    ret = gptlstart (thisfunc)
! perform spherical transform on virt.temp. (t) and specif.hum. (q) field

! Compute geopotential on interfaces (still on Gaussian grid!)

!$OMP PARALLEL DO PRIVATE (i)
    do j=1,jmax
      do i=1,imax
        gauss2d(i,j) = f2(i,j)*grvity  ! f2 = surface height
      end do
    end do
!$OMP END PARALLEL DO

    CALL bilinear_interp (gauss2d, icos2d)
    do ipn=1,nip
      gz_lev(1,ipn) = icos2d(ipn)
    end do

!JR Loop over chunks instead of nvp to save memory. Large arrays need only be dimensioned (nthreads) instead
!JR of (nvp). If memory becomes a problem, cut down the number of threads at the expense of speed.
!JR Also: "k" and "k3d" indices are both needed below, as "k3d" indexes into fully 3D arrays (nvp), and "k"
!JR indexes into arrays dimensioned (maxchunk). maxchunk is probably much smaller than nvp.
    do n=1,nchunks
!$OMP PARALLEL DO PRIVATE (k,ipn,j,i,pid,piu,icos2d)
      do k3d=kbeg(n),kend(n)
        k = k3d - kbeg(n) + 1
        CALL sptez_mod (0, maxwv, idrt, imax, jmax, data%t(:,k3d), vtmp2d(:,:,k), 1)
        CALL bilinear_interp(vtmp2d(1,1,k), icos2d)
        do ipn=1,nip
          t_lyr(k3d,ipn) = icos2d(ipn)
        end do
        do j=1,jmax
          do i=1,imax
            pid = ak(k3d)   + bk(k3d  )*f1(i,j)
            piu = ak(k3d+1) + bk(k3d+1)*f1(i,j)
            if (head%idsl == 2) then
              pl2d(i,j,k) = (pid+piu)/2
            else 
              pl2d(i,j,k) = ((pid**rocp1 - piu**rocp1)/(rocp1*(pid - piu)))**rocpr
            end if
          end do
        end do
      end do
!$OMP END PARALLEL DO

!JR Cannot thread the loop that is a sum over k
      do k3d=kbeg(n),kend(n)
        k = k3d - kbeg(n) + 1
!$OMP PARALLEL DO PRIVATE (i,gauss2d_lo,gauss2d_up)
        do j=1,jmax
          do i=1,imax
            gauss2d_lo = cp*((ak(k3d  ) + bk(k3d  )*f1(i,j))/p1000)**(rd/cp)
            gauss2d_up = cp*((ak(k3d+1) + bk(k3d+1)*f1(i,j))/p1000)**(rd/cp)
            gauss2d(i,j) = gauss2d(i,j) + (gauss2d_lo-gauss2d_up)*vtmp2d(i,j,k)*(p1000/pl2d(i,j,k))**(rd/cp)
          end do
        end do
!$OMP END PARALLEL DO

        CALL bilinear_interp (gauss2d, icos2d)
        do ipn=1,nip
          gz_lev(k3d+1,ipn) = icos2d(ipn)
        end do
      end do
    end do

    write (*,'(a/(5f14.6))') 'ak array:',(head%ak(k),k=1,nvp+1)
    write (*,'(a/(5f14.6))') 'bk array:',(head%bk(k),k=1,nvp+1)

    ret = gptlstop (thisfunc)
    ret = gptlprint_memusage ('end '//thisfunc)
  end subroutine transform_interp_4

  subroutine transform_interp_5 (u_lyr, v_lyr)
    real(4), intent(out) :: u_lyr(nvp,nip)
    real(4), intent(out) :: v_lyr(nvp,nip)

! Local workspace
    character(len=18), parameter :: thisfunc = 'transform_interp_5'
    real(4) :: icos2d1(nip)
    real(4) :: icos2d2(nip)
    real(4) :: gauss2d1(imax,jmax)
    real(4) :: gauss2d2(imax,jmax)
    integer :: k,ipn
    integer :: ret

    ret = gptlstart (thisfunc)
! perform spherical transform on u wind and v wind field
!$OMP PARALLEL DO PRIVATE (ipn, gauss2d1, gauss2d2, icos2d1, icos2d2)
    do k=1,nvp
      call sptezmv_mod (0, maxwv, idrt, imax, jmax, 1, data%d(:,k), data%z(:,k), &
                        gauss2d1, gauss2d2, 1)
!      CALL bilinear_interp_uv(gauss2d1, icos2d1, gauss2d2, icos2d2)
      CALL bilinear_interp_uv2(gauss2d1, icos2d1, gauss2d2, icos2d2)
      do ipn=1,nip
        u_lyr(k,ipn) = icos2d1(ipn)
        v_lyr(k,ipn) = icos2d2(ipn)
      end do
    end do
!$OMP END PARALLEL DO
    
    ret = gptlstop (thisfunc)
    ret = gptlprint_memusage ('end '//thisfunc)
  end subroutine transform_interp_5

  subroutine transform_interp_6 (o3_lyr, qc_lyr)
    real(4),intent(out) :: o3_lyr(nvp,nip)
    real(4),intent(out) :: qc_lyr(nvp,nip)

! Local workspace
    character(len=18), parameter :: thisfunc = 'transform_interp_6'
    integer :: k,ipn
    integer :: ret
    real(4) :: gauss2d(imax,jmax)
    real(4) :: icos2d(nip)                    ! work array on icos grid

    ret = gptlstart (thisfunc)

!$OMP PARALLEL DO PRIVATE (ipn, icos2d, gauss2d)
    do k=1,nvp
! perform spherical transform on ozone (g1) field
      call sptez_mod (0, maxwv, idrt, imax, jmax, data%q(:,k,2), gauss2d, 1)
      CALL bilinear_interp (gauss2d, icos2d)
      do ipn=1,nip
        o3_lyr(k,ipn) = max( icos2d(ipn), 0. )
      end do
      call sptez_mod (0, maxwv, idrt, imax, jmax, data%q(:,k,3), gauss2d, 1)
      CALL bilinear_interp (gauss2d, icos2d)
      do ipn=1,nip
        qc_lyr(k,ipn) = max( icos2d(ipn), 0. )
      end do
    end do
!$OMP END PARALLEL DO

    ret = gptlstop (thisfunc)
    ret = gptlprint_memusage ('end '//thisfunc)
  end subroutine transform_interp_6

  subroutine ss2ggtopo ()

! --- this routine operates on input data consisting of
! ---   (1) a vertical profile of pressure and geopotential;
! ---   (2) a surface elevation that differs from the lowest value in the
!           geopotential profile.

! --- the routine extends/shortens the input profile to bring it into
! --- agreement with the prescribed surface elevation. it does so without
! --- perturbing the above-ground portion of the inut profile.
! --- the hydrostatic equation used in manipulating the vertical profile
! --- is discretized as  d(phi)/d(Pi)=-theta  where Pi is the Exner function.

! Local workspace
    integer :: k,k1,ipn
    real    :: exlo,exup,exsf,th_lyr,pkap,gznew,hmax,hmin,pmax,pmin
    real    :: gzold(nip)
    integer :: kgrnd(nip)
    logical :: vrbos

    print *,'entering ss2ggtopo (hydrostat.reconciliation of srf.hgt/srf.prs)' 

!SMS$PARALLEL (dh,ipn) BEGIN
    do ipn=1,nip			! horiz. loop
      vrbos=perm(ipn).eq.PrintIpnDiag

      if (vrbos) then
       print 102,perm(ipn),'  IN  psrf=',ps_lev(ipn)*.01,		&
        (k,p_lev(k,ipn)*.01,gz_lev(k,ipn)/grvity,t_lyr(k,ipn),k=1,nvp),	&
        nvp+1,p_lev(nvp+1,ipn)*.01,gz_lev(nvp+1,ipn)/grvity
 102   format (i7,' ss2ggtopo',a,f10.2/					&
        '     pressure    geopot      temp'/(i3,3f10.2))
      end if

       do k=1,nvp			! vert. loop
        if (gz_lev(k+1,ipn).gt.hs_lev(ipn)*grvity) then

! --- level k+1 is above ground. integrate hydrostat.eqn down from there.
! --- level k may be either above or below ground.

! --- sequence of operations:
! --- (a) get old midlayer p^kappa from (partial p^(1+kappa))/(partial p)
! --- (b) get old theta from
! ---     partial phi_old / partial pi_old = -theta_old
! --- (c) set theta_new = theta_old (not optimal, but tolerable for now)
! --- (d) get new bottom pressure from
! ---     partial phi_new / partial pi_new = -theta_new
! --- (e) get new surf.temp. from theta_new and new bottom pressure

          exup=cp*(p_lev(k+1,ipn)/p1000)**(rd/cp)
          exlo=cp*(p_lev(k  ,ipn)/p1000)**(rd/cp)
          th_lyr=(gz_lev(k+1,ipn)-gz_lev(k,ipn))/(exlo-exup)
          exsf=exup+(gz_lev(k+1,ipn)-hs_lev(ipn)*grvity)/th_lyr	! new srf.exner
          ps_lev(ipn)=p1000*(exsf/cp)**(cp/rd)			! new srf.pres.
          if (exsf.gt.exlo+.01) then
           pkap=(exsf*ps_lev(ipn)-exup*p_lev(k+1,ipn))/	&
            ((rd+cp)*(ps_lev(ipn)-     p_lev(k+1,ipn)))		! mid-lyr p^kap
          else
           pkap=.5*(exsf+exup)/cp
          end if
          p_lev(1,ipn)=ps_lev(ipn)				! new srf.pres.
          t_lyr(1,ipn)=th_lyr*pkap				! new srf.temp.
          gz_lev(1,ipn)=hs_lev(ipn)*grvity			! new srf.geopot
          do k1=2,k
            p_lev (k1,ipn)=p_lev (1,ipn)
            t_lyr (k1,ipn)=t_lyr (1,ipn)
            gz_lev(k1,ipn)=gz_lev(1,ipn)
          end do
          kgrnd(ipn)=k
          gzold(ipn)=gz_lev(k+1,ipn)
          exit
        end if
      end do			! vert. loop
      if (vrbos) then
       print 101,'ss2ggtopo  ipn,k=',perm(ipn),k,			&
       'exup',exup,'exlo',exlo,'exsf',exsf,'th_lyr',th_lyr,		&
       'p_lev(k+1)',p_lev(k+1,ipn)*.01,'p_lev(k)',p_lev(k,ipn)*.01,	&
       'gz(k+1)',gz_lev(k+1,ipn)/grvity,'gz(k)',gz_lev(k,ipn)/grvity,	&
       'hs_lev',hs_lev(ipn),'cp*pkap',cp*pkap,'t_lyr(1)',t_lyr(1,ipn),	&
       'th_lyr',th_lyr,'gz_lev(1)',gz_lev(1,ipn)/grvity
 101   format (a,i9,i4/(4(a11,f9.3)))
       print 102,perm(ipn),' OUT  psrf=',ps_lev(ipn)*.01,		&
       (k,p_lev(k,ipn)*.01,gz_lev(k,ipn)/grvity,t_lyr(k,ipn),k=1,nvp),	&
        nvp+1,p_lev(nvp+1,ipn)*.01,gz_lev(nvp+1,ipn)/grvity
      end if
    end do			! horiz. loop
!SMS$PARALLEL END

    hmin = minval(hs_lev(1:nip))
    hmax = maxval(hs_lev(1:nip))
    pmin = minval(ps_lev(1:nip))
    pmax = maxval(ps_lev(1:nip))
!SMS$REDUCE(hmax,pmax,max)
!SMS$REDUCE(hmin,pmin,min)

    print 100,'min,max of new srf.height on icos grid:',hmin,hmax
    print 100,'min,max of new srf.press  on icos grid:',pmin,pmax
100 format (a,2f13.2)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! optional: re-compute geopotential on interfaces to check for errors
!SMS$PARALLEL (dh,ipn) BEGIN
    do ipn=1,nip
      vrbos=perm(ipn).eq.PrintIpnDiag
      gznew=gz_lev(1,ipn)
      exup=cp*(p_lev(1,ipn)/p1000)**(rd/cp)
      do k=1,kgrnd(ipn)
        exlo=exup
        exup=cp*(p_lev(k+1,ipn)/p1000)**(rd/cp)
        if (exlo.gt.exup+.01) then
          pkap=(exlo*p_lev(1,ipn)-exup*p_lev(k+1,ipn))/			&
           ((rd+cp)*(p_lev(1,ipn)-     p_lev(k+1,ipn)))
        else
          pkap=.5*(exlo+exup)/cp
        end if
        gznew=gznew+(exlo-exup)*t_lyr(k,ipn)/pkap

        if (vrbos) print 101,'ss2ggtopo  ipn,k=',perm(ipn),k,		&
        'p_lev(k)',p_lev(k,ipn)*.01,'p_lev(k+1)',p_lev(k+1,ipn)*.01,	&
        'gz(k)',gz_lev(k,ipn)/grvity,'gz(k+1)',gz_lev(k+1,ipn)/grvity,	&
        'exup',exup,'exlo',exlo,'cp*pkap',cp*pkap,'t_lyr',t_lyr(k,ipn),	&
        'th_lyr',t_lyr(k,ipn)/pkap,'gzold',gzold(ipn)/grvity,		&
        'gznew',gznew/grvity
      end do
      k=kgrnd(ipn)
      if (abs(gznew-gzold(ipn)) .gt. .05)				&
         print '(a,2i8,f11.1,f9.1,f8.2)',				&
         'height discrepancy at ipn,kgrnd =',ipn,k,gzold(ipn),gznew,	&
         gzold(ipn)-gznew
    end do
!SMS$PARALLEL END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    print *,'...exiting ss2ggtopo'
    return
  end subroutine ss2ggtopo


  subroutine solidrot(deg_lat,deg_lon,u,v,vrbos,ipn)
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! --- generate u,v components representing solid-rotation flow across pole
!
! --- transform lat/lon into spher.coord.system rotated 90 deg
! --- solid rotation will be about k axis in rotated system.

! --- there are 2 options:

! --- (1) rotation about vector i => make vectors i/j/k the new k/i/j vectors

! ---              / cos lat cos lon \
! --- unit vector |  cos lat sin lon  | is represented in the new system as
! ---              \ sin lat         /

! ---  / cos lat sin lon \     / cos latp cos lonp \
! --- |  sin lat          | = |  cos latp sin lonp  |
! ---  \ cos lat cos lon /     \ sin latp          /


! --- (2) rotation about vector j => make vectors i/j/k the new j/k/i vectors

! ---              / cos lat cos lon \
! --- unit vector |  cos lat sin lon  | is represented in the new system as
! ---              \ sin lat         /

! ---  / sin lat         \     / cos latp cos lonp \
! --- |  cos lat cos lon  | = |  cos latp sin lonp  |
! ---  \ cos lat sin lon /     \ sin latp          /

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    real,   intent(IN) :: deg_lat,deg_lon	! latitude/longitude (deg)
    real,  intent(OUT) :: u,v			! velocity components
    integer,intent(IN) :: ipn
    logical,intent(IN) :: vrbos
    real*8 alat,alon,alatp,alonp,dlat,dlon,dlatp,dlonp
    real,parameter :: small=1.e-9, speed=100.

#if ( defined NEED_SINDCOSD )
!JR Define required statement functions for situations where
!JR compiler doesn't support them (e.g. gfortran)
    real*8 :: val, sind, cosd, asind, acosd, x, y, atan2d, tand
    real*8, parameter :: pi = 3.14159265358979
    sind(val) = sin(val*pi/180.)
    cosd(val) = cos(val*pi/180.)
    tand(val) = tan(val*pi/180.)
    
    asind(val) = (180./pi)*asin(val)
    acosd(val) = (180./pi)*acos(val)
    atan2d(x,y) = (180./pi)*atan2(x, y)
#endif

    alat=deg_lat
    alon=deg_lon

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- option 1: rotation about original vector 

! --- forward transform:
    alatp = asind(cosd(alat)*cosd(alon))
    alonp = atan2d(tand(alat),sind(alon))

! --- solid rotation: compute derivatives w.r.t. transformed longitude
! --- (using 'back transform' formulae given below)
    dlat = cosd(alatp)*cosd(alonp)/max(small*1._8,sqrt(1.-(cosd(alatp)*sind(alonp))**2))	! return 0 if 0/0
    dlon =-tand(alatp)*sind(alonp)/max(small*1._8,tand(alatp)**2+cosd(alonp)**2)		! return 0 if 0/0

    u = dlon*cosd(alat) * speed
    v = dlat            * speed

    if (vrbos) then
      print 101,'(solrot1) la/lo,latp/lonp,u,v=',alat,alon,alatp,alonp,u,v

! --- back transform (diagnostic use):
      alat=asind(cosd(alatp)*sind(alonp))
      alon=atan2d(cosd(alonp),tand(alatp))
      alon=mod(alon+360._8,360._8)
      
      print 101,'(solrot1) alat/lon,alatp/lonp=',alat,alon,alatp,alonp

      dlatp=-cosd(alat )*sind(alon )/max(small*1._8,sqrt(1.-(cosd(alat )*cosd(alon ))**2))	! return 0 if 0/0
      dlonp=-tand(alat )*cosd(alon )/max(small*1._8,tand(alat )**2+sind(alon )**2)		! return 0 if 0/0

      print 101,'(solrot1) dlat/lon,dlatp/lonp=',dlat,dlon,dlatp,dlonp

! ---check derivatives by evaluating them in finite difference form
      dlat =(asind (cosd(alatp)*sind(alonp+small))		&
            -asind (cosd(alatp)*sind(alonp-small)))/(2.*small)
      dlon =(atan2d(cosd(alonp+small),tand(alatp))		&
            -atan2d(cosd(alonp-small),tand(alatp)))/(2.*small)
      dlatp=(asind (cosd(alat)*cosd(alon+small))			&
            -asind (cosd(alat)*cosd(alon-small)))/(2.*small)
      dlonp=(atan2d(tand(alat),sind(alon+small))			&
            -atan2d(tand(alat),sind(alon-small)))/(2.*small)
      if (dlon .gt. 179./(2.*small)) dlon =dlon -180./(2.*small)
      if (dlon .lt.-179./(2.*small)) dlon =dlon +180./(2.*small)
      if (dlonp.gt. 179./(2.*small)) dlonp=dlonp-180./(2.*small)
      if (dlonp.lt.-179./(2.*small)) dlonp=dlonp+180./(2.*small)
      
      print 101,'(solrot1) dlat/lon,dlatp/lonp=',dlat,dlon,dlatp,dlonp

! --- print l.h.s., r.h.s. of the 3 equations we are solving
      print 100,cosd(alat)*sind(alon),cosd(alatp)*cosd(alonp)
      print 100,sind(alat)           ,cosd(alatp)*sind(alonp)
      print 100,cosd(alat)*cosd(alon),sind(alatp)
    end if			! vrbos

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! ! --- option 2: rotation about original vector j
! 
! ! --- forward transform:
!    alatp=asind(cosd(alat)*sind(alon))
!    alonp=atan2d(cosd(alon),tand(alat))
! 
! ! --- solid rotation: compute derivatives w.r.t. transformed longitude
! ! --- (using 'back transform' formulae given below)
!    dlat = -cosd(alatp)*sind(alonp)				&
!     /max(small,sqrt(1.-(cosd(alatp)*cosd(alonp))**2))	! return 0 if 0/0
!    dlon =-tand(alatp)*cosd(alonp)				&
!     /max(small,tand(alatp)**2+sind(alonp)**2)		! return 0 if 0/0
! 
!    u=dlon*cosd(alat) * speed
!    v=dlat            * speed
! 
!    if (vrbos) then
!     print 101,'(solrot2) la/lo,latp/lonp,u,v=',alat,alon,alatp,alonp,u,v
! 
! ! --- back transform (diagnostic use):
!     alat=asind(cosd(alatp)*cosd(alonp))
!     alon=atan2d(tand(alatp),sind(alonp))
!     alon=mod(alon+360.,360.)
! 
!     print 101,'(solrot2) alat/lon,alatp/lonp=',alat,alon,alatp,alonp
! 
!     dlatp= cosd(alat )*cosd(alon )				&
!       /max(small,sqrt(1.-(cosd(alat )*sind(alon ))**2))	! return 0 if 0/0
!     dlonp=-tand(alat )*sind(alon )				&
!       /max(small,tand(alat )**2+cosd(alon )**2)		! return 0 if 0/0
! 
!     print 101,'(solrot2) dlat/lon,dlatp/lonp=',dlat,dlon,dlatp,dlonp
! 
! ! ---check derivatives by evaluating them in finite difference form
!     dlat =(asind (cosd(alatp)*cosd(alonp+small))		&
!           -asind (cosd(alatp)*cosd(alonp-small)))/(2.*small)
!     dlon =(atan2d(tand(alatp),sind(alonp+small))		&
!           -atan2d(tand(alatp),sind(alonp-small)))/(2.*small)
!     dlatp=(asind (cosd(alat)*sind(alon+small))			&
!           -asind (cosd(alat)*sind(alon-small)))/(2.*small)
!     dlonp=(atan2d(cosd(alon+small),tand(alat))			&
!           -atan2d(cosd(alon-small),tand(alat)))/(2.*small)
!     if (dlon .gt. 179./(2.*small)) dlon =dlon -180./(2.*small)
!     if (dlon .lt.-179./(2.*small)) dlon =dlon +180./(2.*small)
!     if (dlonp.gt. 179./(2.*small)) dlonp=dlonp-180./(2.*small)
!     if (dlonp.lt.-179./(2.*small)) dlonp=dlonp+180./(2.*small)
! 
!     print 101,'(solrot2) dlat/lon,dlatp/lonp=',dlat,dlon,dlatp,dlonp
! 
! ! --- print l.h.s., r.h.s. of the 3 equations we are solving
!     print 100,sind(alat)           ,cosd(alatp)*cosd(alonp)
!     print 100,cosd(alat)*cosd(alon),cosd(alatp)*sind(alonp)
!     print 100,cosd(alat)*sind(alon),sind(alatp)
!    end if			! vrbos
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
100 format ('equation check (numbers should match)',2f11.5)
101 format (a,4f9.3,2f7.1)

    return
  end subroutine solidrot
end module ss_gfs
