module module_fim_dyn_init

  use findmaxmin2
  use stencilprint
  use stenedgprint
  use module_control  ,only: control,nvlp1,nip,dt,numphr,ArchvStep,ntra,nabl,  &
                             npp,nd,ntrb,ipnDiagLocal,ipnDiagPE,               &
                             TestDiagProgVars,TestDiagNoise,nvlsig,kbl
  use fimnamelist     ,only: nvl,TotalTime,ArchvIntvl,ArchvTimeUnit,itsStart,  &
                             readrestart,fimout,gribout,glvl,PrintDiags,curve, &
                             PrintIpnDiag,PrintMAXMINtimes,FixedGridOrder,     &
                             NumCacheBlocksPerPE,yyyymmddhhmm,nts,             &
                             PrintDiagProgVars, PrintDiagNoise,TimingBarriers, &
                             dt_reducer_numerator,dt_reducer_denominator,      &
                             enkfio_in,PhysicsInterval,biharm_frst,biharm_last,&
                             veldff_bkgnd, veldff_boost,intfc_smooth,          &
                             debugmsg_on,atm_ic
  use module_constants,only: lat,lon,nprox,proxs,prox,area,cs,sn,sidevec_c,pi, &
                             sidevec_e,sideln,rprox_ln,deg_lat,deg_lon,rarea,  &
                             omegx2,raddeg,rsideln,corio,grndspd,thetac,       &
                             dpsig,p1000,cp,rd,qvmin,qwmin,perm,inv_perm,      &
                             nedge,permedge,idxtrl,idxldg,wgttrl,wgtldg,       &
                             trlpt,ldgpt,earthrad,dfflen,dfflen_max,biharw,    &
                             smoo_coeff, corner_xy, hfcoef
  use module_variables,only: u_tdcy,v_tdcy,dp_tdcy,trc_tdcy,trl_tdcy,          &
    massfx,nf,of,vof,us3d,vs3d,ws3d,dp3d,tk3d,			               &
    pr3d,ex3d,tr3d,mp3d,ph3d,rh3d,trdp,relvor,			               &
    sdot,dpl_tdcy,work2d,iwork2d,pw2d,pq2d,psrf,ptdcy,		               &
    dpinit,cumufx,u_edg,v_edg
  use module_sfc_variables,only: rn2d,sn2d,rc2d,qf2d, u10m, v10m, rlut
  use module_core_setup   ,only: iam_compute_root,iam_fim_task,iam_write_task,&
    my_comm,nct,nwt,output_intercomm,use_write_tasks
  use module_dyn_alloc    ,only: dyn_alloc
  use module_hybgen       ,only: hybgen
  use module_hystat       ,only: hystat
  use module_transp3d     ,only: transp0
  use module_wtinfo       ,only: wtinfo
  use global_bounds       ,only: set_global_bounds, set_global_bounds_halo, ims, ime, ips, ipe, ihs, ihe, myrank
  use module_cornwgt      ,only: rdcnfig, cornwgt
  USE module_enkf_io, only: read_enkfio
  use icosio,   only: icosio_setup
  use post,     only: post_init_readnl, post_init_slint, post_init_readnl_called, post_init_slint_called
  use sigio_module
  use units,    only: getunit, returnunit
  use ss_gfs,   only: ss2icos
  use mdul_transects,only: trnsec_init

  implicit none

contains

  subroutine dyn_init (client_server_io_in)
!*********************************************************************
!       Loads the initial variables and constants for the dynamics
!       component.
!       Alexander E. MacDonald  11/27/04
!       J. Lee                  September, 2005
!*********************************************************************

    use headers,only: testcurveheader,testglvlheader

    logical, intent(in), optional :: client_server_io_in

#include <gptl.inc>

! Local variables

!SMS$DISTRIBUTE(dh,1) BEGIN
    logical,allocatable :: in_halo(:)
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE(dh,3) BEGIN
    real,allocatable :: work_edg(:,:,:)  ! array for printing sidevec stencils
!SMS$DISTRIBUTE END

    character*24  :: string
    integer       :: ipn               ! Index for icos point number
    integer       :: isn               ! Index for icos edge number
    integer       :: ivl               ! Index for vertical level
    integer       :: i,its,idx
    integer       :: ipnGlobal
    logical       :: DiagPrint
    logical       :: vrbos

    integer :: iret
    integer :: ndamp
    real    :: wgt
    integer :: lusig
    type(sigio_head)       :: head
    CHARACTER(len=9 )      :: jdate
    CHARACTER(len=2 )      :: hh
    CHARACTER(len=80)      :: sanlFile

! IBM follows the Fortran standard precisely here requiring all
! literal constants to have the same number of characters.  So
! whitespace is significant.
    character(11),dimension(0:3) :: CurveType = (/'IJ         ', &
      'Hilbert    ', &
      'IJ block   ', &
      'SquareBlock'/)
    integer :: ipnx,isnx,ip1,im1,np1,nm1,ierr
    integer :: isncount
    integer :: nprocs = 1
    integer :: ios
!SMS$insert integer              :: HaloSize
!SMS$insert integer,allocatable  :: RegionSize(:)
!SMS$insert CHARACTER(len=80)    :: DecompInfoFile

    integer :: icosio_comm                ! an MPI communicator to pass to icosio
    integer :: idum1, idum2, idum3, idum4 ! dummies
    logical :: client_server_io           ! run icosio in client-server mode?
    logical :: ldum1, ldum2, ldum3        ! dummies
    integer :: unitno                     ! Unit number for I/O
    integer :: ret                        ! return code from GPTL routines
    integer :: k, n, t, k1, k2            ! indices
    real*8  :: tdyn_input                 ! timing for this routine
    real :: taper
    real,parameter :: divby6=1./6.
    real :: utrl,u1trl,u2trl,u3trl,uldg,u1ldg,u2ldg,u3ldg,		&
            vtrl,v1trl,v2trl,v3trl,vldg,v1ldg,v2ldg,v3ldg
    integer ::  i1trl,i2trl,i3trl,i1ldg,i2ldg,i3ldg
    integer :: edg, edgcount

!   taper(k1,k2)=float(k2-k1)**2/(19.+float(k2-k1)**2)  ! range: 0.05...1
    taper(k1,k2)=float(k2-k1)**2/( 9.+float(k2-k1)**2)  ! range: 0.1....1
!   taper(k1,k2)=float(k2-k1)**2/( 4.+float(k2-k1)**2)  ! range: 0.2....1
!   taper(k1,k2)=float(k2-k1)**4/(19.+float(k2-k1)**4)  ! range: 0.05...1
!   taper(k1,k2)=float(k2-k1)**4/( 9.+float(k2-k1)**4)  ! range: 0.1....1
!   taper(k1,k2)=float(k2-k1)**4/( 4.+float(k2-k1)**4)  ! range: 0.2....1

    print *,'entering dyn_init ...'
    ret = gptlstart ('dyn_init')

!sms$comm_size(nprocs)
    if (iam_write_task) nprocs = nct

    call control ()

! Initialize post. If gribout is enabled, ensure there is a max of 1 write task.

    call post_init_readnl(iret)
    if (iret /= 0) then
      write(6,*)'dyn_init: bad return from post_init: stopping'
      stop
    end if

    if (.not. post_init_readnl_called) then
      write(6,*)'dyn_init: logic error: cannot test gribout until post_init has been called'
      stop
    end if

    if (gribout .and. nwt > 1) then
      write(6,*)'dyn_init: Cannot have more than 1 write task when gribout enabled'
      stop
    end if

    if (iam_write_task) then
      call post_init_slint (iret)
      if (iret /= 0) then
        write(6,*)'dyn_init: bad return from post_init_slint: stopping'
        stop
      endif
      if (.not. post_init_slint_called) then
        write(6,*)'dyn_init: logic error: cannot test gribout until post_init_slint has been called'
        stop
      endif
    endif

!TODO:  remove dependence of serial prep code on ComputeTasks!
!SMS$insert ios = 0 ! avoid Lahey complaint on non-root tasks
!SMS$insert allocate(RegionSize(nprocs))
!SMS$insert write(DecompInfoFile,"('DecompInfo_',i0,'.dat')") nprocs
!SMS$insert unitno = getunit ()
!SMS$insert if (unitno < 0) then
!SMS$insert   print*,'dyn_init: getunit failed for file=', trim(DecompInfoFile),'. Stopping'
!SMS$insert   stop
!SMS$insert end if
!SMS$insert open (unitno,file=TRIM(DecompInfoFile),status='old',iostat=ios)
!SMS$insert if (ios /= 0) then
!SMS$insert   print*,'ERROR:  Cannot find ',TRIM(DecompInfoFile),', prep and fim must be run with the same setting for ComputeTasks in FIMnamelist.'
!SMS$insert   stop
!SMS$insert endif
!SMS$insert read (unitno,*) HaloSize
!SMS$insert read (unitno,*) RegionSize
!SMS$insert close(unitno)
!SMS$insert call returnunit (unitno)

!SMS$CREATE_DECOMP(dh,<nip>,<HaloSize>:regionsize=RegionSize)

! Now that SMS$create_decomp has happened, can define some bounds indices.
! CRITICAL!!! Since SMS$UNSTRUCTURED_GRID hasnt happened yet, ihs and ihe
! will not be set correctly. So we need to call set_global_bounds_halo()
! later after the SMS$UNSTRUCTURED_GRID directive below.

    call set_global_bounds ()  ! define ims, ime, ips, ipe, myrank, npes

    if (.not.iam_write_task) then

      call GPUinit(nprocs,myrank,ret)
      if(ret /= 0) then
        print*,'Error in dyn_init.F90 calling GPUinit.cu',ret
        stop
      endif

      if (atm_ic == 1) then			! GFS
        call GetJdate(yyyymmddhhmm,jdate)	! Julian date conversion
        hh = yyyymmddhhmm(9:10)
        sanlFile = jdate // ".gfs.t" // hh // "z.sanl"
      elseif (atm_ic == 2) then			! CFSR
        sanlFile = "siganl.gdas." // yyyymmddhhmm(1:10)
      end if
      PRINT *,' read vertical levels from input data ',sanlFile

!SMS$SERIAL (<nvlsig,OUT> : default=ignore)  BEGIN
! "82" is a magic unit that might have to be byte-swapped
      lusig = getunit(82)
      IF (lusig < 0) THEN
        PRINT*,'dyn_init: getunit failed for unit=82. Stopping'
        STOP
      END IF

      call sigio_sropen(lusig,sanlFile,iret)
      IF (iret .NE. 0) THEN
          PRINT '(a)','dyn_init: error reading '//sanlFile
          STOP
      END IF

      call sigio_srhead(lusig,head,iret)
      call sigio_sclose(lusig,iret)

      CALL returnunit (lusig)
      nvlsig = head%levs                       ! # of layers in input grid
!SMS$SERIAL END

      call dyn_alloc ()
      allocate(in_halo(nip))

! Changes in this block may require changes in output.F90 as well
      ArchvIntvl = min (ArchvIntvl, TotalTime)

      print"('FIM Global Model')"
      print *, " "
      call datetime ()
      print *, " "
      print"(' Curve:    '                           ,A44              )",CurveType(curve)
      print"(' Number of cache blocks per processor:',I8,' blocks'     )",NumCacheBlocksPerPE
      print"(' Grid Level'                           ,I35              )",glvl
      print"(' Number of Processors:'                ,I24,' processors')",nprocs
      print"(' Global size:     '                    ,I28,' points'    )",nip
!SMS$insert print"(' Halo size:       '                    ,I28,' points'    )",HaloSize
      print"(' Forecast duration (',a2,'):'          ,I23)",ArchvTimeUnit,TotalTime
      print"(' Vertical resolution:'                 ,I25,' levels'    )",nvl

      if ((glvl>7 .or. nip==163842) .and.ArchvTimeUnit/='ts'.and.ArchvTimeUnit/='mi'.and.ArchvTimeUnit/='mo') then
        print"(' Default time step:'                   ,I27,' seconds'   )",nint(dt*dt_reducer_denominator/dt_reducer_numerator)
        print"(' Time step reduced to '                ,I24,' seconds'   )",nint(dt)
      elseif ((nip==368642) .and.ArchvTimeUnit/='ts'.and.ArchvTimeUnit/='mi'.and.ArchvTimeUnit/='mo') then
        print"(' Default time step:'                   ,I27,' seconds'   )",60
        print"(' Time step reduced to '                ,I24,' seconds'   )",nint(dt)
      else
        print"(' Length of Time step: '                ,I24,' seconds'   )",nint(dt)
      end if

      print"(' Number of time steps:'                ,I24,' timesteps' )",nts
      print"(' Output every',I33,' timesteps')",ArchvStep
      print"(' Print MAX,MIN routine times',L18)",PrintMAXMINtimes
      print"(' Add timed barriers to measure task skew',L6)",TimingBarriers
      print"(' Output in fixed (IJ) grid ordering',L11)",FixedGridOrder
      print *,'nvlsig =                         ',nvlsig

      if (curve == 0) then !The grid order is already IJ for curve=0.
        FixedGridOrder = .false.
      end if

      if      (PrintDiagProgVars == 0) then
        TestDiagProgVars = 1
      else if (PrintDiagProgVars <  0) then
        TestDiagProgVars = 2*nts
      else
        TestDiagProgVars = PrintDiagProgVars*numphr
      end if

      if    (PrintDiagNoise == 0) then
        TestDiagNoise = 1
      elseif(PrintDiagNoise <  0) then
        TestDiagNoise = 2*nts
      else
        TestDiagNoise = PrintDiagNoise*numphr
      end if
! reduce frequency of diag
      TestDiagProgVars = max(TestDiagProgVars,ArchvStep)
      TestDiagNoise    = max(TestDiagNoise,   ArchvStep)

      print "(' Forecast initial time        ',A16,' YYYYMMDDHHMM')",yyyymmddhhmm
      print "(' Diagnostic prints at ipn           ',I10)",PrintIpnDiag
      print "(' Print diagnostic messages          ',L10)",PrintDiags
      print "(' Print diagnostic prognosis vars    ',I10)",PrintDiagProgVars
      print "(' Test/Print gravity wave noise      ',2I10)",	&
        TestDiagNoise,PrintDiagNoise
      print *,' '

#ifndef DEBUGPRINT
      if (PrintIpnDiag > 0) then
        print*,'Fatal error in dyn_init: PrintIpnDiag > 0 requires DEBUGPRINT=yes in macros.make.*'
        stop
      end if
#endif

! --- 'diffusion length' dfflen = (diffusivity) * (time step) / (mesh size)
! ---                           = (diffusion velocity) x (time step)
      ndamp = 0.1*nvl		! number of layers subjected to enhanced dissip
      do k=1,nvl
        wgt = max (0., real (k-nvl+ndamp) / real (ndamp))
        dfflen(k) = dt*(veldff_bkgnd*(1.-wgt) + max (veldff_bkgnd, veldff_boost)*wgt)
      end do
      dfflen_max = maxval (dfflen)
!
      smoo_coeff(:) = 0.
      do k=kbl+2,nvl	! don't modify lowest -kbl- layers
        smoo_coeff(k) = dt*intfc_smooth*taper(k,kbl+1)*taper(k,nvlp1)
      end do
!
! --- construct weighting function for blending laplacian and biharmonic dissip
! --- biharw = 0 for laplacian mixing, biharw = 1 for biharmonic mixing
      do k=1,nvl
        biharw(k)=0.
        if (k.ge.biharm_frst .and. k.le.biharm_last) biharw(k)=1.
! <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>
! --- optional: gradual transition between laplacian and biharmonic
        if (biharm_frst.gt.1 .and. biharm_frst.lt.nvl)			&
          biharw(k)=max(0.,min(1.,.5+.2*float(k-biharm_frst)))
        if (biharm_last.gt.1 .and. biharm_last.lt.nvl)			&
          biharw(k)=max(0.,min(1.,.5+.2*float(biharm_last-k)))
! <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>

! --- increase dfflen 5-fold in case of biharmonic mixing
        dfflen(k)=dfflen(k)*(1.+4.*biharw(k)**2)
      end do

      if (ips.eq.1) then
        print '(a/(10f8.1))',						&
          'diffusion length (diffusion velocity x time step):',		&
          dfflen
        print '(a/(10f8.3))',						&
          'blending coefficients for lapl/biharm u,v dissipation:',	&
          biharw
        print '(a/(10f8.1))',						&
          'diffusion for interface smoothing',smoo_coeff
      end if
!
! ..................................................................
! Sec. 1. Calculate icosahedral grid description data
! ..................................................................
! The variables to be read in are:

! lat(nip),lon(nip) - lat and lon of icos pts
! nprox(nip) = number of proximity points (6 or 5)
! proxs(npp,nip) = array holding icos side number (1 to 6)
! prox(npp,nip) = array holding index of icos cell across side
! area(nip) = area of cell
! sidevec_c(nd,npp,nip)= side vector projected from cell center
! sidevec_e(nd,npp,nip)= side vector projected from cell edge
! sideln(npp,nip) = length of side (edge)- local projection invarient
! rprox_ln(npp,nip) = reciprocal of length between icos pts
! inv_perm = inverse permutation of the grid

      unitno = getunit ()
      if (unitno < 0) then
        print*,'dyn_init: getunit failed for file=glvl.dat. Stopping'
        stop
      end if

      open (unitno, file='glvl.dat', form='unformatted', action='read', iostat=ios)
      if (ios /= 0) then
        print*,'dyn_init: cannot open glvl.dat for reading. Stopping'
        stop
      end if

! The following serial region is required because the called routines perform
! IO but are not processed by SMS, so that IO actions would otherwise occur on
! all tasks. TODO: These routines may call stop and, since they have not been
! processed by SMS, the stop may hang the run. They should either be processed
! by SMS, or simply return a status that the caller can act on.

!SMS$SERIAL (default=ignore) BEGIN
      call TestGlvlHeader (unitno, 'glvl.dat', 'dyn_init', glvl)
      call TestCurveHeader(unitno, 'glvl.dat', 'dyn_init', curve)
!SMS$SERIAL END

      read (unitno,err=90) lat
      read (unitno, err=90) lon
      read (unitno, err=90) nprox

      do isn=1,size(proxs,1)
        read (unitno, err=90) iwork2d
!sms$ignore begin
!$OMP PARALLEL DO
        do ipn=ips,ipe
          proxs(isn,ipn) = iwork2d(ipn)
        end do
!sms$ignore end
      end do

      do isn=1,size(prox,1)
        read (unitno, err=90) iwork2d
!sms$ignore begin
!$OMP PARALLEL DO
        do ipn=ips,ipe
          prox(isn,ipn) = iwork2d(ipn)
        end do
!sms$ignore end
      end do

      read (unitno, err=90) area

      do isn=1,size(cs,2)
        do idx=1,size(cs,1)
          read (unitno, err=90) work2d
!sms$ignore begin
!$OMP PARALLEL DO
          do ipn=ips,ipe
            cs(idx,isn,ipn) = work2d(ipn)
          end do
!sms$ignore end
        end do
      end do

      do isn=1,size(sn,2)
        do idx=1,size(sn,1)
          read (unitno, err=90) work2d
!sms$ignore begin
!$OMP PARALLEL DO
          do ipn=ips,ipe
            sn(idx,isn,ipn) = work2d(ipn)
          end do
!sms$ignore end
        end do
      end do

      do isn=1,size(sidevec_c,2)
        do idx=1,size(sidevec_c,1)
          read (unitno, err=90) work2d
!sms$ignore begin
!$OMP PARALLEL DO
          do ipn=ips,ipe
            sidevec_c(idx,isn,ipn) = work2d(ipn)
          end do
!sms$ignore end
        end do
      end do

      do isn=1,size(sidevec_e,2)
        do idx=1,size(sidevec_e,1)
          read (unitno, err=90) work2d
!sms$ignore begin
!$OMP PARALLEL DO
          do ipn=ips,ipe
            sidevec_e(idx,isn,ipn) = work2d(ipn)
          end do
!sms$ignore end
        end do
      end do

      do isn=1,size(sideln,1)
        read (unitno, err=90) work2d
!sms$ignore begin
!$OMP PARALLEL DO
        do ipn=ips,ipe
          sideln(isn,ipn) = work2d(ipn)
        end do
!sms$ignore end
      end do

      do isn=1,size(rprox_ln,1)
        read (unitno, err=90) work2d
!sms$ignore begin
!$OMP PARALLEL DO
        do ipn=ips,ipe
          rprox_ln(isn,ipn) = work2d(ipn)
        end do
!sms$ignore end
      end do

!SMS$SERIAL (<inv_perm,perm,printipndiag,out>:default=ignore) BEGIN
      read (unitno,iostat=ierr) inv_perm
      if (ierr.ne.0) then
        write (*,'(a)') 'dyn_init: error reading inv_perm'
        stop
      endif
      do i=1,nip
        perm(inv_perm(i))=i         ! store the inverse of inv_perm in -perm-
      end do

      print 105,'(dyn_init) global-to-local icos index mapping (first 300):',&
        (ipn,inv_perm(ipn),ipn=1,300)
      print 105,'(dyn_init) local-to-global icos index mapping (first 300):',&
        (ipn,perm(ipn),ipn=1,300)
105   format (a/(i7,' =>',i7,i10,' =>',i7,i10,' =>',i7,i10,' =>',i7))

! --- change test point from global to local
      if (PrintIpnDiag.gt.nip) then
        stop '(dyn_init: test point out of range)'
      end if
      if (PrintIpnDiag.gt.0) then
        print '(2(a,i8))','PrintIpnDiag=',PrintIpnDiag,			&
          ' becomes',inv_perm(PrintIpnDiag)
        PrintIpnDiag=inv_perm(PrintIpnDiag)
        if (PrintIpnDiag > nip) then
          print*,'Fatal error in dyn_init: PrintIpnDiag must be <= nip',	&
            PrintIpnDiag,nip
          stop
        end if
      end if
!SMS$SERIAL END

      do idx=1,size(corner_xy,2)
        do isn=1,size(corner_xy,1)
          read(unitno, err=90) work2d
!sms$ignore begin
!$OMP PARALLEL DO
          do ipn=ips,ipe
            corner_xy(isn,idx,ipn) = work2d(ipn)
          end do
!sms$ignore end
        enddo
      enddo

      do isn=1,size(hfcoef,2)
        do idx=1,size(hfcoef,1)
          read (unitno, err=90) work2d
!sms$ignore begin
!$OMP PARALLEL DO
          do ipn=ips,ipe
            hfcoef(idx,isn,ipn) = work2d(ipn)
          end do
!sms$ignore end
        end do
      end do

      close(unitno)
      call returnunit (unitno)

!sms$ignore begin
!$OMP PARALLEL DO
      do ipn=ips,ipe
        if(nprox(ipn)==5)then
          prox(6,ipn) = prox(5,ipn)
        end if
      end do
!sms$ignore end

! Initialize SMS for unstructured grid, including calculating the halo positions which are stored in
! prox.

!SMS$UNSTRUCTURED_GRID(PROX)

! Now that SMS$UNSTRUCTURED_GRID has happened, can properly define ihs, ihe
      call set_global_bounds_halo ()  ! define ihs, ihe

! Update halos of constant arrays.  Cannot rely on automatic halo update from
! file read due to prox not yet being set up prior to !SMS$UNSTRUCTURED_GRID.

!SMS$EXCHANGE(lat,lon,nprox,proxs,area,cs,sn,sidevec_c,sidevec_e,sideln,rprox_ln,corner_xy,hfcoef)

! set up nedge and permedge for computations in the halo (HALO_COMP)
!JR Dont bother to thread because the amount of work is too small
      do ipn=ipe+1,ihe
        in_halo(ipn) = .true.
        nedge(ipn) = 0
      end do

! set up nedge and permedge for interior cells
!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (isn)
      do ipn=ips,ipe
        in_halo(ipn) = .false.
        nedge(ipn) = nprox(ipn)     ! NOOP for interior cells
        do isn=1,nprox(ipn)
          permedge(isn,ipn) = isn   ! NOOP for interior cells
        end do
      end do
!$OMP END PARALLEL DO
!sms$ignore end

! now set up nedge and permedge for halo cells
! NOTE:  *not* owner-computes!
!JR DO NOT THREAD THIS LOOP! DEPENDENCE ON IPN
!sms$ignore begin
      do ipn=ips,ipe
        do isn=1,nprox(ipn)
          ipnx = prox(isn,ipn)
          if (in_halo(ipnx)) then
            nedge(ipnx) = nedge(ipnx) + 1
            permedge(nedge(ipnx),ipnx) = proxs(isn,ipn)
          end if
        end do
      end do
!sms$ignore end

#define DEBUG_HALO_COMP
!TBH:  lots of error checks for debugging
!JR: Dont thread so "stop" behaves correctly if invoked
#ifdef DEBUG_HALO_COMP
!SMS$PARALLEL(dh, ipn) BEGIN
! verify that in_halo is not screwed up in the interior
!$OMP PARALLEL DO
      do ipn=1,nip
        if (in_halo(ipn)) then
          print *,'ERROR C:  in_halo(',ipn,') = ',in_halo(ipn),'  but ipn is an interior cell!'
          stop
        end if
      end do
!$OMP END PARALLEL DO

! verify that prox(isnx,ipnx) == ipn in the interior
!$OMP PARALLEL DO PRIVATE (isn,ipnx,isnx)
      do ipn=1,nip
        do isn=1,nprox(ipn)
          ipnx = prox(isn,ipn)
          isnx = proxs(isn,ipn)
          if (.not.in_halo(ipnx)) then  ! avoid halo points which have not yet been set up
            if (prox(isnx,ipnx) /= ipn) then
              print *,'ERROR C:  prox(',isnx,',',ipnx,') /= ',ipn,'  me = ',myrank,prox(isnx,ipnx)
              stop
            end if
          end if
        end do
      end do
!$OMP END PARALLEL DO

! verify that nedge is OK in the interior
!$OMP PARALLEL DO
      do ipn=1,nip
        if (nedge(ipn) /= nprox(ipn)) then
          print *,'ERROR C:  nedge(',ipn,') /= nprox(',ipn,')  [',nedge(ipn),' /= ',nprox(ipn),']  me = ',myrank
          stop
        end if
      end do
! verify that permedge is OK in the interior
!$OMP PARALLEL DO PRIVATE (isn)
      do ipn=1,nip
        do isn=1,nedge(ipn)
          if (permedge(isn,ipn) /= isn) then
            print *,'ERROR C:  permedge(',isn,',',ipn,') /= ',isn,')  [',permedge(isn,ipn),' /= ',isn,']  me = ',myrank
            stop
          end if
        end do
      end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
#endif
#undef DEBUG_HALO_COMP

!SMS$PARALLEL(dh, ipn) BEGIN
! fix prox for halo cells
! NOTE:  *not* owner-computes!
! TODO:  Move this loop into SMS_UnstructuredGrid along with the exchange of proxs above ??
! TODO:  Would need to pass proxs into !SMS$UNSTRUCTURED_GRID too of course...

!JR Too dangerous to try threading
      do ipn=1,nip
        do isn=1,nprox(ipn)
          ipnx = prox(isn,ipn)
          isnx = proxs(isn,ipn)
          if (in_halo(ipnx)) then
! ipnx is a halo cell
! point prox for halo cell at ipnx back to me (ipn)
            prox(isnx,ipnx) = ipn
! handle pointing prox for halo cell at ipnx to any halo cells that
! are adjacent to itself (ipnx) and I (ipn)
! im1 and ip1 are my edges adjacent to cells that are also adjacent
! to the halo cell at ipnx
            im1 = mod(isn-2+nprox(ipn),nprox(ipn)) + 1
            ip1 = mod(isn             ,nprox(ipn)) + 1
! nm1 and np1 are edges of halo cell at ipnx that are also adjacent
! to me (ipn)
            nm1 = mod(isnx-2+nprox(ipnx),nprox(ipnx)) + 1
            np1 = mod(isnx              ,nprox(ipnx)) + 1
! note that nm1 and ip1 adjoin the same cell
            prox(nm1,ipnx) = prox(ip1,ipn)
! note that np1 and im1 adjoin the same cell
            prox(np1,ipnx) = prox(im1,ipn)
! note that prox(nm1,ipnx) and prox(np1,ipnx) are computed redundantly
          end if
        end do
      end do
!SMS$PARALLEL END

!!SMS$IGNORE BEGIN
!print *,'DEBUG:  setting halo_comp values in first layer of dh__S1 and dh__E1 by brute force'
!print *,'DEBUG:  collapsed_halo_size = ',collapsed_halo_size
!print *,'DEBUG:  dh__NestLevel = ',dh__NestLevel
!print *,'DEBUG:  dh__S1(  1,0,1) = ',dh__S1(  1,0,1)
!print *,'DEBUG:  dh__S1(  1,1,1) = ',dh__S1(  1,1,1)
!print *,'DEBUG:  dh__E1(nip,0,1) = ',dh__E1(nip,0,1)
!print *,'DEBUG:  dh__E1(nip,1,1) = ',dh__E1(nip,1,1)
!!SMS$IGNORE END

#define DEBUG_HALO_COMP
!TBH:  lots of error checks for debugging
#ifdef DEBUG_HALO_COMP
!SMS$PARALLEL(dh, ipn) BEGIN
! verify that in_halo is not screwed up in the interior
!$OMP PARALLEL DO
      do ipn=1,nip
        if (in_halo(ipn)) then
          print *,'ERROR X:  in_halo(',ipn,') = ',in_halo(ipn),'  but ipn is an interior cell!'
          stop
        end if
      end do
!$OMP END PARALLEL DO

! verify that prox(isnx,ipnx) == ipn is still true in the interior
!$OMP PARALLEL DO PRIVATE (isn,ipnx,isnx)
      do ipn=1,nip
        do isn=1,nprox(ipn)
          ipnx = prox(isn,ipn)
          isnx = proxs(isn,ipn)
          if (.not.in_halo(ipnx)) then  ! avoid halo points
            if (prox(isnx,ipnx) /= ipn) then
              print *,'ERROR X:  prox(',isnx,',',ipnx,') /= ',ipn,'  me = ',myrank
              stop
            end if
          end if
        end do
      end do
!$OMP END PARALLEL DO

!SMS$HALO_COMP(<1,1>) BEGIN
! verify that prox(isnx,ipnx) == ipn everywhere that we care
!$OMP PARALLEL DO PRIVATE (isncount,isn,ipnx,isnx)
      do ipn=1,nip
        do isncount=1,nedge(ipn)
          isn = permedge(isncount,ipn)
          ipnx = prox(isn,ipn)
          isnx = proxs(isn,ipn)
          if (prox(isnx,ipnx) /= ipn) then
            print *,'ERROR X:  prox(',isnx,',',ipnx,') /= ',ipn,'  me = ',myrank,'  in_halo(',ipn,') = ', &
              in_halo(ipn),'  in_halo(',ipnx,') = ',in_halo(ipnx)
            stop
          end if
        end do
      end do
!$OMP END PARALLEL DO
!SMS$HALO_COMP END
!SMS$PARALLEL END
#endif
#undef DEBUG_HALO_COMP

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (isn,ipnDiagLocal,ipnDiagPE)
      do ipn=ips,ipe
        corio(ipn) = omegx2*sin(lat(ipn))	! coriolis acceleration
        deg_lat(ipn) = raddeg*lat(ipn)	! latitude in degrees
        deg_lon(ipn) = raddeg*lon(ipn)	! longitude in degrees
        rarea(ipn) = 1./area(ipn)     	! reciprocal of area
        grndspd(ipn) = earthrad*2.*pi*cos(lat(ipn))/86164.	! ground speed
        call GetIpnGlobal (ipn, ipnGlobal, DiagPrint)
        if(DiagPrint) then
          ipnDiagLocal = ipn
          ipnDiagPE    = myrank
        end if
        do isn=1,nprox(ipn)
          rsideln(isn,ipn) = 1./sideln(isn,ipn)
        end do
      end do
!$OMP END PARALLEL DO
!sms$ignore end

!SMS$EXCHANGE(corio,grndspd,deg_lat,deg_lon,rarea,rsideln,perm)

      allocate (trlpt(npp,2,nip),ldgpt(npp,2,nip))
      call rdcnfig ()
      call cornwgt ()

      deallocate (trlpt,ldgpt)

! ...............................................................
! Sec. 2. Load initial state
! ...............................................................

      ret = gptlstart ('dyn_input')
      unitno = getunit ()
      if (unitno < 0) then
        print*,'dyn_init: getunit failed for file=theta_coor.txt. Stopping'
        stop
      end if

      open (unitno, file='theta_coor.txt', form='formatted', action='read', iostat=ios)
      if (ios /= 0) then
        print*,'dyn_init: cannot open theta_coor.txt for reading. Stopping'
        stop
      end if

      read (unitno, *, iostat=ios) thetac
      if (ios /= 0) then
        print*,'dyn_init: error reading theta_coor.txt. Stopping'
        stop
      end if

      close (unitno)
      print '(a/(10f8.1))','thetac (deg K):', thetac
! Re-use the same unit number that was just closed
      open (unitno, file='dpsig.txt', form='formatted', action='read', iostat=ios)
      if (ios /= 0) then
        print*,'dyn_init: cannot open dpsig.txt for reading. Stopping'
        stop
      end if

      read (unitno, *, iostat=ios) dpsig
      if (ios /= 0) then
        print*,'dyn_init: error reading dpsig.txt. Stopping'
        stop
      end if

      close (unitno)
      call returnunit (unitno)
      print '(a/(10f8.1))','dpsig (Pa):',dpsig

      IF (.NOT. readrestart) THEN

        IF (enkfio_in) THEN
          CALL read_enkfio ()
          CALL add_incr ()
        ELSE
          CALL ss2icos ()

!TBH:  ss2icos only sets tr3d(:,:,1:4)
!TBH:  so set remaining values to zero here
!TODO:  This may not be the correct approach.  Need to get
!TODO:  Georg and Rainer together to discuss.  At the moment Georg
!TODO:  believes that tr3d(:,:,5:ntra) are all zero anyway...

        END IF ! enkfio


!SMS$IGNORE BEGIN
! put initialization inside SMS "ignore" so all values are set to
! zero (including halo)

        IF (ntra+ntrb > 4) THEN
          tr3d(:,:,5:ntra+ntrb) = 0.0
        END IF

!SMS$IGNORE END


! Initialize post. If gribout is enabled, ensure there is a max of 1 write task.

!SMS$SERIAL (default=ignore) BEGIN
        CALL post_init_slint (iret)
        IF (iret /= 0) THEN
          WRITE(6,*)'dyn_init: bad return from post_init_slint: stopping'
          STOP
        ENDIF

        IF (.NOT. post_init_slint_called) THEN
          WRITE(6,*)'dyn_init: logic error: cannot test gribout until post_init_slint has been called'
          STOP
        ENDIF
!SMS$SERIAL END

        CALL findmxmn2(ph3d,nvlp1,nip,1,'surf.geopot.')
        ret = gptlstop ('dyn_input')
        ret = gptlget_wallclock ('dyn_input', 0, tdyn_input)  ! The "0" is thread number
        PRINT"(' DYNAMICS INPUT time:',F10.0)", tdyn_input

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (t,k)
        DO ipn=ips,ipe
          DO t=1,ntra+ntrb
            DO k=1,nvl  					! vertical loop
              trdp(k,ipn,t) = tr3d(k,ipn,t)*dp3d(k,ipn)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
!sms$ignore end

        its = itsStart - 1

        PRINT *,'dyn_init calling hybgen ...'

        CALL hybgen(its,         &
          thetac,                 & ! target pot.temperature
          us3d,vs3d,tr3d,         & ! zonal, merid.wind, mass field tracers
          sdot,ex3d,dp3d,pr3d,    & ! intfc displ., exner, lyr thknss, pres
          TimingBarriers)

!JR Can now set initial u10m, v10m since hybgen has modified us3d and vs3d
        print*,'returned from hybgen'

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (k)
        DO ipn=ips,ipe
          u10m(ipn) = us3d(1,ipn)
          v10m(ipn) = vs3d(1,ipn)
          pr3d(1,ipn) = p1000*(ex3d(1,ipn)/cp)**(cp/rd)
          DO k=1,nvl  !  vertical loop
            pr3d(k+1,ipn) = p1000*(ex3d(k+1,ipn)/cp)**(cp/rd)
            dp3d(k  ,ipn) = pr3d(k,ipn) - pr3d(k+1,ipn)
          END DO
          psrf(ipn) = 0.
          ptdcy(ipn,1) = 0.
          ptdcy(ipn,2) = 0.
        END DO
!$OMP END PARALLEL DO
!sms$ignore end

!-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>
      !  Optional: Calculate initial vorticity for use in -pvsurf- as an
      !  alternative to setting relvor=0 at t=0.
      !  Vorticity is calculated as a line integral of the tangential wind
      !  component given by the dot product of the wind vector with sidevec.
      !  Sidevec is a vectorial representation of the edge.
      !  The following are code fragments lifted from edgvar1.F90 and
      !  momtum_compute.F90
 
!SMS$EXCHANGE(us3d,vs3d)

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (vrbos,i1trl,i2trl,i3trl,i1ldg,i2ldg,i3ldg,	&
!$OMP     utrl,u1trl,u2trl,u3trl,uldg,u1ldg,u2ldg,u3ldg,		&
!$OMP     vtrl,v1trl,v2trl,v3trl,vldg,v1ldg,v2ldg,v3ldg,                &
!$OMP     edgcount,edg,k)
  do ipn=ips,ihe
    vrbos = ipn == PrintIpnDiag
    do edgcount=1,nedge(ipn)
  
! --- edge quantities are interpolated from 4 icos cells -- the cells on
! --- either side of the edge plus the 2 immediate neighbors of this pair.
  
      edg = permedge(edgcount,ipn)
      i1trl=idxtrl(edg,1,ipn)		! same as ipn
      i2trl=idxtrl(edg,2,ipn)		! cell across edge 'edg' from ipn
      i3trl=idxtrl(edg,3,ipn)		! cell across edge 'edg-1' from ipn
      i1ldg=idxldg(edg,1,ipn)		! same as ipn
      i2ldg=idxldg(edg,2,ipn)		! cell across edge 'edg' from ipn
      i3ldg=idxldg(edg,3,ipn)		! cell across edge 'edg+1' from ipn
#ifdef DEBUGPRINT
      if (vrbos) then
        print 100,its,perm(ipn),edg,					&
         'glb:', perm(i1trl),perm(i2trl),perm(i3trl),			&
                 perm(i1ldg),perm(i2ldg),perm(i3ldg),			&
         'loc:',        i1trl ,       i2trl ,       i3trl ,		&
                        i1ldg ,       i2ldg ,       i3ldg 
 100    format ('its=',i5,' ipn=',i8,' (edgvar1) edge',i2,		&
          '  interpolation based on cells...',2(/a,2(i14,2i8)))
        print 103,'wgt:',						&
         wgttrl(edg,1,ipn),wgttrl(edg,2,ipn),wgttrl(edg,3,ipn),		&
         wgtldg(edg,1,ipn),wgtldg(edg,2,ipn),wgtldg(edg,3,ipn)
 103    format (a,2(f14.4,2f8.4))
      end if
#endif
      do k=1,nvl

        !    Transform u,v at neighboring icos pt to local coord.system.
        !    cs and sn are coordinate transformation constants.
        !    u_xy,v_xy are values of u and v rotated into local system.
        !    (Unrolled to allow vectorization of the k loop)
 
        !   interpolate rotated wind components to edges
 
        u1trl =  cs(1,edg,ipn)*us3d(k,i1trl) + sn(1,edg,ipn)*vs3d(k,i1trl)
        v1trl = -sn(1,edg,ipn)*us3d(k,i1trl) + cs(1,edg,ipn)*vs3d(k,i1trl)
        u2trl =  cs(2,edg,ipn)*us3d(k,i2trl) + sn(2,edg,ipn)*vs3d(k,i2trl)
        v2trl = -sn(2,edg,ipn)*us3d(k,i2trl) + cs(2,edg,ipn)*vs3d(k,i2trl)
        u3trl =  cs(3,edg,ipn)*us3d(k,i3trl) + sn(3,edg,ipn)*vs3d(k,i3trl)
        v3trl = -sn(3,edg,ipn)*us3d(k,i3trl) + cs(3,edg,ipn)*vs3d(k,i3trl)

        utrl = u1trl*wgttrl(edg,1,ipn) + u2trl*wgttrl(edg,2,ipn)	&
             + u3trl*wgttrl(edg,3,ipn)
        vtrl = v1trl*wgttrl(edg,1,ipn) + v2trl*wgttrl(edg,2,ipn)	&
             + v3trl*wgttrl(edg,3,ipn)

        u1ldg =  cs(1,edg,ipn)*us3d(k,i1ldg) + sn(1,edg,ipn)*vs3d(k,i1ldg)
        v1ldg = -sn(1,edg,ipn)*us3d(k,i1ldg) + cs(1,edg,ipn)*vs3d(k,i1ldg)
        u2ldg =  cs(2,edg,ipn)*us3d(k,i2ldg) + sn(2,edg,ipn)*vs3d(k,i2ldg)
        v2ldg = -sn(2,edg,ipn)*us3d(k,i2ldg) + cs(2,edg,ipn)*vs3d(k,i2ldg)
        u3ldg =  cs(4,edg,ipn)*us3d(k,i3ldg) + sn(4,edg,ipn)*vs3d(k,i3ldg)
        v3ldg = -sn(4,edg,ipn)*us3d(k,i3ldg) + cs(4,edg,ipn)*vs3d(k,i3ldg)

        uldg = u1ldg*wgtldg(edg,1,ipn) + u2ldg*wgtldg(edg,2,ipn)	&
             + u3ldg*wgtldg(edg,3,ipn)
        vldg = v1ldg*wgtldg(edg,1,ipn) + v2ldg*wgtldg(edg,2,ipn)	&
             + v3ldg*wgtldg(edg,3,ipn)

! --- trapezoidal:
!       if (eqwgt) then
          u_edg(k,edg,ipn) = (2.*(u1ldg+u2ldg)+u3trl+u3ldg)*divby6
          v_edg(k,edg,ipn) = (2.*(v1ldg+v2ldg)+v3trl+v3ldg)*divby6
!       else
!         u_edg(k,edg,ipn) = .5*(utrl+uldg)
!         v_edg(k,edg,ipn) = .5*(vtrl+vldg)
!       end if

! --- simpson (assuming -i1,i2- are cells nearest edge midpoint):
!       if (simpson) then
!         if (eqwgt) then
!           u_edg(k,edg,ipn) = (8.*(u1ldg+u2ldg)+u3trl+u3ldg)*divby18
!           v_edg(k,edg,ipn) = (8.*(v1ldg+v2ldg)+v3trl+v3ldg)*divby18
!         else
!           u_edg(k,edg,ipn) = (u_edg(k,edg,ipn)+u1trl+u2trl)*divby3
!           v_edg(k,edg,ipn) = (v_edg(k,edg,ipn)+v1trl+v2trl)*divby3
!         end if
!       end if
      end do			! vert.loop
    end do			! loop through edges
  end do			! horiz.loop
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE (k,edg)
  do ipn=ips,ipe
    do k=1,nvl
      relvor(k,ipn) = 0.
    end do

! loop through edges and compute line integrals of bernoulli function,
! potential temperature, and pressure gradient

    do edg=1,nprox(ipn)
      do k=1,nvl
        relvor(k,ipn) = relvor(k,ipn)				&
            + ((sidevec_e(1,edg,ipn)*u_edg(k,edg,ipn)		&
            +   sidevec_e(2,edg,ipn)*v_edg(k,edg,ipn))) 
      end do
    end do
    do k=1,nvl
      relvor(k,ipn) = relvor(k,ipn)*rarea(ipn)
    end do
  end do
!$OMP END PARALLEL DO
!sms$ignore end
!-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>

        PRINT *,'dyn_init calling hystat ...'
        CALL hystat (its, ph3d, ex3d, mp3d, dp3d, &
          tr3d, trdp, psrf,ptdcy)

! ...............................................................
! Sec. 3. Initialize misc. variables
! ...............................................................

! Initialize tendency arrays
!JR Code as explicit loops for best cache performance
!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (n,k,t)
        DO ipn=ips,ihe
          DO n=1,nabl
            DO k=1,nvl
              u_tdcy(k,ipn,n)   = 0.		! u tendency
              v_tdcy(k,ipn,n)   = 0.		! v tendency
              dp_tdcy(k,ipn,n)  = 0.		! dp tendency
              dpl_tdcy(k,ipn,n) = 0.		! dp tendency, low order
            END DO
          END DO

          DO t=1,ntra+ntrb
            DO n=1,nabl
              DO k=1,nvl
                trc_tdcy(k,ipn,n,t) = 0.		! tracer tendency
                trl_tdcy(k,ipn,n,t) = 0.		! tracer tendency, low order
              END DO
            END DO
          END DO

          DO k=1,nvl
            ws3d(k,ipn) = 0.
            tk3d(k,ipn) = 0.
            rh3d(k,ipn) = 0.
          END DO

          DO t=1,nabl
            DO n=1,npp
              DO k=1,nvl
                massfx(k,n,ipn,t) = 0.	! mass flux (3 time levels)
              END DO
            END DO
          END DO

          pw2d(ipn) = 0.     ! precipitable water
          pq2d(ipn) = 0.     ! vertically integrated hydrometeor condensate
          rn2d(ipn) = 0.     ! accumulated precipitation/rainfall
          sn2d(ipn) = 0.     ! accumulated snowfall
          rc2d(ipn) = 0.
          rlut(ipn) = 0.
        END DO
!$OMP END PARALLEL DO
!sms$ignore end

! initial Adams-Bashforth indices
        nf  = 0	! "new field" index
        of  = 2	! "old field" index
        vof = 1	! "very old field" index

        IF (ntrb > 0) THEN
          CALL transp0(its,cumufx,dp3d,dpinit)	! initialize class B tracer transport
        END IF

      ELSE ! readrestart=.true.

!SMS$SERIAL (default=ignore) BEGIN
        CALL post_init_slint (iret)
        IF (iret /= 0) THEN
          WRITE(6,*)'dyn_init: bad return from post_init_slint: stopping'
          STOP
        ENDIF
        IF (.NOT. post_init_slint_called) THEN
          WRITE(6,*)'dyn_init: logic error: cannot test gribout until post_init_slint has been called'
          STOP
        ENDIF
!SMS$SERIAL END
      END IF ! readrestart: initialize slint

      DEALLOCATE(in_halo)

! --- initialize transects for meridional mass flux diagnostics
      call trnsec_init

! --- exercising stencl routine
      call stencl(deg_lat,1,1.,'latitude')
      call stencl(deg_lon,1,1.,'longitude')

! All the work_edg stuff need only be done in the initial run
      if (.not.readrestart) then
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! --- special diagnostics to figure out what goes on at the poles

!    if (PrintIpnDiag == 1 .or. PrintIpnDiag == nip) then
        allocate (work_edg(1,npp,nip))

        do idx=1,2
!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (isn)
          do ipn=ips,ipe
            do isn=1,nprox(ipn)
              work_edg(1,isn,ipn) = sidevec_e(idx,isn,ipn)
            end do
          end do
!$OMP END PARALLEL DO
!sms$ignore end

!SMS$EXCHANGE(work_edg)

          write (string,'(a,i1,a)') 'latitude, sidevec_e(',idx,')'
          call stenedg(deg_lat,work_edg,1,string)

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (isn)
          do ipn=ips,ipe
            do isn=1,nprox(ipn)
              work_edg(1,isn,ipn) = sidevec_c(idx,isn,ipn)
            end do
          end do
!$OMP END PARALLEL DO
!sms$ignore end

!SMS$EXCHANGE(work_edg)

          write (string,'(a,i1,a)') 'latitude, sidevec_c(',idx,')'
          call stenedg(deg_lat,work_edg,1,string)
        end do

        do idx=1,4
!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (isn)
          do ipn=ips,ipe
            do isn=1,nprox(ipn)
              work_edg(1,isn,ipn) = sn(idx,isn,ipn)
            end do
          end do
!$OMP END PARALLEL DO
!sms$ignore end

!SMS$EXCHANGE(work_edg)

          write (string,'(a,i1,a)') 'latitude, sn(',idx,')'
          call stenedg(deg_lat,work_edg,1,string)

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (isn)
          do ipn=ips,ipe
            do isn=1,nprox(ipn)
              work_edg(1,isn,ipn) = cs(idx,isn,ipn)
            end do
          end do
!$OMP END PARALLEL DO
!sms$ignore end

!SMS$EXCHANGE(work_edg)
          write (string,'(a,i1,a)') 'latitude, cs(',idx,')'
          call stenedg(deg_lat,work_edg,1,string)
        end do

        deallocate (work_edg)
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

        call copytoGPUinit(ims,ime,ips,ipe,ihe,nvl,nvlp1,npp,ntra,PrintIpnDiag,&
          nprox,proxs,prox,cs,sn,nedge,permedge,perm )

      endif ! (not.readrestart)

    endif ! (.not.iam_write_task)

! Set icosio values.

    icosio_comm = my_comm
    if (use_write_tasks) icosio_comm = output_intercomm

! If client/server IO is disabled, write tasks return from the icosio_setup
! call when they receive a shutdown command from the compute root, after the
! last output frame has been processed. Otherwise, they call mpi_finalize,
! followed by Fortran stop.

    client_server_io = .false.
    if (present(client_server_io_in)) client_server_io = client_server_io_in

! If this is a write task, inv_perm is unallocated here. Since Lahey doesn't
! allow null pointer arguments, allocate inv_perm at unit size.

    if (.not.allocated(inv_perm)) allocate(inv_perm(1))

! Tell icosio to use unit 50. Since getunit() with no arguments may return
! different numbers for different MPI tasks, want to be safe.

    unitno = getunit (50)
    if (unitno < 0) then
      print*,'dyn_init: getunit failed for icosio_setup. Stopping'
      stop
    end if

! Send icosio required values.

    if (iam_write_task) then
      call icosio_setup(                      &
        client_server_io_in=client_server_io, &
        debugmsg_on_in=debugmsg_on,           &
        comm_in=icosio_comm,                  &
        i_am_write_task_in=iam_write_task     &
        )
    else
      call icosio_setup(                      &
        binout_in=fimout,                     &
        client_server_io_in=client_server_io, &
        comm_in=icosio_comm,                  &
        debugmsg_on_in=debugmsg_on,           &
        gribout_in=gribout,                   &
        i_am_write_task_in=iam_write_task,    &
        ips_in=ips,                           &
        ipe_in=ipe,                           &
        inv_perm_in=inv_perm,                 &
        lun_in=unitno,                        &
        perm_in=perm,                         &
        permute_in=fixedgridorder,            &
        print_diags_in=printdiags,            &
        using_write_tasks_in=use_write_tasks  &
        )
    endif

    print *,'... exiting dyn_init'
    ret = gptlstop ('dyn_init')
    return

! TODO Replace this generic err= handler with iostat= handlers for each IO
!      operation so that specific, useful feedback is given.

90  write(6,*)'dyn_init: error reading a file'
    call flush(6)
    stop

  end subroutine dyn_init

end module module_fim_dyn_init
