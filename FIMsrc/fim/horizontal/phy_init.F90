module module_fim_phy_init

  implicit none

#include <gptl.inc>

  private
  public :: phy_init, sst_init, sstunit

  integer, save :: sstunit = -1

contains

!*********************************************************************
!       Loads the initial variables and constants for the physics 
!       component.  
!       Alexander E. MacDonald  11/27/04
!       J. Lee                  September, 2005
!*********************************************************************

  subroutine phy_init

    use funcphys ! GFS physics
    use module_sfc_variables

    USE gfs_physics_internal_state_mod, only: gfs_physics_internal_state, gis_phy
    USE gfs_physics_sfc_flx_set_mod,    only: sfcvar_aldata, flxvar_aldata, flx_init
    use gfsphys_nl,                     only: set_gfsphys_vars
    use headers,                        only: testcurveheader,testglvlheader
    use module_control,                 only: nip,dt,numphr,CallPhysics,CallRadiation, &
                                              dospptdt
    USE module_enkf_io,                 only: unitno_enkf
    use namelist_def,                   only: ras
    use resol_def,                      only: nfxr, num_p3d_resol_def=>num_p3d
    use units,                          only: getunit, returnunit

    use fimnamelist,only:  &
      aerosol_file,        &
      alt_topo,            &
      curve,               &
      enkfio_in,           &
      gfsltln_file,        &
      glvl,                &
      mtnvar_file,         &
      num_p3d,             &
      NumCacheBlocksPerPE, &
      nvl,                 &
      PhysicsInterval,     &
      RadiationInterval,   &
      readrestart,         &
      SSTInterval,         &
      UpdateSST,	   &
      ishal_cnv,	   &
      ideep_cnv

! Local variables

    integer       :: mype,idx,ipn,ierr,mylb,myub
    integer       :: unitno
    real*8        :: t0,t1=0.0d0
    logical       :: first_loop
!sms$distribute (dh,1) begin
    real :: work2d(nip)
!sms$distribute end

    integer :: ret
    real*8 :: tphysics_input

!sms$comm_rank(mype)
    ret = gptlstart ('phy_init')

    CallPhysics   = max(1,numphr*PhysicsInterval/3600)
    CallRadiation = max(1,(max(1,numphr*RadiationInterval/3600)/CallPhysics) &
      *CallPhysics)
    dospptdt = CallPhysics ! time step interval to apply SPPT.
 !   dospptdt = numphr ! time step interval to apply SPPT.

    print "(' Call physics every'  ,I27,' timesteps (',I0,' seconds)')", &
      CallPhysics,nint(CallPhysics*dt)
    print "(' Call radiation every',I25,' timesteps (',I0,' seconds)')", &
      CallRadiation,nint(CallRadiation*dt)
    print "(' num_p3d =',I10)",num_p3d

! Verify GFS physics namelist contents
    call set_gfsphys_vars ()

    print "(' ras =',L10)",ras

! On a restart run, these variables will be read in from the restart file so this part can be skipped.

    if (.not.readrestart) then
!SMS$PARALLEL(dh, ipn) BEGIN
      ret = gptlstart ('physics_input')
      unitno = getunit ()
      if (unitno < 0) then
        print*,'phy_init: getunit failed: stopping'
        stop
      end if

      open                (unitno, file="gfsfc.dat", form="unformatted", action='read', err=70)

! The following serial region is required because the called routines perform
! IO but are not processed by SMS, so that IO actions would otherwise occur on
! all tasks. TODO: These routines may call stop and, since they have not been
! processed by SMS, the stop may hang the run. They should either be processed
! by SMS, or simply return a status that the caller can act on.

!SMS$SERIAL (default=ignore) BEGIN
      call TestGlvlHeader (unitno,     "gfsfc.dat",'phy_init',glvl )
      call TestCurveHeader(unitno,     "gfsfc.dat",'phy_init',curve)
!SMS$SERIAL END

      DO idx = 1,SIZE(st3d,1)
        read(unitno, err=90) work2d
        do ipn=1,nip
          st3d(idx,ipn) = work2d(ipn)
        enddo
      enddo
! soil moisture

      do idx = 1,SIZE(sm3d,1)
        read(unitno, err=90) work2d
        do ipn=1,nip
          sm3d(idx,ipn) = work2d(ipn)
        enddo
      enddo
      do idx = 1,SIZE(slc3d,1)
        read(unitno, err=90) work2d
        do ipn=1,nip
          slc3d(idx,ipn) = work2d(ipn)
        enddo
      enddo
!
! skin temperature
!
      read(unitno, err=90) ts2d
!
! snow water equivalent
!
      read(unitno, err=90) sheleg2d
      read(unitno, err=90) tg32d
!
! surface roughness
!
      read(unitno, err=90) zorl2d
!
! maybe conv cloud fraction?
!
      read(unitno, err=90) cv2d
!
! maybe conv cloud bottom pressure?
!
      read(unitno, err=90) cvb2d
!
! maybe conv cloud top pressure?
!
      read(unitno, err=90) cvt2d
!
! mean visible albedo with strong cosz dependence...???...
!
      read(unitno, err=90) alvsf2d
!
! mean vis albedo with weak cosz dependence...???...
!
      read(unitno, err=90) alvwf2d
!
! mean nir albedo with strong cosz dependence...???...
!
      read(unitno, err=90) alnsf2d
!
! mean nir albedo with weak cosz dependence...???...
!
      read(unitno, err=90) alnwf2d
!
! land/sea/ice mask (0:SEA.1:LAND,2:ICE)
      read(unitno, err=90) slmsk2d
!
! assuming veg fraction
!
      read(unitno, err=90) vfrac2d
!
! canopy moisture content in mm
!
      read(unitno, err=90) canopy2d
      read(unitno, err=90) f10m2d
      read(unitno, err=90) t2m2d
      read(unitno, err=90) q2m2d
!
!   vegtype or landuse?
!
      read(unitno, err=90) vtype2d
!
! soilcategory
!
      read(unitno, err=90) stype2d
!
! fractional coverage with strong (facsf) and weak (facwf) cosz dependence
!
      read(unitno, err=90) facsf2d
      read(unitno, err=90) facwf2d
!
! ustar
!
      read(unitno, err=90) uustar2d
!
! looks like surface exchange coeffs for m and h
      read(unitno, err=90) ffmm2d
      read(unitno, err=90) ffhh2d
!
! ice fractions! fice is also used as something different (cloud ice fraction!!)
!
      read(unitno, err=90) hice2d
      read(unitno, err=90) fice2d
      read(unitno, err=90) tprcp2d
      read(unitno, err=90) srflag2d
      read(unitno, err=90) snwdph2d
      read(unitno, err=90) slc2d
      read(unitno, err=90) shdmin2d
      read(unitno, err=90) shdmax2d
      read(unitno, err=90) slope2d
      read(unitno, err=90) snoalb2d
!
! this thing contains such things as origraphuc stand deviation, convexity, asymetry,...for grav wave drag calc
!
      do idx = 1,SIZE(hprm2d,1)
        read(unitno, err=90) work2d
        do ipn=1,nip
          hprm2d(idx,ipn) = work2d(ipn)
        enddo
      enddo
      close(unitno)
      call returnunit (unitno)
      ret = gptlstop ('physics_input')
!SMS$PARALLEL END

      ret = gptlget_wallclock ('physics_input', 0, tphysics_input)  ! The "0" is thread number
      print"(' PHYSICS INPUT time:',F10.0)", tphysics_input
    end if     ! .not. readrestart

    call gfuncphys ()   ! GFS physics

!TODO:  encapsulate this in a new subroutine
! Allocate GFS internal state
! (only on first call, after digital filter it is
!  already allocated)
    if (.not. associated(gis_phy)) then
      allocate( gis_phy )
! set memory bounds from SMS distribution (decomposition)
!sms$ignore begin
      mylb = LBOUND(work2d,1)
      myub = UBOUND(work2d,1)
!sms$ignore end
      gis_phy%ims = mylb
      gis_phy%ime = myub
! set patch (distributed-memory loop) bounds from SMS distribution
!SMS$PARALLEL(dh, ipn) BEGIN
      first_loop = .true.
      do ipn=1,nip
        if (first_loop) then
          gis_phy%ips = ipn
          first_loop = .false.
        endif
        gis_phy%ipe = ipn
      enddo
!SMS$PARALLEL END
      gis_phy%lsoil = 4
!TODO:  fill this in...  
! This will declare surface and flux arrays as (mylb:myub,1), etc.
! Note:  GFS ignores ierr...  
      call sfcvar_aldata(mylb, myub, 1, gis_phy%lsoil, gis_phy%sfc_fld, ierr)
      call flxvar_aldata(mylb, myub, 1, gis_phy%flx_fld, ierr)
      gis_phy%NBLCK = 1
      gis_phy%LEVS = nvl
      gis_phy%num_p2d = 3
      gis_phy%num_p3d = num_p3d
      num_p3d_resol_def = num_p3d
      gis_phy%ras = ras
      gis_phy%NMTVR = 14
      gis_phy%lats_node_r = 1
!TODO:  move nfxr to module resol_def
!JR changing nfxr to 33 (per Bao) fixes bug found by JacquesM--old GFS
!JR physics this value was 14 but there are now hard-wired 33s elsewhere.
      nfxr = 33
      gis_phy%nfxr = nfxr
!JR Set the value in resol_def to match gis_phy. It's critically important
!JR taht these magic numbers only appear ONCE
      allocate( gis_phy%SLAG(mylb:myub,gis_phy%LATS_NODE_R),                   &
        gis_phy%SDEC(mylb:myub,gis_phy%LATS_NODE_R),                   &
        gis_phy%CDEC(mylb:myub,gis_phy%LATS_NODE_R),                   &
        gis_phy%SWH(mylb:myub,gis_phy%LEVS,gis_phy%NBLCK,              &
        gis_phy%LATS_NODE_R),                              &
        gis_phy%HLW(mylb:myub,gis_phy%LEVS,gis_phy%NBLCK,              &
        gis_phy%LATS_NODE_R),                                          &
! jbao added for gfs 2014
        gis_phy%SWHC(mylb:myub,gis_phy%LEVS,gis_phy%NBLCK,             &
        gis_phy%LATS_NODE_R),                                          &
! jbao added for gfs 2014
        gis_phy%HLWC(mylb:myub,gis_phy%LEVS,gis_phy%NBLCK,             &
        gis_phy%LATS_NODE_R),                                          &
! jbao added for gfs 2014
        gis_phy%SFALB(mylb:myub,gis_phy%LATS_NODE_R),                  &
        gis_phy%ACV(mylb:myub,gis_phy%LATS_NODE_R),                    &
        gis_phy%ACVT(mylb:myub,gis_phy%LATS_NODE_R),                   &
        gis_phy%ACVB(mylb:myub,gis_phy%LATS_NODE_R),                   &
        gis_phy%phy_f2d(mylb:myub,gis_phy%LATS_NODE_R,                 &
        gis_phy%num_p2d),                                              &
        gis_phy%XLON(mylb:myub,gis_phy%LATS_NODE_R),                   &
        gis_phy%XLAT(mylb:myub,gis_phy%LATS_NODE_R),                   &
        gis_phy%HPRIME(gis_phy%NMTVR,mylb:myub,gis_phy%LATS_NODE_R),   &
        gis_phy%phy_f3d(mylb:myub,gis_phy%LEVS,gis_phy%NBLCK,          &
        gis_phy%LATS_NODE_R,gis_phy%num_p3d),                          &
        gis_phy%COSZDG(mylb:myub,gis_phy%LATS_NODE_R),                 &
        gis_phy%CLDCOV(gis_phy%LEVS,mylb:myub,gis_phy%LATS_NODE_R),    &
        gis_phy%FLUXR(gis_phy%nfxr,mylb:myub,gis_phy%LATS_NODE_R) )

!JR Set sfalb to something. Otherwise the copy done in cpl_init_phy_to_dyn
! will copy random junk off the heap which may cause a floating point error!
!TODO: Initialize things properly. Are there other arrays allocated above
! that might encounter the same problem? Who knows.

      gis_phy%sfalb(:,:) = -1.e38

      allocate( gis_phy%ps(mylb:myub),               &
        gis_phy%dp(mylb:myub,gis_phy%LEVS),  &
        gis_phy%dpdt(mylb:myub,gis_phy%LEVS),&
        gis_phy%p(mylb:myub,gis_phy%LEVS),   &
        gis_phy%u(mylb:myub,gis_phy%LEVS),   &
        gis_phy%v(mylb:myub,gis_phy%LEVS),   &
        gis_phy%t(mylb:myub,gis_phy%LEVS),   &
        gis_phy%q(mylb:myub,gis_phy%LEVS),   &
        gis_phy%oz(mylb:myub,gis_phy%LEVS),  &
        gis_phy%cld(mylb:myub,gis_phy%LEVS) )

    endif
!TODO:  Connect this properly, at present values set here are overwritten 
!TODO:  in physics().  
    call flx_init(gis_phy%flx_fld, ierr)

    IF (enkfio_in) THEN

      CALL readarr32_r (st3d,gis_phy%lsoil,unitno_enkf)
      CALL readarr32_r (sm3d,gis_phy%lsoil,unitno_enkf)
      CALL readarr32_r (slc3d,gis_phy%lsoil,unitno_enkf)

      CLOSE(unitno_enkf)
      CALL returnunit (unitno_enkf)

    ENDIF


!SMS$PARALLEL(dh, ipn) BEGIN
    do ipn=1,nip
      gis_phy%sfc_fld%CV    (ipn,1) = cv2d    (ipn)
      gis_phy%sfc_fld%CVT   (ipn,1) = cvt2d   (ipn)
      gis_phy%sfc_fld%CVB   (ipn,1) = cvb2d   (ipn)
      gis_phy%sfc_fld%SLMSK (ipn,1) = slmsk2d (ipn)
      gis_phy%flx_fld%SFCDSW(ipn,1) = 0.0
      gis_phy%flx_fld%SFCDLW(ipn,1) = 0.0
      gis_phy%flx_fld%SFCUSW(ipn,1) = 0.0    !hli 09/2014
      gis_phy%flx_fld%SFCULW(ipn,1) = 0.0    !hli 09/2014
      gis_phy%flx_fld%TOPDSW(ipn,1) = 0.0    !hli 09/2014
      gis_phy%flx_fld%TOPUSW(ipn,1) = 0.0    !hli 09/2014
      gis_phy%flx_fld%TOPULW(ipn,1) = 0.0    !hli 09/2014
      gis_phy%sfc_fld%TSEA  (ipn,1) = ts2d    (ipn)
      gis_phy%sfc_fld%STC (:,ipn,1) = st3d  (:,ipn)
      gis_phy%sfc_fld%SHELEG(ipn,1) = sheleg2d(ipn)
      gis_phy%sfc_fld%ZORL  (ipn,1) = zorl2d  (ipn)
      gis_phy%sfc_fld%snoalb(ipn,1) = snoalb2d(ipn)
      gis_phy%HPRIME      (:,ipn,1) = hprm2d(:,ipn)
      gis_phy%sfc_fld%HICE  (ipn,1) = hice2d  (ipn)
      gis_phy%sfc_fld%FICE  (ipn,1) = fice2d  (ipn)
      gis_phy%sfc_fld%TPRCP (ipn,1) = tprcp2d (ipn)
      gis_phy%sfc_fld%SRFLAG(ipn,1) = srflag2d(ipn)
      gis_phy%sfc_fld%SLC (:,ipn,1) = slc3d (:,ipn)
      gis_phy%sfc_fld%SMC (:,ipn,1) = sm3d  (:,ipn)
      gis_phy%sfc_fld%SNWDPH(ipn,1) = snwdph2d(ipn)
      gis_phy%sfc_fld%SLOPE (ipn,1) = slope2d (ipn)
      gis_phy%sfc_fld%SHDMIN(ipn,1) = shdmin2d(ipn)
      gis_phy%sfc_fld%SHDMAX(ipn,1) = shdmax2d(ipn)
      gis_phy%sfc_fld%TG3   (ipn,1) = tg32d   (ipn)
      gis_phy%sfc_fld%VFRAC (ipn,1) = vfrac2d (ipn)
      gis_phy%sfc_fld%CANOPY(ipn,1) = canopy2d(ipn)
      gis_phy%sfc_fld%VTYPE (ipn,1) = vtype2d (ipn)
      gis_phy%sfc_fld%STYPE (ipn,1) = stype2d (ipn)
      gis_phy%sfc_fld%F10M  (ipn,1) = f10m2d  (ipn)
      gis_phy%sfc_fld%FFMM  (ipn,1) = ffmm2d  (ipn)
      gis_phy%sfc_fld%FFHH  (ipn,1) = ffhh2d  (ipn)
      gis_phy%sfc_fld%ALVSF (ipn,1) = alvsf2d (ipn)
      gis_phy%sfc_fld%ALNSF (ipn,1) = alnsf2d (ipn)
      gis_phy%sfc_fld%ALVWF (ipn,1) = alvwf2d (ipn)
      gis_phy%sfc_fld%ALNWF (ipn,1) = alnwf2d (ipn)
      gis_phy%sfc_fld%FACSF (ipn,1) = facsf2d (ipn)
      gis_phy%sfc_fld%FACWF (ipn,1) = facwf2d (ipn)
! Note that T2M and Q2M are overwritten before use in do_physics_one_step()
      gis_phy%sfc_fld%T2M   (ipn,1) = t2m2d   (ipn)
      gis_phy%sfc_fld%Q2M   (ipn,1) = q2m2d   (ipn)
      gis_phy%phy_f3d(ipn,:,1,1,:)  = 0.0
      gis_phy%phy_f2d(ipn,1,:)      = 0.0
      gis_phy%CLDCOV(:,ipn,1)       = 0.0
      gis_phy%flx_fld%HFLX(ipn,1)   = 0.0
      gis_phy%flx_fld%EVAP(ipn,1)   = 0.0
! Note that UUSTAR is overwritten before use here.  uustar2d is not used.  
      gis_phy%sfc_fld%UUSTAR(ipn,1) = 0.01
      gis_phy%oz(ipn,:) = 0.0
    enddo
!SMS$PARALLEL END
    if ( UpdateSST) then
      call sst_init
    endif

    if (ishal_cnv.lt.0 .or. ishal_cnv.gt.3) then
      print *,'ishal_cnv=',ishal_cnv
      stop 'ishal_cnv undefined'
    end if

    if (ideep_cnv.lt.0 .or. ideep_cnv.gt.2) then
      print *,'ideep_cnv=',ideep_cnv
      stop 'ideep_cnv undefined'
    end if

    print *,'... exiting phy_init'

    ret = gptlstop ('phy_init')

    return
70  write(6,*)'phy_init: error opening a file'
    stop
90  write(6,*)'phy_init: error reading a file'
    stop
  end subroutine phy_init

  subroutine sst_init
!*********************************************************************
!       Loads the initial variables and constants for the ocean
!       component.  
!       Alexander E. MacDonald  11/27/04
!       J. Lee                  September, 2005
!*********************************************************************

    use module_control, only: nip,dt,numphr,CallSST
    use fimnamelist   , only: glvl,curve,NumCacheBlocksPerPE, &
        PhysicsInterval,RadiationInterval,SSTInterval
    use units, only: getunit, returnunit

    implicit none

! Local variables

    integer :: mype,idx,ierr,mylb,myub
    real*8  :: t0,t1=0.0d0

! TODO:  Create new decomp "dhp" and use it for all physics declarations 
! TODO:  and loops.  Split modules that contain variables used by both too!  

!sms$comm_rank(mype)

    CallSST = max(1,numphr*SSTInterval/3600)
    print "(' Call SST every'  ,I27,' timesteps (',I0,' seconds)')", &
      CallSST,nint(CallSST*dt)
    call InitSST
    return
  end subroutine sst_init

  subroutine InitSST

    use headers,only: testcurveheader,testglvlheader
    use module_control ,only: nip,prev_date,next_date
    use fimnamelist    ,only: yyyymmddhhmm,glvl,curve,nvl,ptop
    use module_constants, only : lat,lon
    use module_sfc_variables, only : ts2d,sst_prev,sst_next,fice2d,fice2d_prev,fice2d_next,slmsk2d,hice2d
    use gfs_physics_internal_state_mod, only: gfs_physics_internal_state, gis_phy
    use units, only: getunit, returnunit
    USE slint, ONLY: bilinear_init, bl_int

    implicit none

! locals  
    real            , allocatable   :: sst_prev_ll(:,:)
    real            , allocatable   :: fice2d_prev_ll(:,:)
    real            , allocatable   :: sst_next_ll(:,:)
    real            , allocatable   :: fice2d_next_ll(:,:)
    integer :: YEAR,MONTH,DAY  
    integer :: IM_OC
    integer :: JM_OC,ipn

    integer :: MIDMON
    integer :: ID
    integer :: M1
    integer :: M2
    integer :: MIDM
    integer :: MIDP
    integer :: N
    integer :: I
    integer :: NRECS

    integer :: THIS_YEAR,MDATE
    integer :: THIS_MONTH
    integer :: time_header(4),SDATE1,SDATE2,DAYS(12)
    real    :: REAL_VAR
    real    :: DX
    real    :: DY,FAC
    data    DAYS /31,28,31,30,31,30,31,31,30,31,30,31/
    logical :: skip
    integer :: unitno
    integer :: ierr
! BEGIN
    READ(UNIT=yyyymmddhhmm(1:4), FMT='(I4)') YEAR
    READ(UNIT=yyyymmddhhmm(5:6), FMT='(I2)') MONTH
    READ(UNIT=yyyymmddhhmm(7:8), FMT='(I2)') DAY

    IF ( mod(YEAR,4).EQ.0) THEN
      DAYS(2)=29.0
    ENDIF
    IF ( YEAR.EQ.1900) THEN
      DAYS(2)=28.0
    ENDIF
!SMS$SERIAL (<sst_prev,sst_next,fice2d_prev,fice2d_next,OUT> : default=ignore) BEGIN 
    IM_OC=360.0
    JM_OC=180.0
    allocate(sst_prev_ll(IM_OC,JM_OC))
    allocate(fice2d_prev_ll(IM_OC,JM_OC))
    allocate(sst_next_ll(IM_OC,JM_OC))
    allocate(fice2d_next_ll(IM_OC,JM_OC))
    print*,'calling bilinear_init for SST'

    unitno = getunit ()
    if (unitno < 0) then
      print*,'initsst: getunit failed: stopping'
      stop
    end if

    open (unitno, file='glvl.dat', form='unformatted', iostat=ierr)
    if (ierr.ne.0) then
      write (*,'(a)') 'initsst: error opening glvl.dat for reading'
      stop
    endif
    call TestGlvlHeader (unitno, 'glvl.dat', 'sst_init', glvl)
    call TestCurveHeader(unitno, 'glvl.dat', 'sst_init', curve)
    CALL bilinear_init('ocean_bcs_ltln', IM_OC*JM_OC, unitno, nip)
    close(unitno)
    call returnunit (unitno)

    MDATE=YEAR*100+MONTH
    MIDMON = DAYS(MONTH)/2 + 1

    sstunit = getunit ()
    if (sstunit < 0) then
      print*,'initsst: getunit failed for sstunit: stopping'
      stop
    else
      print*,'initsst: sst unit is ', sstunit
    end if

    open (sstunit, file='sst_dat', form='unformatted', iostat=ierr)
    if (ierr.ne.0) then
      write (*,'(a)') 'initsst: error opening sst_dat for reading'
      stop
    endif
    skip=.TRUE.
    do while (skip)
      read (sstunit,iostat=ierr) time_header
      if (ierr.ne.0) then
        write (*,'(a)') 'initsst: error reading time_header'
        stop
      endif
      print*,MDATE,time_header
      read(sstunit,iostat=ierr) sst_prev_ll
      if (ierr.ne.0) then
        write (*,'(a)') 'initsst: error reading sst_prev_ll'
        stop
      endif
      read(sstunit,iostat=ierr) fice2d_prev_ll
      if (ierr.ne.0) then
        write (*,'(a)') 'initsst: error reading fice2d_prev_ll'
        stop
      endif
      SDATE1=time_header(1)*100+time_header(2)
      SDATE2=time_header(3)*100+time_header(4)
!     IF (DAY .LT.MIDMON.AND.MDATE.LE.SDATE2) THEN
      IF (MDATE.LE.SDATE2) THEN
!       SDATE2 needs to equal MDATE for sst_prev_ll
        skip=.FALSE.
      ENDIF
    end do
    IF (DAY .GE.MIDMON) THEN
      read(sstunit,iostat=ierr) time_header
      if (ierr.ne.0) then
        write (*,'(a)') 'initsst: error reading time_header'
        stop
      endif
      read(sstunit,iostat=ierr) sst_prev_ll
      if (ierr.ne.0) then
        write (*,'(a)') 'initsst: error reading sst_prev_ll'
        stop
      endif
      read(sstunit,iostat=ierr) fice2d_prev_ll
      if (ierr.ne.0) then
        write (*,'(a)') 'initsst: error reading fice2d_prev_ll'
        stop
      endif
    ENDIF
    prev_date=time_header
    read(sstunit,iostat=ierr) time_header 
    if (ierr.ne.0) then
      write (*,'(a)') 'initsst: error reading time_header'
      stop
    endif
    read(sstunit,iostat=ierr) sst_next_ll
    if (ierr.ne.0) then
      write (*,'(a)') 'initsst: error reading sst_next_ll'
      stop
    endif
    read(sstunit,iostat=ierr) fice2d_next_ll
    if (ierr.ne.0) then
      write (*,'(a)') 'initsst: error reading fice2d_next_ll'
      stop
    endif
    next_date=time_header

    print*,'FOUND SSTs, calling bl_int',prev_date,next_date

    CALL bl_int (sst_prev_ll(:,:), sst_prev)
    CALL bl_int (sst_next_ll(:,:), sst_next)
    CALL bl_int (fice2d_prev_ll(:,:), fice2d_prev)
    CALL bl_int (fice2d_next_ll(:,:), fice2d_next)

    IF (DAY < MIDMON) THEN

      M1   = MOD(MONTH+10,12) + 1
      M2   = MONTH
      MIDM = DAYS(M1)/2 + 1
      MIDP = DAYS(M1) + MIDMON
      ID = DAY + DAYS(M1)

    ELSE

      M2   = MOD(MONTH,12) + 1
      M1   = MONTH
      MIDM = MIDMON
      MIDP = DAYS(M2)/2 + 1 + DAYS(M1)
      ID = DAY

    ENDIF
!SMS$SERIAL END
!SMS$SERIAL (<FAC,out> : default=ignore) BEGIN
    FAC = (real(ID -   MIDM)*86400)/ &
      (real(MIDP - MIDM)*86400         )
!SMS$SERIAL END
! replace ts2d over ocean points
!SMS$PARALLEL(dh, ipn) BEGIN
    DO ipn=1,nip
!  need logic to keep ice's temperature the same, only update sst and ice fraction
! update ice fraction 1st
! if there is new ice, set it to -1.8 and hice=0.0
! if ice melts, then set ts2d to sst and hice=0.0
      IF (slmsk2d(ipn).NE.1) THEN
        fice2d(ipn)=fice2d_next(ipn)*(FAC)+fice2d_prev(ipn)*(1.0 - FAC)
        if (fice2d(ipn) .GT. 1) fice2d(ipn)=1.0
        if (fice2d(ipn) .LT. 0) fice2d(ipn)=0.0
      ENDIF
      IF (fice2d(ipn) .GE. 0.5 .AND. slmsk2d(ipn) .EQ. 0) THEN ! freeze open ocean
        slmsk2d(ipn)=2.0
        ts2d(ipn)=271.35
        hice2d(ipn)=0.0
      ENDIF
      IF (fice2d(ipn) .LT. 0.5 .AND. slmsk2d(ipn) .EQ. 2) THEN ! melt sea-ice
        slmsk2d(ipn)=0.0
        hice2d(ipn)=0.0
      ENDIF
      IF (slmsk2d(ipn).EQ.0) THEN
        ts2d(ipn)=sst_next(ipn)*(FAC)+sst_prev(ipn)*(1.0 - FAC)
      ENDIF

      gis_phy%sfc_fld%TSEA  (ipn,1) = ts2d    (ipn)
      gis_phy%sfc_fld%HICE  (ipn,1) = hice2d  (ipn)
      gis_phy%sfc_fld%FICE  (ipn,1) = fice2d  (ipn)
      gis_phy%sfc_fld%SLMSK (ipn,1) = slmsk2d (ipn)
    ENDDO
!SMS$PARALLEL END

    RETURN

  end subroutine InitSST

end module module_fim_phy_init
