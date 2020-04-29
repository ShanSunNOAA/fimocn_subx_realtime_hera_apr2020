module hycom_init
use fimnamelist     ,only: PrintIpnDiag,stencl_step
use module_control  ,only: glvl,nip,npp,readrestart,ArchvTimeUnit,	&
                           TotalTime,dt,hrs_in_month,ArchvIntvl
use hycom_diag,only: glob2d, glob3d, glbdia
use stencilprint
use stenedgprint
use findmaxmin1
use findmaxmin2
use edgmaxmin
use units, only: getunit, returnunit
use hycom_control,only: thruflsecs,bdcurrsecs,nsecs,iocnmx,		&
	subcycle_ocn,trcfrq,latdiw,numtr,modesplit,bclin_frq
implicit none

contains

  subroutine ocn_init
   use fimnamelist    ,only: kdm,itest,diag_intvl,test_start,ocnonly,	&
       atmonly,coupled,do_rnoff,sss_rstor,ocn_trnsecdir,inicondir,	&
       ann_core,do_radcor,do_pcpcor,janjic_ocn,UpdateSST,smagcf
   use hycom_variables,only: lm0,ld0,ls0,airtmn
   use hycom_alloc
   use hycom_thermf,   only: opnfor
   use units,          only: getunit, returnunit
   use hycom_constants,only: wet,odepth,land_spval,batrop
   use module_control, only: nts,numphr

   real,external :: its2time
   real          :: TotalRunTime, TotalRunTime_mo
   integer errcon,i,unitno,atm_step,atm_Diag
   character chrglvl*2,flnm*160
   namelist /thrufl/ thruflsecs		! online diag. of flow through passages
   namelist /bdcurr/ bdcurrsecs		! online diag. of bdry.current strength

   print *,'entering ocn_init ...'
   unitno = getunit ()

!   open (unitno,file='./FIMnamelist')
!   read (unitno,OCEANnamelist)
!   write (*,OCEANnamelist)
!   close (unitno)

   call hyc_alloc
! --- synchronize ocean model time step with atmospheric time step -dt-
   batrop=dt/subcycle_ocn	        ! barotropic time step (sec)

   if ((coupled.and.ocnonly) .or. (coupled.and.atmonly) .or.		&
     (ocnonly.and.atmonly)) then
     print *,'STOP: only one can be true: coupled=',coupled,		&
       ' atmonly=',atmonly, 'ocnonly=',ocnonly
     stop 'only one of the 3 can be true (coupled,ocnonly,atmonly)'
   elseif (.not.coupled .and. .not.atmonly .and. .not.ocnonly) then
     print *,'STOP: one has be true: coupled=',coupled,' atmonly=',	&
	atmonly,'ocnonly=',ocnonly
     stop 'one of the 3 has to be true (coupled,ocnonly,atmonly)'
   end if
! --- Calculate run time in months
   TotalRunTime = its2time(nts)

   if (ArchvTimeUnit.eq.'ts') then
     TotalRunTime_mo = TotalRunTime/hrs_in_month/numphr
   else if (ArchvTimeUnit.eq.'mi') then
     TotalRunTime_mo = TotalRunTime/hrs_in_month/60.
   else if (ArchvTimeUnit.eq.'hr') then
     TotalRunTime_mo = TotalRunTime/hrs_in_month
   else if (ArchvTimeUnit.eq.'dy') then
     TotalRunTime_mo = TotalRunTime/30.
   else if (ArchvTimeUnit.eq.'mo') then
     TotalRunTime_mo = TotalRunTime
   else
     write (*,'(a,a)') 'ERROR in ocn_init: unrecognized output time unit: ',ArchvTimeUnit
     stop
   endif

   print *,' Total run time (months) = ',TotalRunTime_mo
! --- If running more than one month, use UpdateSST
   if (TotalRunTime_mo.gt.1.and.atmonly.and. .not.UpdateSST) then
     print *,' atmonly=',atmonly,' UpdateSST=',UpdateSST
     stop 'STOP: UpdateSST=T is required for 1+ mon atmonly runs'
   end if
!
   atm_step=stencl_step
   stencl_step=1
   atm_Diag=PrintIpnDiag
!  PrintIpnDiag=itest

   if (atmonly) then
     print *,' atmonly run, exiting ocn_init'
     return
   end if

   if (len_trim(ocn_trnsecdir).gt.0) then
    thruflsecs(:)=' '
    bdcurrsecs(:)=' '
    write (chrglvl,'(a,i1)') 'g',glvl
    flnm=trim(ocn_trnsecdir)//'transects_'//chrglvl//'.nml'
!SMS$SERIAL BEGIN
    print *,'get ocn.transects for online transport diagnostics from ',	&
      trim(flnm)
    open (unitno,file=flnm,iostat=errcon)
    if (errcon.eq.0) then
     read (unitno,thrufl)
     write (*,thrufl)
     read (unitno,bdcurr)
     write (*,bdcurr)
     close (unitno)
    else
     print *,'no ocn.transects found for online transport diagnostics'
    end if
!SMS$SERIAL END
   end if                       	! ocn_trnsecdir exists

   call geopar			! basin configuration

   if (ocnonly) then
     call opnfor
   else
!    if (.not. ann_core) call rdforf(0.)
   end if			! ocnonly

! moved call inicon after call opnfor, where airtmn is defined
   call inicon			! ocean state initialization

   stencl_step=atm_step
   PrintIpnDiag=atm_Diag
   call returnunit (unitno)
   print *,'... exiting ocn_init'
   return
  end subroutine ocn_init


!*********************************************************************
!   geopar
!     configures the ocean basin and sets misc.geographic parameters
!     R. Bleck				January 2010
!*********************************************************************
  subroutine geopar

  use module_constants,only: nprox,proxs,prox,nedge,permedge,		&
                             deg_lat,deg_lon,area,perm,inv_perm
  use fimnamelist     ,only: kdm,itest,ocn_trnsecdir,diag_intvl,	&
			     inicondir,sss_rstor,ocnonly,atmonly,	&
  			     coupled,do_radcor,do_pcpcor
  use hycom_constants ,only: wet,meshsz,odepth,dpth_edg,assluv,assldp,	&
			     dpmin,land_spval
  use hycom_variables ,only: lgthset,cellset,edgeset,sense,		&
                             ctrlat,ctrlon,transp,ocnarea
  use hycom_embroid
! use mdul_findmean1

!SMS$DISTRIBUTE (dh,1) BEGIN
  real   work (nip)
!SMS$DISTRIBUTE END
  real*4 topo4(nip)
  integer,allocatable :: tmp(:)
  integer i,ipx,edg,edgcount,n,maxlen,len,unitno
  real nfill,npoints
! real,parameter :: maxwest=60.		! latitude of maximum westerlies
  character chrglvl*2,flnm*160
!
  unitno = getunit ()
  if (unitno < 0) then
    print*,'geopar: getunit failed: stopping'
    stop
  end if
!
!SMS$SERIAL BEGIN
  open(unitno, file='ocn_dpmin.txt', form='formatted', status='old')
  read(unitno,*) dpmin
  close(unitno)
!SMS$SERIAL END
!
  print *,'entering hyc_geopar'
  print 102,'icos grid refinement level'          ,glvl
  print 102,'number of ocean layers'              ,kdm
  print 102,'ocean test point (glob)'             ,itest
  print 102,'KPP: full [1] or partial [2] column' ,iocnmx
  print 102,'ocean steps per atmo.step'           ,subcycle_ocn
  print 102,'time steps between tracer tranport'  ,trcfrq
  print 102,'freq. of diagnostic output'          ,diag_intvl
  print 103,'coupled atm & ocean'		  ,coupled
  print 103,'atm only'				  ,atmonly
  print 103,'ocean only'			  ,ocnonly
  print 103,'restore surface salinity'            ,sss_rstor
  print 103,'do radiation correction'             ,do_radcor
  print 103,'do precipitation correction'         ,do_pcpcor
  print 103,'use lat.-dep. intrnl.wave intensity' ,latdiw
  print 104,'time smoothing weights for velocity' ,assluv,1.-2.*assluv,assluv
  print 104,'time smoothing weights for thickness',assldp,1.-2.*assldp,assldp
102 format (a40,' =',i8)
103 format (a40,' =',l8)
104 format (a40,' =',3f8.4)
  print '(a/(5(i8,f6.1)))','minimum layer thickness:',(n,dpmin(n),n=1,kdm)

!..........................................................................
! construct conversion tables for  global <-> local grid indexing
!..........................................................................

! --- change test point from global to local

!SMS$SERIAL (<deg_lat,deg_lon,inv_perm,IN>) BEGIN
   if (itest.gt.nip) stop '(hycom_init: test point out of range)'
   if (itest.gt.0) then
    print '(2(a,i8),a,2f8.2)','test point itest=',itest,		&
      '(glob) becomes', inv_perm(itest),'(locl), lat/lon=',		&
      deg_lat(inv_perm(itest)),deg_lon(inv_perm(itest))
    itest=inv_perm(itest)
    PrintIpnDiag=itest
   end if
!SMS$SERIAL END

!!!!  if (.not. readrestart) then
! odepth read in here will be overwritten by odepth in the restart file
!SMS$SERIAL (<perm,IN>,<odepth,OUT> : default=ignore) BEGIN
!..........................................................................
! read in bathymetry file from "inicondir"/input_gX
!..........................................................................

   write (chrglvl,'(a,i1)') 'g',glvl
   flnm=trim(inicondir)//'input_'//chrglvl//'/topo1_'//chrglvl//'.4bin'
   print *,'get bathymetry from ',trim(flnm)
   open (unitno,file=flnm,form='unformatted',action='read')
   read (unitno) n,topo4
   close (unitno)
   if (n.ne.nip .and. n.ne.glvl) then
    print '(2(a,i7))','bathymetry file valid for',n,' processors, not',nip
    stop '(wrong bathymetry file)'
   end if
   odepth(:)=topo4(perm(:))
!SMS$SERIAL END

!SMS$PARALLEL(dh,i) BEGIN
!$OMP PARALLEL DO
   do i=1,nip
    if (odepth(i).gt.0.) then
!     odepth(i)=5000.				! flat bottom
      odepth(i)=max(odepth(i),100.)		! set minimum depth
      wet(i)=1
    else
      odepth(i)=0.
      wet(i)=0
    end if
   end do
!$OMP END PARALLEL DO
!SMS$EXCHANGE(odepth,wet)
!SMS$PARALLEL END

   work(:)=1.
   call glob2d(work,ocnarea)
   call stencl(odepth,1,1.,'(ocn geopar) raw bottom depth (m)')
   call stencl(deg_lat,1,1.,'(ocn geopar) latitude')
   call stencl(deg_lon,1,1.,'(ocn geopar) longitude')

! --- remove one-cell wide inlets (FIM needs 2 cells to
! --- build up pressure head to balance cross-channel wind stress)

  if (.not. readrestart) then
1  nfill=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(n,ipx) REDUCTION(+:nfill)
   do i=1,nip
    if (wet(i) > 0 ) then
     n=0
     do edg=1,nprox(i)
      ipx=prox(edg,i)
      if (wet(ipx) > 0 ) n=n+1
     end do
     if (n.le.1) then
      nfill=nfill+1.
      print *,'turn ocean cell',perm(i),' into dry land'
      odepth(i)=0.
      wet(i)=0
     end if
    end if
   end do
!$OMP END PARALLEL DO
!SMS$EXCHANGE (wet,odepth)
!SMS$PARALLEL END
!SMS$REDUCE (nfill,SUM)
   if (nfill.gt.0) go to 1

   call stencl(odepth,1,1.,'(ocn geopar) bottom depth (m) prior to embroid')
!  call edgmxmn(dpth_edg,1,'dpth_edg')
!  call stenedg(odepth,dpth_edg,1,'(ocn geopar) bottom depth, cell & edge')

! --- reconcile coastlines between atmospheric and oceanic grid (ocean trumps)
    if (coupled) call embroid

  npoints=0.
!SMS$PARALLEL(dh,i) BEGIN
!$OMP PARALLEL DO REDUCTION(+:npoints)
  do i=1,nip
   if (i.eq.itest) then
    print 100,'depths surrounding i=',perm(i),i,odepth(i)
100 format (/a,i8,' (lcl',i8,')',f8.1)
    do edg=1,nedge(i)
     ipx=prox(permedge(edg,i),i)
     print 105,'prox:',perm(ipx),ipx,odepth(ipx)
105  format (16x,a,i8,' (lcl',i8,')',f8.1)
    end do
    call flush(6)
   end if

   if (wet(i) > 0 ) then
    npoints=npoints+1.
   end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$REDUCE(npoints,SUM)

! --- define surface properties of newly added land points
  write(*,'(a,es11.3,i8)') 'ocean area and number of grid points=',	&
    ocnarea,int(npoints)

  end if    ! if .not. readrestart

! call findmean1(area,nip,'ocn cell size (m^2)',wet)
! call findmean1(area,nip,'atm call size (m^2)')

! --- define water depth on cell edges

!SMS$PARALLEL(dh,i) BEGIN
!$OMP PARALLEL DO
  do i=1,nip					! horizontal loop
   meshsz(i)=sqrt(area(i))			! scale for horiz.mesh size
   dpth_edg(:,i)=land_spval
   do edg=1,nedge(i)				! loop through edges
    dpth_edg(edg,i)=min(odepth(i),odepth(prox(edg,i)))
   end do
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! --- read in edge chains used in mass flux diagnostics

! --- the model is looking for a listing of transects in the namelist
! ---                ocn_trnsecdir//transects_g<glvl>.nml
! --- The namelist contains 2 blocks:
! --- &thrufl lists transects used for diagnosing flow through passages;
! --- &bdcurr lists transects used for western boundary current diagnostics.

! --- Flow through passages is output layer by layer. Thus, vertical-
! --- meridional overturning can be diagnosed by repeating in &thrufl any
! --- zonal cross-basin transects developed for &bdcurr.

  if (len_trim(ocn_trnsecdir).ne.0) then
   print *,'initializing online mass transport diagnostics...'
   write (chrglvl,'(a,i1)') 'G',glvl
   lgthset(:,:)=0

! --- determine space requirements for storing chain info
   do n=1,nsecs			! look for thruflow ("class 1") transects
    flnm=trim(ocn_trnsecdir)//thruflsecs(n)
    if (thruflsecs(n)(1:1).ne.' ') then
     print *,'checking ',trim(flnm),' for length of transect'
     if (index(thruflsecs(n),chrglvl).eq.0) then
      print *,'cannot find ',chrglvl,' in file name ',		&
        trim(thruflsecs(n)),' ... terminating'
      stop
     end if
!
     open (unitno,file=flnm,action='read')
     read (unitno,'(i5)') lgthset(n,1)
     close (unitno)
    end if
   end do			! class 1

   do n=1,nsecs			! look for bdry current ("class 2") transects
    flnm=trim(ocn_trnsecdir)//bdcurrsecs(n)
    if (bdcurrsecs(n)(1:1).ne.' ') then
     print *,'checking ',trim(flnm),' for length of transect'
     if (index(flnm,chrglvl).eq.0) then
      print *,'cannot find ',chrglvl,' in file name ',		&
        trim(bdcurrsecs(n)),' ... terminating'
      stop
     end if
!
     open (unitno,file=flnm,action='read')
     read (unitno,'(i5)') lgthset(n,2)
     close (unitno)
    end if
   end do			! class 2

   maxlen=maxval(lgthset(:,:))
   print *,'max length of transects:',maxlen
   allocate (edgeset(maxlen,nsecs,2),cellset(maxlen,nsecs,2),		&
             sense(maxlen,nsecs,2),ctrlat(maxlen,nsecs,2),		&
             ctrlon(maxlen,nsecs,2),transp(kdm,maxlen,nsecs,2),		&
             tmp(maxlen))
   edgeset(:,:,:)=0
   cellset(:,:,:)=0
   ctrlat (:,:,:)=0.
   ctrlon (:,:,:)=0.
   transp (:,:,:,:)=0.
   sense  (:,:,:)=' '

! --- read the chosen transects (both classes)

   do n=1,nsecs				! loop over class 1 transects
    if (thruflsecs(n)(1:1).ne.' ') then
     flnm=trim(ocn_trnsecdir)//thruflsecs(n)
     print '(/2a)','read cell/edge info for thruflow transect ',	&
      trim(flnm)
     open (unitno,file=flnm,action='read')
!SMS$SERIAL BEGIN
     read (unitno,101) (cellset(len,n,1),edgeset(len,n,1),		&
       sense(len,n,1),ctrlat(len,n,1),ctrlon(len,n,1),len=1,lgthset(n,1))
!SMS$SERIAL END
101  format (/2(i10,i2,1x,a6,2f8.3))
     close (unitno)

     print '(/2a)','(ocn_geopar) class-1 transect (GLOBAL indx) for ',	&
       trim(thruflsecs(n))
     print 101,(cellset(len,n,1),edgeset(len,n,1),sense(len,n,1),	&
      ctrlat(len,n,1),ctrlon(len,n,1),len=1,lgthset(n,1))

!SMS$SERIAL (<inv_perm,IN>) BEGIN
     do len=1,lgthset(n,1)
      tmp(len)=inv_perm(cellset(len,n,1))
!     print '(2(a,i8))','map',cellset(len,n,1),' onto local cell',	&
!       tmp(len)
     end do
     do len=1,lgthset(n,1)
      cellset(len,n,1)=tmp(len)
     end do
!SMS$SERIAL END

     print '(/2a)','(ocn_geopar) class-1 transect (LOCAL indx) for ',	&
       trim(thruflsecs(n))
     print 101,(cellset(len,n,1),edgeset(len,n,1),sense(len,n,1),	&
      ctrlat(len,n,1),ctrlon(len,n,1),len=1,lgthset(n,1))
    end if
   end do			! class 1
   call flush(6)

   do n=1,nsecs			! loop over class 2 transects
    if (bdcurrsecs(n)(1:1).ne.' ') then
     flnm=trim(ocn_trnsecdir)//bdcurrsecs(n)
     print '(/2a)','read cell/edge info for bdry current transect ',	&
      trim(flnm)
     open (unitno,file=flnm,action='read')
!SMS$SERIAL BEGIN
     read (unitno,101) (cellset(len,n,2),edgeset(len,n,2),		&
      sense(len,n,2),ctrlat(len,n,2),ctrlon(len,n,2),len=1,lgthset(n,2))
!SMS$SERIAL END
     close (unitno)

     print '(/2a)','(ocn_geopar) class-2 transect (GLOBAL indx) for ',	&
       trim(bdcurrsecs(n))
     print 101,(cellset(len,n,2),edgeset(len,n,2),sense(len,n,2),	&
      ctrlat(len,n,2),ctrlon(len,n,2),len=1,lgthset(n,2))

!SMS$SERIAL (<inv_perm,IN>) BEGIN
     do len=1,lgthset(n,2)
      tmp(len)=inv_perm(cellset(len,n,2))
!     print '(2(a,i8))','map',cellset(len,n,2),' onto local cell',	&
!       tmp(len)
     end do
     do len=1,lgthset(n,2)
      cellset(len,n,2)=tmp(len)
     end do
!SMS$SERIAL END

     print '(/2a)','(ocn_geopar) class-2 transect (LOCAL indx) for ',	&
       trim(bdcurrsecs(n))
     print 101,(cellset(len,n,2),edgeset(len,n,2),sense(len,n,2),	&
      ctrlat(len,n,2),ctrlon(len,n,2),len=1,lgthset(n,2))
    end if
   end do			! class 2
   call flush(6)
  end if			! ocn_trnsecdir exists


  call returnunit (unitno)
  print *,'... exiting hyc_geopar'
  return
  end subroutine geopar


!*********************************************************************
!    inicon 
!       initializes the ocean state
!       R. Bleck                  		January 2010
!       S. Sun (observed init.cond.)		April   2011
!*********************************************************************
  subroutine inicon
  use module_constants,only: nprox,proxs,prox,nedge,permedge,grvity,	&
                             deg_lat,deg_lon,area,perm
  use module_variables,only: ph3d
  use hycom_constants ,only: wet,odepth,dpth_edg,batrop,thref,onem,	&
                             ArchvStepOcn,TotalStepOcn,land_spval,pref,	&
			     theta,salmin,temmin,nuphill,uphill,dnhill,	&
			     thkmx,tfrez,tmin0,tmax0,smin0,smax0
  use hycom_variables ,only: utrop,vtrop,ptrop,uclin,vclin,dp,pres,	&
                             dpinit,geop,montg,dens,spcifv,passv_tr,	&
                             temp,saln,pbot,gzbot,spvkk,psikk,		&
                             um_edg,vm_edg,un_edg,vn_edg,dp_edg,	&
                             pr_edg,geop_edg,mont_edg,uvsq_edg,		&
                             ubforc,vbforc,curl,massflx,cumuflx,	&
                             btropfx,nstep_ocn,mssfx_ave,curl_ave,	&
                             ssht_ave,srfx_ave,pmne_ave,qf2d_ave,	&
                             hf2d_ave,airt_ave,vpmx_ave,prcp_ave,	&
                             swdn_ave,lwdn_ave,rad_ave,wspd_ave,	&
			     taux_ave,tauy_ave,temice_ave,covice_ave,	&
                             thkice_ave,temp_ave,saln_ave,dens_ave,	&
			     uvel_ave,vvel_ave,dp_ave,			&
                             srfx_bcl,pmne_bcl,ustar_bcl,ustarb_bcl,	&
                             covice,temice,ticeol,thkice,thksno,	&
                             taux,tauy,ustar,ustarb,pmne,srfx,rivflo,	&
                             sw2d_bcl,hmixl,airtmn,hmixl_ave,		&
                             srfxcum,pmnecum,totmass,ocnarea,ssht,	&
                             sstndcy
  use fimnamelist     ,only: kdm,itest,diag_intvl,inicondir,ocnonly,	&
			     atmonly,coupled,do_rnoff,ocn_ic,ocntst
  use hycom_transp3d  ,only: transp0
  use module_diagnoise,only: linout
  use hycom_sigetc
  use hycom_mxkpp
  use hycom_hybgn1
  use hycom_runoff
  use hycom_sigetc

  integer i,k,n,edg,leap,ntr,unitno
  character chrglvl*2,string*30,flnm*160
  logical event
  real q
  real*8 t8,s8,r8,p8,chk_rhor,chk_rhostar,chk_rho
  real, parameter :: thresh=-40.	! Celsus

!SMS$DISTRIBUTE (dh,1) BEGIN
  real :: srfht(nip),work(nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
  real :: pout (kdm+1,nip)
  real :: work1(kdm,nip),work2(kdm,nip)
!SMS$DISTRIBUTE END
! serial for reading data:
  real*4 :: temp4(kdm,nip),saln4(kdm,nip),pout4(kdm+1,nip)
  real*4 :: refp(kdm+1),refmnt(kdm)
  real,parameter :: thkmxnor  = 3.0	! northern hem. max. ice thickness (m)
  real,parameter :: thkmxsou  = 1.0	! southern hem. max. ice thickness (m)
  real,parameter :: tfreznor  =-1.9	! freezing temp,arctic (deg)
  real,parameter :: tfrezsou  =-1.6	! freezing temp,antarctic (deg)

  refmnt(:)=0.
  
! --- specify initial layer structure
!
  unitno = getunit ()
  if (unitno < 0) then
    print*,'inicon: getunit failed: stopping'
    stop
  end if
!
!SMS$SERIAL BEGIN
  open(unitno, file='ocn_target.txt', form='formatted', status='old')
  read(unitno,*) (theta(n),n=1,kdm)
  close(unitno)
!
  t8=temmin
  do k=1,kdm
    r8=theta(k)
    salmin(k)=sofsig(r8,t8)
  end do
!SMS$SERIAL END
!
! initialize arrays which are -- NOT -- in the restart file

!SMS$PARALLEL(dh,i) BEGIN
!SMS$HALO_COMP(<1,1>) BEGIN
!$OMP PARALLEL DO PRIVATE(q)
  do i=1,nip
   q=land_spval
   if (wet(i) > 0 ) q=0.

   do leap=1,3
    utrop(leap,i)=q
    vtrop(leap,i)=q
    ptrop(leap,i)=q
   end do

   srfx_bcl     (i)=q
   sw2d_bcl     (i)=q
   pmne_bcl     (i)=q
   ustar_bcl    (i)=q
   ustarb_bcl   (i)=q

   ticeol       (i)=q
   taux         (i)=q
   tauy         (i)=q
   ustar        (i)=q
   ustarb       (i)=q
   pmne         (i)=q
   srfx         (i)=q
   ssht         (i)=q
   ubforc       (i)=q
   vbforc       (i)=q
   curl		(i)=q

   airt_ave     (i)=q
   vpmx_ave     (i)=q
   prcp_ave     (i)=q
   swdn_ave     (i)=q
   lwdn_ave     (i)=q
    rad_ave     (i)=q
   wspd_ave     (i)=q
   taux_ave     (i)=q
   tauy_ave     (i)=q
   mssfx_ave(:,:,i)=q
   curl_ave     (i)=q
   ssht_ave     (i)=q
   srfx_ave     (i)=q
   pmne_ave     (i)=q
   qf2d_ave     (i)=q
   hf2d_ave     (i)=q
   hmixl_ave    (i)=q
   temice_ave   (i)=q
   covice_ave   (i)=q
   thkice_ave   (i)=q
   montg      (:,i)=q
   temp_ave   (:,i)=q
   saln_ave   (:,i)=q
   dens_ave   (:,i)=q
   uvel_ave   (:,i)=q
   vvel_ave   (:,i)=q
     dp_ave   (:,i)=q
   sstndcy    (:,i)=q

   if (deg_lat(i).gt.0.) then                   ! northern hemi.
     thkmx(i)=thkmxnor
     tfrez(i)=tfreznor
   else                                         ! southern hemi.
     thkmx(i)=thkmxsou
     tfrez(i)=tfrezsou
   end if

  end do			! horiz. loop
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
  do i=1,nip
   do edg=1,npp
    pr_edg (1,edg,i)=0.
    btropfx  (edg,i)=0.
   end do

   do k=1,kdm
     spcifv(k,i)=land_spval
    do edg=1,npp
     um_edg   (k,edg,i)=land_spval
     vm_edg   (k,edg,i)=land_spval
     un_edg   (k,edg,i)=land_spval
     vn_edg   (k,edg,i)=land_spval
     dp_edg   (k,edg,i)=land_spval
     pr_edg (k+1,edg,i)=land_spval
     geop_edg (k,edg,i)=land_spval
     mont_edg (k,edg,i)=land_spval
     uvsq_edg (k,edg,i)=land_spval
     massflx  (k,edg,i)=land_spval
    end do
   end do			! vertical loop
  end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$HALO_COMP END
!SMS$PARALLEL END

! arrays in restart files are initialized only at time=0

  if (.not. readrestart) then		! start from scratch
   print *,'entering hyc_inicon ... ocn_ic (1: Levitus; 2: CFSR) =',ocn_ic
   nstep_ocn=0

!SMS$PARALLEL(dh,i) BEGIN
!SMS$HALO_COMP(<1,1>) BEGIN
!$OMP PARALLEL DO PRIVATE(q)
   do i=1,nip
    q=land_spval
    if (wet(i) > 0 ) q=0.

    covice       (i)=q
    temice       (i)=q
    thkice       (i)=q
    thksno       (i)=q
    pbot         (i)=q
    gzbot        (i)=q
    spvkk        (i)=q
    psikk        (i)=q
    hmixl        (i)=q
    srfxcum      (i)=q
    pmnecum      (i)=q

    pres (1,i)=q
    geop (1,i)=q

    do k=1,kdm
     temp(k,i)=q
     saln(k,i)=q
     do ntr=1,numtr
      passv_tr(k,i,ntr)=q
     end do
     do leap=1,2
      uclin(k,i,leap)=q
      vclin(k,i,leap)=q
      dp   (k,i,leap)=0.		! set dp=0 on land
     end do
     pres(k+1,i)=q
     geop(k+1,i)=q
     montg (k,i)=q
     dens  (k,i)=q
     dpinit(k,i)=q
    end do			! vertical loop

   end do			! horiz. loop
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
   do i=1,nip
    do k=1,kdm
     do edg=1,npp
      cumuflx(k,edg,i)=land_spval
     end do
    end do			! vertical loop
   end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$HALO_COMP END
!SMS$PARALLEL END

!SMS$SERIAL BEGIN
   t8=5.
   s8=36.
   p8=2500.d4			! 2500m
   chk_rho      =39.6858
   if (pref.eq.1.e7) then	! check values for t8,s8,p8
     chk_rhor   =33.0321
     chk_rhostar=32.9379
   elseif (pref.eq.2.e7) then
     chk_rhor   =37.4934
     chk_rhostar=37.4617
   end if
   if (abs(sigloc(t8,s8,p8) - chk_rho).gt.0.001)			&
     print *,'sigloc should be',chk_rho,' now =',sigloc(t8,s8,p8)
   if (abs(sigocn(t8,s8) - chk_rhor).gt.0.001)				&
     print *,'sigocn should be',chk_rhor,' now =',sigocn(t8,s8)
   if (abs(sigma_pref_star(t8,s8,p8) - chk_rhostar).gt.0.001) 		&
     print *,'sigma_pref_star should be',chk_rhostar,' now =',		&
     sigma_pref_star(t8,s8,p8)

   write(*,'(a,es7.1,a,3f10.6)') 'pref= ',pref,' chk rho/rhor/rho*=',	&
     sigloc(t8,s8,p8),sigocn(t8,s8),sigma_pref_star(t8,s8,p8)
   print 101,'target densities:   ',theta
   print 101,'minimum salinities: ',salmin
101 format (a,8f7.3/(20x,8f7.3))
   write(chrglvl,'(a,i1)') 'g',glvl
!SMS$SERIAL END

   if (ocn_ic == 2 ) then 	! CFSR
     call read_ocndata_cfsr
   else if (ocn_ic == 1) then	! Levitus Climatology
!SMS$SERIAL BEGIN
     if (ocntst) then
       write(string,'(a,i1,a,i2,a,i1,a)')'g',glvl,'_k',kdm,'_sig',	&
         int(pref*1.e-7),'b.4bin'
     else
       if (kdm == 26) then
         if (theta(1) == 24.35) then
           write(string,'(a,i1,a,i2,a,i1,a)')'g',glvl,'x',kdm,'dec_sig',	&
             int(pref*1.e-7),'g.4bin'
         elseif (theta(1) == 25.32) then
           write(string,'(a,i1,a,i2,a,i1,a)')'g',glvl,'x',kdm,'dec_sig',	&
             int(pref*1.e-7),'b.4bin'
         end if 
       elseif (kdm == 32) then
         if (theta(1) == 20.82) then
           write(string,'(a,i1,a,i2,a,i1,a)')'g',glvl,'x',kdm,'dec_sig',	&
             int(pref*1.e-7),'a.4bin'
         elseif (theta(1) == 19.52) then
           write(string,'(a,i1,a,i2,a,i1,a)')'g',glvl,'x',kdm,'dec_sig',	&
             int(pref*1.e-7),'b.4bin'
         end if 
       else
         stop 'wrong theta1'
       end if
     end if
!
     flnm=trim(inicondir)//'input_'//chrglvl//'/temp_'//trim(string)
     print *,'get initial temperature from ',trim(flnm)
     open(unitno,file=flnm,form='unformatted',action='read')
     read(unitno) i,temp4	! -9999 on land
     close(unitno)
  
     flnm=trim(inicondir)//'input_'//chrglvl//'/saln_'//trim(string)
     print *,'get initial salinity from    ',trim(flnm)
     open(unitno,file=flnm,form='unformatted',action='read')
     read(unitno) i,saln4	! -9999 on land
     close(unitno)
  
     flnm=trim(inicondir)//'input_'//chrglvl//'/pout_'//trim(string)
     print *,'get initial intfc.prs. from  ',trim(flnm)
     open(unitno,file=flnm,form='unformatted',action='read')
     read(unitno) i,pout4
     close(unitno)
  
     if (i.ne.nip) then
      print '(2(a,i7))','temp file valid for',i,' processors, not',nip
      stop '(wrong temp file)'
     end if

! --- map initial fields onto local index space ("curve")
     do i=1,nip
      if (temp4(1,perm(i)).gt.-99.) then
       temp(:,i)=temp4(:,perm(i))
       saln(:,i)=saln4(:,perm(i))
       pout(:,i)=pout4(:,perm(i))
!
       pres(1,i)=0.
       do k=1,kdm
         saln(k,i)=max(saln(k,i),salmin(k))
         pres(k+1,i)=min(odepth(i),pout(k+1,i))*onem
         if (k.eq.kdm) pres(k+1,i)=odepth(i)*onem
       end do
      end if
     end do

! --- initialize ice-related fields temice/covice/thkice
     do i=1,nip
     if (wet(i) > 0 ) then
       temice(i)=temp(1,i)
       if (temice(i).le.-1.0 .or.	&
        (ocnonly.and.airtmn(i).lt.thresh)) then	! airtmn in Celsus
         covice(i)=1.	! initial ice coverage
         thkice(i)=thkmxnor	! initial ice thickness in NH
         if (deg_lat(i).lt.0.) thkice(i)=thkmxsou	! initial ice thickness in SH
         temice(i)=-1.8 
! --- suppress instant melting by lowering temp in top layers
!ss      temp(1,i)=min(-1.8,temp(1,i))
!ss      temp(2,i)=min(-1.7,temp(2,i))
!ss      temp(3,i)=min(-1.6,temp(3,i))
!ss      temp(4,i)=min(-1.5,temp(4,i))
       else
         covice(i)=0.
         thkice(i)=0.
       end if		! temice <-0.5 .or. ...
     end if		! wet
     end do		! i loop
!SMS$SERIAL END
!
   end if	! ocn_ic=1 or 2: temp/saln/pres/temice/covice/thkice are done & curved right
   
! -- checking temp & saln ranges
!SMS$PARALLEL(dh,i) BEGIN
!$OMP PARALLEL DO
   do i=1,nip
     if (wet(i) > 0 ) then
       do k=1,kdm
         if (abs(temp(k,i)-.5*(tmax0+tmin0)) .ge. .5*(tmax0-tmin0)) then
	   write(*,'(a,es9.2,a,i8,i3,a,2f7.1)') 'chk bad ini temp ',	&
            temp(k,i),' at i=',perm(i),k,				&
		' lat/lon=',deg_lat(i),deg_lon(i)
           temp(k,i)=max(min(temp(k,i),tmax0),tmin0)
         end if
         if (abs(saln(k,i)-.5*(smax0+smin0)) .ge. .5*(smax0-smin0)) then
	   write(*,'(a,es9.2,a,i8,i3,a,2f7.1)') 'chk bad ini saln ',	&
            saln(k,i),' at i=',perm(i),k,				&
		' lat/lon=',deg_lat(i),deg_lon(i)
           saln(k,i)=max(min(saln(k,i),smax0),smin0)
         end if
       end do
     end if

     if (i.eq.itest) then
       if (odepth(i).le.0.) print '(a,i8,a)',			&
        'WARNING: ocean test pt',perm(i),' is land point'
     end if
   end do		! i loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

!SMS$EXCHANGE(temp,saln,pres)
!SMS$PARALLEL(dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(event,t8,s8)
   do i=1,nip
    if (wet(i) > 0 ) then
     if (odepth(i).le.0.) then	! single land point just converted to ocean
      event=.false.
      do edg=1,nprox(i)
       if (odepth(prox(edg,i)).gt.0.) then
        saln(:,i)=saln(:,prox(edg,i))
        temp(:,i)=temp(:,prox(edg,i))
        pres(:,i)=pres(:,prox(edg,i))
        odepth(i)=100.			! ocean depth
        ph3d(1,i)=0.			! terrain height in FIM
        event=.true.
        exit
       end if
      end do
      if (.not. event) then
       print *,'no ocean data at newly created land point',perm(i)
       stop '(point without data)'
      end if
     end if

     do k=1,kdm
      dp(k,i,1)=pres(k+1,i)-pres(k,i)
      t8=temp(k,i)
      s8=saln(k,i)
      dens(k,i)=sigocn(t8,s8)
      if (pres(k,i).gt.(odepth(i)-.01)*onem) then	! massless lyr at bottom
       saln(k,i)=35.
       temp(k,i)=tofsig(theta(k),saln(k,i))
      end if
     end do			! vertical loop
    end if			! ocean point
   end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

   call hybgn1(0,1,theta,dens,dp,pres,temp,saln,uclin,vclin,deg_lat)

   call aperef(dens,pres,theta,refp)

! --- get bottom depth by solving hydrostatic eqn in top-down direction
!SMS$PARALLEL(dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(q)
   do i=1,nip
    rivflo(i)=0.					! start with dry rivers
    srfht(i)=ph3d(1,i)
    if (wet(i) > 0 ) then
     do k=1,kdm
      dp(k,i,2)=dp(k,i,1)
      q=-thref*dens(k,i)			! - dens anomal / ref density
      spcifv(k,i)=thref*q*(1.+q*(1.+q))		! specific volume anomaly
     end do
!-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>
! --- experimental method:
!    refmnt(1)=0.
!    do k=1,kdm-1
!     refmnt(k+1)=refmnt(k)+min(refp(k+1),pres(kdm+1,i))		&
!       *(spcifv(k+1,i)-spcifv(k,i))
!    end do
!    montg(kdm,i)=refmnt(kdm)
!    do k=kdm-1,1,-1
!     montg(k,i)=montg(k+1,i)-pres(k+1,i)*(spcifv(k+1,i)-spcifv(k,i))
!    end do
!-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>
! --- traditional method:
     montg(1,i)=0.
     do k=2,kdm
      montg(k,i)=montg(k-1,i)+pres(k,i)*(spcifv(k,i)-spcifv(k-1,i))
     end do
     geop(1,i)=0.
! --- omit the contribution of reference density to geopotential
! --- (i.e. the dynamically irrelevant purely pressure-dependent part)
     do k=1,kdm
      geop(k+1,i)=geop(k,i)-spcifv(k,i)*(pres(k+1,i)-pres(k,i))
     end do
!-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>

     if (i.eq.itest) then
      print '(i8,a,2f13.2)',perm(i),' (ocn_init) bottom geopot:',	&
        montg(kdm,i)-pres(kdm+1,i)*(thref+spcifv(kdm,i)),		&
        geop(kdm+1,i)-pres(kdm+1,i)*thref
      print '(i8,a/(i3,f9.2,2(f13.1,f9.2)))',perm(i),			&
        '  init.montg.pot. from relaxed and actual state:',		&
        (k,dens(k,i),min(refp(k+1),pres(kdm+1,i))/onem,			&
        refmnt(k)/grvity,pres(k+1,i)/onem,montg(k,i)/grvity,k=1,kdm)
      write(*,'(a,2f8.2)') 'airtmn,SST=',airtmn(i), temp(1,i)
     end if

     psikk(i)=montg (kdm  ,i)
     spvkk(i)=spcifv(kdm  ,i)
     pbot (i)=pres  (kdm+1,i)
     gzbot(i)=geop  (kdm+1,i)

     thksno(i)=0.
     ticeol(i)=temice(i)

     if (i.eq.itest) then
      print '(a,i8,5(a,f7.2))','point',perm(i),' is at lat',deg_lat(i),	&
        ', lon',deg_lon(i),' covice=',covice(i),' thkice=',thkice(i),	&
        ' temice=',temice(i)
      print '(i8,a/(i22,2f8.1,3f8.2,f8.2))',perm(i),			&
       ' column initl.cond  dp    pres    dens    temp    saln   montg',&
       (k,dp(k,i,1)/onem,pres(k+1,i)/onem,dens(k,i),temp(k,i),		&
       saln(k,i),montg(k,i)/grvity,k=1,kdm)
      call flush(6)
     end if
    end if			! ocean point
   end do			! horiz.loop
!$OMP END PARALLEL DO
!SMS$EXCHANGE(pbot,geop)
!SMS$PARALLEL END
   call glob2d(pbot,totmass)
   work1(:,:)=1.
   work2(:,:)=dp(:,:,2)
   call glob3d(work1,work2,totmass)
   totmass=totmass/grvity				! kg
   print '(a,es11.3,7x,a,f9.1)','total ocean mass (1.342e21kg@g5):',	&
     totmass,'avg.water depth (3711.9m@g5):',totmass*grvity/(onem*ocnarea)

!  work(:)=temp(1,:)
!  call stencl(work,1,1.,'(hycom_init) temp0')
!  work(:)=saln(1,:)
!  call stencl(work,1,1.,'(hycom_init) saln0')

   call findmxmn1(spvkk,nip,'(inicon) spvkk',wet)
   call findmxmn1(psikk,nip,'(inicon) psikk',wet)
   do k=1,kdm
    write (string,'(a,i2)') ' k=',k
    call findmxmn2(dens,kdm,nip,k,'(inicon) dens'//string,wet)
    call findmxmn2(temp,kdm,nip,k,'(inicon) temp'//string,wet)
    call findmxmn2(saln,kdm,nip,k,'(inicon) saln'//string,wet)
   end do

! --- specify initial tracer field(s)

   if (numtr.ge.1) then
!SMS$PARALLEL(dh,i) BEGIN
!$OMP PARALLEL DO
    do i=1,nip
     if (wet(i) > 0 ) then
      passv_tr(:,i,1)=deg_lat(i)+90.	! tracer 3: initial latitude + 90
     end if
    end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
   end if

   if (numtr.ge.2) then
!SMS$PARALLEL(dh,i) BEGIN
!$OMP PARALLEL DO
    do i=1,nip
     if (wet(i) > 0 ) then
      do k=1,kdm
       passv_tr(k,i,2)=float(k)		! tracer 4: initial layer index
      end do
     end if
    end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
   end if

! --- determine precipitation pathways (river runoff) to ocean
   if (coupled .and. do_rnoff) then 
    call runoff(srfht,nuphill,uphill,dnhill)
   else
    nuphill(:)=0
    dnhill(:)=0
    uphill(:,:)=0
   end if

! Generate simple wind stress pattern for testing purposes

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!SMS$SERIAL BEGIN
!   print *,'prescribed wind stress (Pa) as function of latitude:'
!   do degr=90.,-90.,-2.
!    if (abs(degr).lt.maxwest) then
!     q=-cosd(degr/maxwest*180.)
!    else
!     q=cosd((abs(degr)-maxwest)/(90.-maxwest)*90.)**2
!    end if
!    call linout(35.+30.*q,'X',int(degr))
!   end do
!!SMS$SERIAL END
!   call flush(6)
!
!!SMS$PARALLEL(dh,i) BEGIN
!   do i=1,nip
!    degr=deg_lat(i)
!    if (abs(degr).lt.maxwest) then
!     taux(i)=-.1*cosd(degr/maxwest*180.)
!    else
!     taux(i)=.1*cosd((abs(degr)-maxwest)/(90.-maxwest)*90.)**2
!    end if
!    tauy(i)=0.
!   end do
!!SMS$PARALLEL END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  print *,'subr.inicon calling get_stress ...'
!  call get_stress(taux,tauy)

!  call stencl(taux,1,100.,'(ocn inicon) eastwd stress (10^-2 Pa)')
!  call findmxmn1(taux,nip,'taux')

! --- initialize transport integrals used for tracer advection
   call transp0(0,1,cumuflx,dp,dpinit)

   print *,'... calling glbdia, nstep=',nstep_ocn
   write(string,'(a,i8)') 'step',nstep_ocn
   call glbdia(nstep_ocn,1,trim(string))	! restart file no ready yet
  end if                      ! .not. readrestart

  if (ArchvTimeUnit.eq.'ts') then
    TotalStepOcn = TotalTime
    ArchvStepOcn = ArchvIntvl
  else if (ArchvTimeUnit.eq.'mi') then
    TotalStepOcn = TotalTime*60./batrop
    ArchvStepOcn = nint(60.*ArchvIntvl/batrop)
  else if (ArchvTimeUnit.eq.'hr') then
    TotalStepOcn = TotalTime*3600./batrop
    ArchvStepOcn = nint(3600.*ArchvIntvl/batrop)
  else if (ArchvTimeUnit.eq.'dy') then 
    TotalStepOcn = TotalTime*24*3600./batrop
    ArchvStepOcn = nint(24.*3600.*ArchvIntvl/batrop)
  else if (ArchvTimeUnit.eq.'mo') then
    TotalStepOcn = TotalTime*hrs_in_month*3600./batrop
    ArchvStepOcn = nint(hrs_in_month*3600.*ArchvIntvl/batrop)  
  else
    write (*,'(a,a)') '(hycom_inicon) unrecognized output time unit: ',	&
      ArchvTimeUnit
    stop
  endif

  print *,'Ocean model batrop/baclin time step (sec):',batrop,batrop*bclin_frq
  print *,'Ocean Total Time, Archvintvl =',TotalTime,ArchvTimeUnit,	&
    ArchvIntvl,ArchvTimeUnit
  print *,'Ocean Total Step, ArchvStep =',TotalStepOcn,ArchvStepOcn
  print *,'Ocean diagnostics interval (time steps) =',diag_intvl

  if (iocnmx.le.2) then 
    call inikpp
  end if

  call returnunit (unitno)

!sms$compare_var(temp, "temp_inicon")
!sms$compare_var(saln, "saln_inicon")

  print *,'... exiting hyc_inicon'
  return
  end subroutine inicon


  subroutine aperef(dnshyb,prshyb,targt,refp)

! --- determine the flat reference state (i.e., the unavailable pot.energy)
! --- for a stack of sloping hybrid-isopycnic layers

  use module_constants,only: area,perm,inv_perm
  use fimnamelist    ,only: kdm,itest
  use hycom_constants ,only: wet,onem
  use hycom_variables ,only: ocnarea

!SMS$DISTRIBUTE (dh,2) BEGIN
  real ,intent(IN) :: prshyb(kdm+1,nip)
  real ,intent(IN) :: dnshyb(kdm  ,nip)
  real             :: pres(kdm+1,nip)
!SMS$DISTRIBUTE END
  real ,intent(IN) :: targt(kdm)
  real,intent(OUT) :: refp(kdm+1)
  integer i,k,n
  logical vrbos
  real   ,parameter :: slithk=1.
  integer,parameter :: nscli=7000./slithk
  real weight(kdm),wgt,slitop,slibot,slisum,sliwgt(nscli),	&
       dnsold(kdm),dnsnew(kdm),prsold(kdm+1),prsnew(kdm+1)

! --- step 1: transform to purely isopycnic layers

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(vrbos,dnsold,prsold,dnsnew,prsnew)
  do i=1,nip
   if (wet(i).gt.0) then
    vrbos=i.eq.itest
    dnsold(:)=dnshyb(:,i)
    prsold(:)=prshyb(:,i)
    call restp_1d(dnsold,prsold,kdm,dnsnew,prsnew,targt,kdm,vrbos,perm(i))
    pres(:,i)=prsnew(:)
   end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! --- step 2: find total mass above each isopycnic interface

  do k=1,kdm
   wgt=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO REDUCTION(+:wgt)
   do i=1,nip
    if (wet(i) > 0 ) wgt=wgt+pres(k+1,i)*area(i)
   end do
!$OMP END PARALLEL DO
!SMS$REDUCE(wgt,SUM)
!SMS$PARALLEL END
   weight(k)=wgt
  end do

! --- step 3: divide basin into stack of thin horizontal slices and find mass
! --- of each slice

  do n=1,nscli
   slitop=float(n-1)*slithk*onem
   slibot=float(n  )*slithk*onem
    wgt=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO REDUCTION(+:wgt)
   do i=1,nip
    if (wet(i) > 0 ) 						&
     wgt=wgt+area(i)*(min(pres(kdm+1,i),slibot)-		&
                      min(pres(kdm+1,i),slitop))
   end do
!$OMP END PARALLEL DO
!SMS$REDUCE(wgt,SUM)
!SMS$PARALLEL END
   sliwgt(n)=wgt
  end do

  do n=1,nscli-99,100
   print '(a,i5,es11.3)','(aperef) mass (m) in slice',n,	&
      sliwgt(n)/(ocnarea*onem)
  end do

! --- step 4: add slices vertically until sum exceeds combined mass of layers
! --- 1...k. this tells us where bottom of layer k will be when flattened

  refp(1)=0.
  refp(kdm+1)=float(nscli)*slithk*onem
  do 5 k=1,kdm-1
   slisum=0.
   do n=1,nscli
    slitop=float(n-1)*slithk*onem
    slibot=float(n  )*slithk*onem
    slisum=slisum+sliwgt(n)
    if (slisum.ge.weight(k)) go to 7
   end do
   print *,'k =',k,'  error: slisum < weight',slisum,weight(k)
   go to 5

! --- interpolate among slices to get precise depth of flattened interface

 7 refp(k+1)=(slibot*(weight(k)-slisum+sliwgt(n))-		&
              slitop*(weight(k)-slisum          ))/sliwgt(n)
   print '(a,i3,f8.1,a,f8.1)','(aperef) mass (m) above intfc',	&
     k+1,weight(k)/(ocnarea*onem),'  flattened intfc depth:',	&
     refp(k+1)/onem
 5 continue

  return
  end subroutine aperef


   subroutine restp_1d(thold,pold,kold,thnew,pnew,targt,knew,vrbos,i)
! --- convert a stairstep (i.e., piecewise constant) theta profile into a
! --- stairstep profile constrained to have prescribed theta ('targt') steps.
       
! --- input   variables: thold,pold,targt,kold,knew,vrbos
! --- output variables: thnew,pnew

   implicit none
   integer,intent(IN)  :: kold,knew,i
   real ,intent(IN)    :: thold(kold),pold(kold+1),targt(knew)
   real ,intent(INOUT) :: thnew(knew),pnew(knew+1)
   logical,intent(IN)  :: vrbos

   integer k,ko
   real oldth(kold)
   real scale,cloutt,colint,pinteg,tha,thb
   real,parameter :: acurcy=1.e-6
       
   if (vrbos)								&
     write (6,101) i,							&
      'restp1 -- input profile:    theta     thknss      press',	&
      pold(1),(k,thold(k),pold(k+1)-pold(k),pold(k+1),k=1,kold)
101 format (i8,2x,a/54x,f11.1/(i32,f11.3,2f11.1))
       
! --- remove theta inversions from input profile
!!     oldth(kold)=thold(kold)
!!     do k=kold,2,-1
!!       oldth(k-1)=min(oldth(k),thold(k-1))
   oldth(1)=thold(1)
   do k=2,kold
     oldth(k)=max(oldth(k-1),thold(k))
   end do
       
   thnew(:)=targt(:)
   thnew(   1)=min(oldth(1),oldth(kold),targt(   1))
   thnew(knew)=max(oldth(1),oldth(kold),targt(knew))
   pnew(     1)=pold(     1)
   pnew(knew+1)=pold(kold+1)
       
! --- column integrals (colin/clout) are computed for diagnostic purposes only
       
   cloutt=0.
   colint=0.
   scale=0.
   do k=1,kold
     colint=colint+oldth(k)*(pold(k)-pold(k+1))
     scale=max(scale,abs(oldth(k)))
   end do
       
! --- find interface pnew(k+1) separating layers k and k+1 by requiring
! --- that integral over prs*d(theta) from thnew(k) to thnew(k+1) be preserved.
       
   ko=1
   do k=1,knew-1
     pinteg=0.
     thb=thnew(k)
 5   tha=thb
     thb=min(thnew(k+1),max(thnew(k),oldth(ko)))

!        if (vrbos) print 103,i,k+1,'  increment pinteg by',
!     .     pold(ko),thb,tha
! 103    format  (i5,i3,a,f11.1,'*(',f6.3,'-',f6.3,')')

   pinteg=pinteg+pold(ko)*(thb-tha)
   if (oldth(ko) < thnew(k+1)) then
     if (ko.lt.kold) then
       ko=ko+1
       go to 5
     end if
     tha=thb
     thb=thnew(k+1)

!    if (vrbos) print 103,i,k+1,'  increment pinteg by',	&
!         pold(kold+1),thb,tha

     pinteg=pinteg+pold(kold+1)*(thb-tha)
   end if
   pnew(k+1)=pnew(k)
   if (thnew(k+1) > thnew(k))					&
      pnew(k+1)=pinteg/(thnew(k+1)-thnew(k))

!        if (vrbos) print 102,i,k+1,'  pinteg,thnew(k/k+1)=',
!     .     pinteg,thnew(k),thnew(k+1),pnew(k+1)
! 102    format (i5,i3,a,f11.1,2f7.3,' => pnew=',f11.1)

     cloutt=cloutt+thnew(k)*(pnew(k)-pnew(k+1))
   enddo
       
   cloutt=cloutt+thnew(knew)*(pnew(knew)-pnew(knew+1))
   if (abs(cloutt-colint).gt.acurcy*scale*pnew(knew+1))		&
     write (6,100) i,'restp1 - column intgl.error',		&
      colint,cloutt,(cloutt-colint)/colint
 100     format (i8,3x,a,2es14.6,es9.1)

   if (vrbos) write (6,101) i,						&
      'restp1 -- outpt profile:    theta     thknss      press',	&
      pnew(1),(k,thold(k),pnew(k+1)-pnew(k),pnew(k+1),k=1,knew)

   return
   end subroutine restp_1d


! !   subroutine get_stress(taux,tauy)
! !   use findmaxmin1
! ! 
! ! ! --- convert surface wind into wind stress
! ! 
! !   use module_variables,only: us3d,vs3d
! !   use module_constants,only: deg_lat
! ! 
! !   implicit none
! ! !SMS$DISTRIBUTE (dh,1) BEGIN
! !   real,intent(OUT) :: taux(nip),tauy(nip)
! !   real usrf(nip),vsrf(nip),wind(nip),zonl(nip)
! ! !SMS$DISTRIBUTE END
! !   real,parameter :: cd=.003		! drag coefficient
! !   real,parameter :: airdns=1.2	! air density (kg/m^3)
! !   integer i
! !   real valmin,valmax
! ! 
! ! !SMS$PARALLEL (dh,i) BEGIN
! !   do i=1,nip
! !    usrf(i)=us3d(1,i)
! !    vsrf(i)=vs3d(1,i)
! !    wind(i)=sqrt(us3d(1,i)**2+vs3d(1,i)**2)
! !   end do
! ! !SMS$PARALLEL END
! !   
! !   call flush(6)
! !   call findmxmn1(usrf,nip,'raw usrf')
! ! ! call findmean1(usrf,nip,'raw usrf')
! !   call findmxmn1(vsrf,nip,'raw vsrf')
! ! ! call findmean1(vsrf,nip,'raw vsrf')
! !   call findmxmn1(wind,nip,'raw wind')
! ! ! call findmean1(wind,nip,'raw wind')
! !   call flush(6)
! ! 
! !   call zonlav(deg_lat,usrf,zonl,'usrf')
! !   usrf(:)=zonl(:)
! !   call zonlav(deg_lat,vsrf,zonl,'vsrf')
! !   vsrf(:)=zonl(:)
! !   call zonlav(deg_lat,wind,zonl,'wind')
! !   wind(:)=zonl(:)
! ! 
! !   call flush(6)
! !   call findmxmn1(usrf,nip,'zon.avg usrf')
! ! ! call findmean1(usrf,nip,'zon.avg usrf')
! !   call findmxmn1(vsrf,nip,'zon.avg vsrf')
! ! ! call findmean1(vsrf,nip,'zon.avg vsrf')
! !   call findmxmn1(wind,nip,'zon.avg wind')
! ! ! call findmean1(wind,nip,'zon.avg wind')
! !   call flush(6)
! ! 
! !   valmin= 1.e33
! !   valmax=-1.e33
! ! !SMS$PARALLEL (dh,i) BEGIN
! !   do i=1,nip
! !    taux(i)=cd*airdns*wind(i)*usrf(i)
! ! !  tauy(i)=cd*airdns*wind(i)*vsrf(i)
! !    tauy(i)=0.
! !   valmin=min(valmin,taux(i))
! !   valmax=max(valmax,taux(i))
! !   end do
! ! !SMS$REDUCE(valmin,MIN)
! ! !SMS$REDUCE(valmax,MAX)
! ! !SMS$PARALLEL END
! !   print '(a,2es15.7)','subr.get_stress - taux min,max =',valmin,valmax
! ! 
! !   return
! !   end subroutine get_stress
! 
! 
!   subroutine zonlav(deg_lat,u,uav,what)
!   use findmaxmin1
! 
! ! --- store zonal average of -u- field in -uav-. overwriting of -u- allowed
! 
!   implicit none
!   character,intent(IN) :: what*(*)
! !SMS$DISTRIBUTE (dh,1) BEGIN
!   real,intent(IN)  :: deg_lat(nip)
!   real,intent(IN)  :: u(nip)
!   real,intent(OUT) :: uav(nip)
! !SMS$DISTRIBUTE END
!   integer,parameter :: nbins=18		! use 10 deg bins
!   integer count(nbins)
!   real bin(nbins+1),binwdth
!   integer i,n
! 
!   bin(:)=0.
!   count(:)=0
!   binwdth=180./nbins
! 
!   call findmxmn1(u,nip,'zonlav input '//what)
! ! call findmean1(u,nip,'zonlav input '//what)
! 
! !SMS$SERIAL (<u,IN>) BEGIN
!   do i=1,nip
!    n=1.+(deg_lat(i)+90.)/binwdth
!    n=max(1,min(nbins,n))
!    count(n)=count(n)+1
!    bin(n)=bin(n)+u(i)
!   end do
! !SMS$SERIAL END
! 
!   print '(a/(2i9))','zonalv -- total number by bin:',(n,count(n),n=1,nbins)
!   do n=1,nbins
!    bin(n)=bin(n)/count(n)
!   end do
!   print '(3a/(i9,f9.1,f9.2))','zonlav -- ',what,' in latitude bins:',	&
!    (n,(n-.5)*binwdth-90.,bin(n),n=1,nbins)
! 
! !SMS$PARALLEL (dh,i) BEGIN
!   do i=1,nip
!    n=1.+(deg_lat(i)+90.)/binwdth
!    n=max(1,min(nbins,n))
!    uav(i)=(bin(n  )*(binwdth*(n+.5)-(deg_lat(i)+90.))		&
!             +bin(n+1)*((deg_lat(i)+90.)-binwdth*(n-.5)))/binwdth
!   end do
! !SMS$PARALLEL END
! 
!   call findmxmn1(uav,nip,'zonlav output '//what)
! ! call findmean1(uav,nip,'zonlav output '//what)
! 
!   return
!   end subroutine zonlav


  subroutine read_ocndata_cfsr	! final output: T/S/P
  use netcdf
  use module_control,only: nip
  USE slint, ONLY: slint_init_ocn, bilinear_interp_ocn,nn_interp_ocn
  use module_constants,only: lat,lon,pi
  use hycom_constants,only: land_spval,odepth,theta,salmin,onem
  use fimnamelist    ,only: kdm,itest,yyyymmddhhmm
  use hycom_variables ,only: temp,saln,pres,thkice,covice,temice
  use hycom_sigetc, only: sigocn
  use hycom_thermf, only: errhandl,readcdf_2d,readcdf_3d

  real :: tout1d(kdm),sout1d(kdm),pout1d(kdm+1),thxtnd(kdm+1),		&
          p_hyb(kdm+1),c,radius0
  real*8 :: t8,s8

  character*132 :: ocn_ts_data, ice_data,ll_data
  character chr10*10
  integer, parameter:: idmr=720,jdmr=410,kdmr40=40,kdmr6=6,kdmr5=5
  real:: tempr(idmr,jdmr,kdmr40), salnr(idmr,jdmr,kdmr40),		&
         iceh (idmr,jdmr,kdmr5), hicer(idmr,jdmr),			&
         icec (idmr,jdmr,kdmr6), cicer(idmr,jdmr),			&
         icet (idmr,jdmr,kdmr5), ticer(idmr,jdmr),			&
         latr (idmr,jdmr),       lonr (idmr,jdmr),			&
         work1(idmr,jdmr),       work2(idmr,jdmr),			&
         tem1d(kdmr40), sal1d(kdmr40), sig1d(kdmr40), z1d(kdmr40),	&
         tempr1d(idmr*jdmr),salnr1d(idmr*jdmr),scal(jdmr),		&
         cice1d(idmr*jdmr),tice1d(idmr*jdmr),hice1d(idmr*jdmr)
  real, allocatable :: grid1(:,:),grid2(:,:),x1dicos(:),tempi(:,:),salni(:,:)
  integer, allocatable :: wet(:)

  integer :: i,j,k,l,kmax1,kr,iz,i0,i1,j1,n,maskr(idmr,jdmr),n2
  integer :: ncid1, ncid2, ncid3 ! netcdf file id returned from nf90_open
  integer :: id_idm, id_jdm, id_kdm40, id_kdm5, id_kdm6  ! dimension id for i,j,k
! Variable ids
  integer :: id_fld1, id_fld2
! Dimension lengths
  integer :: i00=-1, j00=-1, k00_40=-1, k00_5=-1, k00_6=-1   ! init to bad value
! Required for input to nf90_get_var
  integer :: start1d(2),start2d(3),start3d(4)   ! starting indices for the variable in the netcdf file
  integer :: kount1d(2),kount2d(3),kount3d(4)   ! lengths for the variable in the netcdf file
  integer :: itest0, jtest0	! i,j in CFSR corresponding to itest in icos grid
  integer :: dimids(4)
  logical :: vrbos
!
  real, parameter:: zlv(40)=(/ 						&
    5.,   15.,   25.,   35.,   45.,   55.,   65.,   75.,   85.,   95., 	&
  105.,  115.,  125.,  135.,  145.,  155.,  165.,  175.,  185.,  195., 	&
  205.,  215.,  225.,  238.,  262.,  303.,  367.,  459.,  585.,  747., 	&
  950., 1193., 1480., 1807., 2175., 2579., 3017., 3483., 3972., 4478. /)

  print *,'entering read_ocndata_cfsr...'
  allocate (grid2(nip,2),x1dicos(nip),wet(nip),tempi(kdmr40,nip),	&
            salni(kdmr40,nip))
  vrbos=itest > 0

!SMS$SERIAL (<lat,lon,odepth,perm,IN>,<tempi,salni,thkice,covice,temice,OUT> : default=ignore) BEGIN
  ll_data = "grid_spec_05.nc.T382"
  ocn_ts_data = "ocnanl.gdas." // yyyymmddhhmm(1:10) //".ocean_temp_salt.res.nc"
  ice_data    = "ocnanl.gdas." // yyyymmddhhmm(1:10) //".ice_model.res.nc"
  write(*,*)'reading cfsr from ',trim(ocn_ts_data)

  call errhandl(nf90_open (path=trim(ocn_ts_data), mode=NF90_NOWRITE,	&
                ncid=ncid1))   ! netcdf3
  call errhandl(nf90_open (path=trim(   ice_data), mode=NF90_NOWRITE,	&
                ncid=ncid2))   ! netcdf3
  call errhandl(nf90_open (path=trim(    ll_data), mode=NF90_NOWRITE,	&
                ncid=ncid3))   ! netcdf3

  call errhandl (nf90_inq_dimid (ncid1, 'xaxis_1', id_idm))
  call errhandl (nf90_inq_dimid (ncid1, 'yaxis_1', id_jdm))
  call errhandl (nf90_inq_dimid (ncid1, 'zaxis_1', id_kdm40))
  call errhandl (nf90_inq_dimid (ncid2, 'zaxis_1', id_kdm6))
  call errhandl (nf90_inq_dimid (ncid2, 'zaxis_2', id_kdm5))

  call errhandl (nf90_inquire_dimension (ncid1, id_idm,   len=i00))
  call errhandl (nf90_inquire_dimension (ncid1, id_jdm,   len=j00))
  call errhandl (nf90_inquire_dimension (ncid1, id_kdm40, len=k00_40))
  call errhandl (nf90_inquire_dimension (ncid2, id_kdm6,  len=k00_6))
  call errhandl (nf90_inquire_dimension (ncid2, id_kdm5,  len=k00_5))

  if (i00.ne.idmr .or. j00.ne.jdmr .or. k00_40.ne.kdmr40 .or.		&
      k00_6.ne.kdmr6 .or. k00_5.ne.kdmr5) then
    write(*,'(a,5i6)') 'wrong cfsr dimensions: ',i00,j00,k00_40,k00_6,k00_5
    stop 'wrong cfsr dimension'
  else 
    write(*,'(a,5i4)') 'cfsr dimensions=720,410,40,6,5; got ',		&
      i00,j00,k00_40,k00_6,k00_5
  end if

  start2d = (/1,1,1/)
  kount2d = (/idmr,jdmr,1/)
  start3d = (/1,1,1,1/)
  kount3d = (/idmr,jdmr,kdmr40,1/)
  call readcdf_3d(ncid1,'temp', tempr,idmr,jdmr,kdmr40,start3d,kount3d)
  call readcdf_3d(ncid1,'salt', salnr,idmr,jdmr,kdmr40,start3d,kount3d)

  kount3d = (/idmr,jdmr,kdmr6,1/)
  call readcdf_3d (ncid2,'part_size',icec,idmr,jdmr,kdmr6,start3d,kount3d)

  kount3d = (/idmr,jdmr,kdmr5,1/)
  call readcdf_3d (ncid2,'h_ice',iceh,idmr,jdmr,kdmr5,start3d,kount3d)
  call readcdf_3d (ncid2,'t_ice1',icet,idmr,jdmr,kdmr5,start3d,kount3d)

  kount2d = (/idmr,jdmr,1/)
! x_T, y_T: T-cell center or C-cell center
  call readcdf_2d (ncid3,'x_T',lonr,idmr,jdmr,start2d,kount2d)
  call readcdf_2d (ncid3,'y_T',latr,idmr,jdmr,start2d,kount2d)

  do j=1,jdmr
  do i=1,idmr
    lonr(i,j)=mod(lonr(i,j)+360.,360.)  ! to convert lonr to [0:360]
  end do
  end do

  hicer(:,:)=0.
  ticer(:,:)=0.
  do k=1,kdmr5
  hicer(:,:)=hicer(:,:)+icec(:,:,k+1)*iceh(:,:,k)	! mean ice thickness
  ticer(:,:)=ticer(:,:)+icec(:,:,k+1)*icet(:,:,k)	! mean ice temperature
  end do
  cicer(:,:)=1.-icec(:,:,1)

  grid2(:,1)=lat(:)     ! in radians
  grid2(:,2)=lon(:)     ! in radians
  radius0=2.0*pi/180.
  do j=1,jdmr
  scal(j)=cos(latr(360,j)*pi/180.)  ! i=360 is chosen arbitrarily
  end do


  if (vrbos) then
    do j=1,jdmr-1
    do i=1,idmr-1
      if (lon(itest)*180/pi.ge.lonr(i,j) .and.			&
          lon(itest)*180/pi.lt.lonr(i+1,j+1) .and.		&
          lat(itest)*180/pi.ge.latr(i,j) .and.			&
          lat(itest)*180/pi.lt.latr(i+1,j+1)) then
        itest0=i
        jtest0=j
        go to 2
      end if
    end do
    end do
    write(*,*) 'wrong: cannot find i,j in cfsr'
 2  i=itest0
    j=jtest0
    write(*,'(a,i8,a,f6.0,2f7.1,a,2i4,2f7.1)') 'itest=',perm(itest),	&
      ' depth,lat,lon=',odepth(i),lat(itest)*180/pi,lon(itest)*180/pi,	&
      ' found i,j,lat,lon in cfsr=',i,j,latr(i,j),lonr(i,j)
  end if

  n2=2
  if (glvl >= 8) then
    n2=30
  elseif (glvl == 6) then
    n2=6
  elseif (glvl == 7) then
    n2=20
  end if

  print *,'cfsr amend n2=',n2
  do k=1,kdmr40
    do n=1,n2
    do j=1,jdmr
    do i=1,idmr
      if (salnr(i,j,k).gt.0.01) then
        maskr(i,j)= 1   ! points w. good values have mask of 1
      else
        maskr(i,j)=-1   ! land points have mask of -1
      end if
    end do
    end do

    do j=1,jdmr
    do i=1,idmr
      if (salnr(i,j,k).gt.0.01) then    ! nearby +/-1 points are filled with 0 or above
        do i1=i-1,i+1
        i0=mod(i1+idmr-1,idmr)+1
        do j1=max(j-1,1),min(jdmr,j+1)
          maskr(i0,j1)=max(maskr(i0,j1),0)      ! "0" points need to be filled
        end do
        end do
      end if
    end do
    end do

  if (vrbos) then
    i=itest0
    j=jtest0
    write(*,'(a,2i4,3(a,i2))') 'raw i,j=',i,j,' mask=',maskr(i,j),	&
      ' k=',k,' n=',n
    write (chr10,'(a,i2.2)') 'temp k=',k
    call pr_9x9(tempr(1,1,k),idmr,jdmr,itest0,jtest0,0.,1.,chr10)
    write (chr10,'(a,i2.2)') 'saln k=',k
    call pr_9x9(salnr(1,1,k),idmr,jdmr,itest0,jtest0,0.,1.,chr10)
  end if	! vrbos

  call amend2(tempr(1,1,k),salnr(1,1,k),maskr,idmr,jdmr,idmr,jdmr,1,	&
    work1,work2,scal)

  if (vrbos) then
    write(*,'(2(a,i2))') 'after k=',k,' n=',n
    write (chr10,'(a,i2.2)') 'temp k=',k
    call pr_9x9(tempr(1,1,k),idmr,jdmr,itest0,jtest0,0.,1.,chr10)
    write (chr10,'(a,i2.2)') 'saln k=',k
    call pr_9x9(salnr(1,1,k),idmr,jdmr,itest0,jtest0,0.,1.,chr10)
  end if	! vrbos

  end do   ! do n=1:n2

! --- done with amend, start map ocean ice/T/S
    kr=0
    do j=1,jdmr
    do i=1,idmr
      if (salnr(i,j,k).gt.0.01) then     ! ocean point
        kr=kr+1
      end if
    end do
    end do
    write(*,'(a,i2,a,i7)') 'k=',k,' kr=',kr
    allocate (grid1(kr,2))

    tempr1d(:)=0.
    salnr1d(:)=0.
    tice1d(:)=0.
    hice1d(:)=0.
    cice1d(:)=0.
    kr=0
    do j=1,jdmr
    do i=1,idmr
      if (salnr(i,j,k).gt.0.01) then     ! ocean point
        kr=kr+1
        grid1(kr,1)=latr(i,j)*pi/180.
        grid1(kr,2)=lonr(i,j)*pi/180.
        tempr1d(kr) = tempr(i,j,k)
        salnr1d(kr) = salnr(i,j,k)
        if (k==1) then
          tice1d(kr) = ticer(i,j)
          hice1d(kr) = hicer(i,j)
          cice1d(kr) = cicer(i,j)
        end if
!       if (vrbos.and.i==itest0.and.j==jtest0)				&
!         write(*,'(a,3i4,i7,a,2f8.3,es13.4)')				&
!         'input i,j,k,kr=',i,j,k,kr,' lat,lon=',latr(i,j),lonr(i,j)
      end if    ! if ocean pioint
    end do
    end do

    wet(:)=0
    do i=1,nip
      if (odepth(i).ge.zlv(k)) wet(i)=1
!     if (vrbos.and.abs(lat(i)*180/pi-lat0).le.1			&
!              .and.abs(lon(i)*180/pi-lon0).le.1) then
!       write(*,'(a,i6,a,f8.0,a,2f8.2)') 'wrong icos=',i,		&
!         ' topo=',odepth(i),' lat,lon=',lat(i)*180/pi,lon(i)*180/pi
!     end if
    end do

    call slint_init_ocn(grid1,kr,grid2,nip, wet, radius0)
    write(*,'(a,i4)') 'bilinear temp k=',k
    call bilinear_interp_ocn(tempr1d,x1dicos,wet)
    tempi(k,:)=x1dicos(:)
    write(*,'(a,i4)') 'bilinear saln k=',k
    call bilinear_interp_ocn(salnr1d,x1dicos,wet)
    salni(k,:)=x1dicos(:)

    if (k==1) then
      temice(:)=0.
      thkice(:)=0.
      covice(:)=0.
      call nn_interp_ocn(tice1d,temice,wet)
      call nn_interp_ocn(hice1d,thkice,wet)
      call nn_interp_ocn(cice1d,covice,wet)
    end if

    deallocate (grid1)
  end do  ! do k
  write(*,'(a,7i4)') ' done bilinear'

  if (vrbos) then
  i=itest0
  j=jtest0
  l=itest
  write(*,'(a,2i4,a,3f7.1,a,2f8.2)')'cfsr at i,j=',i,j,			&
    ' t,c,h ice=',ticer(i,j),cicer(i,j),hicer(i,j),' ll=',		&
    latr(i,j),lonr(i,j)
  write(*,'(a,i8,a,3f7.1,a,2f8.2)') 'icos at  l=',perm(l),		&
    ' t,c,h ice=',temice(l),covice(l),thkice(l),' ll=',			&
    lat(l)*180/pi,lon(l)*180/pi
  write(*,'(a)')							&
   'hor.interp: k   z     T_bef   S_bef   rho_bef   T_aft    S_aft   rho_aft'
  do k=1,kdmr40
    write(*,'(7x,2i5,6f9.2)') k,int(zlv(k)),tempr(i,j,k),salnr(i,j,k),	&
      sigocn(max(-1.8,tempr(i,j,k))+0.d0,max(5.,salnr(i,j,k))+0.d0),	&
      max(-99.,tempi(k,l)),   max(-99.,salni(k,l)),			&
      sigocn(max(-1.8,tempi(k,l))+0.d0,max(5.,salni(k,l))+0.d0)
  enddo
! call stencl(tempi,kdm,1.,'(bil cfsr) temp')
! call stencl(salni,kdm,1.,'(bil cfsr) saln')

  end if
  write(*,'(a,7i4)') ' start lin2stp'
!SMS$SERIAL END

!SMS$PARALLEL(dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(vrbos,sal1d,tem1d,z1d,sig1d,c,kmax1,thxtnd,	&
!$OMP pout1d,tout1d,sout1d,p_hyb)
    do 955 i=1,nip
    if (odepth(i) > 0.) then
      vrbos=i.eq.itest

      do k=2,kdmr40
        if (odepth(i).lt.zlv(k).and.odepth(i).ge.zlv(k-1)) then   ! fill values at depth level
          salni(k,i)=salni(k-1,i)
          tempi(k,i)=tempi(k-1,i)
          go to 123
        end if
      enddo
 123  continue
!
!--- vertical interpolation
!
      do 961 k=1,kdmr40
      sal1d(k)=salni(k,i)
      tem1d(k)=tempi(k,i)
      z1d(k)=zlv(k)
      if (sal1d(k).lt.0. .or. tem1d(k).lt.-8.) then
        write(*,*) 'wrong cfsr t,s at i,k=',perm(i),k,'t,s=',tem1d(k),sal1d(k)
        stop 'wrong cfsr t,s'
      end if
      sig1d(k)=sigocn(tem1d(k)+0.d0,sal1d(k)+0.d0)
      if(zlv(k).lt.odepth(i)) goto 961
      if(zlv(k).eq.odepth(i)) goto 962
      c=(odepth(i)-zlv(k-1))/(zlv(k)-odepth(i))
      sig1d(k)=(sig1d(k-1)+c*sig1d(k))/(1.+c)
      sal1d(k)=(sal1d(k-1)+c*sal1d(k))/(1.+c)
      tem1d(k)=(tem1d(k-1)+c*tem1d(k))/(1.+c)
      kmax1=k
      z1d(k)=odepth(i)
      goto 963
961   continue
      k=kdmr40
962   kmax1=k
      z1d(kmax1)=odepth(i)
963   continue
!
    do k=kmax1+1,kdmr40
      sig1d(k)=sig1d(kmax1)
      sal1d(k)=sal1d(kmax1)
      tem1d(k)=tem1d(kmax1)
      z1d(k)=z1d(kmax1)
    enddo
!
    if (vrbos) then
      write(*,'(a,10f7.3/(7x,10f7.3))') ' bef sig1d =',(sig1d(k),k=1,kmax1)
    end if

    do k=2,kdmr40	! remove inversions
     sig1d(k)=max(sig1d(k),sig1d(k-1))
    end do
!
    if (vrbos) then
      write(*,'(a,10f7.3/(7x,10f7.3))') 'sig1d =',(sig1d(k),k=1,kmax1)
      write(*,'(a,10f7.2/(7x,10f7.2))') 'sal1d =',(sal1d(k),k=1,kmax1)
      write(*,'(a,10f7.2/(7x,10f7.2))') 'tem1d =',(tem1d(k),k=1,kmax1)
      write(*,'(a,10f7.1/(7x,10f7.1))') 'z1d   =',(z1d  (k),k=1,kmax1)
      do k=1,kdmr40
        write(*,'(i4,4f8.2)')k, z1d(k),sal1d(k),tem1d(k),sig1d(k)
      enddo
    end if
!
!--- turn piecewise linear profiles into stairstep profiles in density space
!
    thxtnd(    1)=sig1d(    1)
    thxtnd(kdm+1)=sig1d(kmax1)
    do k=1,kdm-1
    thxtnd(k+1)=min(thxtnd(kdm+1),max(thxtnd(k),theta(k)))
    end do
    pout1d(kdm+1)=z1d(kmax1)
!
    call lin2stp(sig1d,z1d,kmax1,thxtnd,pout1d,kdm,vrbos,perm(i))
    pout1d(   1)=z1d(    1)
!
    call lin2stp(z1d,tem1d,kmax1,pout1d,tout1d,kdm,vrbos,perm(i))
    call lin2stp(z1d,sal1d,kmax1,pout1d,sout1d,kdm,vrbos,perm(i))
!
    if (vrbos) then
     print *,'after lin2stp:'
     print '(a,i8,7x,a)','i=',perm(i),'pout1d   tem1d   sal1d   sig1d   targt'
     print '(i15,5f8.2)',(k,pout1d(k+1),tout1d(k),sout1d(k),	&
       sigocn(tout1d(k)+0.d0,sout1d(k)+0.d0),theta(k),k=1,kdm)
    endif
    
    do k=1,kdm
      sig1d(k)=sigocn(tout1d(k)+0.d0,sout1d(k)+0.d0)
    end do
    pout1d(1)=0.

    call hybgn2_1d(pout1d,sig1d,tout1d,sout1d,p_hyb,theta,vrbos,perm(i))	!output: tout1d,sout1d,p_hyb

    pres(1,i)=0.
    do k=1,kdm
      temp(k,i)=tout1d(k)
      saln(k,i)=max(sout1d(k),salmin(k))
      pres(k+1,i)=min(odepth(i),p_hyb(k+1))*onem
      if (k.eq.kdm) pres(k+1,i)=odepth(i)*onem
    end do
  end if	! if odepth > 0

 955 continue	! do i=1,nip
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!
  write(*,'(a,7i4)') '...exiting read_ocndata_cfsr'
  return
  end subroutine read_ocndata_cfsr


  subroutine lin2stp(xold,yold,kold,xnew,ynew,knew,vrbos,i)
!
! --- convert piecewise linear curve (xold,yold) into stairstep curve by
! --- integrating -yold- over consecutive -xnew- intervals
!
  implicit none
  integer,intent(in) :: kold,knew,i	!  i=current location in horiz.grid
  real,   intent(in) :: xold(kold),yold(kold),xnew(knew+1)
  real,  intent(out) :: ynew(knew)
  logical,intent(in) :: vrbos		!  if true, print diagnostics
  integer k,ko
  real colin,clout,colmx,yinteg,xlo,xhi,xa,xb,ya,yb,wgt
  real :: xxold(kold),yyold(kold),xxnew(knew+1),yynew(knew)
  logical at_top
! real,parameter :: acurcy=1.e-6
  real,parameter :: acurcy=1.e-3
!
  if (vrbos)								&
    write (*,101) i,'  lin2stp -- old profile:     x           y',	&
    (k,xold(k),yold(k),k=1,kold)
 101  format (i9,a/(i30,f12.3,f13.3))
!
  if (xold(1).lt.xold(kold)) then
    do k=1,kold
      xxold(k)= xold(k)
    enddo
    do k=1,knew+1
      xxnew(k)= xnew(k)
    enddo
  else
    do k=1,kold
      xxold(k)=-xold(k)
    enddo
    do k=1,knew+1
      xxnew(k)=-xnew(k)
    enddo
  end if
  do k=1,kold
    yyold(k)=yold(k)
  enddo
!
! --- column integrals (colin/clout) are computed for diagnostic purposes only
  if (xold(1).ne.xnew(1))						&
     print *,i,'  lin2stp warning - xold(1) and xnew(1) differ',	&
      xold(1),xnew(1)
  if (xold(kold).ne.xnew(knew+1))					&
     print *,i,'  lin2stp warning - xold(kk) and xnew(kk) differ',	&
     xold(kold),xnew(knew+1)
      colin=0.
      clout=0.
      colmx=0.
      do 3 k=1,kold-1
      colmx=max(colmx,abs(yyold(k)))
 3    colin=colin+.5*(yyold(k)+yyold(k+1))*(xxold(k+1)-xxold(k))
!
      at_top=.true.
      do 4 k=1,knew
      xlo=xxnew(k  )
      xhi=xxnew(k+1)
      yynew(k)=yyold(1)
      if (xhi.gt.xlo) then
       at_top=.false.
       yinteg=0.
       do ko=1,kold-1
        xa=max(xlo,min(xhi,xxold(ko  )))
        xb=max(xlo,min(xhi,xxold(ko+1)))
        if (xb.gt.xa) then
         wgt=(xa-xxold(ko))/(xxold(ko+1)-xxold(ko))
         ya=yyold(ko+1)*wgt+yyold(ko)*(1.-wgt)
         wgt=(xb-xxold(ko))/(xxold(ko+1)-xxold(ko))
         yb=yyold(ko+1)*wgt+yyold(ko)*(1.-wgt)
        if (vrbos) print '(2i3,a,4f7.1,2(1p,e11.4))',k,ko,	&
         ' xlo,xhi,xa,xb,ya,yb:',xlo,xhi,xa,xb,ya,yb
         yinteg=yinteg+.5*(ya+yb)*(xb-xa)
        end if
       end do
       yynew(k)=yinteg/(xhi-xlo)
      if (vrbos) print '(i3,a,1p,e12.4)',k,' ynew:',yynew(k)
       clout=clout+yinteg
      else if (at_top) then
       yynew(k)=yyold(   1)
      else
       yynew(k)=yyold(kold)
      end if
 4    continue
      do k=1,knew
      ynew(k)=yynew(k)
      enddo
 
!     if (abs(clout-colin).gt.acurcy*colmx*xold(kold/2))
!    . write (*,100) i,' lin2stp - column intgl.error',
!    .  colin,clout,(clout-colin)/colin
 100  format (i8,a,2(1p,e14.6),(1p,e9.1))
!
  if (vrbos) write (*,101) i,'  lin2stp -- new profile:     x           y',&
    (k,xnew(k),ynew(k),k=1,knew),knew+1,xnew(knew+1)
  return
  end subroutine lin2stp


  subroutine hybgn2_1d(p1d,r1d,t1d,s1d,pnew,targt,vrbos,i)

! --- hycom version 0.9.14

! --- -----------------------------------
! --- single-column hybrid grid generator (based on 'restep' technology)
! --- -----------------------------------

  use fimnamelist    ,only: kdm
  use hycom_constants,only: dpmin
  implicit none

  real,intent(INOUT)  :: p1d(kdm+1),r1d(kdm),t1d(kdm),s1d(kdm)
  real,intent(IN)     :: targt(kdm)
  real,intent(OUT)    :: pnew(kdm+1)
  logical,intent(IN)  :: vrbos
  integer,intent(IN)  :: i

  integer k,k1,ko,last
  real tnew(kdm),snew(kdm),q,dp0,dp1,dpsum,uintg,vintg,old,p_hat,	&
    torho,totem,tosal,tndrho,tndtem,tndsal,arg,arg1,cushn,pbot0
  real,parameter :: tolrnce=1.e-4
  real,parameter :: acurcy=1.e-5
  real c1,c2,c3,c4,c5,c6,c7,c8,c9,r,s,t,pref,sig

! --- for sigma1: Jackett et al. 2006, J. of atm & oceanic technology
  data c1,c2,c3,c4,c5,c6,c7,c8,c9/  & ! T:[-2,30],S:[30,38]
        5.019935E+00, 1.470290E-02, 7.904674E-01,-6.861022E-03, 	&
       -3.051459E-03, 3.099320E-05, 3.050299E-05, 1.046036E-04, 	&
       5.287567E-06/
  data pref/1.e7/

  sig(t,s)=c1+s*(c3+s*c8)+t*(c2+s*c5+t*(c4+s*c7+t*c6)+s*s*c9)

! --- simplified cushion function suitable for nonnegative first argument
  cushn(arg,arg1)=.25*min(arg,2.*arg1)**2/arg1+max(arg,2.*arg1)-arg1
!     cushn(arg,arg1)=arg1                  ! no cushn

  pbot0=p1d(kdm+1)

  if (vrbos) then
    write (*,'(i8,a)') i,' hybgn2  IN:  temp    saln    dens   thkns    dpth'
    do k=1,kdm
      write (*,103) k,t1d(k),s1d(k),r1d(k),p1d(k+1)-p1d(k),p1d(k+1)
    end do
 103   format (i18,3f8.3,f8.2,f8.1)
  end if

  if (vrbos) then
    write (*,99) i,'  hybgn2   o l d   profile'
    do k=1,kdm,10
      write (*,100) (p1d  (k1),k1=k,min(kdm+1,k+10))
      write (*,101) (r1d  (k1),k1=k,min(kdm,k+9))
      write (*,101) (targt(k1),k1=k,min(kdm,k+9))
      write (*,102) (t1d  (k1),k1=k,min(kdm,k+9))
      write (*,102) (s1d  (k1),k1=k,min(kdm,k+9))
    end do
  end if
 99   format (i8,a,'  (5-line groups: dpth,rho,targt,T,S)')
100   format (11f7.1)
101   format (4x,10f7.2)
102   format (4x,10f7.2)

  torho=0.
  totem=0.
  tosal=0.
  do k=1,kdm
    torho=torho+r1d(k)*(p1d(k+1)-p1d(k))
    totem=totem+t1d(k)*(p1d(k+1)-p1d(k))
    tosal=tosal+s1d(k)*(p1d(k+1)-p1d(k))
  end do

  ko=1
  do k=2,kdm
    if (p1d(k).lt.p1d(kdm+1)) ko=k
  end do

! --- step 1: move interfaces to reduce deviations from target.
  if (vrbos)							&
    write (*,'(i8,i3,a/(i8,4f8.3))') i,ko,			&
     ' temp    salt    dens    pres',(k,t1d(k),s1d(k),r1d(k),	&
     p1d(k+1),k=1,kdm)

  call restp_1d(r1d,p1d,kdm,r1d,pnew,targt,kdm,vrbos,i)

! --- step 2: enforce minimum layer thickness constraint

  dpsum=0.
  last=1
  do 7 k=1,kdm-1

! --- set lower limits for layer thknss (dp0) and depth of lower intfc (dpsum)
    dp0=dpmin(k)
!!    if (k.gt.1) then
!!     if (k.gt.last+1) dp0=dp0 * 2.**(.5*(last-k))
!!     if (vrbos .and. k.eq.last+2) print '(2i5,a,i3)',i,
!!   .    '  hybrid layers: k = 1 -',last 
!!    end if

! --- reduce layer thickness in shallow spots, creating sigma coord. effect
    if (4*k.lt.kdm) dp0=dp0*min(1.,pbot0/200.+.4)

    dpsum=dpsum+dp0
    dp1=dp0
    if (k.gt.1) dp1=cushn(max(pnew(k+1),dpsum)-pnew(k),dp0)
    if (vrbos)							&
     print '(i3,a,3f8.3)',k,					&
       ' min.thknss (dp0;dp1) & min.lowr intfc.dpth',		&
        dp0,dp1,dpsum
    dp0=dp1

! --- is lower intfc too close to the surface?
    p_hat=max(dpsum,pnew(k)+dp0)
    if (p_hat.gt.pnew(k+1)) then
      old=pnew(k+1)
      if (k-1.eq.last) last=k			!  still in z domain
      pnew(k+1)=min(pnew(kdm+1),pnew(k)+max(0.001,p_hat-pnew(k)))

      if (vrbos) then
        print '(13x,a/i13,5f11.3)',				&
         '   pold(k+1)    p_hat   pnew(kdm+1)  pnew(k+1)',	&
          k,old,p_hat,pnew(kdm+1),pnew(k+1)
        write (*,107) k,' lower intfc',p1d(k+1),'=>',pnew(k+1)
 107    format (i3,1x,2(1x,a,f9.3))
      end if
    end if
 7  continue

! --- vertical advection

    if (vrbos) then
     print '(i8,a/(i3,2f15.2))',i,				&
        '  old & new pressure used in advection:',		&
         (k,p1d(k),pnew(k),k=1,kdm+1)
     print *,'vertical advection of temperature'
    end if
    call ppmadv(p1d,t1d,pnew,tnew,kdm,vrbos,i)
    if (vrbos) print *,'vertical advection of salinity'
    call ppmadv(p1d,s1d,pnew,snew,kdm,vrbos,i)

    do 6 k=1,kdm
    t1d(k)=tnew(k)
    s1d(k)=snew(k)
 6  r1d(k)=sig(t1d(k),s1d(k))

    tndrho=-torho
    tndtem=-totem
    tndsal=-tosal
    do k=1,kdm
      tndrho=tndrho+r1d(k)*(pnew(k+1)-pnew(k))
      tndtem=tndtem+t1d(k)*(pnew(k+1)-pnew(k))
      tndsal=tndsal+s1d(k)*(pnew(k+1)-pnew(k))
    end do
    if (abs(tndtem).gt.acurcy*30.*p1d(kdm+1))			&
       write (*,104) '  hybgn2 - bad temp.intgl.:',totem,	&
        tndtem,tndtem/(10.*p1d(kdm+1))
    if (abs(tndsal).gt.acurcy*40.*p1d(kdm+1))			&
       write (*,104) '  hybgn2 - bad saln.intgl.:',tosal,	&
        tndsal,tndsal/(35.*p1d(kdm+1))
 104  format (a,2es15.7,es9.1)

    if (vrbos) then
      write (*,99) i,'  hybgn2   n e w   profile'
      do k=1,kdm,10
        write (*,100) (pnew (k1),k1=k,min(kdm+1,k+10))
        write (*,101) (r1d  (k1),k1=k,min(kdm,k+9))
        write (*,101) (targt(k1),k1=k,min(kdm,k+9))
        write (*,102) (t1d  (k1),k1=k,min(kdm,k+9))
        write (*,102) (s1d  (k1),k1=k,min(kdm,k+9))
      end do
    end if

! --- reduce salt content wherever density > targt(kdm)

!     do 31 k=1,kdm
!     if (pnew(k).lt.p1d(kdm+1)-.001) then
!      q=sofsig(targt(kdm),t1d(k))
!      if (q.lt.s1d(k)) then
!       if (vrbos .or. s1d(k)-q.gt..2)
!    .   write (*,'(2i5,i3,4(a,f7.3))') i,j,k,' reduce s=',s1d(k),
!    .   ' ->',q,' to reduce rho=',r1d(k),' ->',targt(kdm)
!       s1d(k)=q
!       r1d(k)=targt(kdm)
!      end if
!     end if
!31   continue

    if (vrbos) then
      write (*,'(i8,a)') i,' hybgn2 OUT:  temp    saln    dens   thkns    dpth'
      do k=1,kdm
        write (*,103) k,t1d(k),s1d(k),r1d(k),pnew(k+1)-pnew(k),pnew(k+1)
      end do
    end if

  return
  end subroutine hybgn2_1d

 
    subroutine ppmadv(xold,fldold,xnew,fldnew,kk,vrbos,i)

!-- PPM-based 1-dim transport routine, extracted from HYCOM's fct3d.f.
!-- Note: |xold-xnew| can exceed cell width (i.e., no CFL constraints)

!-- xold/new	- old/new cell boundaries
!-- fldold/new	- mixing ratio of dep.variable before and after transport
!-- kk		- number of layers
!-- vrbos	- if .true., print diagnostic messages for grid point -i,j-

      implicit none
      integer,intent(IN)    :: kk,i
      real   ,intent(IN)    :: xold(kk+1),xnew(kk+1),fldold(kk)
      real   ,intent(OUT)   :: fldnew(kk)
      logical,intent(IN)    :: vrbos	!  switch for 'verbose' mode

      real zold(kk+1),znew(kk+1),delx(kk+1),delz(kk+1),fco(kk),fcn(kk),	&
           vertfx(kk+1),vertdv(kk)
      real a(kk),b(kk),c(kk),dx,fcdx,yl,yr
      real amount,bfore,after,dpth,scale,slab,dslab
      integer k,lyr
      real,parameter :: athird=1./3.
      real,parameter :: small=1.e-11
      real,parameter :: acurcy=1.e-6

      delx(:)=xnew(:)-xold(:)
!-- make sure -x- increases with -k-. change sign of -xold/new- if necessary
      if (xold(1).lt.xold(kk+1)) then
        zold(:)=xold(:)
        znew(:)=xnew(:)
      else
        zold(:)=-xold(:)
        znew(:)=-xnew(:)
      end if
      delz(:)=znew(:)-zold(:)

      if (vrbos) then
       write (*,100) i,'entering ppmadv: old_x       d(x)   variable',	&
        (k,xold(k),delx(k),fldold(k),k=1,kk),kk+1,xold(kk+1),delx(kk+1)
 100   format (i8,3x,a/(i21,2f11.1,f11.3))
       write (*,102) i,'  list of -xnew- values:',xnew
 102   format (i8,a/(10f8.2))
      end if

!-- deduce old and new cell width from -zold,znew-
      do 15 k=1,kk
      fco(k)=max(0.,zold(k+1)-zold(k))
15    fcn(k)=max(0.,znew(k+1)-znew(k))

      bfore=0.
      scale=0.
      dpth=0.
      do k=1,kk
        bfore=bfore+fldold(k)*fco(k)
        dpth=dpth+fco(k)
        scale=max(scale,abs(fldold(k)))
      end do
      fldnew(:)=fldold(:)

!-- start by filling zero-width cells with data from neighboring cells

      do 17 k=kk-1,1,-1
 17   fldnew(k)=(fldnew(k)*fco(k)+fldnew(k+1)*small)	&
               /(          fco(k)+            small)
      do 18 k=2,kk
 18   fldnew(k)=(fldnew(k)*fco(k)+fldnew(k-1)*small)	&
               /(          fco(k)+            small)

!-- fit 0th, 1st, or 2nd deg. polynomial to -fldnew- in each cell
      a(1 )=fldnew(1 )
      b(1 )=0.
      c(1 )=0.
      a(kk)=fldnew(kk)
      b(kk)=0.
      c(kk)=0.

      do 16 k=2,kk-1
!-- uncomment one of the following 3 options to activate pcm,plm,ppm resp.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- piecewise constant method:
!!!      a(k)=fldnew(k)
!!!      b(k)=0.
!!!      c(k)=0.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- piecewise linear method:
!-- fit linear function a+bx to tracer in each cell (-.5 < x < +.5)
!!!      a(k)=fldnew(k)
!!!      b(k)=0.
!!!      if (fldnew(k).le.min(fldnew(k-1),fldnew(k+1)) .or.		&
!!!          fldnew(k).ge.max(fldnew(k-1),fldnew(k+1))) then
!!!        b(k)=0.
!!!      else if ((fldnew(k+1)-fldnew(k-1))*(fldnew(k-1)+fldnew(k+1)	&
!!!        -2.*fldnew(k)).gt.0.) then
!!!        b(k)=fldnew(k)-fldnew(k-1)
!!!      else
!!!        b(k)=fldnew(k+1)-fldnew(k)
!!!      end if
!!!      c(k)=0.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- piecewise parabolic method:
!-- construct parabola  a+bx+cx^2  whose integral over [-.5,+.5] equals
!-- fldnew(k) and which passes though points yl,yr at [-.5,+.5] resp.
!!      yl=.5*(fldnew(k-1)+fldnew(k))
!!      yr=.5*(fldnew(k+1)+fldnew(k))
      yl=(max(small,fco(k-1))*fldnew(k)+max(small,fco(k))*fldnew(k-1))/	&
         (max(small,fco(k-1))          +max(small,fco(k)))
      yr=(max(small,fco(k+1))*fldnew(k)+max(small,fco(k))*fldnew(k+1))/	&
         (max(small,fco(k+1))          +max(small,fco(k)))
      a(k)=1.5*fldnew(k)-.25*(yl+yr)
      b(k)=yr-yl
      c(k)=6.*(.5*(yl+yr)-fldnew(k))
      if (abs(yr-yl) .lt. 6.*abs(.5*(yl+yr)-fldnew(k))) then
!-- apex of parabola lies inside interval [-.5,+.5], implying an over-
!-- or undershoot situation. change curve to prevent over/undershoots.
        if (abs(yr-yl) .gt. 2.*abs(.5*(yl+yr)-fldnew(k))) then
!-- put apex of parabola on edge of interval [-.5,+.5]
          if ((yr-yl)*(.5*(yl+yr)-fldnew(k)) .gt. 0.) then
!-- apex at x=-.5
            a(k)=.25*(3.*fldnew(k)+yl)
            c(k)=3.*(fldnew(k)-yl)
            b(k)=c(k)
          else
!-- apex at x=+.5
            a(k)=.25*(3.*fldnew(k)+yr)
            c(k)=3.*(fldnew(k)-yr)
            b(k)=-c(k)
          end if
        else			!  -1/6 < x < +1/6
!-- moving apex won't help. replace parabola by constant.
          a(k)=fldnew(k)
          b(k)=0.
          c(k)=0.
        end if
      end if
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 16   continue

      if (vrbos) print '(a,i12,a/(i20,2x,3es12.3))','i=',i,	&
         '  ppm coeff.  a           b           c',		&
          (k,a(k),b(k),c(k),k=1,kk)

!-- get flux by summing -fldnew- over upstream slab of thickness -delz-

      do 22 k=2,kk
      slab=0.
      amount=0.
      vertfx(k)=0.
      if (delz(k).gt.0.) then			! interface moves in +k dir.
        lyr=k-1
24      lyr=lyr+1
        if (slab.ge.delz(k)) goto 23
        if (fco(lyr).gt.0.) then
          dslab=min(slab+fco(lyr), delz(k))	&
               -min(slab         , delz(k))
          dx=dslab/fco(lyr)
          fcdx=a(lyr)			&
              +b(lyr)*.5*(dx-1.)	&	!  not needed in pcm
              +c(lyr)*(.25-dx*(.5-dx*athird))	!  not needed in pcm,plm
          amount=amount+fcdx*dslab
          slab=slab+dslab
        end if
        if (lyr.lt.kk) go to 24
      else if (delz(k).lt.0.) then		! interface moves in -k dir.
        lyr=k
25      lyr=lyr-1
        if (slab.ge.-delz(k)) goto 23
        if (fco(lyr).gt.0.) then
          dslab=min(slab+fco(lyr),-delz(k))	&
               -min(slab         ,-delz(k))
          dx=dslab/fco(lyr)
          fcdx=a(lyr)			&
              +b(lyr)*.5*(1.-dx)	&	!  not needed in pcm
              +c(lyr)*(.25-dx*(.5-dx*athird))	!  not needed in pcm,plm
          amount=amount+fcdx*dslab
          slab=slab+dslab
        end if					! delz < or > 0
        if (lyr.gt.1) go to 25
      end if
23    if (slab.ne.0.) vertfx(k)=-delz(k)*amount/slab
22    continue

      vertfx(   1)=0.			!  don't allow flux through lower bdry
      vertfx(kk+1)=0.			!  don't allow flux through upper bdry
      do 26 k=1,kk
26    vertdv(k)=vertfx(k+1)-vertfx(k)

      	if (vrbos)							&
         write (*,'(a/(i3,4es12.3))')					&
          'ppmadv:   flux  flx.div/thk    old_thk     new_thk',		&
          (k,vertfx(k),vertdv(k)/max(small,fcn(k)),fco(k),fcn(k),	&
          k=1,kk),kk+1,vertfx(kk+1)

      do 4 k=1,kk
      if (fcn(k).gt.0.) then
       amount=fldnew(k)*fco(k)-vertdv(k)
!      if (abs(amount).lt.100.*small*abs(vertdv(k))) amount=0.
       fldnew(k)=(fldnew(k)*small+amount)/(small+fcn(k))
      end if
 4    continue

      after=0.
      do k=1,kk
        after=after+fldnew(k)*fcn(k)
      end do

      if (abs(bfore-after).gt.acurcy*scale*dpth)			&
       write (*,104) i,'ppmadv - bad column intgl.:',bfore,after
104    format (i8,2x,a,2es15.7)

      if (vrbos)							&
       write (*,100) i,'exiting ppmadv:   d(x)      new_x   variable',	&
        (k,delx(k),xnew(k),fldnew(k),k=1,kk),kk+1,delx(kk+1),xnew(kk+1)

      return
      end subroutine ppmadv


      subroutine pr_9x9(array,idm,jdm,iz,jz,offset,scale,what)
!
! --- (sames as prt9x9 except that array dimensions are passed as arguments)
!
! --- write 9 x 9 point cluster of 'array' values centered on (iz,jz).
! --- the printed numbers actually represent (array(i,j) + offset) * scale
!
      implicit none
      integer,intent(IN) :: idm,jdm,iz,jz
      real   ,intent(IN) :: array(idm,jdm),scale,offset
      character what*(*),text*12
      integer i,j,jwrap
      jwrap(j)=mod(j-1+jdm,jdm)+1               !  for use in cyclic domain
!
      text='            '
      text(1:min(12,len_trim(what)))=what(1:min(12,len_trim(what)))
      write (*,100) text,(jwrap(j),j=jz-4,jz+4)
      write (*,101) (i,(scale*(array(i,jwrap(j))+offset),	&
         j=jz-4,jz+4),i=max(1,iz-4),min(idm,iz+4))
 100  format(a12,9i7)
 101  format(i10,3x,9f7.1)
!
      return
      end subroutine pr_9x9

      subroutine amend2(data1,data2,nreps,idm,jdm,ii,jj,	&
        minum,work1,work2,scal)
!
! --- given an array incompletely filled with grid point values representing
! --- grid-square averages of raw observations, use an objective analysis
! --- approach to (a) fill in missing grid point values and (b) spatially
! --- average grid point values that are based on too few observations.
!
! --- i n p u t   variables:
!
! --- data  : the array of grid point values to be amended and smoothed.
! --- nreps : the number of raw observations individual grid point values are
! ---         based on; a zero indicates that the corresponding value in 'data'
! ---         is missing; points to be ignored are identified by neg. values.
! --- idm   : 1st dimension of 'data' and 'nreps' in the calling program.
! --- jdm   : 2nd dimension of the input data fields.
! --- ii,jj : the number of grid points in i,j direction to be worked on.
! --- minum : the minimum number of observations needed to obtain a 'reliable' 
! ---         grid point average.
! --- work  : work space, same dimensions as 'data' and 'nreps'.
! --- scal  : mesh size in N-S direction
!
      implicit none
! --- o u t p u t  : the amended and smoothed field is returned in 'data'
      integer,parameter :: nsrch=16
! --- nsrch = radius of influence (in grid units) searched for data to be used
! --- in averaging process.
      integer,intent(IN) :: minum,idm,jdm,ii,jj,nreps(idm,jdm)
      real,   intent(IN) :: scal(jdm)
      real,intent(INOUT) :: data1(idm,jdm),data2(idm,jdm)
      real,  intent(OUT) :: work1(idm,jdm),work2(idm,jdm)
      integer jcyc,i,j,i1,j1,i2,j2,imin,jmin,m,n,n4,nrp
      real dist(nsrch+1,nsrch+1),wfunc,d,scx,scy,value1,value2,	&
        sum,dmin,wgt,weight1,scalx(idm,jdm),scaly(idm,jdm)
      real,parameter :: huge=1.e33
!
! --- cyclic function in zonal direction
!
      jcyc(n)=mod(n+jdm-1,jdm)+1
!
! --- objective analysis weight function
      wfunc(d)=exp(-d/9.)
!
! --- map scale factors in i,j direction for lateral distance calculation
!
!$OMP PARALLEL DO
      do 8 j=1,jj
      do 8 i=1,ii
      scalx(i,j)=1.
      scaly(i,j)=scal(j)
      work1(i,j)=data1(i,j)
      work2(i,j)=data2(i,j)
 8    continue
!$OMP END PARALLEL DO
!
! --- main analysis loop
!
!$OMP PARALLEL DO PRIVATE(scx,scy,value1,value2,sum,dmin,imin,jmin,n,n4,wgt,weight1,nrp,dist)
      do  1 j=1,jj
      do 11 i=1,ii
!
! --- don't modify data points based on more than -minum- observations
      if (nreps(i,j).lt.0.or.nreps(i,j).ge.minum) go to 11
!
! --- sort neighboring points by distance
      scx=scalx(i,j)
      scy=scaly(i,j)
!
      do 2 j1=0,nsrch
      do 2 i1=0,nsrch
 2    dist(i1+1,j1+1)=(scx*i1)**2+(scy*j1)**2
!
      nrp=0
      value1=0.
      value2=0.
      sum=0.
!
      do 4 m=1,(nsrch+1)**2
!
      dmin=huge
      do 3 j1=0,nsrch
      do 3 i1=0,nsrch
      if (dmin.gt.dist(i1+1,j1+1)) then
        imin=i1
        jmin=j1
        dmin=dist(i1+1,j1+1)
      end if
 3    continue
!
! --- (imin,jmin) is the closest point not yet used
      dist(imin+1,jmin+1)=huge          
!
! --- each (imin,jmin) combination yields up to 4 data points.
      n4=0
      do 51 i2=i-imin,i+imin,max(1,2*imin)
      if (i2.lt.1.or.i2.gt.ii) go to 51
      do 52 j2=j-jmin,j+jmin,max(1,2*jmin)
! --- jj>jdm indicates that both user and input grids are global in j direction
      if (jj.le.jdm.and.(j2.lt.1.or.j2.gt.jj)) go to 52
      if (nreps(i2,jcyc(j2)).gt.0) n4=n4+nreps(i2,jcyc(j2))
 52   continue
 51   continue
      if (n4.le.0) go to 4
      wgt=wfunc(dmin/(scx*scx+scy*scy))
!
      do 61 i2=i-imin,i+imin,max(1,2*imin)
      if (i2.lt.1.or.i2.gt.ii) go to 61
      do 62 j2=j-jmin,j+jmin,max(1,2*jmin)
! --- jj>jdm indicates that both user and input grids are global in j direction
      if (jj.le.jdm.and.(j2.lt.1.or.j2.gt.jj)) go to 62
      n=nreps(i2,jcyc(j2))
      if (n.le.0) go to 62
!
! --- compute weighted average of original grid point values
! --- weights depend on distance and number of reports at each grid point
!cc      weight=wgt*sqrt(float(n))
      weight1=wgt*float(n)
      if (nrp+n4.gt.minum) weight1=weight1*float(minum-nrp)/float(n4)
      sum=sum+weight1
!diag if (iabs(i-itest)+iabs(j-jtest).le.1) print '('' i,j='',2i3,
!diag.   ''  i2,j2='',2i3,''  nreps,wgt,value='',i5,f9.3,1pe13.3)',
!diag.   i,j,i2,j2,nreps(i2,j2),weight1,work1(i2,j2)
      value1=value1+work1(i2,jcyc(j2))*weight1
      value2=value2+work2(i2,jcyc(j2))*weight1
 62   continue
 61   continue
!
! --- stop once enough reports have been found
      nrp=nrp+n4
      if (nrp.ge.minum) go to 7
 4    continue
!c      if (nrp.le.0)
!c     . print '('' search around i,j='',2i3,'' yields only'',i3,
!c     . '' reports. increase search radius.'')',i,j,nrp
!c      if (nrp.le.0) stop 'error in amend'
!
 7    if (sum.ne.0.) data1(i,j)=value1/sum
      if (sum.ne.0.) data2(i,j)=value2/sum
!diag if (iabs(i-itest)+iabs(j-jtest).le.1) print '('' i,j='',2i3,
!diag.   ''   final value:''f9.3)',i,j,data1(i,j)
!
 11   continue
 1    continue
!$OMP END PARALLEL DO
!
      return
      end subroutine amend2

end module hycom_init
