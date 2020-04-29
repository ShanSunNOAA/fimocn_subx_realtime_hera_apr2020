module hycom_run
use stencilprint
use findmaxmin1
implicit none
contains

  subroutine hycom(its)
  use module_control,  only: nip,dt,numphr,ArchvStep0
  use fimnamelist,     only: PrintIpnDiag,ArchvTimeUnit,stencl_frst,	&
			     stencl_last,stencl_step,			&
  			     kdm,itest,test_start,diag_intvl,ann_core,	&
			     ocnonly,atmonly,coupled,do_rnoff,do_radcor,&
			     do_pcpcor
  use hycom_control   ,only: diag_sst
  use module_constants,only: area,rarea,grvity,deg_lat,deg_lon,corio
  use module_variables,only: us3d,vs3d
  use module_sfc_variables,only:	&	! air-sea coupling data
     ts2d	&! skin temperature (K, blend of SST and ice temp)
    ,hf2d	&! sensible heat flux (W/m^2)
    ,qf2d	&! latent heat flux (W/m^2)
    ,rsds	&! downward shortwave radiation (W/m^2)
    ,rsus	&! upward shortwave radiation (W/m^2)
    ,rlds	&! downward longwave radiation (W/m^2)
    ,rlus	&! upward longwave radiation (W/m^2)
    ,rain	&! precip in last time step (m/timestep)
    ,runoff	&! runoff in last time step (m/timestep)
    ,strs2d	&! surface stress (N/m^2)
    ,us2d	&! u-star in atmosphere (m/s)
    ,t2m2d	&! 2m temperature (K)
    ,q2m2d	&! 2m vapor mixing ratio (kg/kg)
    ,slmsk2d	&! sea-land-ice mask (0-1-2)
    ,hice2d	&! ice thickness (m, grid cell average)
    ,fice2d	 ! fractional ice coverage (rel.units)
  use gfs_physics_internal_state_mod, only: gfs_physics_internal_state,	&
                                            gis_phy
  use hycom_variables, only: nstep_ocn,pbot,gzbot,utrop,vtrop,ptrop,    &
                             uclin,vclin,pres,dp,geop,montg,		&
                             temp,saln,dens,spcifv,			&
                             um_edg,vm_edg,un_edg,vn_edg,pr_edg, 	&
                             dp_edg,uvsq_edg,geop_edg,mont_edg,		&
                             uvsq_edg,spcv_edg,passv_tr,ubforc,vbforc,	&
                             psikk,spvkk,massflx,btropfx,		&
                             wspd,srfx,ssht,pmne,taux,tauy,		&
                             ustar,ustarb,cumuflx,dpinit,dpfinl,	&
                             srfx_bcl,sw2d_bcl,pmne_bcl,ustar_bcl,	&
                             ustarb_bcl,rivflo,hmixl,jerlov,		&
                             covice,thkice,thksno,temice,curl,prcp,	&
                             curl_ave,mssfx_ave,srfxcum,pmnecum,	&
                             radcor,pcpcor
  use hycom_constants, only: odepth,batrop,dtleap,onem,thref,wet,	&
			     ArchvStepOcn,TotalStepOcn,theta,nuphill,	&
			     uphill,dnhill,spcifh,evaplh,KelvC,airdns,	&
                             rho,rhoice,tmin0,tmax0,smin0,smax0
  use hycom_control,   only: subcycle_ocn,modesplit,buoyfor,iocnmx,	&
                             bclin_frq

  use hycom_diag         ,only: transec,glbdia, slabsum
  use hycom_hystat       ,only: hystat
  use hycom_edgvar       ,only: bclin_edgvar,uvedge
  use hycom_output       ,only: output
  use hycom_cnuity       ,only: cnuity
  use hycom_hybgn1       ,only: hybgn1
  use hycom_momtum       ,only: momtum
! use hycom_mxlayr			! kraus-turner
! use hycom_mxlyr1			! kraus-turner
! use hycom_mxlyr2			! kraus-turner
! use hycom_maintn			! kraus-turner
  use hycom_mxkpp        ,only: mxkpp
  use hycom_diamix       ,only: diamix
  use hycom_enloan       ,only: enloan
  use hycom_transp3d     ,only: transp0,transp1,transp2
  use hycom_convec       ,only: convec
  use hycom_thermf       ,only: thermf_core,cd_coare
  use hycom_sverbal      ,only: getcurl
  use hycom_sigetc       ,only: qsatur
  use hycom_barotp       ,only: barotp
! use ersatzpipe

#include <gptl.inc>

  integer,intent(IN) :: its		! atmospheric time step
!SMS$DISTRIBUTE  (dh,1) BEGIN
  real outflo(nip)			! river water sent downstream
!SMS$DISTRIBUTE  END
!SMS$DISTRIBUTE  (dh,2) BEGIN
  real tsrflx(kdm,nip)			! temp subjected only to srf heat flux
!SMS$DISTRIBUTE  END
  integer,parameter :: ncalls=20
  real*8  valmin,valmax
  real*8  rivinp,rivout,rvchng,rivtot
  integer i,leapm,leapn,ix,k,nsubr,nsub,edg,itest_sav,			&
          atm_frst,atm_last,atm_step,atm_Diag
  character name(ncalls)*12,string*80
  real           :: udiff,vdiff,c_d,time,sqrtdns
  real,external  :: its2time

  logical :: vrbos,do_baclin
  integer :: ret				! return code from GPTL routines
  logical,parameter :: chksum=.false.		! evaluate global T/S integrals
  integer :: start_pipe				! time step for starting pipe

!  call StartTimer(t0)
  ret = gptlstart ('iHycom') 

! start_pipe=numphr*24+1	! start pipe output after 24 hrs
  start_pipe=huge(0)            ! no pipe output
  sqrtdns=sqrt(airdns*.001)	! lahey compiler doesn't allow sqrt in initlzn
  time=its2time(its)
  
! --- switch to ocean stencilprint settings (save atmospheric settings)
  itest_sav=itest
  if (itest .gt. 0 .and. its > test_start) then
   atm_frst=stencl_frst
   atm_last=stencl_last
   atm_step=stencl_step
   atm_Diag=PrintIpnDiag
! --- the next 4 lines activate stencl diagnostics at ocean test point
   stencl_frst=1
   stencl_last=kdm
   stencl_step=5
   PrintIpnDiag=itest
   print '(i9,2(a,i8),i8,a/18x,2(a,3i4))',nstep_ocn,                    &
     ' (hycom) replacing PrintIpnDiag=',atm_Diag,'  by:', itest,        &
     perm(itest),' (locl/glob)',                                        &
     'replacing stencl_frst/last/step=',atm_frst,atm_last,atm_step,	&
     '   by:',stencl_frst,stencl_last,stencl_step
  else
   itest=-1
  end if

  do 1 nsub=1,subcycle_ocn

  nstep_ocn=nstep_ocn+1				! main time loop, RH nstep_ocn=0 
  leapm=mod(nstep_ocn  ,2)+1
  leapn=mod(nstep_ocn+1,2)+1
  dtleap=batrop*min(nstep_ocn,2)		! leapfrog time step

  if (nstep_ocn.eq.1)							&
   print '(2(a,i8,4x),a,2i2,4x,a,f9.3)','starting ocn step',		&
   nstep_ocn,'its=',its,'leapm,leapn =',leapm,leapn,'d a y',		&
   nstep_ocn*batrop/(86400.*subcycle_ocn)
  write(string,'(a,i9)') 'step',nstep_ocn

  nsubr=0
! tstop(:)=0.
! call StartTimer(tstart)

  if (mod(nstep_ocn,ArchvStepOcn).eq.0)					&
    call glbdia(nstep_ocn,leapn,trim(string))

  if (ocnonly) then
    if (.not. ann_core) call thermf_core(its)
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
    do i=1,nip
      if (wet(i) > 0 ) then
        rain(i)=prcp(i)*dt
        ssht(i)=(montg(1,i)+thref*ptrop(1,i))/grvity      ! surf.hgt. (m)
      end if
    end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  else if (coupled) then	! coupled run

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (k,c_d,udiff,vdiff)
   do i=1,nip
    if (wet(i) > 0 ) then
     if (do_pcpcor) rain(i)=rain(i)*(1.d0+pcpcor)
     if (do_radcor) then
       srfx(i)=(rsds(i)+rsus(i)+rlds(i)+rlus(i))*(1.d0+radcor)		&
        +hf2d(i)+qf2d(i)					! pos.down W/m^2
     else
       srfx(i)=rsds(i)+rsus(i)+rlds(i)+rlus(i)+hf2d(i)+qf2d(i)	! pos.down W/m^2
     end if
     wspd(i)=sqrt(us3d(1,i)**2+vs3d(1,i)**2)		! surf.wind (m/s)
     ustar(i)=us2d(i)*sqrtdns*(1.-.9*covice(i))		! friction velocity m/s
     ssht(i)=(montg(1,i)+thref*ptrop(1,i))/grvity	! srf.hgt. (m)

! --- stress calculation:
!!   udiff=us3d(1,i)-uclin(1,i,leapm)-utrop(leapm,i)	! atm minus ocn velo
!!   vdiff=vs3d(1,i)-vclin(1,i,leapm)-vtrop(leapm,i)	! atm minus ocn velo
     udiff=us3d(1,i)					! disregard ocean velo
     vdiff=vs3d(1,i)					! disregard ocean velo
     wspd(i)=sqrt(udiff*udiff+vdiff*vdiff)
     if (wspd(i).gt.1.e-5) then				! scaled by ice coverage
       taux(i)=strs2d(i)*us3d(1,i)/wspd(i)*(1.-.9*covice(i))	! eastwd stress
       tauy(i)=strs2d(i)*vs3d(1,i)/wspd(i)*(1.-.9*covice(i))	! northwd stress
     else 
       taux(i)=0.
       tauy(i)=0.
     end if
!!   ustar(i)=sqrt(thref*strs2d(i))			! friction velo
!! --- discard stress supplied by FIM, use hycom's bulk formula instead
!!    c_d = max(0.,cd_coare(wspd(i),.8*qsatur(ts2d(i)-KelvC),		&
!!    ts2d(i)-KelvC,temp(1,i)))
!!    taux(i)=airdns*c_d*wspd(i)*udiff*(1.-.9*covice(i))
!!    tauy(i)=airdns*c_d*wspd(i)*vdiff*(1.-.9*covice(i))
!!    ustar(i)=sqrt(thref*sqrt(taux(i)**2+tauy(i)**2))
! 
! --- allow seaice to intercept precip
     if (covice(i) > 0.) then
       thkice(i)=thkice(i)+rain(i)*covice(i)*rho/rhoice	! 'rain' in m/timestep
       rain(i)=rain(i)*(1.-covice(i))
     end if

     do k=1,kdm
      if (abs(temp(k,i)-.5*(tmax0+tmin0)).ge..5*(tmax0-tmin0)) then
        write(*,'(i7,a,es9.2,a,i8,i3,a,2f7.1)')	nstep_ocn,		&
        ' chk bad temp ',temp(k,i),' at i=',perm(i),k,' lat/lon=',	&
        deg_lat(i),deg_lon(i)
        temp(k,i)=min(max(temp(k,i),tmin0+1.),tmax0-1.)
      end if
      if (abs(saln(k,i)-.5*(smax0+smin0)).ge..5*(smax0-smin0)) then
        write(*,'(i7,a,es9.2,a,i8,i3,a,2f7.1)')	nstep_ocn,		&
        ' chk bad saln ',saln(k,i),' at i=',perm(i),k,' lat/lon=',	&
        deg_lat(i),deg_lon(i)
        saln(k,i)=min(max(saln(k,i),smin0+1.),smax0-1.)
      end if
     end do
    end if			! ocean point
    pmne(i)=rain(i)/dt+qf2d(i)*thref/evaplh	! pos. down m/s; on land & ocean
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! --- transport river water downstream

   if (do_rnoff) then
    nsubr=min(ncalls,nsubr+1)
    name(nsubr)='hycom_rivers'
    ret = gptlstart (name(nsubr))
!   call stencl(rivflo,1,-1.,'rivflo')
    if (buoyfor) then
     rivinp=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO REDUCTION(+:rivinp)
     do i=1,nip
      outflo(i)=0.
      if (wet(i) == 0 ) then
! --- pmne on land or runoff goes to rivflo: use pmne for better conservation
       rivflo(i)=rivflo(i)+max(0.,pmne(i)*dt*area(i))		! m^3/timestep
       rivinp   =rivinp   +max(0.,pmne(i)*dt*area(i))
!!     rivflo(i)=rivflo(i)+        runoff(i)*area(i)		! m^3/timestep
!!     rivinp   =rivinp   +        runoff(i)*area(i) 
       if (dnhill(i).gt.0) then
        outflo(i)=rivflo(i)		! water to be flushed downstream
        rivflo(i)=0.
       end if				! downstream neighbor exists
      end if				! land point
     end do
!$OMP END PARALLEL DO
!SMS$EXCHANGE (outflo)

     rivout=0.
!$OMP PARALLEL DO PRIVATE(ix) REDUCTION(+:rivout)
     do i=1,nip
      if (nuphill(i).gt.0) then		! one or more upstream neighbors exist
       do edg=1,nuphill(i)
        ix=uphill(edg,i)
        rivflo(i)=rivflo(i)+outflo(ix)	! water arriving from upstream
       end do
      end if				! upstream neighbors exist

      if (wet(i) > 0 ) then				! river has reached ocn
       pmne(i)=pmne(i)+rivflo(i)*rarea(i)/batrop		! m/sec
       rivout=rivout+rivflo(i)
       rivflo(i)=0.
      end if				! ocean point
     end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
     if (mod(nstep_ocn,ArchvStepOcn).eq.0) then
!SMS$REDUCE(rivinp,rivout,SUM)
      print '(a,2f8.4)','river input & output (Sv):',rivinp*1.e-6/dt,	&
        rivout*1.e-6/dt
     end if
    end if			! buoyfor
    ret = gptlstop (name(nsubr))
   end if			! do_rnoff
  end if			! ocnonly or coupled

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos)
  do i=1,nip
   if (wet(i) > 0 ) then
! --- build up long-term averages
    srfxcum(i)=srfxcum(i)+srfx(i)
    pmnecum(i)=pmnecum(i)+pmne(i)

    vrbos=i.eq.itest .and. mod(nstep_ocn,diag_intvl).eq.0
    if (vrbos) then
     if (ocnonly) then
      print 101,'i=',perm(i),ArchvTimeUnit,time,'  surface fluxes:',	&
        'heatfx',srfx(i),'swflx',rsds(i),'lwflx',rlds(i),		&
	'P-E.E8',pmne(i)*1.e8,'taux',taux(i)*100,'tauy',tauy(i)*100,	&
	'ustar',ustar(i)*100,'SST',temp(1,i),'temice',temice(i),	&
	'covice',covice(i),'thkice',thkice(i),'hmixl',hmixl(i),		&
        'lat',deg_lat(i),'lon',deg_lon(i)
101  format ("(hycom_run) ",a,i8,2x,a,f8.1,a/(5(a7,"=",f8.3)))
     else
      print 102,'i=',perm(i),ArchvTimeUnit,time,'  surface fluxes:',	&
        'heatfx',srfx(i),'rsds',rsds(i),'rsus',rsus(i),			&
        'rlds',rlds(i),'rlus',rlus(i),'hf2d',hf2d(i),'qf2d',qf2d(i),	&
        'P-E.E8',pmne(i)*1.e8,'taux',taux(i)*100,'tauy',tauy(i)*100,	&
        'ustar',ustar(i)*100,'tair',ts2d(i)-KelvC,'SST',temp(1,i),	&
        'temice',temice(i),'covice',covice(i),'thkice',thkice(i)
102  format ("(hycom_run) ",a,i8,2x,a,f8.1,a/(5(a7,"=",f8.3)))
     endif	! if ocnonly
    end if	! if vrbos
   end if	! ocean point
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  if (mod(nstep_ocn,ArchvStepOcn).eq.0) then
   call findmxmn1(taux,nip,'(hycom_run) taux',wet)
   call findmxmn1(tauy,nip,'(hycom_run) tauy',wet)
   call findmxmn1(srfx,nip,'(hycom_run) srfx',wet)
  end if

! call IncrementTimer(tstart,tstop(nsubr))

! if (mod(nstep_ocn,diag_intvl).eq.0 .and. itest.gt.0)	then
!  call findmxmn1(taux,nip,'taux',wet)
!  call findmxmn1(tauy,nip,'tauy',wet)
!  call findmxmn1(srfx,nip,'srfx',wet)
!  call findmxmn1(rsds,nip,'rsds',wet)
!  call findmxmn1(rlds,nip,'rlds',wet)
!  call findmxmn1(pmne,nip,'pmne',wet)
!  call findmxmn1(us2d,nip,'ustar',wet)
!  call findmxmn1(temice,nip,'temice',wet)
!  call findmxmn1(covice,nip,'covice',wet)
!  call findmxmn1(thkice,nip,'thkice',wet)
!  if (ocnonly) then
!   call findmxmn1(rivflo,nip,'rivflo')
!  else
!   call findmxmn1(t2m2d,nip,'t2m2d',wet)
!   call findmxmn1(q2m2d,nip,'q2m2d',wet)
!   call findmxmn1(hf2d,nip,'hf2d',wet)
!   call findmxmn1(qf2d,nip,'qf2d',wet)
!   call findmxmn1(rain,nip,'rain',wet)
!  endif
! end if

  if (nstep_ocn .eq. ArchvStep0) call output(0,1)

! call IncrementTimer(tstart,tstop(0))

! --- solve hydrostatic equation
  nsubr=min(ncalls,nsubr+1)
  name(nsubr)='hycom_hystat'
  ret = gptlstart (name(nsubr))
  call hystat(nstep_ocn,leapm,						&
    geop,montg,dens,spcifv,temp,saln,dp,pres,ptrop,psikk,spvkk,pbot,gzbot)
  ret = gptlstop (name(nsubr))
! call IncrementTimer(tstart,tstop(nsubr))

! --- determine edge variables
  nsubr=min(ncalls,nsubr+1)
  name(nsubr)='hycom_edgvar'
  ret = gptlstart (name(nsubr))
  call bclin_edgvar(nstep_ocn,leapm,leapn,				&
    uclin,vclin,utrop,vtrop,dp,pres,geop,montg,spcifv,um_edg,vm_edg,	&
    un_edg,vn_edg,dp_edg,pr_edg,geop_edg,mont_edg,uvsq_edg,spcv_edg)
  ret = gptlstop (name(nsubr))
  if (mod(nstep_ocn,diag_intvl).eq.0 .and. itest.gt.0)			&
     call glbdia(nstep_ocn,leapn,name(nsubr))

! --- solve continuity equation
  nsubr=min(ncalls,nsubr+1)
  name(nsubr)='hycom_cnuity'
  ret = gptlstart (name(nsubr))
! call prt_column(nstep_ocn,'bef cnu')
  call cnuity(nstep_ocn,leapm,leapn,					&
    um_edg,vm_edg,un_edg,vn_edg,dp_edg,dp,pres,massflx,mssfx_ave)
  ret = gptlstop (name(nsubr))
! call IncrementTimer(tstart,tstop(nsubr))
  if (mod(nstep_ocn,diag_intvl).eq.0 .and. itest.gt.0)			&
     call glbdia(nstep_ocn,leapn,name(nsubr))

  do_baclin=.false.
  if (mod(nstep_ocn,bclin_frq).eq.0) do_baclin=.true.

! --- build up time integral of mass fluxes
  if (mod(nstep_ocn,2).eq.0)						&
    call transp1(nstep_ocn,cumuflx,massflx)

  if (do_baclin) then

   if (diag_sst) call sstbud(0,'  initialization',temp,dp,leapn)

! --- solve tracer transport equation on long time step
   nsubr=min(ncalls,nsubr+1)
   name(nsubr)='hyc_transp3d'
   ret = gptlstart (name(nsubr))
   dpfinl(:,:)=dp(:,:,leapn)
   if (chksum) call slabsum(nstep_ocn,leapn,'bfore transp3d')
   call transp2(nstep_ocn,temp,saln,passv_tr,cumuflx,dpinit,dpfinl)
   if (chksum) call slabsum(nstep_ocn,leapn,'after transp3d')
   ret = gptlstop (name(nsubr))
!  call IncrementTimer(tstart,tstop(nsubr))

   if (diag_sst) call sstbud(1,' horiz.advection',temp,dp,leapn)

! --- convective adjustment
   nsubr=min(ncalls,nsubr+1)
   name(nsubr)='hycom_convec'
   ret = gptlstart (name(nsubr))
   call convec(nstep_ocn,leapn,uclin,vclin,dens,dp,pres,temp,saln)
   if (chksum) call slabsum(nstep_ocn,leapn,'after convec')
   ret = gptlstop (name(nsubr))
!  call IncrementTimer(tstart,tstop(nsubr))

   if (diag_sst) call sstbud(2,' convec.adjustmt',temp,dp,leapn)

  end if				! do_baclin

! --- solve baroclinic momentum equation
! call findmxmn1(temp(1,:),nip,'sst bef mom',wet)
! call findmxmn1(saln(1,:),nip,'sss bef mom',wet)
! call findmxmn1(uclin(1,:,leapn),nip,'u bef mom',wet)
! call findmxmn1(vclin(1,:,leapn),nip,'v bef mom',wet)

  nsubr=min(ncalls,nsubr+1)
  name(nsubr)='hycom_momtum'
  ret = gptlstart (name(nsubr))
  call momtum(nstep_ocn,leapm,leapn,					&
    utrop,vtrop,uclin,vclin,um_edg,vm_edg,geop_edg,mont_edg,		&
    uvsq_edg,spcv_edg,taux,tauy,ubforc,vbforc,dp,pres,ustarb)
  ret = gptlstop (name(nsubr))
! call IncrementTimer(tstart,tstop(nsubr))

  if (modesplit) then
! --- solve barotropic continuity and momentum equations
   nsubr=min(ncalls,nsubr+1)
   name(nsubr)='hycom_barotp'
   ret = gptlstart (name(nsubr))
   call barotp(nstep_ocn,leapn,utrop,vtrop,ptrop,btropfx,ubforc,vbforc)
   ret = gptlstop (name(nsubr))
!  call IncrementTimer(tstart,tstop(nsubr))
  end if		! modesplit

  if (buoyfor) then
! --- build up time integrals of surface forcing on long time steps
!SMS$PARALLEL (dh,i) BEGIN
   do i=1,nip
    if (wet(i) > 0 ) then
     srfx_bcl  (i)=srfx_bcl  (i)+        srfx  (i)/bclin_frq
     sw2d_bcl  (i)=sw2d_bcl  (i)+(rsds(i)+rsus(i))/bclin_frq
     pmne_bcl  (i)=pmne_bcl  (i)+        pmne  (i)/bclin_frq
     ustar_bcl (i)=ustar_bcl (i)+        ustar (i)/bclin_frq
     ustarb_bcl(i)=ustarb_bcl(i)+        ustarb(i)/bclin_frq
    end if
   end do
!SMS$PARALLEL END

   if (do_baclin) then

! --- thermodynamic ice model
    nsubr=min(ncalls,nsubr+1)
    name(nsubr)='hycom_enloan'
    ret = gptlstart (name(nsubr))
    call enloan(nstep_ocn,leapn,covice,thkice,thksno,temice,srfx_bcl,	&
      pmne_bcl,temp,saln,dp,pres)
    if (chksum) call slabsum(nstep_ocn,leapn,'after enloan')
    ret = gptlstop (name(nsubr))
!   call IncrementTimer(tstart,tstop(nsubr))

    if (diag_sst) call sstbud(3,' ice melt/freeze',temp,dp,leapn)

    if (diag_sst) then
! --- isolate effect of air-sea fluxes on SST (diagnostic use only)
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos)
     do i=1,nip
      if (wet(i) > 0 ) then
       vrbos=i.eq.itest
       tsrflx(:,i)=temp(:,i)
       tsrflx(1,i)=tsrflx(1,i)+srfx_bcl(i)*batrop*bclin_frq		&
         *grvity/(spcifh*dp(1,i,leapn))
       if (vrbos) print '(i7,a,es11.3,2f9.1,f9.3)',perm(i),		&
         ' heatflx effect on SST:',					&
         srfx_bcl(i)*batrop*bclin_frq*grvity/(spcifh*dp(1,i,leapn)),	&
         srfx_bcl(i),batrop*bclin_frq,dp(1,i,leapn)/onem
      end if
     end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
     call sstbud(4,'  air-sea fluxes',tsrflx,dp,leapn)
    end if		! diag_sst

    nsubr=min(ncalls,nsubr+1)
    name(nsubr)='hycom_mxlayr'
    ret = gptlstart (name(nsubr))
    call mxkpp(nstep_ocn,leapn,theta,dens,dp,pres,temp,saln,		&
      uclin,vclin,srfx_bcl,sw2d_bcl,pmne_bcl,ustar_bcl,ustarb_bcl,	&
      jerlov,deg_lat,corio,hmixl)
    if (chksum) call slabsum(nstep_ocn,leapn,'after mxkpp')
    ret = gptlstop (name(nsubr))
!   call IncrementTimer(tstart,tstop(nsubr))
    if (mod(nstep_ocn,diag_intvl).eq.0 .and. itest.gt.0)		&
       call glbdia(nstep_ocn,leapn,name(nsubr))

   if (diag_sst) call sstbud(5,'turb.ML exchange',temp,dp,leapn)


   end if			! do_baclin
  end if			! buoyfor

  if (do_baclin) then
   nsubr=min(ncalls,nsubr+1)
   name(nsubr)='hycom_hybgn1'
   ret = gptlstart (name(nsubr))
   call hybgn1(nstep_ocn,leapn,theta,dens,dp,pres,temp,saln,		&
               uclin,vclin,deg_lat)
!  call prt_column(nstep_ocn,'aft hybgen')
   if (chksum) call slabsum(nstep_ocn,leapn,'after hybgen')
   ret = gptlstop (name(nsubr))
!  call IncrementTimer(tstart,tstop(nsubr))
   if (mod(nstep_ocn,diag_intvl).eq.0 .and. itest.gt.0)			&
     call glbdia(nstep_ocn,leapn,name(nsubr))

   if (diag_sst) call sstbud(6,' vert.regridding',temp,dp,leapn)

!  call pipe2(nstep_ocn,temp,kdm,'Thybgn1')
!  call pipe2(nstep_ocn,saln,kdm,'Shybgn1')
!  call pipe2(nstep_ocn,uclin(:,:,leapm),kdm,'UMhybgn1')
!  call pipe2(nstep_ocn,uclin(:,:,leapn),kdm,'UNhybgn1')
!  call pipe2(nstep_ocn,vclin(:,:,leapm),kdm,'VMhybgn1')
!  call pipe2(nstep_ocn,vclin(:,:,leapn),kdm,'VNhybgn1')

! --- diapycnal column mixing
   if (iocnmx.eq.2) then
    nsubr=min(ncalls,nsubr+1)
    name(nsubr)='hycom_diamix'
    ret = gptlstart (name(nsubr))
    call diamix(nstep_ocn,leapn,dp,dens,pres,temp,saln,uclin,vclin)
!   call prt_column(nstep_ocn,'aft diapfl')
    if (chksum) call slabsum(nstep_ocn,leapn,'after diamix')
    ret = gptlstop (name(nsubr))
!   call IncrementTimer(tstart,tstop(nsubr))
   else
     stop 'wrong iocnmx'
   end if

   if (diag_sst) call sstbud(7,'   diapyc.mixing',temp,dp,leapn)

! --- re-initialize time integrals for tracer transport
   name(nsubr)='hycom_transp0'
   ret = gptlstart (name(nsubr))
   call transp0(nstep_ocn,leapn,cumuflx,dp,dpinit)
   ret = gptlstop (name(nsubr))

   if (buoyfor) then

!SMS$PARALLEL (dh,i) BEGIN
    do i=1,nip
     if (wet(i) > 0 ) then
      srfx_bcl  (i)=0.
      sw2d_bcl  (i)=0.
      pmne_bcl  (i)=0.
      ustar_bcl (i)=0.
      ustarb_bcl(i)=0.
     end if 
    end do 
!SMS$PARALLEL END
   end if			! buoyfor
  end if			! do_baclin

  if (coupled) then

! --- feed SST back to atmosphere

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
   do i=1,nip
    if (wet(i) > 0 ) then
     if (covice(i) .le. .50) then
      slmsk2d(i)=0.	! water point
     else
      slmsk2d(i)=2.	! ice point
     end if
     ts2d(i)=temp(1,i)*(1.-covice(i))+temice(i)*covice(i) + KelvC	! ts2d = skin temp
     hice2d(i)=thkice(i)
     fice2d(i)=covice(i)
     gis_phy%sfc_fld%TSEA  (i,1) = ts2d    (i)
     gis_phy%sfc_fld%HICE  (i,1) = hice2d  (i)
     gis_phy%sfc_fld%FICE  (i,1) = fice2d  (i)
     gis_phy%sfc_fld%SLMSK (i,1) = slmsk2d (i)
    end if
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

   if (mod(nstep_ocn,diag_intvl).eq.0 .and. itest.gt.0)		&
      call findmxmn1(ts2d,nip,'(hyc_run) ts2d')

  end if			! coupled run

! --- diagnose wind stress curl for sverdrup transport calculation
  call getcurl(nstep_ocn,taux,tauy,curl,curl_ave)
!
! --- archive results at regular intervals

  if (mod(nstep_ocn,ArchvStepOcn).eq.0) then
! --- misc.diagnostics
    nsubr=min(ncalls,nsubr+1)
    name(nsubr)='hycom_glbdia'
    ret = gptlstart (name(nsubr))
    call glbdia(nstep_ocn,leapn,trim(string))
    call transec(nstep_ocn,leapn)
    ret = gptlstop (name(nsubr))

!   call IncrementTimer(tstart,tstop(nsubr))
  end if         ! if (mod(nstep_ocn,ArchvStepOcn).eq.0)

    nsubr=min(ncalls,nsubr+1)
    name(nsubr)='hycom_output'
    ret = gptlstart (name(nsubr))
    if (nstep_ocn .ge. ArchvStep0) call output(nstep_ocn,leapn)
    ret = gptlstop (name(nsubr))
!   print *,'back from hycom_output, step',nstep_ocn

!   call IncrementTimer(tstart,tstop(nsubr))

!! --- timing diagnostics
!!  if (mod(nstep_ocn,diag_intvl).eq.0 .or.				&
!!      mod(nstep_ocn,10*bclin_frq).eq.0) then
!!   do k=1,nsubr
!!   valmin=(tstop(k)-tstop(k-1))*1.e6
!!   valmax=valmin
!!SMS$REDUCE(valmin,MIN)
!!SMS$REDUCE(valmax,MAX)
!!    print '(2a,2(i8,a))','time spent in subr ',name(k),		&
!!      nint(valmin),' -',nint(valmax),'  usecs'
!!   end do
!!  end if
!!  valmin=(tstop(nsubr)-tstop(0))*1.e6
!!  valmax=valmin
!!!SMS$REDUCE(valmin,MIN)
!!!SMS$REDUCE(valmax,MAX)
!!  print '(3(a,i9),a)','ocn step',nstep_ocn,' completed in ',		&
!!    nint(valmin),' -',nint(valmax),'  usecs'

! if (nstep_ocn.ge.start_pipe) then
!   call pipe2(nstep_ocn,temp            ,kdm, 'temp_end')
!   call pipe2(nstep_ocn,saln            ,kdm, 'salt_end')
!   call pipe2(nstep_ocn,   dp(:,:,leapm),kdm,  'dpm_end')
!   call pipe2(nstep_ocn,   dp(:,:,leapn),kdm,  'dpn_end')
!   call pipe2(nstep_ocn,uclin(:,:,leapm),kdm,'uvelm_end')
!   call pipe2(nstep_ocn,vclin(:,:,leapm),kdm,'vvelm_end')
!   call pipe2(nstep_ocn,uclin(:,:,leapn),kdm,'uveln_end')
!   call pipe2(nstep_ocn,vclin(:,:,leapn),kdm,'vveln_end')
! end if

1 continue				! time loop

! --- restore atmospheric stencil settings
  itest=itest_sav
  if (itest > 0 .and. its > test_start) then
   stencl_frst=atm_frst
   stencl_last=atm_last
   stencl_step=atm_step
   PrintIpnDiag=atm_Diag
   print '(i9,a,i8/18x,a,3i4)',nstep_ocn,				&
     ' (hycom) restoring PrintIpnDiag=',atm_Diag,			&
     'restoring stencl_frst/last/step =',stencl_frst,stencl_last,	&
     stencl_step
  end if

! call IncrementTimer(t0,iHycomLoopTime)
  ret = gptlstop ('iHycom')

  return
  end subroutine hycom
end module hycom_run
