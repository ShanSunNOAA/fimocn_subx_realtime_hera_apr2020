module hycom_enloan
contains

!**************************************************************************
!  enloan
!
!  'energy loan' ice model with advection (free drift) but no ice dynamics.
!  ice amount represents energy 'loaned' to water column to prevent
!  wintertime cooling below freezing level. 'loan' is paid back in summer.
!  ice surface temperature calculation based on assumption of zero heat
!  flux divergence (i.e. air-ice flux = ice-water flux).
!
!  S. Sun, R.Bleck	May   2011
!  R.Bleck		March 2013	ice advection by surface flow
!**************************************************************************

  subroutine enloan(nstep,leapn,					&
                    covice,thkice,thksno,temice,surflx,pmevap,		&
                    temp,saln,dp,pres)

  use module_control  ,only: nip,npp
  use module_constants,only: grvity,area,rarea,deg_lat,nprox,prox,	&
                             proxs,sidevec_e,perm
  use fimnamelist    ,only: kdm,itest,diag_intvl,ocnonly,sss_rstor
  use hycom_variables ,only: salnow,ticeol,um_edg,vm_edg,ocnarea,	&
                             qf2d_ave
  use hycom_constants ,only: wet,epsil,tenm,onem,onecm,onemm,onemu,     &
                             dpthmn,saldif,rhoice,spcifh,thref,huuge,	&
			     batrop,land_spval,fusion,thkmx,tfrez,	&
			     KelvC
  use hycom_control   ,only: bclin_frq
  use hycom_dffusn,only: dffusn_lev
  use hycom_diag,only: glob2d
  use findmaxmin1
  use stencilprint
  use stenedgprint

  implicit none
  integer,intent(IN) :: nstep,leapn
!SMS$DISTRIBUTE (dh,1) BEGIN
  real,intent(INOUT) :: covice(nip)	! ice coverage (rel.units)
  real,intent(INOUT) :: thkice(nip) 	! grid-cell averaged ice thickness (m)
  real,intent(INOUT) :: thksno(nip) 	! grid-cell averaged snow thickness (m)
  real,intent(INOUT) :: temice(nip)	! ice surface temperature
  real,intent(INOUT) :: surflx(nip)	! surface net heat flux (W/m^2)
  real,intent(INOUT) :: pmevap(nip)	! freshwater flux, P-E (m/s)
  real    :: secnd(nip),icex(nip),origfx(nip),sxchg(nip),mxchg(nip)
  integer :: mask(nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
  real,intent(INOUT) :: temp(kdm,nip)	! temperature
  real,intent(INOUT) :: saln(kdm,nip)	! salinity
  real,intent(IN)    :: dp  (kdm,nip,2)	! layer thickness
  real,intent(IN)    :: pres(kdm+1,nip)	! interface pressure
  real    :: flux(npp,nip)
!SMS$DISTRIBUTE END
  real*8 tmxl,dpth,paybak,borrow,availb,totice,pp,qq,sold,		&
       twrk,old,tnew,q,top,bot,thkinv,dlt,secmn,secmx,			&
       valu1,valu2,outgo,sum,xcess,dpgt0,prorat,colini,colout,		&
       colinii,coliniw,colouti,coloutw
  real dfflen
  integer i,ix,k,iter,edg,edgx
  logical vrbos
  logical, parameter :: brine_pump=.true.
  character title*24
  real,parameter :: thresh=4.			! diagnostic use
  real,parameter :: acurcy=1.e-5
  real,parameter :: piston=12./(86400.*300.)    ! 300 days over 12m depth scale

! real :: thin   = 0.05		! min. ice thickness (m)
  real :: thin   = 0.15		! min. ice thickness (m)
  real :: kice   = 2.04		! heat conductivity in ice (W/m/deg)
  real :: rate   = 2./86400.	! max. ice melting rate (m/sec)
  real :: stfbol = 5.67e-8	! stefan-boltzmann const. (W/m^2/deg^4)
! real :: brntop = 0.		! top of dpth.intvl for brine deposition
! real :: brnbot = 500.   	! bottom of dpth.intvl brine deposition
! to match hycom
  real :: brntop = 50.		! top of dpth.intvl for brine deposition
  real :: brnbot = 300.   	! bottom of dpth.intvl brine deposition

  dlt=batrop*bclin_frq		! routine-specific time step
  if (mod(nstep,diag_intvl).eq.0) then
   write (title,'(a,i9)') '(enloan) step',nstep
   call stencl(thkice,1,100.,title//'ice thickness (cm) IN')
  end if

! --- allow ice to drift with surface flow

!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE (thkice)
!$OMP PARALLEL DO PRIVATE (vrbos,ix,edgx,valu1,valu2,outgo)
  do i=1,nip
   vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
   flux(:,i)=0.
   if (wet(i) > 0) then
    if (vrbos) then
      write (*,'(i8,a,i8/(i27,3f10.2))') nstep,			&
      ' entering enloan t,s,p @i=',perm(i),			&
      (k,temp(k,i),saln(k,i),pres(k+1,i)/onem,k=1,kdm)
      write (*,'(a,2f12.2)') 'ice thickness(m)/tem=',thkice(i),temice(i)
    end if

    do edg=1,nprox(i)			! loop through edges
     ix=prox(edg,i)
     if (wet(ix) > 0) then
      edgx=proxs(edg,i)
! --- outgo = outging velocity component on edge times edge length
      valu1=sidevec_e(2,edg ,i )*um_edg(1,edg ,i )		&
          - sidevec_e(1,edg ,i )*vm_edg(1,edg ,i )		! m^2/sec
      valu2=sidevec_e(2,edgx,ix)*um_edg(1,edgx,ix)		&
          - sidevec_e(1,edgx,ix)*vm_edg(1,edgx,ix)		! m^2/sec
      outgo=.5*(valu1-valu2)		! avg.to prevent roundoff errors
      if (outgo.ge.0.) then
       flux(edg,i)=thkice(i )*outgo
      else
       flux(edg,i)=thkice(ix)*outgo
      end if
     end if
    end do
   end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,sum,colinii,coliniw,		&
!$OMP dpgt0,sold,colouti,coloutw,colini,colout)
  do i=1,nip
   vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
   if (wet(i) > 0) then
    sum=0.
    do edg=1,nprox(i)			! loop through edges
     sum=sum+flux(edg,i)
    end do
    thkice(i)=thkice(i)-sum*rarea(i)*dlt

! --- convert snow to ice

    if (thksno(i).gt.0.) then
     if (vrbos) then
      colinii=-saldif*thkice(i)*rhoice*grvity			! g m/sec^2
      coliniw=0.
      do k=1,kdm
       coliniw=coliniw+saln(k,i)*dp(k,i,leapn)			! g m/sec^2
      end do
     end if

! --- enloan doesn't discriminate between snow and ice. as part of converting
! --- snow to ice, transfer some salt from ocean to snow to match ice salinity.
! --- (enloan treats ice as neg.salt reservoir => a sign change is involved)
 
! --- virtual salt flux into ocean if snow were to melt right away:
! ---      snomelt=0.
! --- virtual salt flux into ocean if snow turns into ice and melts later:
! ---      icemelt=-saldif*thksno(i)*rhoice
! --- net salt flux accompanying snow => ice conversion: snomelt-icemelt
     sxchg(i)=thksno(i)*saldif*rhoice				! g/m^2
     dpgt0=max(epsil,dp(1,i,leapn))
     sold=saln(1,i)
     saln(1,i)=saln(1,i)+sxchg(i)*grvity/dpgt0			! g/kg
     thkice(i)=thkice(i)+thksno(i)
     thksno(i)=0.

     if (vrbos) then
      print '(a,2i8,a,es11.3)','step',nstep,perm(i),			&
        ' salt needed for snow-to-ice conv (g/m^2):',sxchg(i)
      print 100,nstep,perm(i),'thkold',old,'thknew',thkice(i),		&
        'sold',sold,'snew',saln(1,i),'dplyr1',dp(1,i,leapn)/onem,	&
        'dplyr2',dp(2,i,leapn)/onem
     end if

     if (vrbos) then
      colouti=-saldif*thkice(i)*rhoice*grvity			! g m/sec^2
      coloutw=0.
      do k=1,kdm
       coloutw=coloutw+saln(k,i)*dp(k,i,leapn)			! g m/sec^2
      end do
      colini=colinii+coliniw
      colout=colouti+coloutw
      print '(i8,a,2es15.7)',perm(i),'  (enloan-1) colini,colout=',	&
       colini,colout
     end if
    end if		! thksno > 0
   end if		! ocean point
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  if (mod(nstep,diag_intvl).eq.0) then
   call stenedg(thkice,flux,1,title//'  ice thickness, ice flux')
   call stencl(thkice,1,100.,title//'  ice thickness (cm) after transport')
  end if

! --- energy loan: add extra energy to the ocean to keep SST from dropping
! --- below tfrez in winter. return this borrowed energy to the 'energy bank'
! --- in summer as quickly as surflx > 0 allows.
! --- salt loan: analogous to energy loan.
  
!SMS$PARALLEL (dh,i) BEGIN
!OMP PARALLEL DO
  do i=1,nip
  icex(i)=land_spval
  mask(i)=0
   if (thkice(i) > 0.) mask(i)=1
  end do
!OMP END PARALLEL DO
!SMS$PARALLEL END

  dfflen=0.1*sqrt(5.e13/float(nip))	! SQRT[earth surface/(10*cell count)]
  call dffusn_lev(thkice,1,1,1,dfflen,mod(nstep,diag_intvl).eq.0,mask)

  xcess=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(vrbos,colinii,coliniw,dpth,tmxl,borrow,	&
!$OMP old,bot,top,thkinv,sum,dpgt0,prorat,availb,paybak,		&
!$OMP colouti,coloutw,colini,colout) REDUCTION(+:xcess)
  do 10 i=1,nip
  vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
  if (wet(i) > 0 ) then
!
   origfx(i)=surflx(i)
   colinii=-saldif*thkice(i)*rhoice*grvity			! g m/sec^2
   coliniw=0.
   do k=1,kdm
    coliniw=coliniw+saln(k,i)*dp(k,i,leapn)			! g m/sec^2
   end do

! --- distribute heat flux over shallow layer to estimate its effect on SST
   dpth=max(dp(1,i,leapn),dpthmn*onem)
   tmxl=temp(1,i)+surflx(i)*dlt*grvity/(spcifh*dpth)
   if (tmxl.lt.tfrez(i)) then
!
! --- water  f r e e z e s. borrow energy to raise tmxl<tfrez back to tfrez

    borrow=(tfrez(i)-tmxl)*spcifh*dpth/grvity			! > 0, J/m^2
    old=thkice(i)
    thkice(i)=thkice(i)+borrow/(fusion*rhoice)
    surflx(i)=surflx(i)+borrow/dlt				! W/m^2
    mxchg(i)=-grvity*(borrow/fusion)				! kg/m/sec^2
    sxchg(i)=saldif*(borrow/fusion)				! g/m^2

    qf2d_ave(i)=qf2d_ave(i)+borrow/batrop

    if (vrbos) then
     print '(i8,a,i7,a,f9.4,a)',nstep,' (enloan) i=',perm(i),		&
      '  freezing:',-mxchg(i)/onecm,' cm'
     if (thkice(i).eq.0.) print '(30x,a)','(new ice)'
    end if

    if (brine_pump .and. deg_lat(i).le.0.) then
! --- don't leave brine at the surface. dump it in deeper layers
     bot=min(brnbot*onem,    pres(kdm+1,i))
     top=min(brntop*onem,.75*pres(kdm+1,i))
     thkinv=1./(bot-top)
     sum=0.
     do k=1,kdm
      dpgt0=max(epsil,pres(k+1,i)-pres(k,i))
      prorat=max(0.,min(pres(k+1,i),bot)-max(pres(k,i),top))*thkinv
      sum=sum+prorat
      old=saln(k,i)
      saln(k,i)=saln(k,i)+sxchg(i)*grvity/dpgt0*prorat	! g/kg
 
      if (vrbos) then
       print 102,perm(i),k,'sold',old,'snew',saln(k,i),	&
         'dpgt0',dpgt0/onem,'prorat',prorat
      end if
 
      if (pres(k,i).gt.bot) exit
     end do
     if (abs(sum-1.).gt.acurcy) then
      print '(a,i8,a,f9.6,a)','i=',perm(i),			&
        ' enloan warning: sum(prorat)=',sum,', should be 1'
!     stop '(error)'
     end if
    else	! no brine_pump
     dpgt0=max(epsil,dp(1,i,leapn))
     old=saln(1,i)
     saln(1,i)=saln(1,i)+sxchg(i)*grvity/dpgt0	! g/kg

     if (vrbos) then
      print 102,perm(i),1,'sold',old,'snew',saln(1,i),	&
        'sxchg',sxchg(i),'dpgt0',dpgt0/onem
 102  format (i8,i3,2x,4(a8,'=',f7.3))
     end if
    end if	! brine_pump
 
    if (vrbos) then
     if (old.eq.0. .and. borrow.gt.0.) print '(i9,a,i8)',	&
      nstep,'  new ice forms at i=',i
     print 100,nstep,perm(i),'sst',temp(1,i),'sss',saln(1,i),	&
      'thkice',thkice(i),'loan-6',borrow*1.e-6
 100  format (i8," (enloan) i=",i7/5(a8,'=',f7.3))
    end if
 
   else if (tmxl.gt.tfrez(i) .and. thkice(i).gt.0.) then
 
! --- ice   m e l t s. return the borrowed energy whenever tmxl > tmelt
 
!   availb=(tmxl-tfrez(i))*spcifh*dpth/grvity			!  < 0, J/m^2
    availb=(tmxl-tfrez(i))*spcifh*dp(1,i,leapn)/grvity		!  < 0, J/m^2
    paybak=min(thkice(i)*fusion*rhoice,availb,		&
	       dlt*rate*fusion*rhoice)				!  > 0, J/m^2
 
    thkice(i)=thkice(i)-paybak/(fusion*rhoice)
    surflx(i)=surflx(i)-paybak/dlt				! W/m^2
    mxchg(i)=grvity *(paybak/fusion)				! kg/m/sec^2
    sxchg(i)=-saldif*(paybak/fusion)				! g/m^2

    qf2d_ave(i)=qf2d_ave(i)-paybak/batrop

    old=saln(1,i)
    dpgt0=max(epsil,dp(1,i,leapn))
    saln(1,i)=saln(1,i)+sxchg(i)*grvity/dpgt0			! g/kg

    if (vrbos) then
     print '(i8,a,i7,a,f9.4,a)',nstep,' (enloan) i=',perm(i),		&
      '  melting:',mxchg(i)/onecm,' cm'
     print 100,nstep,perm(i),						&
      'thkold',thkice(i),'avail-6',thkice(i)*fusion*rhoice*1.e-6,	&
      'availb',availb*1.e-6,'paybak',paybak*1.e-6
     print 100,nstep,perm(i),'sold',old,'snew',saln(1,i),		&
      'sxchg-3',sxchg(i)*1.e-3,'dpgt0',dpgt0/onem
    end if

    if (saln(1,i).lt.0.) then
      print '(a,i8,a)','i=',perm(i),' enloan error: salinity < 0'
      stop '(error)'
    end if

    if (vrbos) then
     if (thkice(i).le.0.) print '(i8,a,i7)',nstep,			&
      '  ice depleted at  i=',perm(i)
     print 100,nstep,perm(i),'sst',temp(1,i),'sss',saln(1,i),		&
      'thkice',thkice(i),'loan-6',-paybak*1.e-6
    end if

   end if                    ! melting or freezing

! --- check column conservation
   colouti=-saldif*thkice(i)*rhoice*grvity
   coloutw=0.
   do k=1,kdm
    coloutw=coloutw+saln(k,i)*dp(k,i,leapn)
   end do
   colini=colinii+coliniw
   colout=colouti+coloutw
   if (vrbos)								&
    print '(i8,a,2es15.7)',perm(i),'  (enloan-2) colini,colout=',	&
     colini,colout
   if (abs(colini-colout).gt.acurcy*max(abs(colini),abs(colout)))	&
     print '(i8,a,2es14.6,es9.1)',perm(i),				&
      ' (enloan) bad column intgl',					&
      colini,colout,(colout-colini)/max(abs(colini),abs(colout))

! --- restore to observed surface salinity
   if (ocnonly .and. sss_rstor) then
    pmevap(i)=pmevap(i)							&
      -(salnow(i)-saln(1,i))*piston*bclin_frq/(saln(1,i)*dlt)
    if (vrbos) write(*,'(i8,a,2f8.2)')					&
      nstep,' current and observed SSS=',saln(1,i),salnow(i)
   end if

   icex(i)=max(thkice(i)-thkmx(i),0.)	! ice exceeding thkmx
   if (icex(i).gt.0.) xcess=xcess+1.
  end if		! ocean point
10 continue
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$REDUCE (xcess,SUM)

! call findmxmn1(icex,nip,'icex',wet)

  call glob2d(thkice,totice)
  if (totice.gt.0.) write (*,'(i9,a,2pf9.3)') nstep,			&
    '  globally averaged ice thickness (cm):',totice/ocnarea

! --- spread out portion of ice thicker than thkmx
  if (xcess.gt.0.) then
   dfflen=sqrt(5.e13/float(nip))	! SQRT[earth surface/(10*cell count)]
   call dffusn_lev(icex,1,1,1,dfflen,mod(nstep,diag_intvl).eq.0)
  end if

  secmn= huuge
  secmx=-huuge
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,q,old,tnew,twrk,pp,qq)			&
!$OMP REDUCTION(min:secmn) REDUCTION(max:secmx)
  do 20 i=1,nip
  secnd(i)=0.
  vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
  if (wet(i) > 0 ) then
   thkice(i)=icex(i)+min(thkice(i),thkmx(i))

   if (thkice(i).gt.1.e-5*thin) then
! --- compute fractional ice coverage for energy flux calculation
    covice(i)=min(1.,thkice(i)/thin)

! --- compute ice surface temperature

    q=kice/max(thkice(i),thin)					! W/deg
    old=temice(i)
    tnew=old
    do iter=1,4
     twrk=tnew
     pp=.5*(q+4.*stfbol*(twrk+KelvC)**3)/(6.*stfbol*(twrk+KelvC)**2)
     qq=(origfx(i)+q*(tfrez(i)-twrk))/   (6.*stfbol*(twrk+KelvC)**2)
     tnew=twrk+(sign(sqrt(max(0.,pp*pp+qq)),pp)-pp)
     if (vrbos) then
      print '(i8,a,i1,5(a,f5.1),2(a,f8.1))',perm(i),' enloan it',	&
       iter,' atm.flx=',origfx(i),' ice flx=',q*(twrk-tfrez(i)),' ->',	&
        q*(tnew-tfrez(i)),' t=',twrk,' ->',tnew,' pp=',pp,' qq=',qq
     end if
     tnew=max(-50.,min(tfrez(i),tnew))
    end do					! iter
    temice(i)=temice(i)+.05*(tnew-temice(i))	! under-relax

    if (old.lt.0.) then
      secnd(i)=temice(i)+ticeol(i)-2.*old
      secmn=min(secmn,secnd(i))
      secmx=max(secmx,secnd(i))
    end if
    ticeol(i)=old

    if (temice(i).lt.-70.) then
     write (*,'(2i8,a,f9.1)') nstep,perm(i),'  ice surf.tmp.',temice(i)
     temice(i)=-70.
    end if
   else				! no ice
    covice(i)=0.
    temice(i)=temp(1,i)
    thkice(i)=max(thkice(i),0.)
   end if
   if (vrbos) then
    print 100,nstep,perm(i),'atmflx',origfx(i),'covice',covice(i),	&
      'thkice',thkice(i),'oldT',old,'newT',temice(i)
   end if

!  ziceav(i)=ziceav(i)+thkice(i)
!  ticeav(i)=ticeav(i)+temice(i)

     if (vrbos) then
       write (*,'(i8,a,i8/(i27,3f10.2))') nstep,			&
        ' exiting enloan t,s,p @i=',perm(i),				&
        (k,temp(k,i),saln(k,i),pres(k+1,i)/onem,k=1,kdm)
       write (*,'(a,2f12.2)') 'ice thickness(m)/tem=',thkice(i),temice(i)
     end if

  end if				! ocean point
20 continue
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- optional: find places where ice temperature oscillates wildly
!SMS$REDUCE(secmn,MIN)
!SMS$REDUCE(secmx,MAX)
  if (max(-secmn,secmx).gt.thresh) then

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
   do 23 i=1,nip
   if (secnd(i).eq.secmn .and. secmn.lt.-thresh) then
    write (*,101) nstep,perm(i),secmn
   end if
   if (secnd(i).eq.secmx .and. secmx.gt. thresh) then
    write (*,101) nstep,perm(i),secmx
   end if
23 continue
!$OMP END PARALLEL DO
!SMS$PARALLEL END

101 format ('step',2i8,'  (enloan) warning: large ice temp fluctuation',f9.1)
  end if
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (mod(nstep,diag_intvl).eq.0) then
   call findmxmn1(temice,nip,'ice  temp',wet)
   call findmxmn1(thkice,nip,'ice thkns',wet)
   call findmxmn1(secnd,nip,'(enloan) 2.deriv',wet)
   call stencl(thkice,1,100.,title//'  ice thickness (cm) OUT')
   call stencl(temice,1,1.,title//'  ice surface temperature')
   call stencl(secnd,1,1.,title//'  ice temp. 2nd derivative')
  end if

  return
  end subroutine enloan

end module hycom_enloan
