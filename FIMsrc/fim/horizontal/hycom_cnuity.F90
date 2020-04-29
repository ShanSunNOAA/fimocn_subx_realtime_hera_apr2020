module hycom_cnuity
use stencilprint
use stenedgprint
use findmaxmin1
use findmaxmin2
use findmaxmin3
use edgmaxmin
contains 
 
!*********************************************************************
!  cnuity
!    Solves continuity equation in ocean model
!    R. Bleck						Dec 2009
!    R. Bleck	added biharm.intfc smoothing option	Jan 2014
!*********************************************************************

  subroutine cnuity(nstep,leapm,leapn, 					&
    um_edg,vm_edg,un_edg,vn_edg,dpm_edg,dp,pres,			&
    massflx,mssfx_avg)

use hycom_constants ,only: wet,onem,tencm,onecm,onemm,r_onem,epsil,	&
                           dtleap,sigjmp,assldp,grvity,land_spval
use hycom_variables ,only: ocnarea,massglb0,totmass,masscor
use fimnamelist     ,only: kdm,itest,diag_intvl,thkdff,biharm
use hycom_control   ,only: kgap
use module_control  ,only: nip,npp
use module_constants,only: nedge,permedge,area,rarea,nprox,prox,	&
                           proxs,sidevec_e,rprox_ln,sideln,rsideln,perm
use hycom_diag      ,only: glob2d

   implicit none
   integer,intent (IN) :: nstep				! model time step
   integer,intent (IN) :: leapm,leapn			! leapfrog time slots

!SMS$DISTRIBUTE (dh,1) BEGIN
   real rpbot(nip),aux1(nip)
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE (dh,2) BEGIN
   real,intent(INOUT) :: dp   (kdm  ,nip,2)		! layer thickness
   real,intent(INOUT) :: pres (kdm+1,nip  )		! interface pressure
   real rmnus(kdm,nip),rplus(kdm,nip),flxdv(kdm,nip),			&
        dpmax(kdm,nip),dpmin(kdm,nip),dptmp(kdm,nip),			&
        oldp(kdm+1,nip),pb_edg(npp,nip),clippd(npp,nip),		&
        flxplus(kdm,nip),flxmnus(kdm,nip),				&
        del2p(kdm,nip),dp_lo(kdm,nip),dpold(kdm,nip)
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE (dh,3) BEGIN
   real,intent(IN)    :: um_edg    (kdm,npp,nip)	! u on edge, mid-time
   real,intent(IN)    :: vm_edg    (kdm,npp,nip)	! v on edge, mid-time
   real,intent(IN)    :: un_edg    (kdm,npp,nip)	! u on edge, old
   real,intent(IN)    :: vn_edg    (kdm,npp,nip)	! v on edge, old
   real,intent(INOUT) :: dpm_edg   (kdm,npp,nip)	! dp on edge, mid-time
   real,intent(INOUT) :: massflx   (kdm,npp,nip)	! baroclinic mass flx
   real,intent(INOUT) :: mssfx_avg (kdm,npp,nip)	! time-avg'd mass flx
   real bolusfx(kdm,npp,nip),antiflx(kdm,npp,nip)
!SMS$DISTRIBUTE END

   integer edg,k,k1,k2,i,ix,edgx,edgcount,ncount
   real outgom,outgon,dpedg,clip,q,old,pbot1,pbot2,p1,p2,		&
        lo_ord(kdm),hi_ord(kdm),flxmx,flxmn,dtinv,valu1,valu2,		&
        difmax,boost,boomx,booz,sum(kdm),deltap,day
   real*8 real8
   character string*24,text*24
   logical vrbos,succes

! --- boost = 1 at all depths except within the lowest 'booz' meters where
! --- it increases linearly to 'boomx'
   boost(pbot1,pbot2,p1,p2,boomx,booz)=max(1.,boomx+(1.-boomx)		&
     *min(pbot1-p1,pbot2-p2)/(booz*onem))

   dtinv=1./dtleap
   if (mod(nstep,diag_intvl).eq.0) then
    print *,'(ocn cnuity) time step',nstep
    do k=1,kdm,kgap
     write (string,'(a,i2)') 'ocnuity old dp k=',k
     call findmxmn3(dp,kdm,nip,2,k,leapn,string,wet)
     write (string,'(a,i2)') 'ocnuity mid dp k=',k
     call findmxmn3(dp,kdm,nip,2,k,leapm,string,wet)
    end do
   end if

   write (text,'(a,i8)') 'ocn_cnuity  step',nstep
   if (mod(nstep,diag_intvl).eq.0) then
    dptmp(:,:)=dp(:,:,leapn)
    call stencl(dptmp,kdm,r_onem,text//'  dp IN')
   end if

! --- lo_ord = low-order flux, piecewise constant (PCM), at old time level.
! --- hi_ord = high-order flux, piecewise linear (PLM), at mid-time.
! --- antiflx = antidiffusive flux (high- minus low-order).
!
!  difmax=epsil
!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE (dpm_edg,dp(:,:,leapn))
!SMS$HALO_COMP(<1,1>) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,edg,edgx,ix,valu1,valu2,		&
!$OMP outgom,q,dpedg,hi_ord,outgon,lo_ord)
   do i=1,nip
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    antiflx(:,:,i)=land_spval
    massflx(:,:,i)=land_spval
    pb_edg (  :,i)=0.
    if (wet(i) > 0 ) then
     do edgcount=1,nedge(i)		! loop through edges
      edg=permedge(edgcount,i)
      ix=prox(edg,i)
      if (wet(ix) > 0 ) then
       edgx=proxs(edg,i)
       do k=1,kdm
        pb_edg(edg,i)=pb_edg(edg,i)+dpm_edg(k,edg,i)

! --- outgo = outging velocity component on edge times edge length
        valu1=sidevec_e(2,edg ,i )*um_edg(k,edg ,i )		&
            - sidevec_e(1,edg ,i )*vm_edg(k,edg ,i )		! m^2/sec
        valu2=sidevec_e(2,edgx,ix)*um_edg(k,edgx,ix)		&
            - sidevec_e(1,edgx,ix)*vm_edg(k,edgx,ix)		! m^2/sec
        outgom=.5*(valu1-valu2)		! avg.to prevent roundoff errors

! --- q = fraction of cell diameter traveled during current time step
! --- rprox_ln = (distance between icos pts)^-1 ; sideln = edge length
        valu1=rsideln(edg ,i )*dtleap*rprox_ln(edg ,i )		! sec/m^2
        valu2=rsideln(edgx,ix)*dtleap*rprox_ln(edgx,ix)		! sec/m^2
        q=min(abs(outgom)*.5*(valu1+valu2),1.)
        dpedg=.5*(dpm_edg(k,edg,i)+dpm_edg(k,edgx,ix))
        if (outgom.ge.0.) then					! outgoing
         hi_ord(k)=dpedg*(1.-q)+dp(k,i ,leapm)*q	! Pa
        else							! incoming
         hi_ord(k)=dpedg*(1.-q)+dp(k,ix,leapm)*q	! Pa
        end if
        hi_ord(k)=hi_ord(k)*outgom				! N/sec

        valu1=sidevec_e(2,edg ,i )*un_edg(k,edg ,i )		&
            - sidevec_e(1,edg ,i )*vn_edg(k,edg ,i )		! m^2/sec
        valu2=sidevec_e(2,edgx,ix)*un_edg(k,edgx,ix)		&
            - sidevec_e(1,edgx,ix)*vn_edg(k,edgx,ix)		! m^2/sec
        outgon=.5*(valu1-valu2)

        if (outgon.ge.0.) then					! outgoing
         lo_ord(k)=outgon*dp(k,i ,leapn)			! kg m/sec^3
        else							! incoming
         lo_ord(k)=outgon*dp(k,ix,leapn)			! kg m/sec^3
        end if
       end do			! vert.loop

       do k=1,kdm
        antiflx(k,edg,i)=hi_ord(k)-lo_ord(k)			! N/sec
        massflx(k,edg,i)=lo_ord(k)				! N/sec
       end do			! vert.loop
      else			! neighbor is land point
       massflx(:,edg,i)=0.
       antiflx(:,edg,i)=0.
      end if			! neighbor is ocean/land point
     end do			! loop through edges
    end if			! ocean point
   end do			! horiz.loop
!$OMP END PARALLEL DO
!SMS$HALO_COMP END
!!SMS$REDUCE (difmax,MAX)
!SMS$PARALLEL END

! --- check for bit-identity of incoming/outgoing fluxes
!  call chkdup(massflx,kdm,'massflx 0',succes)
!  if (.not.succes) stop '(error massflx 0)'
!  call chkdup(antiflx,kdm,'antiflx 0',succes)
!  if (.not.succes) stop '(error antiflx 0)'

! --- advance -dp- field using low-order (diffusive) flux values
 
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos)
   do i=1,nip
    dp_lo(:,i)=0.
    flxdv(:,i)=0.
    aux1(i)=0.
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    if (wet(i) > 0 ) then
     do edg=1,nprox(i)
      do k=1,kdm
       flxdv(k,i)=flxdv(k,i)+massflx(k,edg,i)			! N/sec
      end do
     end do
     do k=1,kdm
      dpold(k,i)=dp(k,i,leapn)
      pres(k+1,i)=pres(k,i)+dp(k,i,leapn)
      dp_lo(k,i)=dp(k,i,leapn)-flxdv(k,i)*rarea(i)*dtleap	! Pa
      aux1(i)=aux1(i)+(dp_lo(k,i)-dpold(k,i))			! diagnostic
! --- don't let roundoff errors cause dp<0
      if (dp_lo(k,i)      .lt.0. .and.				&
          dp_lo(k,i)+onecm.gt.0.) dp_lo(k,i)=0.

      if (dp_lo(k,i).lt.-onemm) then
       print 104, '(ocn cnuity) error -- neg. dp_lo, i,k=',		&
        perm(i),k,'old,low_ord dp:',dp(k,i,leapn),dp_lo(k,i)
       print 103,'low-order fluxes: ',					&
          (massflx(k,edg,i),edg=1,nprox(i))
       print 103,'old/new dp, flxdv:',					&
          dp(k,i,leapn),dp_lo(k,i),-flxdv(k,i)*rarea(i)*dtleap
 103   format (a,6es10.2)
      end if

     end do
    end if			! ocean point
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
 
   if (mod(nstep,diag_intvl).eq.0) then
!   call glob2d(aux1,real8)
!   print '(a,i8,f9.2)','1:avg bottom prs change (mm) during step',	&
!     nstep,real8/(ocnarea*onemm)

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (k2)
    do i=1,nip
     if (i.eq.itest) then
      print 100,'ocn_cnuity  step',nstep,'i=',perm(i),		&
       ' old & low-ord dp, old pres'
 100  format (2(a,i8,2x),a)
      do k1=1,kdm,kgap
       k2=min(k1+5,kdm)
!      print 101,(dp(k,i,leapn)*r_onem,k=k1,k2)
!101  format (6x,6f11.3)
!      print 101,(dp_lo(k,i)*r_onem,k=k1,k2)
!      print 102,(pres(k,i)*r_onem,k=k1,k2+1)
!102  format (7f11.3)
       print 101,(dp(k,i,leapn),k=k1,k2)
 101  format (6x,6es11.3)
       print 101,(dp_lo(k,i),k=k1,k2)
       print 102,(pres(k,i),k=k1,k2+1)
 102  format (7es11.3)
       print *
      end do
     end if
    end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
     
    dptmp(:,:)=dp(:,:,leapn)
    call stenedg(dptmp,massflx,kdm,text//'  old dp, low-ord flx')
    call stenedg(dp_lo,massflx,kdm,text//'  low-ord dp, low-ord flx')
    call stencl(dp_lo,kdm,r_onem,text//'  low-ord dp')
   end if

! --- segregate incoming and outgoing antidiffusive fluxes

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos)
   do i=1,nip
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    flxplus(:,i)= epsil
    flxmnus(:,i)=-epsil
    if (wet(i) > 0 ) then
     do edg=1,nprox(i)
      do k=1,kdm
       if (antiflx(k,edg,i).ge.0.) then			! outgoing
        flxplus(k,i)=flxplus(k,i)+antiflx(k,edg,i)	! N/sec
       else						! incoming
        flxmnus(k,i)=flxmnus(k,i)+antiflx(k,edg,i)	! N/sec
       end if
      end do
     end do
    end if			! ocean point
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! --- at each grid point, determine the ratio of the largest permissible
! --- pos. (neg.) change in -dp- to the sum of all incoming (outgoing) fluxes
 
!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE (dp_lo)
!$OMP PARALLEL DO PRIVATE (edgx)
   do i=1,nip
    rplus(:,i)=0.		! limit on  i n c o m i n g  fluxes
    rmnus(:,i)=0.		! limit on  o u t g o i n g  fluxes
    dpmax(:,i)=0.		! only needed if passed to stencl
    dpmin(:,i)=0.		! only needed if passed to stencl
    if (wet(i) > 0 ) then
     do k=1,kdm
      dpmax(k,i)=max(dp_lo(k,i),dp(k,i,leapn))
      dpmin(k,i)=min(dp_lo(k,i),dp(k,i,leapn))
     end do
     do edg=1,nprox(i)
      edgx=prox(edg,i)
      if (wet(edgx) > 0 ) then
       do k=1,kdm
        dpmax(k,i)=max(dpmax(k,i),dp_lo(k,edgx),dp(k,edgx,leapn))
        dpmin(k,i)=min(dpmin(k,i),dp_lo(k,edgx),dp(k,edgx,leapn))
       end do
      end if			! neighbor is ocean point
     end do			! loop through edges
     do k=1,kdm
      dpmin(k,i)=max(0.,dpmin(k,i))
      rmnus(k,i)=(dp_lo(k,i)-dpmin(k,i))*area(i)/(dtleap*flxplus(k,i))
      rplus(k,i)=(dp_lo(k,i)-dpmax(k,i))*area(i)/(dtleap*flxmnus(k,i))
     end do
    end if			! ocean point
   end do			! horiz.loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

!    call stencl(dpmax,kdm,1./onecm,text//'  dpmax (cm)')
!    call stencl(dpmin,kdm,1./onecm,text//'  dpmin (cm)')
!    call stenedg(rplus,antiflx,kdm,text//'  rplus, antiflx')
!    call stenedg(rmnus,antiflx,kdm,text//'  rmnus, antiflx')

!    do k=1,kdm,kgap
!     write (string,'(a,i2)') 'ocnuity k=',k
!     call findmxmn2(rmnus,kdm,nip,k,trim(string)//' rmnus',wet)
!     call findmxmn2(rplus,kdm,nip,k,trim(string)//' rplus',wet)
!    end do

! --- limit antidiffusive fluxes
! --- (keep track in -clippd- of discrepancy between high-order
! --- fluxes and the sum of low-order and clipped antidiffusive fluxes.
! --- this will be used later to restore nondivergence of barotropic flow)

!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE (rmnus,rplus)
!SMS$HALO_COMP(<1,1>) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,edg,ix,clip)
   do i=1,nip
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    clippd(:,i)=0.
    if (wet(i) > 0 ) then
     do edgcount=1,nedge(i)		! loop through edges
      edg=permedge(edgcount,i)
      ix=prox(edg,i)
      if (wet(ix) > 0 ) then
       do k=1,kdm

!        if (vrbos) then
!         print '(2i5,2i2,a,4es9.1)',perm(i),perm(ix),edg,k,	&
!          ' rplus(i/ix),rmnus(i/ix)',rplus(k,i),rplus(k,ix),	&
!          rmnus(k,i),rmnus(k,ix)
!        end if

        if (antiflx(k,edg,i).ge.0.) then
         clip=min(1.,rmnus(k,i),rplus(k,ix))		! outgoing
        else
         clip=min(1.,rplus(k,i),rmnus(k,ix))		! incoming
        end if
        clippd(edg,i)=clippd(edg,i)+antiflx(k,edg,i)*(1.-clip)	! N/sec
        antiflx(k,edg,i)=antiflx(k,edg,i)*clip			! N/sec
        massflx(k,edg,i)=massflx(k,edg,i)+antiflx(k,edg,i)	! N/sec
       end do			! vert.loop
      end if			! neighbor is ocean point
     end do			! loop through edges
    end if			! ocean point
   end do			! horiz.loop
!$OMP END PARALLEL DO
!SMS$HALO_COMP END
!SMS$PARALLEL END
 
! --- evaluate effect of antidiffusive fluxes on -dp- field
 
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos)
   do i=1,nip
    flxdv(:,i)=0.
    aux1(i)=0.
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    if (wet(i) > 0 ) then
     do edg=1,nprox(i)
      do k=1,kdm
       flxdv(k,i)=flxdv(k,i)+antiflx(k,edg,i)			! N/sec
      end do
     end do
     do k=1,kdm
      dp(k,i,leapn)=dp_lo(k,i)-flxdv(k,i)*rarea(i)*dtleap 	! Pa
      aux1(i)=aux1(i)+(dp(k,i,leapn)-dpold(k,i))		! diagnostic
! --- don't let roundoff errors cause dp<0
      if (dp(k,i,leapn)      .lt.0. .and.			&
          dp(k,i,leapn)+onecm.gt.0.) dp(k,i,leapn)=0.

      if (dp(k,i,leapn).lt.-onemm) then
       print 104,'(ocn cnuity) error -- neg. dp  stage 1, i,k=',	&
        perm(i),k,'old,mid,new dp:',dpold(k,i),dp(k,i,leapm),		&
        dp(k,i,leapn)
104    format (a,i7,i3/9x,a,3es11.3)
       print 103,'antidiff fluxes:  ',					&
          (antiflx(k,edg,i),edg=1,nprox(i))
       print 103,'old/new dp, flxdv:',					&
          dp_lo(k,i),dp(k,i,leapn),-flxdv(k,i)*rarea(i)*dtleap
      end if

      pres(k+1,i)=pres(k,i)+dp(k,i,leapn)
     end do			! vert.loop
     rpbot(i)=1./pres(kdm+1,i)
    end if			! ocean point
   end do			! horiz.loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

!  call edgmxmn(massflx,kdm,'(ocnuity-1) massflx',wet)

! --- check for bit-identity of incoming/outgoing fluxes
!  call chkdup(massflx,kdm,'massflx 1',succes)
!  if (.not.succes) stop '(error massflx 1)'
!  call chkdup(antiflx,kdm,'antiflx 1',succes)
!  if (.not.succes) stop '(error massflx 1)'

    if (mod(nstep,diag_intvl).eq.0) then
!    call glob2d(aux1,real8)
!    print '(a,i8,f9.2)','2:avg bottom prs change (mm) during step',	&
!      nstep,real8/(ocnarea*onemm)
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (k2)
     do i=1,nip
      if (i.eq.itest) then
       print 100,'ocn_cnuity',nstep,'i=',perm(i),			&
        ' low-ord dp, new dp & pres'
        do k1=1,kdm,kgap
         k2=min(k1+5,kdm)
!        print 101,(dp_lo(k,i)*r_onem,k=k1,k2)
!        print 101,(dp(k,i,leapn)*r_onem,k=k1,k2)
!        print 102,(pres(k,i)*r_onem,k=k1,k2+1)
         print 101,(dp_lo(k,i),k=k1,k2)
         print 101,(dp(k,i,leapn),k=k1,k2)
         print 102,(pres(k,i),k=k1,k2+1)
         print *
        end do
      end if
     end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
      
     dptmp(:,:)=dp(:,:,leapn)
     call stenedg(dptmp,antiflx,kdm,text//'  new dp, antidiff flx')
    end if

! --- recover fluxes lost in the flux limiting process.
! --- treat these fluxes as an 'upstream' barotropic correction to
! --- the sum of diffusive and antidiffusive fluxes obtained so far.
 
!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE (rpbot,dp(:,:,leapn))
!SMS$HALO_COMP(<1,1>) BEGIN
!$OMP PARALLEL DO PRIVATE (edg,ix,q)
   do i=1,nip
    flxdv(:,i)=0.
    if (wet(i) > 0 ) then
     do edgcount=1,nedge(i)		! loop through edges
      edg=permedge(edgcount,i)
      ix=prox(edg,i)
      if (wet(ix) > 0 ) then
       do k=1,kdm
        q=0.
        if      (clippd(edg,i).gt.0.) then
         q=clippd(edg,i)*(dp(k,i ,leapn)*rpbot(i ))
        else if (clippd(edg,i).lt.0.) then
         q=clippd(edg,i)*(dp(k,ix,leapn)*rpbot(ix))
        end if
        flxdv(k,i)=flxdv(k,i)+q					! N/sec
        massflx(k,edg,i)=massflx(k,edg,i)+q			! N/sec
       end do
      end if			! neighbor is ocean point
     end do
    end if			! ocean point
   end do
!$OMP END PARALLEL DO
!SMS$HALO_COMP END
!SMS$PARALLEL END

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (old)
   do i=1,nip
    aux1(i)=0.
    if (wet(i) > 0 ) then
     do k=1,kdm
      old=dp(k,i,leapn)
      dp(k,i,leapn)=dp(k,i,leapn)-flxdv(k,i)*rarea(i)*dtleap	! Pa
      aux1(i)=aux1(i)+(dp(k,i,leapn)-dpold(k,i))		! diagnostic
! --- don't let roundoff errors cause dp<0
      if (dp(k,i,leapn)      .lt.0. .and.			&
          dp(k,i,leapn)+onecm.gt.0.) dp(k,i,leapn)=0.
      pres(k+1,i)=pres(k,i)+dp(k,i,leapn)

      if (dp(k,i,leapn).lt.-onemm) then
       print 104,'(ocn cnuity) error -- neg. dp  stage 2, i,k=',	&
        perm(i),k,'old,mid,new dp:',dpold(k,i),dp(k,i,leapm),		&
        dp(k,i,leapn)
       print 103,'old/new dp, flxdv:',					&
          old,dp(k,i,leapn),-flxdv(k,i)*rarea(i)*dtleap
      end if

     end do
    end if			! ocean point
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!   if (mod(nstep,diag_intvl).eq.0) then
!    call glob2d(aux1,real8)
!    print '(a,i8,f9.2)','3:avg bottom prs change (mm) during step',	&
!     nstep,real8/(ocnarea*onemm)
!   end if

!  call edgmxmn(massflx,kdm,'(ocnuity-2) massflx',wet)

! --- check for bit-identity of incoming/outgoing fluxes
!  call chkdup(massflx,kdm,'massflx 2',succes)
!  if (.not.succes) stop '(error massflx 2)'
!  call chkdup(antiflx,kdm,'antiflx 2',succes)
!  if (.not.succes) stop '(error massflx 2)'

! --- -----------------------------------
! --- interface smoothing (=> bolus flux)
! --- -----------------------------------

   if (thkdff.eq.0.)  go to 1

!sms$compare_var(massflx, "hycom_cnuity.F90 - massflx1 ")

   if (mod(nstep,diag_intvl).eq.0) then
    dptmp(:,:)=dp(:,:,leapn)
    call stencl(dptmp,kdm,r_onem,text//'  dp bfor intfc_smoo')
    call stencl(pres(2:kdm+1,:),kdm,r_onem,text//' pres bfor intfc_smoo')
   end if

! --- to avoid roundoff errors in flxmx,flxmn calc., deduce -dp- from -p-.
!SMS$EXCHANGE (pres)
!SMS$PARALLEL (dh,i) BEGIN
!SMS$HALO_COMP(<1,1>) BEGIN
!$OMP PARALLEL DO
   do i=1,nip
    oldp(:,i)=land_spval
    if (wet(i) > 0 ) then
     do k=1,kdm
      dp(k,i,leapn)=pres(k+1,i)-pres(k,i)
      oldp(k+1,i)=pres(k+1,i)		! diag.use only
     end do
    end if
   end do
!$OMP END PARALLEL DO

! --- compute 'pressure fluxes' related to interface slope

!$OMP PARALLEL DO
   do i=1,nip
    del2p(:,i)=land_spval
    if (biharm.eq.0.) then
     do k=2,kdm
      del2p(k,i)=pres(k,i)
     end do
    end if
   end do
!$OMP END PARALLEL DO
!SMS$HALO_COMP END
!SMS$PARALLEL END

   if (biharm.gt.0.) then

! --- compute a quantity proportional to the (negative) laplacian of -p- by
! --- subtracting from each -p- value the average of its neighbors.

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (sum,ncount,ix,deltap)
    do i=1,nip
     if (wet(i) > 0) then
      sum(:)=0.
      ncount=0
      do edg=1,nprox(i)			! loop through edges
       ix=prox(edg,i)			! neighbor across shared edge
       if (wet(ix) > 0) then
        ncount=ncount+1
        do k=2,kdm			! verical loop
         if (pres(kdm+1,ix).lt.pres(k,i)) then

! --- low bottom pressure at neighboring point suggests that neighboring point
! --- is forced up in column by shallow bathymetry. don't use it for smoothing.
          deltap=pres(kdm+1,ix)-pres(k,ix)

         else if (pres(kdm+1,i).lt.pres(k,ix)) then

! --- low bottom pressure at central point suggests that central point
! --- is forced up in column by shallow bathymetry. don't use it for smoothing.
          deltap=pres(k,i)-pres(kdm+1,i)

         else				! "normal" case
          deltap=pres(k,i)-pres(k,ix)
         end if
         sum(k)=sum(k)+deltap
        end do				! vert.loop
       end if				! neighbor is ocean point
      end do				! loop through edges
      do k=2,kdm
       del2p(k,i)=biharm*sum(k)/float(ncount)+(1.-biharm)*pres(k,i)
      end do
     end if				! ocean point
    end do				! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!SMS$EXCHANGE (del2p)

    if (mod(nstep,diag_intvl).eq.0)				&
     call stencl(del2p(2:kdm,:),kdm-1,r_onem,text//' del2p bfor intfc_smoo')

   end if 			! biharm > 0

!   call chkdup(sideln,1,'sideln',succes)		! should be done in hycom_init
!   if (.not.succes) call cpydup(sideln,1,.false.)
!   if (.not.succes) stop '(error sideln)'

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (ix,q,flxmx,flxmn)
   do i=1,nip
!   vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    bolusfx(:,:,i)=land_spval
    if (wet(i) > 0 ) then
     bolusfx(:,:,i)=0.
     do edg=1,nprox(i)		! loop through edges
      ix=prox(edg,i)
      if (wet(ix) > 0 ) then
       do k=2,kdm

! --- pressure flux is proportional to slope of -del2p- field
        q=dtleap*thkdff*(del2p(k,i)-del2p(k,ix))*sideln(edg,i)
!          *boost(pres(kdm+1,i),pres(kdm+1,ix),pres(k,i),pres(k,ix),	&
!               3.,30.)							! N

! --- keep interfaces from going underground
        flxmx= (pres(kdm+1,ix)-pres(k,ix))*area(ix)			&
          /float(nprox(ix))						! N
        flxmn=-(pres(kdm+1,i )-pres(k,i ))*area(i )			&
          /float(nprox(i))						! N
        bolusfx(k,edg,i)=min(flxmx*.999,max(flxmn*.999,q))		! N

! --- difference of 2 interface 'pressure fluxes' becomes thickness flux
! --- (outgoing > 0)
        bolusfx(k-1,edg,i)=bolusfx(k,edg,i)-bolusfx(k-1,edg,i)		! N
       end do			! vert.loop
       bolusfx(kdm,edg,i)=-bolusfx(kdm,edg,i)				! N
      end if			! neighbor is ocean point
     end do			! loop through edges
    end if			! ocean point
   end do			! horiz.loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!sms$compare_var(bolusfx, "hycom_cnuity.F90 - bolusfx1 ")

   if (mod(nstep,diag_intvl).eq.0) then
    dptmp(:,:)=dp(:,:,leapn)
    call stenedg(dptmp,bolusfx,kdm,text//'  dp bfor intfc_smoo, bolus flx')
    call stencl(pres(2:kdm+1,:),kdm,r_onem,text//' pres bfor intfc_smoo')
   end if

! --- check for bit-identity of incoming/outgoing fluxes
!  call chkdup(bolusfx,kdm,'bolusfx 1',succes)
!  if (.not.succes) stop '(error bolusfx 1)'
!  call cpydup(bolusfx,kdm,.true.)

   if (biharm.gt.0.) then

! --- segregate outgoing bolus fluxes

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
    do i=1,nip
!    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
     flxplus(:,i)=epsil
     if (wet(i) > 0 ) then
      do edg=1,nprox(i)		! loop through edges
       flxplus(:,i)=flxplus(:,i)+max(0.,bolusfx(:,edg,i))		! N
      end do

! --- at each grid point, determine the ratio of the largest permissible
! --- mass loss to the sum of all outgoing bolus fluxes

      rmnus(:,i)=dp(:,i,leapn)*area(i)/flxplus(:,i)		! dim.less
     end if			! ocean point
    end do			! horiz.loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! --- limit bolus fluxes

!SMS$EXCHANGE (rmnus)
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (ix,clip)
    do i=1,nip
!    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
     clippd(:,i)=0.
     if (wet(i) > 0 ) then
      do edg=1,nprox(i)		! loop through edges
       ix=prox(edg,i)
       if (wet(ix) > 0) then
        do k=1,kdm

!       if (wet(ix).eq.0 .and. bolusfx(k,edg,i).ne.0.)			&
!        print '(2(a,i8),es11.3)',					&
!          'error - nonzero flux across coastline between',perm(i),	&
!          ' and',perm(ix),bolusfx(k,edg,i)

         
         if (bolusfx(k,edg,i).ne.0.) then
          if      (bolusfx(k,edg,i).gt.0.) then
           clip=min(1.,rmnus(k,i ))		! outgoing
          else if (bolusfx(k,edg,i).lt.0.) then
           clip=min(1.,rmnus(k,ix))		! incoming
          end if

! --- clipped part of bolus flux is kept to restore zero column integral later
          clippd(edg,i)=clippd(edg,i)+bolusfx(k,edg,i)*(1.-clip)	! N
          bolusfx(k,edg,i)=bolusfx(k,edg,i)*clip			! N
         end if			! bolusfx .ne. 0
        end do			! vert.loop
       end if			! neighbor is ocean point
      end do			! loop through edges
     end if			! ocean point
    end do			! horiz.loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!sms$compare_var(rmnus, "hycom_cnuity.F90 - rmnus ")
!sms$compare_var(flxplus, "hycom_cnuity.F90 - flxplus ")
!sms$compare_var(bolusfx, "hycom_cnuity.F90 - bolusfx2 ")

! --- check for bit-identity of incoming/outgoing fluxes
!   call chkdup(bolusfx,kdm,'bolusfx 2',succes)
!   if (.not.succes) stop '(error bolusfx 2)'
!   call chkdup(massflx,kdm,'massflx 2',succes)
!   if (.not.succes) stop '(error massflx 2)'
!   call chkdup(clippd,1,'clippd   ',succes)
!   if (.not.succes) stop '(error clippd)'
!   call cpydup(bolusfx,kdm,.true.)

   end if			! biharm > 0

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,old)
   do i=1,nip
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    flxdv(:,i)=0.
    aux1(i)=0.
    if (wet(i) > 0 ) then
     do edg=1,nprox(i)
      flxdv(:,i)=flxdv(:,i)+bolusfx(:,edg,i)			! N
     end do
     do k=1,kdm
      old=dp(k,i,leapn)
      dp(k,i,leapn)=dp(k,i,leapn)-flxdv(k,i)*rarea(i)		! Pa
      aux1(i)=aux1(i)+(dp(k,i,leapn)-dpold(k,i))		! diagnostic
! --- don't let roundoff errors cause dp<0
      if (dp(k,i,leapn)      .lt.0. .and.			&
          dp(k,i,leapn)+onecm.gt.0.) dp(k,i,leapn)=0.
      dp(k,i,leapn)=dp(k,i,leapn)*(1.d0+masscor)
      pres(k+1,i)=pres(k,i)+dp(k,i,leapn)
      dptmp(k,i)=dp(k,i,leapn)

      if (dp(k,i,leapn).lt.-onemm) then
       print 104,'(ocn cnuity) error -- neg. dp  stage 3, i,k=',	&
        perm(i),k,'old,mid,new dp:',dpold(k,i),dp(k,i,leapm),		&
        dp(k,i,leapn)
       print 103,'bolus fluxes:    ',					&
          (bolusfx(k,edg,i),edg=1,nprox(i))
       print 103,'old/new dp, flxdv:',					&
          old,dp(k,i,leapn),-flxdv(k,i)*rarea(i)*dtleap
      end if

     end do			! vert.loop

     if (vrbos) then
      print '(2(a,i8),a/a,6i11)','(ocn_cnuity) step',nstep,	&
         '  i=',i,'  bolus fluxes','neighbors',			&
         (perm(prox(edg,i)),edg=1,nprox(i))
      do k=1,kdm
       print '(i3,7x,6es11.3)',k,(bolusfx(k,edg,i),edg=1,nprox(i))
      end do
      print '(2(a,i8),a)','step',nstep,'  i=',i,		&
         ' (ocn_cnuity) old/new lower intfc_pres across edges'
      do k=2,kdm+1
       print '(i2,7f11.2)',k,oldp(k,i)*r_onem,			&
          (oldp(k,prox(edg,i))*r_onem,edg=1,nprox(i))
       print '(i2,7f11.2)',k,pres(k,i)*r_onem,			&
          (pres(k,prox(edg,i))*r_onem,edg=1,nprox(i))
      end do
     end if

     rpbot(i)=1./pres(kdm+1,i)
    end if			! ocean point
   end do			! horiz.loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!   if (mod(nstep,diag_intvl).eq.0) then
!    call glob2d(aux1,real8)
!    print '(a,i8,f9.2)','4:avg bottom prs change (mm) during step',	&
!     nstep,real8/(ocnarea*onemm)
!   end if

!  call edgmxmn(massflx,kdm,'(ocnuity-3) massflx',wet)

   if (biharm.gt.0.) then

! --- distribute clipped part of bolusflx over column as upstream correction

!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE(dptmp,rpbot)
!$OMP PARALLEL DO PRIVATE (ix,q)
    do i=1,nip
     flxdv(:,i)=0.
     aux1(i)=0.
     if (wet(i) > 0 ) then
      do edg=1,nprox(i)
       ix=prox(edg,i)

       if (wet(ix).eq.0 .and. clippd(edg,i).ne.0.)			&
          print '(a,2i8,es11.3)',					&
          '(ocn_cnuity) error - nonzero flux across coastline',		&
          perm(i),perm(ix),clippd(edg,i)

       q=0.
       do k=1,kdm
        if      (clippd(edg,i).gt.0.) then
         q=clippd(edg,i)*dptmp(k,i )*rpbot(i )		! N
        else if (clippd(edg,i).lt.0.) then
         q=clippd(edg,i)*dptmp(k,ix)*rpbot(ix)		! N
        end if
        flxdv(k,i)=flxdv(k,i)+q
        massflx(k,edg,i)=massflx(k,edg,i)+q*dtinv		! N/sec
       end do			! vert.loop
      end do			! loop through edges
      do k=1,kdm
       dp(k,i,leapn)=dp(k,i,leapn)-flxdv(k,i)*rarea(i)		! Pa
       aux1(i)=aux1(i)+(dp(k,i,leapn)-dpold(k,i))		! diagnostic
       pres(k+1,i)=pres(k,i)+dp(k,i,leapn)
      end do
     end if			! ocean point
    end do			! horiz.loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!   if (mod(nstep,diag_intvl).eq.0) then
!    call glob2d(aux1,real8)
!    print '(a,i8,f9.2)','5:avg bottom prs change (mm) during step',	&
!     nstep,real8/(ocnarea*onemm)
!   end if

   end if			! biharm > 0

! --- add instantaneous mass flux to time average
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
   do i=1,nip
    if (wet(i) > 0 ) then
     do edg=1,nprox(i)
      mssfx_avg(:,edg,i)=mssfx_avg(:,edg,i)+massflx(:,edg,i)
     end do
    end if			! ocean point
   end do			! horiz.loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! --- check for bit-identity of incoming/outgoing fluxes
!  call chkdup(massflx,kdm,'massflx 3',succes)
!  if (.not.succes) stop '(error massflx 3)'
!  call cpydup(massflx,kdm,.true.)

   if (mod(nstep,diag_intvl).eq.0) then
    dptmp(:,:)=dp(:,:,leapn)
    call stenedg(dptmp,bolusfx,kdm,text//'  dp aftr intfc_smoo, bolus flx')
    call stencl(pres(2:kdm+1,:),kdm,r_onem,text//' pres aftr intfc_smoo')
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (k2)
    do i=1,nip
     if (i.eq.itest) then
      print 100,'(ocn_cnuity) step',nstep,'i=',perm(i),			&
       ' old/new dp, pres aftr intfc_smoo'
      do k1=1,kdm,kgap
       k2=min(k1+5,kdm)
       print 101,(dpold(k,i)*r_onem,k=k1,k2)
       print 101,(dp(k,i,leapn)*r_onem,k=k1,k2)
       print 102,(pres(k,i)*r_onem,k=k1,k2+1)
       print *
      end do
     end if
    end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
   end if
 1 continue

!sms$compare_var(massflx, "hycom_cnuity.F90 - massflx2 ")

! --- time smoothing

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos)
   do i=1,nip
    aux1(i)=0.
    if (wet(i) > 0 ) then
     vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
     do k=1,kdm
      dp(k,i,leapm)=assldp*(dpold(k,i)+dp(k,i,leapn))			&
        +(1.-assldp-assldp)*           dp(k,i,leapm)

      if (min(dp(k,i,leapm),dp(k,i,leapn)).lt.-onemm) then
       print 104,'(ocn cnuity) error -- neg. dp  stage 4, i,k=',	&
        perm(i),k,'old,mid,new dp:',dpold(k,i),dp(k,i,leapm),		&
        dp(k,i,leapn)
      end if

      aux1(i)=aux1(i)+(dp(k,i,leapn)-dpold(k,i))		! diagnostic
     end do
    end if			! ocean point
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
!  call edgmxmn(massflx,kdm,'(ocnuity-4) massflx',wet)

   if (mod(nstep,diag_intvl).eq.0) then
    print *,'(ocn cnuity) time step',nstep
    do k=1,kdm,kgap
     write (string,'(a,i2)') 'ocnuity mid dp k=',k
     call findmxmn3(dp,kdm,nip,2,k,leapm,string,wet)
     write (string,'(a,i2)') 'ocnuity new dp k=',k
     call findmxmn3(dp,kdm,nip,2,k,leapn,string,wet)
    end do
    if (nstep.gt.1)							&
    call findmxmn1(aux1,nip,'ocnuity botm.pr.chg',wet)
   end if

!  if (mod(nstep,diag_intvl).eq.0) then
!   call glob2d(aux1,real8)
!   print '(a,i8,f9.2)','6:avg bottom prs change (mm) during step',	&
!     nstep,real8/(ocnarea*onemm)
!   call stencl(aux1,1,1./onecm,text//'  bottom pres.tendency (cm)')
!   dptmp(:,:)=dp(:,:,leapn)
!   call stencl(dptmp,kdm,r_onem,text//'  dp OUT')
!  end if

! --- check for bit-identity of incoming/outgoing fluxes
   call chkdup(massflx,kdm,'massflx 4',succes)
   if (.not.succes) stop '(error massflx 4)'

! --- determine drift correction factor for total mass
   day=float(nstep)*.5*dtleap/86400.
!  if (nstep==1 .or. mod(day+.0001,30.).lt..0002) then		!  once a month
   if (nstep==1 .or. mod(day+.0001,10.).lt..0002) then		!  every 10 days

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
   do i=1,nip
    aux1(i)=0.
    if (wet(i) > 0 ) aux1(i)=pres(kdm+1,i)
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

    call glob2d(aux1,totmass)
    totmass=totmass/grvity				! N => kg
    if (nstep==1) massglb0=totmass
!   masscor=dtleap*(massglb0-totmass)/(massglb0*91.*86400.)		! 3mo
!   masscor=dtleap*(massglb0-totmass)/(massglb0*30.*86400.)		! 30dy
    masscor=dtleap*(massglb0-totmass)/(massglb0*10.*86400.)		! 10dy
    write (*,'(i9,a,f11.7,a,es9.2)') nstep, '  overall rel.mass gain:',	&
      (totmass-massglb0)/massglb0,'  => thknss factor',masscor
   end if

   return
   end subroutine cnuity


   subroutine cpydup(flux,kdm,oppsign)

! --- to enforce bit-identity of duplicate edge fluxes, copy flux
! --- across a given edge from lower- to higher-indexed cell

   use module_control  ,only: nip,npp
   use module_constants,only: nprox,prox,proxs,perm
   use fimnamelist     ,only: itest
   use hycom_constants ,only: wet

   implicit none
   integer,intent(IN) :: kdm
   logical,intent(IN) :: oppsign	! if true, vrbls have opposite sign
!SMS$DISTRIBUTE  (dh,3) BEGIN
   real,intent(INOUT) :: flux(kdm,npp,nip)
!SMS$DISTRIBUTE END
   integer i,ix,edg,edgx,k
   logical vrbos

   print *,'entering cpydup with itest=',perm(itest)
!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE (flux)
!$OMP PARALLEL DO PRIVATE (vrbos,ix,edgx)
   do i=1,nip 
    vrbos=i.eq.itest
    if (wet(i) > 0 ) then
     do edg=1,nprox(i)          	! loop through edges
      ix=prox(edg,i)
      if (wet(ix) > 0 ) then
       edgx=proxs(edg,i)
       if (perm(ix).lt.i) then
        do k=1,kdm
         if (oppsign .and.						&
           flux(k,edg,i)+flux(k,edgx,ix).ne.0.) then

!         if (vrbos) then
           print 100,'(cpydup) replace edge vrbl',flux(k,edg,i),	&
           ' at k,i,edg=',k,i,edg,'by',-flux(k,edgx,ix),		&
           ' at k,i,edg=',k,perm(ix),edgx
 100       format (a,es16.8,a,i3,i8,i3/a26,es16.8,a,i3,i8,i3)
!         end if

          flux(k,edg,i)=-flux(k,edgx,ix)
         else if (.not.oppsign .and.					&
           flux(k,edg,i)-flux(k,edgx,ix).ne.0.) then

!         if (vrbos) then
           print 100,'(cpydup) replace edge vrbl',flux(k,edg,i),	&
           ' at k,i,edg=',k,i,edg,'by', flux(k,edgx,ix),		&
           ' at k,i,edg=',k,perm(ix),edgx
!         end if

          flux(k,edg,i)=flux(k,edgx,ix)
         end if
        end do			! vert.loop
       end if			! ix < i
      end if			! neighbor is ocean point
     end do			! loop though edges
    end if			! ocean point
   end do			! horiz.loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

   return
   end subroutine cpydup


   subroutine chkdup(flux,kdm,info,succes)

! --- check duplicate edge fluxes for bit-identity

   use module_control  ,only: nip,npp
   use module_constants,only: nprox,prox,proxs,perm
   use fimnamelist     ,only: itest
   use hycom_constants ,only: wet

   implicit none
   integer  ,intent(IN) :: kdm
   character,intent(IN) :: info*(*)
   logical ,intent(OUT) :: succes
!SMS$DISTRIBUTE  (dh,3) BEGIN
   real,intent(INOUT) :: flux(kdm,npp,nip)
!SMS$DISTRIBUTE END
   integer i,ix,edg,edgx,k,ncount

!  print *,'entering chkdup with variable ',trim(info)
   ncount = 0
!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE (flux)
!$OMP PARALLEL DO PRIVATE (ix,edgx) REDUCTION(+:ncount)
   do i=1,nip 
    if (wet(i) > 0 ) then
     do edg=1,nprox(i)          	! loop through edges
      ix=prox(edg,i)
      if (wet(ix) > 0 ) then
       edgx=proxs(edg,i)
       do k=1,kdm
        if (flux(k,edg,i)+flux(k,edgx,ix).ne.0. .and.			&
            flux(k,edg,i)-flux(k,edgx,ix).ne.0.) then
         print '(2a,i3,2(i8,i3)/2es16.8)',info,				&
          ' -- edge variables differ at k,i,edg,ix,edgx=',k,		&
          perm(i),edg,perm(ix),edgx,flux(k,edg,i),flux(k,edgx,ix)
         ncount = 1
        end if
       end do
      end if			! neighbor is ocean point
     end do			! loop though edges
    end if			! ocean point
   end do			! horiz.loop
!SMS$PARALLEL END
!SMS$REDUCE (ncount,SUM)
   if (ncount > 0) then
     succes=.false.
   else
     succes=.true.
   end if

   return
   end subroutine chkdup
end module hycom_cnuity
