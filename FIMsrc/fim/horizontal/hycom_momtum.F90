module hycom_momtum
use stencilprint
use stenedgprint
use findmaxmin2
use findmaxmin3
contains

!*********************************************************************
!     momtum
!       Solves momentum equations in HYCOM
!       R. Bleck					Jan 2010
!       R. Bleck   switch to biharmonic dissipation	Jul 2013
!*********************************************************************

  subroutine momtum(nstep,leapm,leapn,					&
    utrop,vtrop,uclin,vclin,um_edg,vm_edg,geop_edg,mont_edg,uvsq_edg,	&
    spcv_edg,taux,tauy,ubforc,vbforc,dp,pres,ustarb)

use module_control  ,only: npp,nip
use module_constants,only: prox,nprox,rarea,sideln,sidevec_c,		&
                           sidevec_e,sn,cs,corio,grvity,deg_lat,perm,	&
                           grvity
use hycom_constants ,only: wet,meshsz,thref,dpth_edg,batrop,drcoef,	&
			   tenm,onem,tencm,onemm,epsil,dtleap,assluv,	&
			   land_spval
use hycom_control   ,only: modesplit,bclin_frq
use fimnamelist     ,only: kdm,itest,diag_intvl,janjic_ocn,veldff,smagcf
use hycom_variables ,only: covice
use hycom_edgvar    ,only: blokf
use hycom_dissip

  implicit none
! External variables

  integer,intent (IN)  :: nstep			! model time step
  integer,intent (IN)  :: leapm,leapn		! leapfrog time slots
!SMS$DISTRIBUTE (dh,1) BEGIN
  real,intent (IN)    :: taux       (nip)	! eastw'd wind stress
  real,intent (IN)    :: tauy       (nip)	! northw'd wind stress
  real,intent (OUT)   :: ubforc     (nip)	! forcing of ba'trop u
  real,intent (OUT)   :: vbforc     (nip)	! forcing of ba'trop v
  real,intent (INOUT) :: ustarb     (nip)
! Local variables
  real drag(nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
  real,intent (IN)    :: dp     (kdm,nip,2)	! layer thickness
  real,intent (INOUT) :: pres (kdm+1,nip)	! interface pressure
  real,intent (INOUT) :: uclin  (kdm,nip,2)	! eastw'd baroclinic velocity
  real,intent (INOUT) :: vclin  (kdm,nip,2)	! northw'd baroclinic velocity
  real,intent (IN)    :: utrop    (3,nip)	! eastw'd barotropic velocity
  real,intent (IN)    :: vtrop    (3,nip)	! northw'd barotropic velocity
! Local variables
  real pgfx(kdm,nip),pgfy(kdm,nip),gradx(kdm,nip),grady(kdm,nip),	&
       utotn(kdm,nip),vtotn(kdm,nip),uzeta(kdm,nip),vzeta(kdm,nip),	&
       udissp(kdm,nip),vdissp(kdm,nip),ustres(kdm,nip),vstres(kdm,nip),	&
       dfflxu(kdm,nip),dfflxv(kdm,nip),wgt(kdm,nip),worku(kdm,nip),	&
       workv(kdm,nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,3) BEGIN
  real,intent (IN)    :: um_edg  (kdm,npp,nip)	! u velocity on edges, mid-time
  real,intent (IN)    :: vm_edg  (kdm,npp,nip)	! v velocity on edges, mid-time
  real,intent (IN)    :: geop_edg(kdm,npp,nip)	! geopotential on edges
  real,intent (IN)    :: mont_edg(kdm,npp,nip)	! Montgomery potential on edges
  real,intent (IN)    :: uvsq_edg(kdm,npp,nip)	! velocity-squared on edges
  real,intent (IN)    :: spcv_edg(kdm,npp,nip)	! specific volume on edges
!SMS$DISTRIBUTE END
  integer   :: i,ix			! icos point index
  integer   :: k
  integer   :: edg			! icos edge index
  logical   :: vrbos
  character :: text*32,string*24
  real relvor(kdm),sum(kdm),uold(kdm),vold(kdm),thktop,thkbot,		&
       dudx(kdm),dudy(kdm),dvdx(kdm),dvdy(kdm),deform(kdm),smag(kdm),	&
       ubot,vbot,botvel,q,slab,pbot1,p1,boost,hxtpol,arg1,arg2,		&
       valmin,valmax
  real,parameter :: cbar=.5		! tidal flow speed for bottom drag, m/s
  logical,parameter :: lateral=.false.	! if true, extrapolate PGF sideways

  real harmon				! harmonic mean
  harmon(arg1,arg2)=max(onemm,arg1)*max(onemm,arg2)/			&
                   (max(onemm,arg1)+max(onemm,arg2))*2.
  boost(pbot1,p1)=max(1.,1.5-(pbot1-p1)/min(3.*tenm,.125*pbot1))

  thktop=50.*onem		! depth of wind stress penetration
  thkbot=10.*onem		! depth of bottom drag penetration
! hxtpol=10.*onem		! depth threshold for sideways PGF extrapol.
  hxtpol=5.*onem		! depth threshold for sideways PGF extrapol.

  write (text,'(a,i8)') '(ocn_momtum)  step',nstep

!SMS$PARALLEL (dh,i) BEGIN

!sms$compare_var(uclin, "hycom_momtum.F90 - uclin1 ")
!sms$compare_var(vclin, "hycom_momtum.F90 - vclin1 ")
!sms$compare_var(utrop, "hycom_momtum.F90 - utrop1 ")
!sms$compare_var(vtrop, "hycom_momtum.F90 - vtrop1 ")
!sms$compare_var(um_edg, "hycom_momtum.F90 - um_edg1 ")
!sms$compare_var(vm_edg, "hycom_momtum.F90 - vm_edg1 ")
!sms$compare_var(mont_edg, "hycom_momtum.F90 - mont_edg1 ")
!sms$compare_var(dp, "hycom_momtum.F90 - dp1 ")

!$OMP PARALLEL DO PRIVATE (ubot,vbot,q,botvel)
  do i=1,nip			! horizontal loop
   utotn(:,i)=0.
   vtotn(:,i)=0.
   pgfx (:,i)=0.
   pgfy (:,i)=0.
   gradx(:,i)=0.
   grady(:,i)=0.
   wgt  (:,i)=0.
   worku(:,i)=land_spval
   workv(:,i)=land_spval
   if (wet(i) > 0 ) then
    do k=1,kdm
     pres(k+1,i)=pres(k,i)+dp(k,i,leapn)
    end do

! --- determine near-bottom velocity for bottom drag calculation

    ubot=0.
    vbot=0.
    do k=1,kdm
     utotn(k,i)=uclin(k,i,leapn)+utrop(leapn,i)
     vtotn(k,i)=vclin(k,i,leapn)+vtrop(leapn,i)
     q=min(thkbot,pres(kdm+1,i))
     ubot=ubot+utotn(k,i)*(min(pres(k+1,i),pres(kdm+1,i)-q)	&
                          -min(pres(k  ,i),pres(kdm+1,i)-q))
     vbot=vbot+vtotn(k,i)*(min(pres(k+1,i),pres(kdm+1,i)-q)	&
                          -min(pres(k  ,i),pres(kdm+1,i)-q))
    end do
    botvel=sqrt(ubot*ubot+vbot*vbot)/thkbot
    drag(i)=min(drcoef*(botvel+cbar)*onem/thkbot,1./batrop)	! units: 1/s
    ustarb(i)=sqrt(drcoef)*botvel

    do edg=1,nprox(i)		! loop through edges
     do k=1,kdm

      if (janjic_ocn) then
! --- integrate mid-lyr geopotential around cell edge to get pressure force
       pgfx(k,i) = pgfx(k,i)-(geop_edg(k,edg,i)+uvsq_edg(k,edg,i))	&
              *sidevec_c(2,edg,i)
       pgfy(k,i) = pgfy(k,i)+(geop_edg(k,edg,i)+uvsq_edg(k,edg,i))	&
              *sidevec_c(1,edg,i)
!      pgfx(k,i) = pgfx(k,i)-geop_edg(k,edg,i)*sidevec_c(2,edg,i)
!      pgfy(k,i) = pgfy(k,i)+geop_edg(k,edg,i)*sidevec_c(1,edg,i)

      else
! --- integrate bernoulli fcn and density around cell edge to get pressure force
       pres(k+1,i)=pres(k,i)+dp(k,i,leapm)
! --- pressure gradient, term 1:  grad (montg + v^2/2)
       pgfx(k,i)=pgfx(k,i)-(mont_edg(k,edg,i)+uvsq_edg(k,edg,i))	&
              *sidevec_c(2,edg,i)				! m^3/sec^2
       pgfy(k,i)=pgfy(k,i)+(mont_edg(k,edg,i)+uvsq_edg(k,edg,i))	&
              *sidevec_c(1,edg,i)				! m^3/sec^2

! --- pressure gradient, term 2:  p grad alpha
       pgfx(k,i)=pgfx(k,i)+.5*(pres(k,i)+pres(k+1,i))			&
                 *spcv_edg(k,edg,i)*sidevec_c(2,edg,i)		! m^3/sec^2
       pgfy(k,i)=pgfy(k,i)-.5*(pres(k,i)+pres(k+1,i))			&
                 *spcv_edg(k,edg,i)*sidevec_c(1,edg,i)		! m^3/sec^2
      end if			! janjic
     end do
    end do
    do k=1,kdm
     gradx(k,i)=pgfx(k,i)*rarea(i)				! m/sec^2
     grady(k,i)=pgfy(k,i)*rarea(i)				! m/sec^2

! --- in preparation for sideways extrapolation, weight PGF by layer thickness
     if (lateral) then
      wgt(k,i)=max(0.,min(dp(k,i,leapm),hxtpol))
      pgfx(k,i)=gradx(k,i)*wgt(k,i)
      pgfy(k,i)=grady(k,i)*wgt(k,i)
     end if
     worku(k,i)=uclin(k,i,leapn)
     workv(k,i)=vclin(k,i,leapn)
    end do
   end if			! ocean point
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  if (mod(nstep,diag_intvl).eq.0) then
   call stencl(utotn,kdm,1.e3,trim(text)//'  utotn (mm/s)')
   call stencl(vtotn,kdm,1.e3,trim(text)//'  vtotn (mm/s)')
   call stencl(gradx,kdm,dtleap*1.e4,trim(text)//'  pgfx x 10^4')
   call stencl(grady,kdm,dtleap*1.e4,trim(text)//'  pgfy x 10^4')
   if (lateral)				&
   call stencl(wgt,kdm,100./onem,trim(text)//'  PGF extrapol.weight x 100')
  end if

! --- dissipation
  call del4prep(nstep,worku,workv)	! comment out for Laplacian dissip
  call dissip(nstep,leapn,worku,workv,dp)

!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE(utotn,vtotn,pgfx,pgfy,dp,pres,wgt)
!$OMP PARALLEL DO PRIVATE (sum,dudx,dudy,dvdx,dvdy,vrbos,ix,q,		&
!$OMP relvor,deform,smag,uold,vold,slab)
  do i=1,nip			! horizontal loop
   ustres(:,i)=0.
   vstres(:,i)=0.
   udissp(:,i)=0.
   vdissp(:,i)=0.
   dfflxu(:,i)=0.
   dfflxv(:,i)=0.
   uzeta(:,i)=0.
   vzeta(:,i)=0.
!  aux1(:,:,i)=0.
!  aux2(:,:,i)=0.
!  aux3(:,:,i)=0.
!  aux4(:,:,i)=0.
   if (lateral) then
    gradx(:,i)=0.
    grady(:,i)=0.
    sum(:)=0.
   end if
   dudx(:)=0.
   dudy(:)=0.
   dvdx(:)=0.
   dvdy(:)=0.
   if (wet(i) > 0 ) then
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0

    if (vrbos) then
     write (*,100) nstep,perm(i)
 100 format ('ocn_momtum     time step',i8,'  i =',i8/		&
      ' k  uold  unew gradp corio dissp stres',2x,		&
      '    vold  vnew gradp corio dissp stres')
    end if

    do edg=1,nprox(i)		! loop through edges
     ix=prox(edg,i)
     do k=1,kdm

! --- integrate tangential velocity along edge to find velocity derivatives
      dudx(k)=dudx(k)+sidevec_e(2,edg,i)*um_edg(k,edg,i)	! m^2/sec
      dudy(k)=dudy(k)-sidevec_e(1,edg,i)*um_edg(k,edg,i)	! m^2/sec
      dvdx(k)=dvdx(k)+sidevec_e(2,edg,i)*vm_edg(k,edg,i)	! m^2/sec
      dvdy(k)=dvdy(k)-sidevec_e(1,edg,i)*vm_edg(k,edg,i)	! m^2/sec

! --- sideways extrapolation of pressure force to (near-)massless grid points
      if (lateral) then
       sum(k)=sum(k)+wgt(k,ix)
       gradx(k,i)=gradx(k,i)+pgfx(k,ix)
       grady(k,i)=grady(k,i)+pgfy(k,ix)
      end if
     end do			! vertical loop
    end do			! loop through edges

    do k=1,kdm
     if (lateral) then
      sum(k)=max(sum(k),onem)
      q=hxtpol-wgt(k,i)
      gradx(k,i)=(pgfx(k,i)+q*gradx(k,i)/sum(k))/hxtpol		! m^3/sec^2
      grady(k,i)=(pgfy(k,i)+q*grady(k,i)/sum(k))/hxtpol		! m^3/sec^2
     end if
     dudx(k)=dudx(k)*rarea(i)					! 1/sec
     dudy(k)=dudy(k)*rarea(i)					! 1/sec
     dvdx(k)=dvdx(k)*rarea(i)					! 1/sec
     dvdy(k)=dvdy(k)*rarea(i)					! 1/sec
    end do
  
! --- compute vorticity and total deformation
    do k=1,kdm
     relvor(k)=dvdx(k)-dudy(k)					! 1/sec
     deform(k)=sqrt((dudx(k)-dvdy(k))**2+(dudy(k)+dvdx(k))**2)	! 1/sec
     smag(k)=deform(k)*meshsz(i)*smagcf		! smagorinsky diffusion "velo." 

! --- lateral momentum dissipation fluxes
     udissp(k,i)=worku(k,i)*max(veldff,smag(k))			! m/sec^2
     vdissp(k,i)=workv(k,i)*max(veldff,smag(k))			! m/sec^2

! --- coriolis term

     uzeta(k,i)=uclin(k,i,leapm)*(relvor(k)+corio(i))		&
                   +utrop(leapm,i)* relvor(k)
     vzeta(k,i)=vclin(k,i,leapm)*(relvor(k)+corio(i))		&
                   +vtrop(leapm,i)* relvor(k)

! --- add wind stress

     ustres(k,i)=taux(i)*grvity/thktop				&
      *(min(pres(k+1,i)+onemm,thktop)				&
       -min(pres(k  ,i)      ,thktop))/max(dp(k,i,leapn),onemm)
     vstres(k,i)=tauy(i)*grvity/thktop				&
      *(min(pres(k+1,i)+onemm,thktop)				&
       -min(pres(k  ,i)      ,thktop))/max(dp(k,i,leapn),onemm)

! --- add bottom drag

     q=(max(pres(kdm+1,i)-thkbot,    pres(k+1,i))		&
       -max(pres(kdm+1,i)-thkbot,min(pres(k+1,i)-onemm,		&
            pres(k,i))))/max(dp(k,i,leapn),onemm)
     ustres(k,i)=ustres(k,i)-drag(i)*utotn(k,i)*q
     vstres(k,i)=vstres(k,i)-drag(i)*vtotn(k,i)*q

! --- advance velocity field in time
! --- u/v tendcy is the sum of prs.grad., coriolis term, and stresses

     uold(k)=uclin(k,i,leapn)
     vold(k)=vclin(k,i,leapn)
     uclin(k,i,leapn)=uclin(k,i,leapn)+dtleap*			&
        (gradx(k,i)+vzeta(k,i)+udissp(k,i)+ustres(k,i))
     vclin(k,i,leapn)=vclin(k,i,leapn)+dtleap*			&
        (grady(k,i)-uzeta(k,i)+vdissp(k,i)+vstres(k,i))

! --- downward momentum mixing into thin bottom layers
     if (k.gt.1) then
      slab=max(tencm,pres(k+1,i)-pres(kdm+1,i)+tenm)
      q=min(dp(k,i,leapn),slab)
      uclin(k,i,leapn)=(uclin(k  ,i,leapn)*      q		&
                       +uclin(k-1,i,leapn)*(slab-q))/slab
      vclin(k,i,leapn)=(vclin(k  ,i,leapn)*      q		&
                       +vclin(k-1,i,leapn)*(slab-q))/slab
     end if

     if (vrbos) then
      write (*,101) k,						&
      uold(k),uclin(k,i,leapn),gradx(k,i)*dtleap,		&
      vzeta(k,i)*dtleap,udissp(k,i)*dtleap,ustres(k,i)*dtleap,	&
      vold(k),vclin(k,i,leapn),grady(k,i)*dtleap,		&
      -uzeta(k,i)*dtleap,vdissp(k,i)*dtleap,vstres(k,i)*dtleap
      call flush(6)
 101 format (i2,2p,2f6.1,4f6.2,4x,2f6.1,4f6.2)
     end if

    end do			! vertical loop
    if (vrbos) print '(i8,a/(10f8.3))',perm(i),' (ocn_momtum) smagcoeff:',smag

! --- subtract barotropic component from velocity field; use it
! --- to force barotropic flow

    if (modesplit) then
     ubforc(i)=0.
     vbforc(i)=0.
     do k=1,kdm
      ubforc(i)=ubforc(i)+uclin(k,i,leapn)*dp(k,i,leapn)
      vbforc(i)=vbforc(i)+vclin(k,i,leapn)*dp(k,i,leapn)
     end do			! vertical loop

     ubforc(i)=ubforc(i)/pres(kdm+1,i)
     vbforc(i)=vbforc(i)/pres(kdm+1,i)
     do k=1,kdm
      uclin(k,i,leapn)=uclin(k,i,leapn)-ubforc(i)
      vclin(k,i,leapn)=vclin(k,i,leapn)-vbforc(i)
     end do
    end if

! --- time smoothing

    do k=1,kdm
     uclin(k,i,leapm)=assluv*(uold(k)+uclin(k,i,leapn))	&
          +(1.-assluv-assluv)*        uclin(k,i,leapm)
     vclin(k,i,leapm)=assluv*(vold(k)+vclin(k,i,leapn))	&
          +(1.-assluv-assluv)*        vclin(k,i,leapm)
    end do			! vertical loop

    if (uclin(kdm,i,leapn)**2+vclin(kdm,i,leapn)**2.gt.9.)		&
     print '(2i8,a,2f9.2)',nstep,perm(i),' excessive bottom velocity',	&
      uclin(kdm,i,leapn),vclin(kdm,i,leapn)
   end if			! ocean point

  end do			! horizontal loop
!$OMP END PARALLEL DO

!sms$compare_var(uzeta, "hycom_momtum.F90 - uzeta2 ")
!sms$compare_var(vzeta, "hycom_momtum.F90 - vzeta2 ")
!sms$compare_var(dfflxu, "hycom_momtum.F90 - dfflxu2 ")
!sms$compare_var(dfflxv, "hycom_momtum.F90 - dfflxv2 ")
!sms$compare_var(udissp, "hycom_momtum.F90 - udissp2 ")
!sms$compare_var(vdissp, "hycom_momtum.F90 - vdissp2 ")
!sms$compare_var(ustres, "hycom_momtum.F90 - ustres2 ")
!sms$compare_var(vstres, "hycom_momtum.F90 - vstres2 ")
!sms$compare_var(gradx, "hycom_momtum.F90 - gradx2 ")
!sms$compare_var(grady, "hycom_momtum.F90 - grady2 ")
!sms$compare_var(uclin, "hycom_momtum.F90 - uclin2 ")
!sms$compare_var(vclin, "hycom_momtum.F90 - vclin2 ")
!SMS$PARALLEL END

!   valmin=minval(uclin)
!   valmax=maxval(uclin)
! !SMS$REDUCE(valmin,MIN)
! !SMS$REDUCE(valmax,MAX)
!   print 102,'uclin',valmin,valmax
!   valmin=minval(vclin)
!   valmax=maxval(vclin)
! !SMS$REDUCE(valmin,MIN)
! !SMS$REDUCE(valmax,MAX)
!   print 102,'vclin',valmin,valmax
!   valmin=minval(dp)
!   valmax=maxval(dp)
! !SMS$REDUCE(valmin,MIN)
! !SMS$REDUCE(valmax,MAX)
!   print 102,'thkns',valmin,valmax
102 format ('(ocn_momtum) min,max of ',a,2es15.7)

! if (mod(nstep,diag_intvl).eq.0) then
!  if (lateral) then
!   call stencl(gradx,kdm,dtleap*1.e4,trim(text)//'  weighted pgfx x 10^4')
!   call stencl(grady,kdm,dtleap*1.e4,trim(text)//'  weighted pgfy x 10^4')
!  end if

!  do k=1,kdm
!   write (string,'(a,i2)') 'um velocity k =',k
!   call findmxmn3(uclin,kdm,nip,2,k,leapm,string,wet)
!   write (string,'(a,i2)') 'vm velocity k =',k
!   call findmxmn3(vclin,kdm,nip,2,k,leapm,string,wet)
!   write (string,'(a,i2)') 'un velocity k =',k
!   call findmxmn3(uclin,kdm,nip,2,k,leapn,string,wet)
!   write (string,'(a,i2)') 'vn velocity k =',k
!   call findmxmn3(vclin,kdm,nip,2,k,leapn,string,wet)
!   write (string,'(a,i2)') 'lyrm thknss k =',k
!   call findmxmn3(dp,kdm,nip,2,k,leapm,string,wet)
!   write (string,'(a,i2)') 'lyrn thknss k =',k
!   call findmxmn3(dp,kdm,nip,2,k,leapn,string,wet)
!   write (string,'(a,i2)') 'intfc pres k =',k+1
!   call findmxmn2(pres,kdm+1,nip,k+1,string,wet)
!  end do
!  call flush(6)
!  call stencl(uzeta,1,1.e5,trim(text)//'  uzeta x 10^5')
!  call stencl(vzeta,1,1.e5,trim(text)//'  vzeta x 10^5')
!  call stencl(ustres,1,1.e3,trim(text)//'  ustres x 10^3')
!  call stencl(vstres,1,1.e3,trim(text)//'  vstres x 10^3')
! end if

  return
  end subroutine momtum
end module hycom_momtum
