module hycom_hybgn2
contains
!
! --- ---------------------
! --- hybrid grid generator (based on 'restep' technology)
! --- ---------------------
!
  subroutine hybgn2(nstep,leapn,theta,dens,dp,pres,temp,saln,uvel,vvel,	&
                    latij,passv_tr)

  use module_control  ,only: nip
  use module_constants,only: perm
  use fimnamelist    ,only: kdm,numtr,itest,diag_intvl
  use hycom_constants ,only: wet,epsil,onem,onecm,onemm,salmin,dpmin
  use hycom_sigetc

  implicit none
  integer,intent(IN) :: nstep                   ! model time step
  integer,intent(IN) :: leapn                   ! leapfrog time slot (n='new')
  real   ,intent(IN) :: theta(kdm)              ! target densities
!SMS$DISTRIBUTE  (dh,1) BEGIN
  real,intent(IN)    :: latij     (nip)         ! latitude in degree
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE  (dh,2) BEGIN
  real,intent(INOUT) :: uvel (kdm,nip,2)     ! u velocity
  real,intent(INOUT) :: vvel (kdm,nip,2)     ! v velocity
  real,intent(INOUT) :: dens (kdm,nip  )     ! density
  real,intent(INOUT) :: dp   (kdm,nip,2)     ! layer thickness
  real,intent(INOUT) :: pres (kdm+1,nip)     ! interface pressure
  real,intent(INOUT) :: temp (kdm,nip)       ! temperature
  real,intent(INOUT) :: saln (kdm,nip)       ! salinity
  real,intent(INOUT),optional:: passv_tr(kdm,nip,numtr)
!SMS$DISTRIBUTE END
  logical vrbos,dotrcr
  integer i,k,k1,ko,nt,last,ntot,nwrk,ntotd,nwrkd
  real*8 dpcol(kdm),pcol(kdm+1),delp,tem,sal,dns,			&
         colint,colins,cloutt,clouts,q,q1,q2,abv,dtwarm,dtcool,		&
         tdp,sdp,told,sold,told1,sold1,rold,rold1
  real*8,parameter :: tolrnce=1.e-4, acurcy=1.e-11

  real*8 dp0,dp1,dpsum,uintg,vintg,phi,plo,pa,pb,old,p_hat,		&
       torho,totem,tosal,totrc,totu,totv,tndrho,tndtem,tndsal,tndtrc,	&
       tdcyu,tdcyv,scale,arg,arg1,cushn,pbot0,latij0
  real*8 targt(kdm+1),r1d(kdm),p1d(kdm+1),dp1d(kdm),t1d(kdm),s1d(kdm),	&
       pnew(kdm+1),trac(kdm,numtr),u1d(kdm),v1d(kdm),dpold(kdm),	&
       uold(kdm),vold(kdm)
  real*8,parameter :: uvscl=0.02			!  2 cm/s
!
!cc   parameter (slak=.5/86400.)	! intfc nudging time scale: 2 days
!cc   parameter (slak=1./86400.)	! intfc nudging time scale: 1 day
      real*8, parameter :: slak=2./86400.	! intfc nudging time scale: 12 h
!cc   parameter (slak=4./86400.)	! intfc nudging time scale: 6 hrs
!
! --- linear taper functions (latitude and depth-dependent) for slak
      real*8 slakf
!     real tapr,wgtf
!     tapr(q)=1.+9.*max(0.,1.-.02e-4*q)			! q = pressure (Pa)
!     wgtf(q)=max(0.,min(1.,(abs(q)-50.)*.1))		! 0->1 for q=50->60
!     slakf(q)=0.7*wgtf(q)+tapr(p_hat)*slak*delt1*(1.-wgtf(q))	! q = lat.(deg)
!     slakf(q)=tapr(p_hat)*slak*delt1			! no lat.dependence
      slakf(q)=1.					! no slak whatsoever
!
! --- simplified cushion function suitable for nonnegative first argument
      cushn(arg,arg1)=.25*min(arg,2.*arg1)**2/arg1+max(arg,2.*arg1)-arg1
!     cushn(arg,arg1)=arg1                  ! no cushn

  dotrcr=.false.
! print *,'entering ocn_hybgen, itest =',itest

!SMS$PARALLEL (dh,i) BEGIN
  do i=1,nip
   if (wet(i) > 0 ) then
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
 
    latij0=latij(i)
! --- extract t,s,rho column from 3-d grid 
    do k=1,kdm
     t1d(k)=temp(k,i)
     s1d(k)=saln(k,i)
     r1d(k)=dens(k,i)
     u1d(k)=uvel(k,i,leapn)
     v1d(k)=vvel(k,i,leapn)
     dp1d(k)= dp(k,i,leapn)
     p1d(k)=pres(k,i)
     targt(k)=theta(k)
     if (present(passv_tr)) trac(k,:)=passv_tr(k,i,:)
    end do
    p1d(kdm+1)=pres(kdm+1,i)
    pbot0=p1d(kdm+1)

    if (vrbos) then
     write (*,103) nstep,perm(i),					&
     ' ocn_hybgen  IN: temp   saln   dens  thkns   dpth     u      v',	&
      (k,t1d(k),s1d(k),r1d(k),dp1d(k)/onem,p1d(k+1)/onem,u1d(k)*100.,	&
       v1d(k)*100,k=1,kdm)
 103 format (2i7,a/(i28,3f7.3,f7.2,f7.1,2f7.2))
    end if
!
    ntot=0
    nwrk=0
!
    if (vrbos) then
     write (*,99) nstep,perm(i),'  ocnhyb   o l d   profile'
     do k=1,kdm,10
      write (*,100) (p1d  (k1)/onem,k1=k,min(kdm+1,k+10))
      write (*,101) (r1d  (k1),k1=k,min(kdm,k+9))
      write (*,101) (targt(k1),k1=k,min(kdm,k+9))
      write (*,102) (t1d  (k1),k1=k,min(kdm,k+9))
      write (*,102) (s1d  (k1),k1=k,min(kdm,k+9))
     end do
    end if
 99 format (2i7,a,'  (5-line groups: dpth,dens,targt,T,S)')
100 format (11f7.1)
101 format (4x,10f7.2)
102 format (4x,10f7.2)
!
    torho=0.
    totem=0.
    tosal=0.
    do k=1,kdm
     torho=torho+r1d(k)*(p1d(k+1)-p1d(k))
     totem=totem+t1d(k)*(p1d(k+1)-p1d(k))
     tosal=tosal+s1d(k)*(p1d(k+1)-p1d(k))
    end do
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- absorb near-massless layers on sea floor in layer above
    ko=1
    do 4 k=2,kdm
    if (p1d(k).lt.p1d(kdm+1)-onecm) then
     ko=k		! this layer absorbs subsequent layers
    else
     q1=max(epsil,p1d(k  )-p1d(ko))
     q2=max(   0.,p1d(k+1)-p1d(k ))
     q=q1/(q1+q2)
     if (q.lt.0. .or. q.gt.1.) then
      write (*,*) 'error hybgn2 - i,q1,q2,q=',perm(i),q1,q2,q
     end if
     t1d(ko)=t1d(ko)*q+t1d(k)*(1.-q)
     s1d(ko)=s1d(ko)*q+s1d(k)*(1.-q)
     t1d(k)=t1d(ko)
     s1d(k)=max(s1d(ko),salmin(k))
     if (present(passv_tr)) then
      if (dotrcr) then
       trac(ko,:)=trac(ko,:)*q+trac(k,:)*(1.-q)
       trac(k,:)=trac(ko,:)
      end if
     end if
     if (vrbos) then
      write (*,'(2i7,a,i3,a,i3,5x,a,3f7.3)') nstep,perm(i),		&
       '  absorb lyr',k,' in',ko,'new t,s,rho:',t1d(ko),s1d(ko),	&
        sigocn(t1d(ko),s1d(ko))
     end if
    end if
 4  continue
!
    r1d(ko)=sigocn(t1d(ko),s1d(ko))
    do 11 k=ko+1,kdm
    r1d(k)=r1d(ko)
 11 p1d(k)=p1d(kdm+1)
!
    do 23 k=1,kdm
 23 dpold(k)=p1d(k+1)-p1d(k)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! --- step 1: move interfaces to reduce deviations from target.
    if (vrbos) then
     write (*,'(i7,i8,i3,a/(i15,3f8.3,f8.1))') nstep,perm(i),ko,	&
      ' temp    salt    dens    pres',(k,t1d(k),s1d(k),r1d(k),		&
       p1d(k+1)/onem,k=1,kdm)
    end if

    call restp_1d(r1d,p1d,kdm,r1d,pnew,targt,kdm,vrbos,perm(i))
!
! --- expand thin layers at the expense of thick layers above and below.
    call equilb(r1d,pnew,kdm,vrbos,perm(i))

! --- step 2: enforce minimum layer thickness constraint
!
    dpsum=0.
    last=1
    do 7 k=1,kdm-1
    ntot=ntot+1
!
! --- set lower limits for layer thknss (dp0) and depth of lower intfc (dpsum)
    dp0=dpmin(k)*onem
    if (k.gt.1) then
     if (k.gt.last+1) dp0=dp0 * 2.**(.5*(last-k))
     if (vrbos .and. k.eq.last+2) print '(i5,a,i3)',perm(i),		&
      '  hybrid layers: k = 1 -',last 
!
! --- optional: reduce spacing of z layers near equator, but hide transition
! --- in a subtropical latitude band where z layers are least likely to exist
     dp0=dp0*max(.6,min(1.,(abs(latij0)+5.)*.04))
    end if
!
! --- reduce layer thickness in shallow spots, creating sigma coord. effect
    if (4*k.lt.kdm) dp0=dp0*min(1.,pbot0/(200.*onem)+.4)
!
    dpsum=dpsum+dp0
    dp1=dp0
    if (k.gt.1) dp1=cushn(max(pnew(k+1),dpsum)-pnew(k),dp0)
    if (vrbos) print '(i3,a,3f8.3)',k,				&
     ' min.thknss (dp0;dp1) & min.lowr intfc.dpth',		&
      dp0/onem,dp1/onem,dpsum/onem
    dp0=dp1
!
! --- is lower intfc too close to the surface?
    p_hat=max(dpsum,pnew(k)+dp0)
    if (p_hat.gt.pnew(k+1)) then
     nwrk=nwrk+1
     old=pnew(k+1)
     if (k-1.eq.last) last=k			!  still in z domain
     pnew(k+1)=min(pnew(kdm+1),						&
     pnew(k)+max(onemm,slakf(latij0)*(p_hat-pnew(k))))

     if (vrbos) then
      print '(10x,a/i7,i3,5f11.3)',					&
       '  pold(k+1)     p_hat   pnew(kdm+1)    slakf   pnew(k+1)',	&
        perm(i),k,old/onem,p_hat/onem,pnew(kdm+1)/onem,slakf(latij0),	&
         pnew(k+1)/onem
      write (*,107) perm(i),k,' lower intfc',p1d(k+1)/onem,'=>',	&
        pnew(k+1)/onem
 107  format (i7,i3,1x,2(1x,a,f9.3))
     end if
    end if
 7  dp1d(k  )=pnew(k  +1)-pnew(k  )
    dp1d(kdm)=pnew(kdm+1)-pnew(kdm)
!
! --- vertical advection
!
    if (vrbos) then
     print *,'vertical advection of temperature'
    end if
    call ppmadv(p1d,t1d,pnew,t1d,kdm,vrbos,perm(i))
    if (vrbos) then
     print *,'vertical advection of salinity'
    end if
    call ppmadv(p1d,s1d,pnew,s1d,kdm,vrbos,perm(i))
!
    do 5 k=1,kdm
 5  r1d(k)=sigocn(t1d(k),s1d(k))

    if (present(passv_tr)) then
     if (dotrcr) then
!
! --- evaluate effect of regridding on tracer field(s)
!
      do nt=1,numtr
       scale=1.e-33
       totrc=0.
       do 22 k=1,kdm
       totrc=totrc+trac(k,nt)*(p1d(k+1)-p1d(k))
 22    scale=max(scale,abs(trac(k,nt)))
!
       if (vrbos) then
        print *,'vertical advection of tracer',nt
       end if
       call ppmadv(p1d,trac(1,nt),pnew,trac(1,nt),kdm,vrbos,perm(i))
!
       tndtrc=-totrc
       do 20 k=1,kdm
 20    tndtrc=tndtrc+trac(k,nt)*dp1d(k)
!
       if (abs(tndtrc).gt.acurcy*scale*pnew(kdm+1)) then
        write (*,104) perm(i),'  ocn_hybgen - bad trcr.intgl.:',	&
          totrc,tndtrc,tndtrc/(scale*pnew(kdm+1))
       end if
      end do				!  numtr
     end if				!  dotrcr
    end if				!  dotrcr
!
    tndrho=-torho
    tndtem=-totem
    tndsal=-tosal
    do k=1,kdm
     tndrho=tndrho+r1d(k)*dp1d(k)
     tndtem=tndtem+t1d(k)*dp1d(k)
     tndsal=tndsal+s1d(k)*dp1d(k)
    end do
    if (abs(tndtem).gt.acurcy*30.*p1d(kdm+1))				&
     write (*,104) perm(i),'  ocn_hybgen - bad temp.intgl.:',totem,	&
      tndtem,tndtem/(10.*p1d(kdm+1))
    if (abs(tndsal).gt.acurcy*40.*p1d(kdm+1))				&
     write (*,104) perm(i),'  ocn_hybgen - bad saln.intgl.:',tosal,	&
      tndsal,tndsal/(35.*p1d(kdm+1))
 104 format (i7,a,2es15.7,es9.1)
!
    if (vrbos) then
     write (*,99) nstep,perm(i),'  ocnhyb   n e w   profile'
     do k=1,kdm,10
      write (*,100) (pnew (k1)/onem,k1=k,min(kdm+1,k+10))
      write (*,101) (r1d  (k1),k1=k,min(kdm,k+9))
      write (*,101) (targt(k1),k1=k,min(kdm,k+9))
      write (*,102) (t1d  (k1),k1=k,min(kdm,k+9))
      write (*,102) (s1d  (k1),k1=k,min(kdm,k+9))
     end do
    end if
!
! --- reduce salt content wherever density > targt(kdm)
!
!   do 31 k=1,kdm
!   if (pnew(k).lt.p1d(kdm+1)) then
!    q=sofsig(targt(kdm),t1d(k))
!    if (q.lt.s1d(k)) then
!     if (vrbos .or. s1d(k)-q.gt..2)				&
!      write (*,'(i9,i4,i3,4(a,f7.3))') nstep,i,k,		&
!       ' reduce s=',s1d(k),' ->',q,' to reduce rho=',		&
!        r1d(k),' ->',targt(kdm)
!     s1d(k)=q
!     r1d(k)=targt(kdm)
!    end if
!   end if
!31 continue
!
!   ntotd(j)=ntot
!   nwrkd(j)=nwrk
!
! --- temporarily store old layer thickness at -u,v- points in -pu,pv-
!
    do 1 k=1,kdm
 1  p1d(k+1)=p1d(k)+dpold(k)		!  old pressure
!
    totu=0.
    totv=0.
    uold(:)=u1d(:)
    vold(:)=v1d(:)
    do 15 k=1,kdm
    totu=totu+uold(k)*(p1d(k+1)-p1d(k))
 15 totv=totv+vold(k)*(p1d(k+1)-p1d(k))
    tdcyu=-totu
    tdcyv=-totv
!
    do 18 k=1,kdm
    phi=pnew(k+1)
    plo=pnew(k  )
    if (phi.gt.plo) then
     uintg=0.
     vintg=0.
     pb=plo
     do 16 ko=1,kdm
     if (p1d(ko+1).le.plo) go to 16
     pa=pb
     pb=min(phi,p1d(ko+1))
     uintg=uintg+uold(ko)*(pb-pa)
     vintg=vintg+vold(ko)*(pb-pa)
     if (pa.ge.phi) go to 17
 16  continue
 17  tdcyu=tdcyu+uintg
     tdcyv=tdcyv+vintg
     u1d(k)=uintg/(phi-plo)
     v1d(k)=vintg/(phi-plo)
    end if
 18 continue
!
    if (abs(tdcyu).gt.acurcy*uvscl*p1d(kdm+1))				&
     write (*,104) perm(i),'  hybgen - bad u intgl.',totu,		&
      tdcyu,tdcyu/(uvscl*p1d(kdm+1))
    if (abs(tdcyv).gt.acurcy*uvscl*p1d(kdm+1))				&
     write (*,104) perm(i),'  hybgen - bad v intgl.',totv,		&
      tdcyv,tdcyv/(uvscl*p1d(kdm+1))
 14 continue
!
    if (vrbos) then
     write (*,103) nstep,perm(i),					&
     ' ocn_hybgen OUT: temp   saln   dens  thkns   dpth     u      v',	&
      (k,t1d(k),s1d(k),r1d(k),dp1d(k)/onem,pnew(k+1)/onem,u1d(k)*100.,	&
       v1d(k)*100,k=1,kdm)
    end if
!
! need work
!  if (mod(time+.0001,1.).lt..0002) then
!   nwrk=0
!   ntot=0
!   nwrk=nwrk+nwrkd(j)
!   ntot=ntot+ntotd(j)
!   end do
!   write (*,'(a,f6.1,a,i9,a)') 'hybgen - grid restoration at',		&
!    100.*float(nwrk)/float(ntot),' per cent of',ntot,' points'
!  end if
!
! --- put 1-d column back into 3-d grid
    do k=1,kdm
     temp(k,i)   = t1d(k)
     saln(k,i)   = s1d(k)
     dens(k,i)      = r1d(k)
     uvel(k,i,leapn)= u1d(k)
     vvel(k,i,leapn)= v1d(k)
     dp  (k,i,leapn)=pnew(k+1)-pnew(k)
     pres(k,i)      =pnew(k)
!    diaflx(k)=diaflx(k)+(dp1d(k)-dpold(k))	!  diapyc.flux
     if (present(passv_tr)) passv_tr(k,i,:)=trac(k,:)
    end do
    pres(kdm+1,i)=pnew(kdm+1)
   end if                        ! ocean point
  end do                        ! horiz. loop
!SMS$PARALLEL END

! print *,'...exiting ocn_hybgen'
  return
  end subroutine hybgn2
 
 
   subroutine restp_1d(thold,pold,kold,thnew,pnew,targt,knew,vrbos,i)
! --- convert a stairstep (i.e., piecewise constant) theta profile into a
! --- stairstep profile constrained to have prescribed theta ('targt') steps.
       
! --- input   variables: thold,pold,targt,kold,knew,vrbos
! --- output variables: thnew,pnew

   implicit none
   integer,intent(IN)    :: kold,knew,i
   real*8 ,intent(IN)    :: thold(kold),pold(kold+1),targt(knew)
   real*8 ,intent(INOUT) :: thnew(knew),pnew(knew+1)
   logical,intent(IN)    :: vrbos

   integer k,ko
   real*8 oldth(kold)
   real*8 scale,cloutt,colint,pinteg,tha,thb
   real*8, parameter :: acurcy=1.e-11
       
   if (vrbos) then
     write (6,101) i,							&
      'restp1 -- input profile:    theta     thknss      press',	&
      pold(1),(k,thold(k),pold(k+1)-pold(k),pold(k+1),k=1,kold)
101 format (i7,2x,a/54x,f11.1/(i32,f11.3,2f11.1))
   end if
       
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
   if (abs(cloutt-colint).gt.acurcy*scale*pnew(knew+1)) then
     write (6,100) i,'restp1 - column intgl.error',		&
      colint,cloutt,(cloutt-colint)/colint
 100     format (i7,3x,a,2es14.6,es9.1)
   end if     

   if (vrbos) then
     write (6,101) i,							&
      'restp1 -- outpt profile:    theta     thknss      press',	&
      pnew(1),(k,thold(k),pnew(k+1)-pnew(k),pnew(k+1),k=1,knew)
   end if

   return
   end subroutine restp_1d


   subroutine ppmadv(xold,fldold,xnew,fldnew,kk,vrbos,i)

!-- PPM-based 1-dim transport routine, extracted from HYCOM's fct3d.f.
!-- Note: |xold-xnew| can exceed cell width (i.e., no CFL constraints)

!-- xold/new	- old/new cell boundaries
!-- fldold/new	- mixing ratio of dep.variable before and after transport
!-- kk		- number of layers
!-- vrbos	- if .true., print diagnostic messages for grid point -i-

   implicit none
   integer,intent(IN)    :: kk,i
   real*8 ,intent(IN)    :: xold(kk+1),xnew(kk+1),fldold(kk)
   real*8 ,intent(INOUT) :: fldnew(kk)
   logical,intent(IN)    :: vrbos	!  switch for 'verbose' mode

   real*8 zold(kk+1),znew(kk+1),delx(kk+1),delz(kk+1),fco(kk),fcn(kk),	&
        vertfx(kk+1),vertdv(kk)
   real*8 a(kk),b(kk),c(kk),dx,fcdx,yl,yr
   real*8 amount,bfore,after,dpth,scale,slab,dslab
   integer k,lyr
   real*8,parameter :: athird=1./3.
   real*8,parameter :: small=1.e-11
   real*8,parameter :: acurcy=1.e-11

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
100 format (i8,3x,a/(20x,i3,2f11.1,es11.3))
   end if

!-- deduce old and new cell width from -zold,znew-
   do 15 k=1,kk
   fco(k)=max(0.,zold(k+1)-zold(k))
15 fcn(k)=max(0.,znew(k+1)-znew(k))

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
17 fldnew(k)=(fldnew(k)*fco(k)+fldnew(k+1)*small)		&
            /(          fco(k)+            small)
   do 18 k=2,kk
18 fldnew(k)=(fldnew(k)*fco(k)+fldnew(k-1)*small)		&
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
!cc   a(k)=fldnew(k)
!cc   b(k)=0.
!cc   c(k)=0.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- piecewise linear method:
!-- fit linear function a+bx to tracer in each cell (-.5 < x < +.5)
!cc   a(k)=fldnew(k)
!cc   b(k)=0.
!cc   if (fldnew(k).le.min(fldnew(k-1),fldnew(k+1)) .or.		&
!cc       fldnew(k).ge.max(fldnew(k-1),fldnew(k+1))) then
!cc     b(k)=0.
!cc   else if ((fldnew(k+1)-fldnew(k-1))*(fldnew(k-1)+fldnew(k+1)	&
!cc     -2.*fldnew(k)).gt.0.) then
!cc     b(k)=fldnew(k)-fldnew(k-1)
!cc   else
!cc     b(k)=fldnew(k+1)-fldnew(k)
!cc   end if
!cc   c(k)=0.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!-- piecewise parabolic method:
!-- construct parabola  a+bx+cx^2  whose integral over [-.5,+.5] equals
!-- fldnew(k) and which passes though points yl,yr at [-.5,+.5] resp.
!!   yl=.5*(fldnew(k-1)+fldnew(k))
!!   yr=.5*(fldnew(k+1)+fldnew(k))
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
16 continue

    if (vrbos) then
     print '(a,i7,a/(i17,2x,3es12.3))','ipn=',i,		&
      '  ppm coeff.  a           b           c',		&
      (k,a(k),b(k),c(k),k=1,kk)
    end if

!-- get flux by summing -fldnew- over upstream slab of thickness -delz-

   do 22 k=2,kk
   slab=0.
   amount=0.
   vertfx(k)=0.
   if (delz(k).gt.0.) then			! interface moves in +k dir.
     lyr=k-1
24   lyr=lyr+1
     if (slab.ge.delz(k)) goto 23
     if (fco(lyr).gt.0.) then
       dslab=min(slab+fco(lyr), delz(k))	&
            -min(slab         , delz(k))
       dx=dslab/fco(lyr)
       fcdx=a(lyr)				&
           +b(lyr)*.5*(dx-1.)			&!  not needed in pcm
           +c(lyr)*(.25-dx*(.5-dx*athird))	 !  not needed in pcm,plm
       amount=amount+fcdx*dslab
       slab=slab+dslab
     end if
     if (lyr.lt.kk) go to 24
   else if (delz(k).lt.0.) then			! interface moves in -k dir.
     lyr=k
25   lyr=lyr-1
     if (slab.ge.-delz(k)) goto 23
     if (fco(lyr).gt.0.) then
       dslab=min(slab+fco(lyr),-delz(k))	&
            -min(slab         ,-delz(k))
       dx=dslab/fco(lyr)
       fcdx=a(lyr)				&
           +b(lyr)*.5*(1.-dx)			&!  not needed in pcm
           +c(lyr)*(.25-dx*(.5-dx*athird))	 !  not needed in pcm,plm
       amount=amount+fcdx*dslab
       slab=slab+dslab
     end if					! delz < or > 0
     if (lyr.gt.1) go to 25
   end if
23 if (slab.ne.0.) vertfx(k)=-delz(k)*amount/slab
22 continue

   vertfx(   1)=0.			!  don't allow flux through lower bdry
   vertfx(kk+1)=0.			!  don't allow flux through upper bdry
   do 26 k=1,kk
26 vertdv(k)=vertfx(k+1)-vertfx(k)

   if (vrbos) then
    write (*,'(a/(i3,4es12.3))')					&
    'ppmadv:   flux  flx.div/thk    old_thk     new_thk',		&
     (k,vertfx(k),vertdv(k)/max(small,fcn(k)),fco(k),fcn(k),k=1,kk),	&
      kk+1,vertfx(kk+1)
   end if

   do 4 k=1,kk
   if (fcn(k).gt.0.) then
    amount=fldnew(k)*fco(k)-vertdv(k)
!   if (abs(amount).lt.100.*small*abs(vertdv(k))) amount=0.
    fldnew(k)=(fldnew(k)*small+amount)/(small+fcn(k))
   end if
 4 continue

   after=0.
   do k=1,kk
     after=after+fldnew(k)*fcn(k)
   end do

   if (abs(bfore-after).gt.acurcy*scale*dpth) then
     write (*,104) i,'ppmadv - bad column intgl.:',bfore,after
104 format (i8,2x,a,2es15.7)
   end if

   if (vrbos) then
    write (*,100) i,'exiting ppmadv:   d(x)      new_x   variable',	&
     (k,delx(k),xnew(k),fldnew(k),k=1,kk),kk+1,delx(kk+1),xnew(kk+1)
   end if

   return
   end subroutine ppmadv
 
 
   subroutine equilb(thet,pres,kk,vrbos,i)
! --- expand thin layers at the expense of thick layers above and below.
! --- do this without changing theta
   use hycom_constants,only: onem,onemm,onemu

   implicit none
   integer,intent(IN) :: kk,i	! no. of layers, test point location
   logical,intent(IN) :: vrbos	! if true, print results at test point
   real*8,intent(IN)    :: thet(kk)	! pot.temp. in layers
   real*8,intent(INOUT) :: pres(kk+1)	! Exner fcn on interfaces
   integer k,k1,k2,ncount,iter
   real*8 dp1,dp2,dp3,dp4,dp5,th1,th2,th3,th4,th5,dis2,dis3,dis4,	&
      ratio,goal,bfore,after,pnew(kk+1),pwrk(kk+1),heatfx(kk),cnv
   logical event
   real*8,parameter :: slak=.5		! retardation coefficient
   real*8,parameter :: dffudt=.1	! therm.diffu.coeff x time step [m^2]

   if (vrbos) then
    write (6,99) i,'  ocn_equilb input profile:'
    do k2=1,kk,10
     write (6,100) (pres(k1),k1=k2,min(kk+1,k2+10) )
     write (6,101) (thet(k1),k1=k2,min(kk  ,k2+9 ) )
     write (6,100)
    end do
99  format ('i=',i8,a)
100 format (-4p,11f7.1)
102 format (11i7)
101 format (5x,10f7.2)
   end if

   ncount=0
   pnew(:)=pres(:)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- invoke heat diffusion (McDougall-Dewar) to inflate thin layers

!  heatfx(kk)=0.
!  heatfx( 1)=0.
!  do iter=1,3				! apply the scheme repeatedly
!   pwrk(:)=pnew(:)
!   do k=2,kk-1
!    heatfx(k)=0.
!    if (pwrk(k  ).gt.pwrk(   1)+onem .and.
!  .     pwrk(k+1).lt.pwrk(kk+1)-onem) then
!     cnv=onem				! pressure units per meter
!     heatfx(k)=dffudt*cnv*cnv*.5*(thet(k+1)-thet(k-1))
!  .   /max(.01*cnv,pwrk(k+1)-pwrk(k))
!    end if
!   end do
!   do k=1,kk-1
!    if (pwrk(k+1).gt.pwrk(   1)+onem .and.
!  .     pwrk(k+1).lt.pwrk(kk+1)-onem)
!  .  pnew(k+1)=min(.5*(pwrk(k+2)+pwrk(k+1)),
!  .            max(.5*(pwrk(k  )+pwrk(k+1)),
!  .    pwrk(k+1)+(heatfx(k)-heatfx(k+1))/(thet(k+1)-thet(k))))
!   end do
!   event=.false.
!   do k=2,kk+1
!    if (pnew(k).lt.pnew(k-1)-onemu) then
!     event=.true.
!     print '(i5,i3,a)',i,k,' dp<0 due to heat diffusion'
!    end if
!   end do

!   if (vrbos .or. event) then
!    print '(i5,a,i2,a)',i,' heat diffusion, iter',iter,
!  .  '  pres (3-line groups: old,new,dif x 10^4)'
!    do k2=1,kk,10
!     write (6,108) (pwrk(k1),k1=k2,min(kk+1,k2+9) )
!     write (6,108) (pnew(k1),k1=k2,min(kk+1,k2+9) )
!     write (6,109) (int(pnew(k1)-pwrk(k1)),k1=k2,
!  .    min(kk+1,k2+9) )
!     write (6,108)
!    end do
!108  format (-4p,10f8.2)
!109  format (10i8)
!   end if
!   if (.not.event) exit
!  end do                        ! iter
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   bfore=0.
   do k=1,kk
    bfore=bfore+(pnew(k+1)-pnew(k))*thet(k)
   end do

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- scenario 1: sequence of 5 thin-thick-thin-thick-thin layers

   do 1 k=4,kk-3
   if (pnew(k-2).lt.pnew(   1)+onemm .or.			&
       pnew(k+3).gt.pnew(kk+1)-onemm) go to 1
   dp1=pnew(k-1)-pnew(k-2)
   dp2=pnew(k  )-pnew(k-1)
   dp3=pnew(k+1)-pnew(k  )
   dp4=pnew(k+2)-pnew(k+1)
   dp5=pnew(k+3)-pnew(k+2)
   th1=thet(k-2)
   th2=thet(k-1)
   th3=thet(k  )
   th4=thet(k+1)
   th5=thet(k+2)
! --- look for small dp1,dp3,dp5 in combination with large dp2,dp4
   if (dp2.gt.dp1 .and. dp4.gt.dp5) then
    goal=.5*(dp3+min(dp2,dp4))		!  desired thknss of lyr 3
    if (dp2.gt.goal .and. dp4.gt.goal) then
! --- thin-thick-thin-thick-thin combination found -> inflate lyr 3
     dis3=min(dp2-goal,goal-dp3) * slak
     dis4=min(dp4-goal,goal-dp3) * slak
     if (th3.gt.th2 .and. th4.gt.th3) then
! --- theta conservation requires  dis3*(th3-th2)+dis4*(th4-th3)=0
      ratio=(th4-th3)/(th3-th2)
      if (ratio.gt.1.) then
       dis4=dis3/ratio
      else
       dis3=dis4*ratio
      end if
     end if
! --- ready to expand middle layer
     pnew(k  )=pnew(k  )-dis3
     pnew(k+1)=pnew(k+1)+dis4

     if (vrbos) then
      write (6,'(a,i3,a,-4p,5f9.3/a,5f9.3)') 'k=',k,		&
       ' thknss quintuplet',dp1,dp2,dp3,dp4,dp5,		&
        '                becomes',dp1,pnew(k)-pnew(k-1),	&
         pnew(k+1)-pnew(k),pnew(k+2)-pnew(k+1),dp5
     end if
     ncount=ncount+1
    end if
   end if
 1 continue

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- scenario 2: sequence of 3 thin-thick-thin layers

   do 2 k=3,kk-2
   if (pnew(k-1).lt.pnew(   1)+onemm .or.			&
       pnew(k+2).gt.pnew(kk+1)-onemm) go to 2
   dp2=pnew(k  )-pnew(k-1)
   dp3=pnew(k+1)-pnew(k  )
   dp4=pnew(k+2)-pnew(k+1)
   th2=thet(k-1)
   th3=thet(k  )
   th4=thet(k+1)
! --- look for small dp2,dp4 in combination with large dp3
   if (dp3.gt.dp2 .and. dp3.gt.dp4) then
    goal=.5*(dp3+max(dp2,dp4))		!  desired thknss of lyr 3
    if (dp2.lt.goal .and. dp4.lt.goal) then
! --- thin-thick-thin combination found -> deflate lyr 3
     dis3=.5*(goal-dp2) * slak
     dis4=.5*(goal-dp4) * slak
     if (th3.gt.th2 .and. th4.gt.th3) then
! --- theta conservation requires  dis3*(th3-th2)+dis4*(th4-th3)=0
      ratio=(th4-th3)/(th3-th2)
      if (ratio.gt.1.) then
       dis4=dis3/ratio
      else
       dis3=dis4*ratio
      end if
     end if
! --- ready to shrink middle layer
     pnew(k  )=pnew(k  )+dis3
     pnew(k+1)=pnew(k+1)-dis4

     if (vrbos) then
      write (6,'(a,i3,a,-4p,3f9.3/a,3f9.3)') 'k=',k,' thknss triple',	&
       dp2,dp3,dp4,'            becomes',pnew(k)-pnew(k-1),		&
        pnew(k+1)-pnew(k),pnew(k+2)-pnew(k+1)
     end if
     ncount=ncount+1
    end if
   end if
 2 continue

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- scenario 3: sequence of 3 thick-thin-thick layers

   do 3 k=3,kk-2
   if (pnew(k-1).lt.pnew(   1)+onemm .or.			&
       pnew(k+2).gt.pnew(kk+1)-onemm) go to 3
   dp2=pnew(k  )-pnew(k-1)
   dp3=pnew(k+1)-pnew(k  )
   dp4=pnew(k+2)-pnew(k+1)
   th2=thet(k-1)
   th3=thet(k  )
   th4=thet(k+1)
! --- look for large dp2,dp4 in combination with small dp3
   if (dp3.lt.dp2 .and. dp3.lt.dp4) then
! --- thick-thin-thick combination found -> inflate lyr 3
    goal=.5*(dp3+min(dp2,dp4))		!  desired thknss of lyr 3
    dis3=.5*(goal-dp3) * slak
    dis4=dis3
    if (th3.gt.th2 .and. th4.gt.th3) then
! --- theta conservation requires  dis3*(th3-th2)+dis4*(th4-th3)=0
     ratio=(th4-th3)/(th3-th2)
     if (ratio.gt.1.) then
      dis4=dis3/ratio
     else
      dis3=dis4*ratio
     end if
    end if
! --- ready to inflate middle layer
    pnew(k  )=pnew(k  )-dis3
    pnew(k+1)=pnew(k+1)+dis4

    if (vrbos) then
     write (6,'(a,i3,a,-4p,3f9.3/a,3f9.3)') 'k=',k,' thknss triple',	&
      dp2,dp3,dp4,'            becomes',pnew(k)-pnew(k-1),		&
       pnew(k+1)-pnew(k),pnew(k+2)-pnew(k+1)
    end if
    ncount=ncount+1
   end if
 3 continue

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- scenario 4: sequence of thick-(thin-thin)-thick layers

   do 4 k=4,kk-2
   if (pnew(k-2).lt.pnew(   1)+onemm .or.			&
       pnew(k+2).gt.pnew(kk+1)-onemm) go to 4
   dp1=pnew(k-1)-pnew(k-2)
   dp2=pnew(k  )-pnew(k-1)
   dp3=pnew(k+1)-pnew(k  )
   dp4=pnew(k+2)-pnew(k+1)
   th1=thet(k-2)
   th2=thet(k-1)
   th3=thet(k  )
   th4=thet(k+1)
! --- look for small dp1,dp3,dp5 in combination with large dp2,dp4
   if (dp1.gt.dp2+dp3 .and. dp4.gt.dp2+dp3) then
! --- thick-thin-thin-thick combination found -> inflate lyrs 2,3
    goal=.5*(dp2+dp3+min(dp1,dp4))	!  desired thknss of lyr 2+3
    dis2=(goal-dp2-dp3) * slak
    dis4=dis2
    if (th2.gt.th1 .and. th3.gt.th2 .and. th4.gt.th3) then
! --- theta conservation requires  dis2*(th2-th1)+dis4*(th3-th4)=0
     ratio=(th4-th3)/(th2-th1)
     if (ratio.gt.1.) then
      dis4=dis2/ratio
     else
      dis2=dis4*ratio
     end if
    end if
! --- ready to expand middle layers
    pnew(k-1)=pnew(k-1)-dis2
    pnew(k+1)=pnew(k+1)+dis4

    if (vrbos) then
     write (6,'(a,i3,a,-4p,4f9.3/a,4f9.3)') 'k=',k,		&
      ' thknss quadruplet',dp1,dp2,dp3,dp4,			&
       '                becomes',pnew(k-1)-pnew(k-2),		&
        pnew(k)-pnew(k-1),pnew(k+1)-pnew(k),pnew(k+2)-pnew(k+1)
    end if
    ncount=ncount+1
   end if
 4 continue

!  event = ncount>17			!  find interesting cases
   event=.false.

   do k=1,kk
    if (pnew(k+1).lt.pnew(k)-onemm) then
     event=.true.
     write (*,'(i5,i8,a)') i,k,                         	&
     ' error: nonmonotonic pressure on return from ocn_equilb'
    end if
   end do

   if (event) then
    write (6,99) i,'  ocn_equilb input profile:'
    do k2=1,kk,10
     write (6,100) (pres(k1),k1=k2,min(kk+1,k2+10) )
     write (6,101) (thet(k1),k1=k2,min(kk  ,k2+9 ) )
     write (6,100)
    end do
   end if
   if (event .or. vrbos) then
    write (6,99) i,'  ocn_equilb output profile:'
    write (6,*) ncount,' inflations'
    do k2=1,kk,10
     write (6,100) (pnew(k1),k1=k2,min(kk+1,k2+10) )
     write (6,102) (int(.01*(pnew(k1)-pres(k1))),k1=k2,		&
       min(kk+1,k2+10) )
     write (6,101) (thet(k1),k1=k2,min(kk  ,k2+9 ) )
     write (6,100)
    end do
   end if

   pres(:)=pnew(:)

   after=0.
   do k=1,kk
    after=after+(pres(k+1)-pres(k))*thet(k)
   end do
   if (abs(after-bfore).gt.1.e-15*abs(bfore))				&
    print '(i5,a,es14.6,2es11.3)',i,' (ocn_equilb) bad column intgl.',	&
     bfore,after-bfore,(after-bfore)/bfore
!  if (vrbos) stop					!  take out

   return
   end subroutine equilb

end module hycom_hybgn2
