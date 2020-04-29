module hycom_diapfl
contains
!!*********************************************************************
!!  diapfl
!!    single-column version of HYCOM's diapycnal mixing routine
!!
!!    R. Bleck                                   June 2011
!!*********************************************************************

   subroutine diapfl(nstep,delt,kk,dp,rho,tem,sal,vrbos,i)

   use module_constants,only: grvity
   use fimnamelist     ,only: diapyc,diapyn
   use hycom_constants ,only: thref,onem,tencm,onecm,onemm
   use hycom_sigetc

   implicit none
   integer,intent(IN)    :: kk,nstep,i
   real   ,intent(IN)    :: delt
   real*8 ,intent(INOUT) :: dp(kk)
   real*8 ,intent(INOUT) :: rho(kk)
   real*8 ,intent(INOUT) :: tem(kk)
   real*8 ,intent(INOUT) :: sal(kk)
   logical,intent(IN)    :: vrbos

   real*8 flxu(kk),flxl(kk),pcol(kk+1),pdot(kk+1),flngth(kk),clip(kk),	&
          ennsq(kk),xchcoef(kk),alfa,beta,qmin,qmax,amount,		&
          tflxu(0:kk+1),tflxl(0:kk+1),clipt,				&
          sflxu(0:kk+1),sflxl(0:kk+1),clips,				&
          dpold(kk),told(kk),sold(kk),scale,				&
          delp,q,totem,tosal,totra,tndcyt,tndcys,tndtra,		&
          boost,pbot1,p1,boomx,booz
   integer k,ka,kb,iter,kfrst,klast,incrm
   character text*18
   logical event,abort
   real,parameter :: small=1.e-22, acurcy=1.e-11
 
! --- boost = 1 at all depths except within the lowest 'booz' meters where
! --- it increases linearly to 'boomx'
   boost(pbot1,p1,boomx,booz)=max(1.,boomx+(1.-boomx)			&
     *(pbot1-p1)/(booz*onem))

   if (kk.lt.3) return				! need at least 3 layers

   if (vrbos) then
    print 103,nstep,i,							&
     '  entering diapfl:  temp    saln    dens    thkns    tracer',	&
      (k,tem(k),sal(k),rho(k),dp(k)/onem,k=1,kk)
103 format(2i8,a/(30x,i3,3f8.3,f9.3))
   end if

!! --- T/S conservation diagnostics (optional):
   totem=0.
   tosal=0.
   totra=0.
   scale=1.e-33
   do k=1,kk
     totem=totem+tem(k)*dp(k)
     tosal=tosal+sal(k)*dp(k)
!    if (dotrcr) then
!      totra=totra+tracer(k,1)*dp(k)
!      scale=scale+abs(tracer(k,1))
!    end if
   end do

   pcol(1)=0.
   do k=1,kk
    pcol(k+1)=pcol(k)+dp(k)
   end do

! --- find buoyancy frequency for each layer
 
   do 43 k=2,kk-1
! --- ennsq = buoy.freq.^2 / g^2 [s^2/m^2]
     ennsq(k)=max(0.,min(rho(k+1)-rho(k  ),				&
                         rho(k  )-rho(k-1))/max(dp(k),tencm))

! --- store (exch.coeff x buoy.freq.^2 / g x time step) in -flngth-
! --- (dimensions of flngth: length in pressure units)
! -----------------------------------------------------------------------
! --- use the following if exch.coeff. = diapyc / buoyancy frequency
!cc     flngth(k)=diapyn*sqrt(ennsq(k)) * delt * onem
! -----------------------------------------------------------------------
! --- use the following if exch.coeff. = diapyc
!cc     flngth(k)=diapyc*ennsq(k)*g * delt * onem
! -----------------------------------------------------------------------
! --- use the smaller of the above
!cc  flngth(k)=min(diapyn*sqrt(ennsq(k)),				&
!cc                diapyc*ennsq(k)*grvity) * delt * onem
 
! --- use the larger of the above
     flngth(k)=max(diapyn*sqrt(ennsq(k)),				&
                   diapyc*ennsq(k)*grvity) * delt * onem
! -----------------------------------------------------------------------
! --- boost diffusion near bottom
!        flngth(k)=flngth(k) * boost(pcol(kk+1),pcol(k),3.,300.)

! --- optional: save effective exchange coefficient in -xchcoef- array
        xchcoef(k)=flngth(k)/						&
         (max(1.e-11,ennsq(k)*grvity)*delt*onem)
43 continue
   ennsq( 1)=0.
   ennsq(kk)=0.
   flngth( 1)=0.
   flngth(kk)=0.
   xchcoef( 1)=0.
   xchcoef(kk)=0.

   if (vrbos) then
    print '(2i8,a/(i11,3es13.3))',nstep,i,				&
    '  ennsq     flngth[m]   diffu[m^2/s]',(k,ennsq(k),			&
      flngth(k)/onem,xchcoef(k),k=1,kk)
   end if

! --- find fluxes at the upper and lower interface of each layer
! --- (compute only the part common to T/S/mass fluxes)

   flxu( 1)=0.
   flxl( 1)=0.
   flxu(kk)=0.
   flxl(kk)=0.
   do 37 k=2,kk-1
   alfa=thref*dsigdt(tem(k),sal(k))
   beta=thref*dsigds(tem(k),sal(k))

   flxu(k)=flngth(k)/				&
     max(beta*(sal(k)-sal(k-1))			&
        +alfa*(tem(k)-tem(k-1)),small)
   flxl(k)=flngth(k)/				&
     max(beta*(sal(k+1)-sal(k))			&
        +alfa*(tem(k+1)-tem(k)),small)
37 continue

! --- determine mass flux -pdot- that will leave layer densities unchanged
! --- despite nonzero T/S fluxes.

   pdot(   1)=0.
   pdot(kk+1)=0.
   flxl(kk-1)=min(flxl(kk-1),dp(kk))
   clip(1)=1.
   do 40 k=2,kk
   clip(k)=1.
   if (k.gt.1) pdot(k)=flxu(k)-flxl(k-1)
40 continue

   if (vrbos) then
    print '(i8,a/(i6,f11.5,5es11.3))',i,				&
    '   pres      thkns       pdot     flngth       flxu       flxl',	&
     (k,pcol(k+1)/onem,dp(k)/onem,pdot(k+1)/onem,			&
      flngth(k)/onem,flxu(k)/onem,flxl(k)/onem,k=1,kk)
   end if

! --- now clip mass fluxes to prevent dp < 0.
! --- since pdot,flxu,flxl are inter-related, this is an iterative process

   do iter=1,5		!  go up and down the column repeatedly
   event=.false.

   if (mod(iter,2).eq.1) then
    kfrst=2
    klast=kk
    incrm=1
   else
    kfrst=kk
    klast=2
    incrm=-1
   end if
   do 42 k=kfrst,klast,incrm
   if (pdot(k).eq.0.) then
     event=.true.
     clip(k )=0.
   else if (pdot(k).gt.0.) then
     q=max(0.,.5*dp(k-1))/flxu(k  )
     if (q.lt.clip(k  )) then
       event=.true.
       clip(k  )=q
     end if

     if (vrbos) then
      print 101,i,k,							&
       '  pdot(k)   dp(k-1)  flxu(k)  flxl(k-1)  clip(k  )   iter',	&
       pdot(k)/onem,dp(k-1)/onem,flxu(k)/onem,				&
       flxl(k-1)/onem,clip(k  ),iter
 101  format (i8,i3,a/10x,5es10.2,i6)
     end if

   else if (pdot(k).lt.0.) then
     q=max(0.,.5*dp(k  ))/flxl(k-1)
     if (q.lt.clip(k-1)) then
       event=.true.
       clip(k-1)=q
     end if

     if (vrbos) then
      print 101,i,k,							&
       '  pdot(k)   dp(k  )  flxu(k)  flxl(k-1)  clip(k-1)   iter',	&
       pdot(k)/onem,dp(k  )/onem,flxu(k)/onem,				&
       flxl(k-1)/onem,clip(k-1),iter
     end if

   end if
42 continue

   do 44 k=1,kk

   if (vrbos .and. clip(k  ).lt.1. .and. k.lt.kk) then
    print '(2i8)',nstep,i
    print 102,k-1,flxu(k-1)/onem,flxl(k-1)/onem,			&
       flxu(k-1)*clip(k-1)/onem,flxl(k-1)*clip(k-1)/onem,clip(k-1)
    print 102,k  ,flxu(k  )/onem,flxl(k  )/onem,			&
       flxu(k  )*clip(k  )/onem,flxl(k  )*clip(k  )/onem,clip(k  )
    print 102,k+1,flxu(k+1)/onem,flxl(k+1)/onem,			&
      flxu(k+1)*clip(k+1)/onem,flxl(k+1)*clip(k+1)/onem,clip(k+1)
102 format (i3,4x,							&
     'flxu(k)   flxl(k)   flxu(k)*clip   flxl(k)*clip   clip(k)'/	&
      4x,2es10.2,es14.2,es15.2,es11.2)
   end if

   flxu(k)=flxu(k)*clip(k)
   flxl(k)=flxl(k)*clip(k)
   if (k.gt.1) pdot(k)=flxu(k)-flxl(k-1)
   clip(k)=1.
44 continue

   if (.not.event) exit
   end do			!  iter

! --- convert flxu,flxl into actual T/S fluxes

   tflxu( 1)=0.
   tflxl( 1)=0.
   sflxu( 1)=0.
   sflxl( 1)=0.
   tflxu(kk)=0.
   tflxl(kk)=0.
   sflxu(kk)=0.
   sflxl(kk)=0.
   tflxl( 0)=0.
   tflxu(kk+1)=0.
   sflxl(   0)=0.
   sflxu(kk+1)=0.
!  if (dotrcr) then
!    trflxu( 1,:)=0.
!    trflxl( 1,:)=0.
!    trflxu(kk,:)=0.
!    trflxl(kk,:)=0.
!    trflxl(   0,:)=0.
!    trflxu(kk+1,:)=0.
!  end if

   do 35 k=2,kk-1
   tflxu(k)=flxu(k)*tem(k-1)
   tflxl(k)=flxl(k)*tem(k+1)
   sflxu(k)=flxu(k)*sal(k-1)
   sflxl(k)=flxl(k)*sal(k+1)
!  if (dotrcr) then
!    trflxu(k,:)=flxu(k)*tracer(k-1,:)
!    trflxl(k,:)=flxl(k)*tracer(k+1,:)
!  end if
35 continue

   do 34 k=1,kk
   sold(k)=sal(k)
   told(k)=tem(k)
!  if (dotrcr) then
!    trold(k,:)=tracer(k,:)
!  end if
34 continue

   clipt=0.
   clips=0.
!  if (dotrcr) cliptr(:)=0.

! --- update interface pressure and layer temperature/salinity
   abort=.false.
   do 39 k=kk,1,-1
   ka=max( 1,k-1)
   kb=min(kk,k+1)
   dpold(k)=dp(k)
   pcol(k)=pcol(k)-pdot(k)
   dp(k)=pcol(k+1)-pcol(k)
   if (dp(k).lt.-onemm) then
     print '(a,i8,i4,es11.3)','diapfl: dp<0 at',i,k,dp(k)/onem
     abort=.true.
   end if

   delp=dp(k)
   if (delp.gt.0.) then
     amount=tem(k)*dpold(k)-(tflxu(k+1)-tflxu(k)+tflxl(k-1)-tflxl(k))
     q=amount
     qmax=max(told(k),told(ka),told(kb))
     qmin=min(told(k),told(ka),told(kb))
     amount=max(qmin*delp,min(amount,qmax*delp))
     clipt=clipt+(q-amount)
     tem(k)=amount/delp

     amount=sal(k)*dpold(k)-(sflxu(k+1)-sflxu(k)+sflxl(k-1)-sflxl(k))
     q=amount
     qmax=max(sold(k),sold(ka),sold(kb))
     qmin=min(sold(k),sold(ka),sold(kb))
     amount=max(qmin*delp,min(amount,qmax*delp))
     clips=clips+(q-amount)
     sal(k)=amount/delp

!    if (dotrcr) then
!      do nt=1,ntrcr
!        amount=tracer(k,nt)*dpold(k)			&
!            -(trflxu(k+1,nt)-trflxu(k,nt)		&
!             +trflxl(k-1,nt)-trflxl(k,nt))		&
!        q=amount
!        qmax=max(trold(k,nt),trold(ka,nt),trold(kb,nt))
!        qmin=min(trold(k,nt),trold(ka,nt),trold(kb,nt))
!        amount=max(qmin*delp,min(amount,qmax*delp))
!        cliptr(nt)=cliptr(nt)+(q-amount)
!        tracer(k,nt)=amount/delp
!      end do
!    end if				!  dotrcr
   end if				!  delp > 0
 39 continue
!
!  if (dotrcr) cliptr(:)=cliptr(:)/p(kk+1)
   clipt=clipt/pcol(kk+1)
   clips=clips/pcol(kk+1)

   do 41 k=1,kk

! --- restore 'clipped' T/S amount to column
   tem(k)=tem(k)+clipt
   sal(k)=sal(k)+clips
   rho(k)=sigocn(tem(k),sal(k))
!  if (dotrcr) tracer(k,:)=tracer(k,:)+cliptr(:)

!  diaflx(k)=diaflx(k)+(dp(k)-dpold(k))	! diapyc.flux
! --- make sure p is computed from dp, not the other way around (roundoff!)
41 pcol(k+1)=pcol(k)+dp(k)

! --- T/S conservation diagnostics (optional):
   tndcys=-tosal
   tndcyt=-totem
   tndtra=-totra
   do k=1,kk
     tndcys=tndcys+sal(k)*dp(k)
     tndcyt=tndcyt+tem(k)*dp(k)
!    if (dotrcr) tndtra=tndtra+tracer(k,1)*dp(k)
   end do
   if (abs(tndcyt).gt.acurcy*10.*pcol(kk+1)) then
    print 100,i,' (diapfl) bad temp.intgl.',totem,tndcyt,clipt
   end if
   if (abs(tndcys).gt.acurcy*35.*pcol(kk+1)) then
    print 100,i,' (diapfl) bad saln.intgl.',tosal,tndcys,clips
   end if
!  if (dotrcr) then
!   if (abs(tndtra)*kk.gt.acurcy*scale*pcol(kk+1)) then
!    print 100,i,' (diapfl) bad trcr.intgl.',totra,tndtra,cliptr(1)
!   end if
!  end if
100 format(i8,a,es16.8,es12.4,' clip=',es12.4)

   if (vrbos) then
    print 103,nstep,i,							&
     '  exiting  diapfl:  temp    saln    dens    thkns   tracer',	&
      (k,tem(k),sal(k),rho(k),dp(k)/onem,k=1,kk)
   end if

31 continue
   if (abort) stop

!  if (dotrcr) print '(a)','tracer diapycnal mixing done'
   return
   end subroutine diapfl
end module hycom_diapfl
