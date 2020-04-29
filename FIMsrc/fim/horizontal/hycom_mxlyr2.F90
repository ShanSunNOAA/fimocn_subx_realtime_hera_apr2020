module hycom_mxlyr2
!*********************************************************************
!  mxlyr2
!   Kraus-Turner mixed layer for isopycnic-coordinate ocean model
!
!   R. Bleck					~1980 (MICOM)
!   S. Sun, R. Bleck      			March 2011
!   R. Bleck (addition of fossil layer)		May 2011
!
!   Reference: Bleck,Hanson,Hu,Kraus  J.Phys.Oceanogr. 1989
!*********************************************************************

  use findmaxmin1
  use findmaxmin2
  contains

  subroutine mxlyr2(nstep,leapn,theta,dens,dp,pres,temp,saln,		&
                    uvel,vvel,surflx,pmevap,ustar)

  use module_control  ,only: nip
  use module_constants,only: cp,grvity,area,corio,perm
  use fimnamelist     ,only: kdm,itest,diag_intvl
  use hycom_control   ,only: numtr,bclin_frq
  use hycom_constants ,only: wet,batrop,spcifh,thref,epsil,sigjmp,	&
                             salmin,dpthmn,onem,onecm,onemm,		&
                             c1,c2,c3,c4,c5,c6,c7       ! state eqn. coeff.
  use hycom_sigetc

  implicit none
  integer,intent(IN) :: nstep			! model time step
  integer,intent(IN) :: leapn			! leapfrog time slot (n='new')
  real   ,intent(IN) :: theta(kdm)		! target densities
!SMS$DISTRIBUTE  (dh,1) BEGIN
  real,intent(IN)    :: surflx(nip)		! rad.+sensib+latnt heat flux
  real,intent(IN)    :: pmevap(nip)		! downwd freshwater flux (p-e)
  real,intent(IN)    :: ustar (nip)		! friction velocity
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE  (dh,2) BEGIN
  real,intent(INOUT) :: uvel  (kdm,nip,2)	! u velocity
  real,intent(INOUT) :: vvel  (kdm,nip,2)	! v velocity
  real,intent(INOUT) :: dens  (kdm,nip  )	! density
  real,intent(INOUT) :: dp    (kdm,nip,2)	! layer thickness
  real,intent(INOUT) :: pres  (kdm+1,nip)	! interface pressure
  real,intent(INOUT) :: temp  (kdm,nip)		! temperature
  real,intent(INOUT) :: saln  (kdm,nip)		! salinity
! real   :: dpchg(kdm,nip)
!SMS$DISTRIBUTE END

  real*8 :: ustar3,thkold,thknew,thermg,buoyfl,turgen,fosold,		&
          salflx,em,en,sum1,sum2,sum3,tdp,sdp,z,q,q1,q2,delp,pnew,  	&
          tcol(kdm),scol(kdm),dpcol(kdm),prcol(kdm+1),puv(kdm+1),	&
          delt,obukhv,belo,heat,salt,colint,colins,cloutt,clouts,	&
          dtuppr,dsuppr,dtlowr,dslowr,tmax,tmin,smax,smin,		&
          tlowr,slowr,t1,s1,t2,s2,t3,s3
  real    :: sal4
  integer :: i,k,ka,kdetr

  logical :: vrbos,convec
  character :: string*16
  real,parameter :: acurcy = 1.e-11	! integral conservation tolerance

! --- -----------------------------------
! --- mixed layer entrainment/detrainment
! --- -----------------------------------

  delt=batrop*bclin_frq			! routine-specific time step

!SMS$PARALLEL (dh,i) BEGIN
  do 8 i=1,nip				! horiz. loop
! dpchg(:,i)=dp(:,i,leapn)		! diagnostic use only
  if (wet(i) > 0 ) then			! ocean point
   vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0

   do k=1,kdm
    tcol(k)=temp(k,i)
    scol(k)=saln(k,i)
    dens(k,i)=sigocn(tcol(k),scol(k))
    dpcol(k)=dp(k,i,leapn)
   end do
   salflx=-scol(1)*pmevap(i)/thref				! g/m^2/s
   heat=surflx(i)*delt*grvity/spcifh				! deg Pa
   salt=salflx   *delt*grvity					! (g/kg) Pa

   if (vrbos) then
    print 104,nstep,perm(i),						&
      'MXLAYR  in:  thkns    temp    saln    dens      u       v'
    print 105,(k,dpcol(k)/onem,tcol(k),scol(k),dens(k,i),		&
      uvel(k,i,leapn),vvel(k,i,leapn),k=1,kdm)
104 format (i9,i7,3x,a)
105 format (26x,i3,f8.1,5f8.2)
   end if

! --- column integrals are computed for diagnostic purposes only
   colint=heat
   colins=salt
   do k=1,kdm
    colint=colint+tcol(k)*dpcol(k)
    colins=colins+scol(k)*dpcol(k)
   end do

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- remove static instabilities (convective adjustment)
   convec=.false.
   do ka=1,2			! look for layers lighter than layer -ka-
5   q1=dpcol(ka)
    do k=ka+1,kdm
     if (dens(k,i).gt.dens(ka,i)) exit
     q2=dpcol(k)
     if (q2.gt.0.) then

      if (vrbos) then
       print '(i9,i7,a,2f7.3,2(a,i3))',nstep,perm(i),			&
         ' density inversion:',dens(1,i),dens(k,i),' -> entrain lyr',	&
         k,' into lyr',ka
      end if

      tcol(ka)=(tcol(ka)*q1+tcol(k)*q2)/(q1+q2)
      scol(ka)=(scol(ka)*q1+scol(k)*q2)/(q1+q2)
      uvel(ka,i,leapn)=(uvel(ka,i,leapn)*q1+uvel(k,i,leapn)*q2)/(q1+q2)
      vvel(ka,i,leapn)=(vvel(ka,i,leapn)*q1+vvel(k,i,leapn)*q2)/(q1+q2)
      dens(ka,i)=sigocn(tcol(ka),scol(ka))
      if (k.eq.2) dens(k,i)=dens(ka,i)
      dpcol(ka)=q1+q2
      dpcol(k)=0.
      convec=.true.
      if (k.lt.kdm) go to 5
     end if
    end do		! vertical loop
   end do		! ka loop
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ustar3 =ustar(i)**3
   thkold =dpcol(1)
   prcol(1)=pres(1,i)
   do k=1,kdm
    prcol(k+1)=prcol(k)+dpcol(k)
   end do

   if (vrbos .and. convec) then
    print 104,nstep,perm(i),						&
      'aft.convec:  thkns    temp    saln    dens      u       v'
    print 105,(k,dpcol(k)/onem,tcol(k),scol(k),dens(k,i),		&
      uvel(k,i,leapn),vvel(k,i,leapn),k=1,kdm)
   end if

! --- surflx = downward energy flux (W/m^2)
! --- salflx = downward salt flux (g/m^2/s)
! --- buoyfl = downward buoyancy flux, w_prime_buoyancy_prime_bar (m^2/s^3)
! --- surface density increases (column is destabilized) if buoyfl < 0

   buoyfl=-grvity*thref**2*(dsigds(tcol(1),scol(1))*salflx		&
                           +dsigdt(tcol(1),scol(1))*surflx(i)/spcifh)

! --- determine turb.kin.energy generation due to wind stirring and
! --- diabatic forcing
   em=1.25				!   kraus-turner 'm' value
   en=0.4				!   kraus-turner 'n' value
!  em=0.8*exp(-prcol(2)/(50.*onem))	!   hadley ctr. choice
!  en=0.15				!   hadley ctr. choice
   thermg=-.5*((en+1.)*buoyfl+(en-1.)*abs(buoyfl))		! m^2/s^3
   turgen=delt*(2.*em*grvity*ustar3/thref+thkold*thermg)/thref**2

   if (vrbos) then
    print '(i9,i7,a/16x,6es10.2)',nstep,perm(i),			&
      '   surflx    pmevap     ustar    buoyfl    thermg    turgen',	&
       surflx(i),pmevap(i),ustar(i),buoyfl,thermg,turgen
   end if

   thknew=thkold
   if (thkold.lt.dpthmn*onem) turgen=max(0.,turgen)

   if (turgen.ge.0.) then

! --- (mixed layer  d e e p e n s)

! --- find pnew in case of mixed layer deepening (turgen > 0).
! --- entrain as many layers as needed to deplete -turgen-.

    sum1=dens(1,i)*thkold
    sum2=dens(1,i)*thkold**2

    do k=2,kdm
     if (dpcol(k).gt.0.) then
      belo=max(dens(1,i)+sigjmp,dens(k,i))
      pnew=(2.*turgen+belo*prcol(k)**2-sum2)/				&
            max(epsil,belo*prcol(k)   -sum1)

!      if (vrbos) then
!       print '(i9,i7,i3,a,3es11.3)',nstep,perm(i),k,			&
!         ' p(k),p(k+1),pnew=',prcol(k)/onem,prcol(k+1)/onem,pnew/onem
!      end if

! --- stop iterating for 'pnew' as soon as pnew < k-th interface pressure
      if (pnew.lt.prcol(k)) exit
      thknew=pnew
      sum1=sum1+belo*dpcol(k)
      sum2=sum2+belo*(prcol(k+1)**2-prcol(k)**2)
     end if
    end do

! --- keep mixed layer depth within proper range
    thknew=max(dpthmn*onem,min(thknew,prcol(kdm+1)))

! --- 'tdp,sdp' will be needed for temp./salin. mixing during entrainment
    sdp=scol(1)*thkold
    tdp=tcol(1)*thkold
    do k=2,kdm
     tdp=tdp+tcol(k)*(min(thknew,prcol(k+1))-min(thknew,prcol(k)))
     sdp=sdp+scol(k)*(min(thknew,prcol(k+1))-min(thknew,prcol(k)))
     dpcol(k)=max(prcol(k+1),thknew)-max(prcol(k),thknew)
     if (k.gt.2 .and. dpcol(k).eq.0.) then
      sal4=scol(1)
      sal4=max(sal4,salmin(k))
      scol(k)=sal4
      tcol(k)=tofsig(theta(k),sal4)
     end if
    end do

! --- redistribute temp. and salin. during entrainment

    if (vrbos) then
     print '(i9,i7,'' entrain'',es11.3,'' m of water'')',		&
       nstep,perm(i),(thknew-thkold)/onem
    end if

    scol(1)=sdp/thknew
    tcol(1)=tdp/thknew
    dens(1,i)=sigocn(tcol(1),scol(1))
    dpcol(1)=thknew
    if (dpcol(2).eq.0.) then
! --- if fossil ML is empty, assign ML values to it
     scol(2)=scol(1)
     tcol(2)=tcol(1)
     dens(2,i)=dens(1,i)
    end if

    if (vrbos) then
     print 102,nstep,perm(i),1,' t,s,dp after entrainment',tcol(1),	&
       scol(1),thknew/onem
102  format (i9,i7,i3,a,2f8.3,f8.2)
    end if

   else					! turgen < 0

! --- (mixed layer  r e c e d e s)

! --- new mixed layer depth is given by the obukhov length which
! --- is found by setting the l.h.s. in the turgen expression to zero

    obukhv=2.*em*grvity*ustar3/max(epsil,-thref*thermg)

! --- don't allow mixed layer to get too shallow
    thknew=max(dpthmn*onem,obukhv)

    if (vrbos) then
     print '(i9,i7,a,es11.3,a,f7.1)',nstep,perm(i),' detrain',		&
       (thkold-thknew)/onem,' m of water. obukhov lgth=',obukhv/onem
    end if

    tcol(2)=(tcol(2)*dpcol(2)+tcol(1)*(thkold-thknew))/			&
            (        dpcol(2)+        (thkold-thknew))
    scol(2)=(scol(2)*dpcol(2)+scol(1)*(thkold-thknew))/			&
            (        dpcol(2)+        (thkold-thknew))
    dens(2,i)=sigocn(tcol(2),scol(2))

    dpcol(1)=thknew
    dpcol(2)=dpcol(2)+(thkold-thknew)

   end if			! turgen > 0 or < 0

! --- apply surface forcing

! --- temperature change due to surflx = downward energy flux (W/m^2)
   tcol(1)=tcol(1)+heat/dpcol(1)

! --- salinity change due to pmevap = precipitation minus evaporation (m/s)
   scol(1)=scol(1)+salt/dpcol(1)

   dens(1,i)=sigocn(tcol(1),scol(1))

! --- now try to detrain fossil mixed layer (layer 2) into interior. to
! --- facilitate this, 'unmix' fossil mixed layer into an upper sublayer whose
! --- density matches the active mixed layer and a lower sublayer (of depth z)
! --- having a density close to that of layer -kdetr-. detrain lower sublayer.

   kdetr=3
   do k=3,kdm
    if (dens(2,i).lt.theta(k)-.0001) exit
     kdetr=k+1
   end do
   if (kdetr.gt.kdm) then
    print '(i9,i7,a,i4)',nstep,perm(i),					&
      ' (mxlyr2) fossil layer denser than bottom layer'
   end if

   if (vrbos) then
    print 104,nstep,perm(i),						&
      'aft.forcng:  thkns    temp    saln    dens      u       v'
    print 105,(k,dpcol(k)/onem,tcol(k),scol(k),dens(k,i),		&
      uvel(k,i,leapn),vvel(k,i,leapn),k=1,min(kdm,kdetr))
   end if

   if (kdetr.le.kdm .and. dpcol(2).gt.onem) then
    t1=tcol(1)
    s1=scol(1)
    t2=tcol(2)
    s2=scol(2)
! --- find representative water mass properties below fossil layer
    tdp =0.
    sdp =0.
    delp=0.
    do k=1,kdm
     q=min(prcol(3)+onem ,prcol(k+1))			&
      -max(prcol(3)-onemm,prcol(k  ))
     if (q.gt.0.) then
      tdp=tdp+tcol(k)*q
      sdp=sdp+scol(k)*q
      delp=delp+q
     end if
     if (prcol(k+1).ge.prcol(3)+onem) exit
    end do
    t3=tdp/delp
    s3=sdp/delp
    k=kdetr
! --- make upper sublayer slightly denser than mixed layer
    if (t2.lt.t1) t1=t1+max(t2-t1,-.20)
    if (s2.gt.s1) s1=s1+min(s2-s1, .02)

! --- establish bounding box for T/S unmixing in the fossil layer.
! --- this is to avoid generating unrealistic T,S combinations during unmixing
    tmin=min(t1,t2,t3)
    smin=min(s1,s2,s3)
    tmax=max(t1,t2,t3)
    smax=max(s1,s2,s3)

    if (vrbos) then
     print '(i9,i7,2(a,3f7.3))',nstep,perm(i),'  t1,2,3=',t1,t2,t3,	&
                                              '  s1,2,3=',s1,s2,s3
     print '(i9,i7,a,2f8.3/30x,a,2f8.3)',nstep,perm(i),			&
      ' bounding box: tmin,tmax=',tmin,tmax,' smin,smax=',smin,smax
    end if

    if (t2.lt.tmin+.001 .or. t2.gt.tmax-.001) then	! S unmixing only
     q=dsigds(t2,s2)						! > 0
     dslowr=min(smax-s2,(theta(k) -dens(2,i))/q)		! > 0
     dsuppr=max(smin-s2,(dens(1,i)-dens(2,i))/q)		! < 0
     z=0.
     if (dsuppr*dslowr.lt.0.) then
      z=dpcol(2)*dsuppr/(dsuppr-dslowr)
      if (z.gt.0.) then
       scol(2)=s2+dsuppr
       scol(k)=(scol(k)*dpcol(k)+(s2+dslowr)*z)/(dpcol(k)+z)
       tcol(k)=(tcol(k)*dpcol(k)+ t2        *z)/(dpcol(k)+z)
      end if
     end if

    if (vrbos) then
     print '(i9,i7,a,4f8.3)',nstep,perm(i),				&
       ' S unmixing: dslo,dsup,dp,z=',dslowr,dsuppr,dpcol(2)/onem,z/onem
    end if

    else	&
    if (s2.lt.smin+.001 .or. s2.gt.smax-.001) then	! T unmixing only
     q=dsigdt(t2,s2)						! < 0
     dtlowr=max(tmin-t2,(theta(k) -dens(2,i))/q)		! < 0
     dtuppr=min(tmax-t2,(dens(1,i)-dens(2,i))/q)		! > 0
     z=0.
     if (dtuppr*dtlowr.lt.0.) then
      z=dpcol(2)*dtuppr/(dtuppr-dtlowr)
      if (z.gt.0.) then
       tcol(2)=t2+dtuppr
       tcol(k)=(tcol(k)*dpcol(k)+(t2+dtlowr)*z)/(dpcol(k)+z)
       scol(k)=(scol(k)*dpcol(k)+ s2        *z)/(dpcol(k)+z)
      end if
     end if

    if (vrbos) then
     print '(i9,i7,a,4f8.3)',nstep,perm(i),				&
       ' T unmixing: dtlo,dtup,dp,z=',dtlowr,dtuppr,dpcol(2)/onem,z/onem
    end if

    else						! T and S unmixing 
     tlowr=t1+(t2-t1)*(s3-s1)/(s2-s1)
     slowr=s1+(s2-s1)*(t3-t1)/(t2-t1)
     if (tlowr.ge.tmin) then
      z=dpcol(2)*(s1-s2)/(s1-s3)
      slowr=s3
     else 	&
     if (slowr.le.smax) then
      z=dpcol(2)*(t1-t2)/(t1-t3)
      tlowr=t3
     else
      print '(a,i7/a,8f8.3)','(mxlyr2) unmix error at ipn=',perm(i),	&
        '  t1/2/3, s1/2/3, tlowr,slowr=',t1,t2,t3,s1,s2,s3,tlowr,slowr
      stop
     end if
     tcol(2)=t1
     scol(2)=s1
     tcol(k)=(tcol(k)*dpcol(k)+(tlowr)*z)/(dpcol(k)+z)
     scol(k)=(scol(k)*dpcol(k)+(slowr)*z)/(dpcol(k)+z)

     if (vrbos) then
      print '(i9,i7,a,4f8.3)',nstep,perm(i),				&
        ' T/S unmixing: dtlo,dslo,dp,z=',tlowr-t2,slowr-s2,		&
        dpcol(2)/onem,z/onem
     end if

    end if
   
    fosold=dpcol(2)
    dpcol(2)=dpcol(2)-z
    dpcol(k)=dpcol(k)+z

    dens(2,i)=sigocn(tcol(2),scol(2))
    dens(k,i)=sigocn(tcol(k),scol(k))

    if (vrbos) then
     print '(i9,i7,3(a,f9.3))',nstep,perm(i),				&
       ' reducing fossil ML depth',fosold/onem,'  by',z/onem,' m'
     print 102,nstep,perm(i),k,' t,s,dp after detrainment',tcol(k),	&
       scol(k),dpcol(k)/onem
    end if
   end if			! dpcol(2) > 0

   temp(:,i)=tcol(:)
   saln(:,i)=scol(:)
   dp(:,i,leapn)=dpcol(:)

! --- -----------------------
! --- momentum redistribution
! --- -----------------------

! --- store pre-entrainment/detrainment interface pressures in -puv-
! --- note that -dp- contains post-entrainment/detrainment layer thknss

   puv(:)=pres(:,i)
   do k=1,kdm
    pres(k+1,i)=pres(k,i)+dpcol(k)		! update -pres-
   end do

! --- redistribute momentum in the vertical -- entrainment case

   if (thknew.gt.thkold) then
    sum1=0.
    sum2=0.
    sum3=0.
    do k=1,kdm
     delp=max(0.,min(dpcol(k),puv(k+1))			&
                -min(dpcol(k),puv(k  )))
     sum1=sum1                +delp
     sum2=sum2+uvel(k,i,leapn)*delp
     sum3=sum3+vvel(k,i,leapn)*delp
    end do
    if (sum1.gt.0.) then
     uvel(1,i,leapn)=sum2/sum1
     vvel(1,i,leapn)=sum3/sum1
    end if

! --- redistribute momentum in the vertical -- detrainment case

   else			! thknew < thkold

    sum1=puv(2)
    do k=2,kdm
     puv(k)=puv(k-1)+dpcol(k-1)
     q=max(0.,min(1.,(sum1-puv(k))/max(dpcol(k),epsil)))
     uvel(k,i,leapn)=uvel(1,i,leapn)*q+uvel(k,i,leapn)*(1.-q)
     vvel(k,i,leapn)=vvel(1,i,leapn)*q+vvel(k,i,leapn)*(1.-q)
    end do
   end if			! entrainment/detrainment

! --- check conservation of column integrals
   cloutt=0.
   clouts=0.
   do k=1,kdm
    cloutt=cloutt+tcol(k)*dpcol(k)
    clouts=clouts+scol(k)*dpcol(k)
   end do
   if (abs(colint-cloutt).gt.acurcy*10.*prcol(kdm+1)) then
    print 103,perm(i),'(mxlyr2) T column intgl.error',			&
          colint,cloutt,(cloutt-colint)/(10.*prcol(kdm+1))
   end if
   if (abs(colins-clouts).gt.acurcy*35.*prcol(kdm+1)) then
    print 103,perm(i),'(mxlyr2) S column intgl.error',			&
          colins,clouts,(clouts-colins)/(35.*prcol(kdm+1))
   end if
103 format (i7,3x,a,2es14.6,es9.1)

   if (vrbos) then
    print 104,nstep,perm(i),						&
      'MXLAYR out:  thkns    temp    saln    dens      u       v'
    print 105,(k,dpcol(k)/onem,tcol(k),scol(k),dens(k,i),		&
      uvel(k,i,leapn),vvel(k,i,leapn),k=1,kdm)
   end if

!  dpchg(:,i)=dp(:,i,leapn)-dpchg(:,i)
  end if				! ocean point
8 continue				! horiz. loop
!SMS$PARALLEL END

  if (mod(nstep,diag_intvl).eq.0) then
   call findmxmn1(surflx,nip,'(mxlayr) surflx',wet)
   call findmxmn1(pmevap,nip,'(mxlayr) pmevap',wet)
   do k=1,kdm
    write (string,'(a,i2)') ' k=',k
    call findmxmn2(dens,kdm,nip,k,'(mxlayr) dens'//string,wet)
    call findmxmn2(temp,kdm,nip,k,'(mxlayr) temp'//string,wet)
    call findmxmn2(saln,kdm,nip,k,'(mxlayr) saln'//string,wet)
!   call findmxmn2(dpchg,kdm,nip,k,'(mxlayr) dpchg'//string,wet)
   end do
  end if

  print *,'ocn mixed layer updated,  time step',nstep
  return
  end subroutine mxlyr2
end module hycom_mxlyr2
