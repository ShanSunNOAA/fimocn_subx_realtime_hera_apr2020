module hycom_mxlayr
!*********************************************************************
!  mxlayr
!   Kraus-Turner mixed layer for isopycnic-coordinate ocean model
!   R. Bleck                   ~1980 (MICOM)
!   S. Sun, R. Bleck      March 2011
!
!   Reference: Bleck,Hanson,Hu,Kraus  J.Phys.Oceanogr. 1989
!*********************************************************************

  use findmaxmin1
  use findmaxmin2
  contains

  subroutine mxlayr(nstep,leapn,theta,dens,dp,pres,temp,saln,		&
                    uvel,vvel,surflx,pmevap,ustar)

  use module_control  ,only: nip
  use module_constants,only: cp,grvity,area,perm
  use fimnamelist     ,only: kdm,itest,diag_intvl
  use hycom_control   ,only: bclin_frq
  use hycom_constants ,only: wet,batrop,spcifh,thref,epsil,sigjmp,	&
			     dpthmn,onem,onecm,onemm,salmin,		&
                             c1,c2,c3,c4,c5,c6,c7	! state eqn. coeff.
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
  real,intent(INOUT) :: uvel (kdm,nip,2)	! u velocity
  real,intent(INOUT) :: vvel (kdm,nip,2)	! v velocity
  real,intent(INOUT) :: dens (kdm,nip  )	! density
  real,intent(INOUT) :: dp   (kdm,nip,2)	! layer thickness
  real,intent(INOUT) :: pres (kdm+1,nip)	! interface pressure
  real,intent(INOUT) :: temp (kdm,nip)		! temperature
  real,intent(INOUT) :: saln (kdm,nip)		! salinity
! real   :: dpchg(kdm,nip)
!SMS$DISTRIBUTE END

  real*8 :: ustar3,thkold,thknew,thermg,buoyfl,turgen,dtem,dsal,dsg,	&
          salflx,em,en,sum1,sum2,sum3,tdp,sdp,z,q,q1,q2,delp,pnew,  	&
          tem1d(kdm),sal1d(kdm),dpcol(kdm),prcol(kdm+1),puv(kdm+1),	&
          delt,obukhv,belo,heat,salt,colint,colins,cloutt,clouts
  real    :: sal4
  integer :: i,k,kdetr

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

    if (vrbos) then
    print 104,nstep,perm(i),						&
      'MXLAYR  in:  thkns    temp    saln    dens      u       v'
    print 105,(k,dpcol(k)/onem,temp(k,i),saln(k,i),		&
      dens(k,i),uvel(k,i,leapn),vvel(k,i,leapn),k=1,kdm)
104 format (i9,i7,3x,a)
105 format (26x,i3,f8.1,5f8.2)
    end if

   tem1d(:)=temp(:,i)
   sal1d(:)=saln(:,i)
   dpcol(:)=dp(:,i,leapn)
   salflx=-sal1d(1)*pmevap(i)/thref				! g/m^2/s
   heat=surflx(i)*delt*grvity/spcifh				! deg Pa
   salt=salflx   *delt*grvity					! (g/kg) Pa

! --- column integrals are computed for diagnostic purposes only
   colint=heat
   colins=salt
   do k=1,kdm
    colint=colint+tem1d(k)*dpcol(k)
    colins=colins+sal1d(k)*dpcol(k)
   end do

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- remove static instabilities (convective adjustment)
   convec=.false.
5  q1=dpcol(1)
   do k=2,kdm
    if (dens(k,i).gt.dens(1,i)) exit
    q2=dpcol(k)
    if (q2.gt.0.) then

     if (vrbos) then
      print '(i9,i7,a,2f8.4,a,i3,a)',nstep,perm(i),			&
        ' density inversion:',dens(1,i),dens(k,i),' -> entrain lyr',	&
        k,' into ML'
     end if

     tem1d(1)=(tem1d(1)*q1+tem1d(k)*q2)/(q1+q2)
     sal1d(1)=(sal1d(1)*q1+sal1d(k)*q2)/(q1+q2)
     uvel(1,i,leapn)=(uvel(1,i,leapn)*q1+uvel(k,i,leapn)*q2)/(q1+q2)
     vvel(1,i,leapn)=(vvel(1,i,leapn)*q1+vvel(k,i,leapn)*q2)/(q1+q2)
     dens(1,i)=sigocn(tem1d(1),sal1d(1))
     dpcol(k)=0.
     dpcol(1)=q1+q2
     convec=.true.
     if (k.lt.kdm) go to 5
    end if
   end do		! vertical loop
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! --- column integrals are computed for diagnostic purposes only
   colint=heat
   colins=salt
   do k=1,kdm
    colint=colint+tem1d(k)*dpcol(k)
    colins=colins+sal1d(k)*dpcol(k)
   end do

   ustar3 =ustar(i)**3
   thkold =dpcol(1)
   prcol(1)=pres(1,i)
   do k=1,kdm
    prcol(k+1)=prcol(k)+dpcol(k)
   end do

   if (vrbos .and. convec) then
    print 104,nstep,perm(i),						&
      'aft.convec:  thkns    temp    saln    dens      u       v'
    print 105,(k,dpcol(k)/onem,temp(k,i),saln(k,i),		&
      dens(k,i),uvel(k,i,leapn),vvel(k,i,leapn),k=1,kdm)
   end if

! --- surflx = downward energy flux (W/m^2)
! --- salflx = downward salt flux (g/m^2/s)
! --- buoyfl = downward buoyancy flux, w_prime_buoyancy_prime_bar (m^2/s^3)
! --- surface density increases (column is destabilized) if buoyfl < 0

   buoyfl=-grvity*thref**2*(dsigds(tem1d(1),sal1d(1))*salflx		&
                           +dsigdt(tem1d(1),sal1d(1))*surflx(i)/spcifh)

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
      pnew=(2.*turgen+belo*prcol(k)**2-sum2)/			&
            max(epsil,belo*prcol(k)   -sum1)

      if (vrbos) then
       print '(i9,i7,i3,a,3es11.3)',nstep,perm(i),k,		&
         ' p(k),p(k+1),pnew=',prcol(k)/onem,prcol(k+1)/onem,pnew/onem
      end if

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
    sdp=sal1d(1)*thkold
    tdp=tem1d(1)*thkold
    do k=2,kdm
     tdp=tdp+tem1d(k)*(min(thknew,prcol(k+1))-min(thknew,prcol(k)))
     sdp=sdp+sal1d(k)*(min(thknew,prcol(k+1))-min(thknew,prcol(k)))
     dpcol(k)=max(prcol(k+1),thknew)-max(prcol(k),thknew)
     if (dpcol(k).eq.0.) then
      sal4=sal1d(1)
      sal4=max(sal4,salmin(k))
      sal1d(k)=sal4
      tem1d(k)=tofsig(theta(k),sal4)
     end if
    end do

! --- redistribute temp. and salin. during entrainment

    if (vrbos) then
     print '(i9,i7,'' entrain'',es11.3,'' m of water'')',		&
       nstep,i,(thknew-thkold)/onem
    end if

    sal1d(1)=sdp/thknew
    tem1d(1)=tdp/thknew
    dpcol(1)=thknew
    dens(1,i)=sigocn(tem1d(1),sal1d(1))

    if (vrbos) then
     print 102,nstep,perm(i),1,' t,s,dp after entrainment',tem1d(1),	&
       sal1d(1),thknew/onem
    end if

   else					! turgen < 0

! --- (mixed layer  r e c e d e s)

! --- new mixed layer depth is given by the obukhov length which
! --- is found by setting the l.h.s. in the turgen expression to zero

    obukhv=2.*em*grvity*ustar3/max(epsil,-thref*thermg)

! --- don't allow mixed layer to get too shallow
    thknew=max(dpthmn*onem,obukhv)

    if (vrbos) then
     print '(i9,i7,a,es11.3,a,f7.1)',nstep,perm(i),' detrain up to',	&
       (thkold-thknew)/onem,' m of water. obukhov lgth=',obukhv/onem
    end if

! --- mixed layer water will be detrained into layer defined in -kdetr-.
    kdetr=2
    do k=2,kdm
     if (dens(1,i).lt.dens(k,i)-.0001) exit
      kdetr=min(k+1,kdm)
    end do

! --- unmix fossil mixed layer into an upper sublayer whose T/S properties
! --- match the actual mixed layer and a lower sublayer (of depth z)
! --- having the same density as layer -kdetr-. merge upper sublayer with
! --- the actual mixed layer and lower sublayer with layer -kdetr-.

    dtem=max(0.,heat/thknew)
    dsal=min(0.,salt/thknew)
    dsg=dtem*dsigdt(tem1d(1),sal1d(1))+dsal*dsigds(tem1d(1),sal1d(1))
    if (dsg.lt.0.) then
     k=kdetr
     q=max(0.,dens(k,i)-dens(1,i))/dsg			! note: q < 0
     z=(thkold-thknew)/(1.-q)
     if (z.gt.0) then

      if (vrbos) then
       print '(a,2es9.1,2f8.2,es9.1)','dtem,dsal,dns(1),dns(k),q=',	&
         dtem,dsal,dens(1,i),dens(k,i),q
       print 102,nstep,perm(i),k,' t,s,dens in upper sublyr',		&
         tem1d(1)+dtem  ,sal1d(1)+dsal  ,sigocn(tem1d(1)+dtem  ,sal1d(1)+dsal  )
       print 102,nstep,perm(i),k,' t,s,dens in lower sublyr',		&
         tem1d(1)+dtem*q,sal1d(1)+dsal*q,sigocn(tem1d(1)+dtem*q,sal1d(1)+dsal*q)
       print '(20x,a,es9.1)','net T/S gain after sublyr formation:',	&
         (thkold-thknew-z)+q*z
      end if

      tem1d(k)=(tem1d(k)*dpcol(k)+(tem1d(1)+dtem*q)*z)/(dpcol(k)+z)
      sal1d(k)=(sal1d(k)*dpcol(k)+(sal1d(1)+dsal*q)*z)/(dpcol(k)+z)
      tem1d(1)=tem1d(1)+dtem*(thkold-thknew-z)/(thkold-z)
      sal1d(1)=sal1d(1)+dsal*(thkold-thknew-z)/(thkold-z)
      dpcol(k)=dpcol(k)+z
      dpcol(1)=thkold  -z
     end if			! z   > 0
    end if			! dsg < 0
    dens(1,i)=sigocn(tem1d(1),sal1d(1))
!   dens(k,i)=sigocn(tem1d(k),sal1d(k))

    if (vrbos) then
     print '(i9,i7,3(a,f9.3))',nstep,perm(i),' reducing ML depth',	&
       thkold/onem,' by',z/onem,' m'
     print 102,nstep,perm(i),k,' t,s,dp after detrainment',tem1d(k),	&
       sal1d(k),dpcol(k)/onem
102  format (i9,i7,i3,a,2f8.3,f8.2)
    end if

   end if			! turgen > 0 or < 0

! --- apply surface forcing

! --- temperature change due to surflx = downward energy flux (W/m^2)
   tem1d(1)=tem1d(1)+heat/dpcol(1)

! --- salinity change due to pmevap = precipitation minus evaporation (m/s)
   sal1d(1)=sal1d(1)+salt/dpcol(1)
   dens(1,i)=sigocn(tem1d(1),sal1d(1))

   temp(:,i)=tem1d(:)
   saln(:,i)=sal1d(:)
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
    cloutt=cloutt+tem1d(k)*dpcol(k)
    clouts=clouts+sal1d(k)*dpcol(k)
   end do
   if (abs(colint-cloutt).gt.acurcy*10.*prcol(kdm+1)) then
    print 103,perm(i),'mxlayr - T column intgl.error',			&
          colint,cloutt,(cloutt-colint)/(10.*prcol(kdm+1))
   end if
   if (abs(colins-clouts).gt.acurcy*35.*prcol(kdm+1)) then
    print 103,perm(i),'mxlayr - S column intgl.error',			&
          colins,clouts,(clouts-colins)/(35.*prcol(kdm+1))
   end if
103 format (i7,3x,a,2es14.6,es9.1)

   if (vrbos) then
    print 104,nstep,perm(i),						&
      'MXLAYR out:  thkns    temp    saln    dens      u       v'
    print 105,(k,dpcol(k)/onem,temp(k,i),saln(k,i),			&
      dens(k,i),uvel(k,i,leapn),vvel(k,i,leapn),k=1,kdm)
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
  end subroutine mxlayr
end module hycom_mxlayr
