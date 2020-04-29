module  hycom_edgvar
use stencilprint
use stenedgprint
use edgmaxmin
contains
  subroutine uvedge(nstep,kdm,uvel,vvel,u_edg,v_edg)

! --- interpolate velocity components -uvel,vvel- to edges of icos cells

  use module_control  ,only: npp,nip
  use module_constants,only: prox,nprox,cs,sn,nedge,permedge,perm
  use hycom_constants ,only: wet,slip
  use fimnamelist     ,only: itest,diag_intvl

  implicit none
  integer,intent(IN) :: nstep,kdm
!SMS$DISTRIBUTE (dh,2) BEGIN
  real,intent(IN)  :: uvel (kdm,    nip),vvel (kdm,    nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,3) BEGIN
  real,intent(OUT) :: u_edg(kdm,npp,nip),v_edg(kdm,npp,nip)
!SMS$DISTRIBUTE END
  integer i,ipx,im1,ip1,edg,edgcount,k
  logical vrbos
! real,parameter :: divby6=1./6.	! use in trapezoidal rule
  real,parameter :: divby18=1./18.	! use in simpson's rule

  !  The following are u and v at neighboring icos points, NOT on edge
  real    :: u_xy1,u_xy2,u_xy3,u_xy4 ! u on the xy local grid (m/s), at prox pt
  real    :: v_xy1,v_xy2,v_xy3,v_xy4 ! v on the xy local grid (m/s), at prox pt

!SMS$PARALLEL (dh,i) BEGIN

!sms$compare_var(uvel, "uvedge.F90 - uvel1 ")
!sms$compare_var(vvel, "uvedge.F90 - vvel1 ")

!SMS$EXCHANGE(uvel,vvel)

!SMS$HALO_COMP(<1,1>) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,edg,ipx,im1,ip1,u_xy1,v_xy1,		&
!$OMP u_xy2,v_xy2,u_xy3,v_xy3,u_xy4,v_xy4)
  do i=1,nip				! horizontal loop
   vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
   u_edg(:,:,i)=0.
   v_edg(:,:,i)=0.
   do edgcount=1,nedge(i)		! loop through edges
    edg=permedge(edgcount,i)
  
! --- edge quantities are interpolated from 4 icos cells -- the cells on
! --- either side of the edge plus the 2 joint neighbors of this pair.
  
    ipx=prox(edg,i)
    im1=mod(edg-2+nprox(i),nprox(i))+1
    ip1=mod(edg           ,nprox(i))+1
    im1=prox(im1,i)
    ip1=prox(ip1,i)
  
    if (vrbos) then
     print 100,perm(i),edg,perm(i),wet(i),perm(ipx),wet(ipx),perm(im1),	&
       wet(im1),perm(ip1),wet(ip1)
 100 format (i8,' edg=',i1,' intpol based on icos cells',		&
       2(i8,i2),' (dbl wgt),',/41x,2(i8,i2),' (sngl wgt)')
    end if
  
!DIR$ vector always
    do k=1,kdm				! loop over layers
  
  !  Transform u,v at neighboring icos pt to local coord.system.
  !  cs and sn are coordinate transformation constants.
  !  u_xy,v_xy are values of u and v rotated into local system.
  !  (Unrolled to allow vectorization of the k loop)
  
     u_xy1= cs(1,edg,i)*uvel(k,i)+sn(1,edg,i)*vvel(k,i)
     v_xy1=-sn(1,edg,i)*uvel(k,i)+cs(1,edg,i)*vvel(k,i)
     if (wet(ipx) > 0 ) then
      u_xy2= cs(2,edg,i)*uvel(k,ipx)+sn(2,edg,i)*vvel(k,ipx)
      v_xy2=-sn(2,edg,i)*uvel(k,ipx)+cs(2,edg,i)*vvel(k,ipx)
      if (wet(i) == 0 ) then		! ipx is ocean point but i is not
       u_xy1=u_xy2*slip
       v_xy1=v_xy2*slip
      end if
     else				! ipx is land point
      u_xy2=u_xy1*slip
      v_xy2=v_xy1*slip
     end if
     if (wet(im1) > 0 ) then
      u_xy3= cs(3,edg,i)*uvel(k,im1)+sn(3,edg,i)*vvel(k,im1)
      v_xy3=-sn(3,edg,i)*uvel(k,im1)+cs(3,edg,i)*vvel(k,im1)
     else				! im1 is land point
!     u_xy3=.5*(u_xy1+u_xy2)*slip
!     v_xy3=.5*(v_xy1+v_xy2)*slip
      u_xy3=0.
      v_xy3=0.
     end if
     if (wet(ip1) > 0 ) then
      u_xy4= cs(4,edg,i)*uvel(k,ip1)+sn(4,edg,i)*vvel(k,ip1)
      v_xy4=-sn(4,edg,i)*uvel(k,ip1)+cs(4,edg,i)*vvel(k,ip1)
     else				! ip1 is land point
!     u_xy4=.5*(u_xy1+u_xy2)*slip
!     v_xy4=.5*(v_xy1+v_xy2)*slip
      u_xy4=0.
      v_xy4=0.
     end if
  
! --- interpolate rotated velocity components to edges
  
     u_edg(k,edg,i)=(8.*(u_xy1+u_xy2)+u_xy3+u_xy4)*divby18
     v_edg(k,edg,i)=(8.*(v_xy1+v_xy2)+v_xy3+v_xy4)*divby18
     if (vrbos) then
      print '(a,2i2,a,5es11.2/10x,a,5es11.2)','k,edg=',k,edg,		&
       '  u_edg:',u_edg(k,edg,i),u_xy1,u_xy2,u_xy3,u_xy4,		&
       '  v_edg:',v_edg(k,edg,i),v_xy1,v_xy2,v_xy3,v_xy4
     end if
    end do				! vertical loop
   end do				! loop over edges
  end do				! horizontal loop
!$OMP END PARALLEL DO
!SMS$HALO_COMP END

!sms$compare_var(u_edg, "uvedge.F90 - u_edg2 ")
!sms$compare_var(v_edg, "uvedge.F90 - v_edg2 ")
!SMS$PARALLEL END

  return
  end subroutine uvedge


!*********************************************************************
!      bclin_edgvar
!	Interpolates  data to cell edges in the local stereographic grid
!	J.  Lee                    Sep 2005
!	A.  E. MacDonald           Nov 2005  fim conversion
!       R. Bleck                   Oct 2009  adapted for ocean model use
!*********************************************************************

  subroutine bclin_edgvar (nstep,leapm,leapn,				&
   uclin,vclin,utrop,vtrop,dp,pres,geop,montg,spcifv,um_edg,vm_edg,	&
   un_edg,vn_edg,dp_edg,pr_edg,geop_edg,mont_edg,uvsq_edg,spcv_edg)
  
  use module_control  ,only: npp,nip
  use module_constants,only: prox,nprox,grvity
  use hycom_constants ,only: wet,onem,tenm,odepth,thref,land_spval
  use fimnamelist     ,only: kdm,itest,diag_intvl,janjic_ocn

  implicit none
  
! Type and dimension external variables:
  
  integer,intent(IN)    :: nstep		! model time step
  integer,intent(IN)    :: leapm,leapn		! leapfrog time slots
!SMS$DISTRIBUTE (dh,2) BEGIN
  real,intent(IN)    :: utrop (3  ,nip)		! eastward ba'trop velocity
  real,intent(IN)    :: vtrop (3  ,nip)		! northward ba'trop velocity
  real,intent(IN)    :: uclin (kdm,nip,2)	! eastward ba'clin velocity
  real,intent(IN)    :: vclin (kdm,nip,2)	! northward ba'clin velocity
  real,intent(IN)    :: dp    (kdm,nip,2)	! layer thickness
  real,intent(INOUT) :: pres(kdm+1,nip)   	! pressure on interfaces
  real,intent(IN)    :: geop(kdm+1,nip)		! geopotential anomaly
  real,intent(IN)    :: montg (kdm,nip)		! montgomery potential
  real,intent(IN)    :: spcifv(kdm,nip)		! specific volume (sigma units)
! Type and dimension of local variables:
  real utotm(kdm,nip),utotn(kdm,nip),			&
       vtotm(kdm,nip),vtotn(kdm,nip),geop_lyr(kdm,nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,3) BEGIN
  real,intent(INOUT) :: um_edg  (kdm,npp,nip)	! totl.u on edges, mid-time
  real,intent(INOUT) :: vm_edg  (kdm,npp,nip)	! totl.v on edges, mid-time
  real,intent(INOUT) :: un_edg  (kdm,npp,nip)	! totl.u on edges, old
  real,intent(INOUT) :: vn_edg  (kdm,npp,nip)	! totl.v on edges, old
  real,intent(INOUT) :: dp_edg  (kdm,npp,nip)	! lyr thknss on edges, mid-time
  real,intent(INOUT) :: pr_edg(kdm+1,npp,nip)	! pressure on edges, mid-time
  real,intent(INOUT) :: geop_edg(kdm,npp,nip)	! geopotential on edges, mid-lyr
  real,intent(INOUT) :: mont_edg(kdm,npp,nip)	! Montgomery potential on edges
  real,intent(INOUT) :: uvsq_edg(kdm,npp,nip)	! velocity-squared on edges
  real,intent(INOUT) :: spcv_edg(kdm,npp,nip)	! specific volume on edges
  real geopx(kdm,npp,nip)
!SMS$DISTRIBUTE END
  integer :: im1,ipx,ip1		! indices of 3 consec neighbors of i
  integer :: i				! icos index
  integer :: k,kx,ka,kb			! layer Index
  integer :: edg,em1,ep1		! icos edge index
!  The following are u and v at neighboring icos points, NOT on edge
! real    :: u_xy1,u_xy2,u_xy3,u_xy4 ! u on the xy local grid (m/s), at prox pt
! real    :: v_xy1,v_xy2,v_xy3,v_xy4 ! v on the xy local grid (m/s), at prox pt
  real    :: q,blok,var,thbot,thtop,pbot,ptop,thav,p_lyr(kdm),		&
             silhou,tapr,plyrk
  real    :: ngborm1,ngborm2,ngborm3,ngborv1,ngborv2,ngborv3
  logical :: vrbos
  character      :: text*25
! real,parameter :: divby6=1./6.	! use in trapezoidal rule
  real,parameter :: divby18=1./18.	! use in simpson's rule
  real           :: hfharm,a,b
  hfharm(a,b) = a*b/(a+b)	! harmonic average x 0.5

  write (text,'(a,i8)') '(ocn_edgvar) step',nstep

  if (janjic_ocn) then

! --- integrate hydrostatic eqn at neighboring points to common pressure level

!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE(pres,geop,spcifv)
!$OMP PARALLEL DO PRIVATE (vrbos,p_lyr,ipx,kx,thav,silhou,tapr,plyrk)
   do i=1,nip			! horizontal loop
    do k=1,kdm
     geop_lyr(k,i)=land_spval
     geopx(k,:,i)=land_spval
    end do
    if (wet(i) > 0) then
     vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
     do k=1,kdm
      p_lyr(k)=.5*(pres(k,i)+pres(k+1,i))
      geop_lyr(k,i)=.5*(geop(k,i)+geop(k+1,i))
     end do			! vertical loop

! --- find 'silhoutte depth' = min.depth at surrounding points
     silhou=1.e33
     do edg=1,nprox(i)		! loop through edges
      geopx(:,edg,i)=geop_lyr(:,i)
      ipx=prox(edg,i)
      if (wet(ipx) > 0) then
       silhou=min(silhou,pres(kdm+1,ipx))
      end if
     end do

! --- hydrostatic integration to -p_lyr- in neighboring columns
! --- (don't extend downward integration past -plyrk-)

     do edg=1,nprox(i)		! loop through edges
      ipx=prox(edg,i)
      if (wet(ipx) > 0) then

       do k=1,kdm		! vertical loop
        plyrk=min(p_lyr(k),silhou)
!       plyrk=    p_lyr(k)
        kx=1
        do while (kx.lt.kdm .and. pres(kx+1,ipx).lt.plyrk)
         kx=kx+1
        end do
        if (pres(kx+1,ipx).ge.plyrk) then
! --- option 1: keep theta constant within each layer
         thav=spcifv(kx,ipx)
! --- option 2: allow theta to vary linearly (use PLM limiters)
!        ka=max(  1,kx-1)
!        kb=min(kdm,kx+1)
!        var=.5*min(spcifv(kb,ipx)-spcifv(kx,ipx),		&
!                   spcifv(kx,ipx)-spcifv(ka,ipx))
!        thbot=spcifv(kx,ipx)-var
!        thtop=spcifv(kx,ipx)+var
!        pbot=pres(kx+1,ipx)
!        ptop=pres(kx  ,ipx)
!        thav=.5*(thbot+(thbot*(plyrk-ptop)			&
!                       +thtop*(pbot-plyrk))			&
!                             /(pbot -ptop))
         geopx(k,edg,i)=geop(kx,ipx)-thav*(plyrk-pres(kx,ipx))
        end if			! neighboring intfc deeper than plyrk
       end do			! vertical loop
      end if			! neighbor is ocean point
     end do			! loop over edges

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- nudge pressure force to zero at points below silhouette depth
     if (vrbos) print '(a,i8,f11.3)',				&
       '(ocn_edgvar) silhou at i=',perm(i),silhou/onem
     do k=1,kdm
      tapr=max(0.,min(1.,(p_lyr(k)-silhou)/(100.*onem)))
      do edg=1,nprox(i)
       ipx=prox(edg,i)
       geopx(k,edg,i)=(1.-tapr)*geopx(k,edg,i)+tapr*geop_lyr(k,i)
      end do
     end do
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     if (vrbos) then
      print '(2(a,i8),a/(i9,2(f15.2,f11.2),f13.2))',			&
       '(ocn_edgvar)  neighb column',perm(ipx),				&
       '    central column',perm(i),'    interpol geop',		&
       (k,pres(k,ipx)/onem,geop(k,ipx),p_lyr(k)/onem,geop_lyr(k,i),	&
       geopx(k,edg,i),k=1,kdm)
     end if
    end if			! ocean point
   end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  end if			! janjic = true
!..........................................................................
!  Sec. 1.  Interpolate Edge Variables u_edg,v_edg
!..........................................................................
  
!SMS$PARALLEL (dh,i) BEGIN
  
!sms$compare_var(uclin, "edgvar.F90 - uclin1 ")
!sms$compare_var(vclin, "edgvar.F90 - vclin1 ")
!sms$compare_var(montg, "edgvar.F90 - montg1 ")

!$OMP PARALLEL DO
  do i=1,nip
   utotm(:,i)=land_spval
   utotn(:,i)=land_spval
   vtotm(:,i)=land_spval
   vtotn(:,i)=land_spval
   if (wet(i) > 0 ) then
    do k=1,kdm
     utotm(k,i)=uclin(k,i,leapm)+utrop(leapm,i)
     utotn(k,i)=uclin(k,i,leapn)+utrop(leapn,i)
     vtotm(k,i)=vclin(k,i,leapm)+vtrop(leapm,i)
     vtotn(k,i)=vclin(k,i,leapn)+vtrop(leapn,i)
     pres(k+1,i)=pres(k,i)+dp(k,i,leapm)
    end do
   end if			! ocean point
  end do
!$OMP END PARALLEL DO
!SMS$EXCHANGE(dp,montg,pres,spcifv)
!SMS$PARALLEL END

  call uvedge(nstep,kdm,utotm,vtotm,um_edg,vm_edg)
  call uvedge(nstep,kdm,utotn,vtotn,un_edg,vn_edg)

!..............................................................................
!  Sec. 2.  Interpolate Edge Variables dp_edg,pr_edg,mont_edg,uvsq_edg,spcv_edg
!..............................................................................
  
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,ipx,im1,ip1,blok,ngborm1,ngborv1,	&
!$OMP ngborm2,ngborv2,ngborm3,ngborv3)
  do i=1,nip				! horizontal loop
   vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
   dp_edg  (:,:,i)=0.
   pr_edg  (:,:,i)=0.
   mont_edg(:,:,i)=0.
   uvsq_edg(:,:,i)=0.
   spcv_edg(:,:,i)=0.
   if (wet(i) > 0 ) then
    do edg=1,nprox(i)			! loop through edges
     ipx=prox(edg,i)
     if (wet(ipx) > 0 ) then
      pr_edg(kdm+1,edg,i)=min(pres(kdm+1,i),pres(kdm+1,ipx))
!     pr_edg(kdm+1,edg,i)=2.*hfharm(pres(kdm+1,i),pres(kdm+1,ipx))
     else
      ipx=i				! avoid stepping on land
     end if
     im1=mod(edg-2+nprox(i),nprox(i))+1
     ip1=mod(edg           ,nprox(i))+1
     im1=prox(im1,i)
     if (wet(im1) == 0 ) im1=i		! avoid stepping on land
     ip1=prox(ip1,i)
     if (wet(ip1) == 0 ) ip1=i		! avoid stepping on land

     do k=kdm,2,-1
      pr_edg(k,edg,i)=pr_edg(k+1,edg,i)					&
        -.5*min(tenm,dp(k,i,leapm),dp(k,ipx,leapm))
     end do

     do k=1,kdm

! --- interpolate layer thickness to edges

      if (k.lt.kdm) pr_edg(k+1,edg,i)=min(pr_edg(k+1,edg,i),            &
       pr_edg(k,edg,i)+(8.*(dp(k,i  ,leapm)+dp(k,ipx,leapm))		&
                          + dp(k,im1,leapm)+dp(k,ip1,leapm))*divby18)

      dp_edg(k,edg,i)=max(0.,pr_edg(k+1,edg,i)-pr_edg(k,edg,i))

      if (vrbos .and. k.eq.kdm/2) then
       print 101,k,edg,dp_edg(k,edg,i),perm(i),dp(k,i,leapm),		&
        perm(ipx),dp(k,ipx,leapm),perm(im1),dp(k,im1,leapm),		&
        perm(ip1),dp(k,ip1,leapm)
 101    format ('edgvar k =',i4,'  edge',i2,'  dp_edg =',f10.1,		&
         '  based on...'/4(i8,f10.1))
      end if
     end do			! loop over layers

! --- interpolate montg.pot. to edges
   
     ipx=prox(edg,i)
     im1=mod(edg-2+nprox(i),nprox(i))+1
     ip1=mod(edg           ,nprox(i))+1
     im1=prox(im1,i)
     ip1=prox(ip1,i)

     do k=1,kdm

! --- determine the extent (rel.units) to which a sidewall is present
      if (wet(ipx) > 0 ) then
       blok=blokf(dp(k,i,leapm),pres(k+1,i),pres(kdm+1,ipx),pres(k,ipx))
      else
       blok=1.
      end if
      ngborm1=montg (k,ipx)*(1.-blok)+montg (k,i)*blok
      ngborv1=spcifv(k,ipx)*(1.-blok)+spcifv(k,i)*blok

      if (wet(im1) > 0 ) then
       blok=blokf(dp(k,i,leapm),pres(k+1,i),pres(kdm+1,im1),pres(k,im1))
      else
       blok=1.
      end if
      ngborm2=montg (k,im1)*(1.-blok)+montg (k,i)*blok
      ngborv2=spcifv(k,im1)*(1.-blok)+spcifv(k,i)*blok

      if (wet(ip1) > 0 ) then
       blok=blokf(dp(k,i,leapm),pres(k+1,i),pres(kdm+1,ip1),pres(k,ip1))
      else
       blok=1.
      end if
      ngborm3=montg (k,ip1)*(1.-blok)+montg (k,i)*blok
      ngborv3=spcifv(k,ip1)*(1.-blok)+spcifv(k,i)*blok

! --- mont.pot. + kinetic energy
      mont_edg(k,edg,i)=					&
        (8.*(montg(k,i)+ngborm1)+ngborm2+ngborm3)*divby18
      uvsq_edg(k,edg,i)=					&
        .5*(um_edg(k,edg,i)**2+vm_edg(k,edg,i)**2)

! --- specific volume anomaly
      spcv_edg(k,edg,i)=					&
        (8.*(spcifv(k,i)+ngborv1)+ngborv2+ngborv3)*divby18
     end do			! loop over layers
    end do			! loop over edges
   end if			! ocean point
  end do			! horizontal loop
!$OMP END PARALLEL DO

  if (mod(nstep,diag_intvl).eq.0) then
   call stenedg(pres,pr_edg,kdm+1,text//'  intfc pressure, cell & edge')
   call stenedg(dp(:,:,leapm),dp_edg,kdm,text//'  layer thknss, cell & edge')
   call stenedg(montg,mont_edg,kdm,text//'  montg.pot., cell & edge')
   call stenedg(spcifv,spcv_edg,kdm,text//'  spcifv in cell and on edge')
   call stenedg(utotm,um_edg,kdm,text//'  utotm & um_edg (cm/s)')
   call stenedg(vtotm,vm_edg,kdm,text//'  vtotm & vm_edg (cm/s)')
  end if
  
!sms$compare_var(um_edg, "hycom_edgvar.F90 - um_edg2 ")
!sms$compare_var(pr_edg, "hycom_edgvar.F90 - pr_edg2 ")
!SMS$PARALLEL END
  
  if (mod(nstep,diag_intvl).eq.0) then
   call edgmxmn(um_edg,  kdm,'(ocn) um_edg  ')
   call edgmxmn(vm_edg,  kdm,'(ocn) vm_edg  ')
   call edgmxmn(dp_edg,  kdm,'(ocn) dp_edg  ')
   call edgmxmn(mont_edg,kdm,'(ocn) mont_edg')
   call edgmxmn(spcv_edg,kdm,'(ocn) spcv_edg')
  end if

  if (janjic_ocn) then

! --- interpolate geopotential to edges
  
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (ipx,em1,ep1,im1,ip1)
   do i=1,nip
    if (wet(i) > 0) then
     do edg = 1,nprox(i)	! loop through edges
!     vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
      ipx = prox(edg,i)
      em1 = mod(edg-2+nprox(i),nprox(i))+1
      ep1 = mod(edg           ,nprox(i))+1
      im1 = prox(em1,i)
      ip1 = prox(ep1,i)
      do k=1,kdm
       geop_edg(k,edg,i)=					&
         (8.*(geop_lyr(k,i)+geopx(k,edg,i))			&
            +geopx(k,em1,i)+geopx(k,ep1,i))*divby18
      end do			! vertical loop
     end do			! loop over edges
    end if			! ocean point
   end do			! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

   if (mod(nstep,diag_intvl).eq.0) then
    call stenedg(geop_lyr,geop_edg,kdm,text//'  geopot, cell & edge')
   end if
  end if			! janjic = true
  return
  end subroutine bclin_edgvar


  real function blokf(dp,plo,pbx,pupx)
  use hycom_constants,only: epsil

! --- blokf determines the extent to which a sidewall is present
! --- dp   = layer thickness
! --- plo  = lower interface pressure
! --- pbx  = neighboring bottom pressure
! --- pupx = pressure at neighboring upper interface

  real,intent(IN) :: dp,plo,pupx,pbx
! blokf=0.
! if (plo.gt.pbx)					&
!   blokf=1.-max(0.,min(1.,(plo-min(pbx,pupx+dp))/max(epsil,dp)))
  blokf=max(0.,min(1.,(plo-pbx)/max(epsil,dp)))
  return
  end function blokf
end module hycom_edgvar
