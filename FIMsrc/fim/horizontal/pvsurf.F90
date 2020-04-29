module mdul_pvsurf
contains
  subroutine pvsurf(its,relvor,delp,tracr,potvor,th_pvsrf,pr_pvsrf,	&
                    us_pvsrf,vs_pvsrf,pvsnap)
  use global_bounds,   only: ims, ime, ips, ipe, ihe
  use module_control,  only: nip, nvlp1, ntra, ntrb
  use fimnamelist,     only: nvl, PrintIpnDiag,	&
                       PVparam, TimingBarriers, ArchvTimeUnit
  use module_constants,only: grvity,corio,deg_lat,rd,cp,perm
  use module_variables,only: worka,us3d,vs3d
  use stencilprint
  use findmaxmin1
  use mdul_smoofld    ,only: smolight

! --- find pot.temp. on surface(s) of constant potential vorticity

  implicit none
  integer,intent(IN)  :: its			! model time step
  real   ,intent(IN)  :: relvor(nvl,ims:ime)	! relative vorticity
  real   ,intent(IN)  :: delp(nvl,ims:ime)	! layer thickness
  real   ,intent(IN)  :: tracr(nvl,ims:ime,ntra+ntrb)	! pot.temperature
  real   ,intent(OUT) :: potvor(nvl,ims:ime)	! pot.vorticity
  real   ,intent(OUT) :: th_pvsrf(ims:ime)	! pot.temp. on PV iso-surface
  real   ,intent(OUT) :: pr_pvsrf(ims:ime)	! pressure on PV iso-surface
  real   ,intent(OUT) :: us_pvsrf(ims:ime)	! u component on PV iso-surface
  real   ,intent(OUT) :: vs_pvsrf(ims:ime)	! v component on PV iso-surface
  real   ,intent(OUT),optional :: pvsnap(5,ims:ime)	! PV at SNAP levels

  integer ico,k,k0,k1,k2,n
  real    arg1,arg2,arg3,z,z1,z2,sumpos,sumneg,pvmin,prs,exner(nvlp1),q
  logical vrbos
  character string*8
  real,parameter :: small=1.e-5
  real, external :: its2time
  real, parameter:: thsnap(5)=(/1850., 1450., 840., 650., 340./)

  integer :: ret

#include <gptl.inc>

  print *,'(pvsurf) diagnosing PV and theta on PV surface',PVparam

! --- smooth rel.vorticity field

!sms$ignore begin
!$OMP PARALLEL DO schedule (static)
   do ico=ips,ipe
    do k=1,nvl
     worka(k,ico)=relvor(k,ico)
    end do
   end do
!$OMP END PARALLEL DO
!sms$ignore end

   call stencl(worka,nvl,1.e7,'relvort x 10^7 before smoothing')
   call smolight(worka,nvl,2,'rel.vorticity')
   call stencl(worka,nvl,1.e7,'relvort x 10^7 after smoothing')

!sms$ignore begin

! --- compute potential vorticity

!$OMP PARALLEL DO PRIVATE(k0,k1,k2,z,z1,z2,pvmin,sumpos,sumneg,		&
!$OMP             prs,exner,arg1,arg2,arg3,vrbos),			&
!$OMP             SCHEDULE (static)
   do ico=ips,ipe
    vrbos=ico.eq.PrintIpnDiag
    prs=0.
    do k=nvl,1,-1
     prs=prs+delp(k,ico)
     exner(k)=prs**(rd/cp)
    end do
    exner(nvlp1)=0.
    do k=1,nvl
     exner(k)=.5*(exner(k)+exner(k+1))
    end do
    do k=1,nvl
     arg1=0.
     arg2=0.
     if (k.eq.1) then
      arg3=max(tracr(  2,ico,1)-tracr(    1,ico,1),small)
     else if (k.eq.nvl) then
!     arg3=max(tracr(nvl,ico,1)-tracr(nvl-1,ico,1),small)
      arg3=max(tracr(nvl,ico,1)-tracr(nvl-1,ico,1),100.)	! experimental
     else
      arg1=max(tracr(k  ,ico,1)-tracr(k-1,ico,1),small)
      arg2=max(tracr(k+1,ico,1)-tracr(k  ,ico,1),small)
      arg3=2.*arg1*arg2/(arg1+arg2)                     ! harmonic average
     end if
     potvor(k,ico)=1.e6*grvity*(worka(k,ico)+corio(ico))*arg3		&
       /max(1.,delp(k,ico))
     if (corio(ico).lt.0.) potvor(k,ico)=-potvor(k,ico)	! sou.hem. sign change

!<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>
     if (present(pvsnap)) then
      if (k.gt.1) then
       do k1=1,5
        if (thsnap(k1).le.tracr(k  ,ico,1) .and.			&
            thsnap(k1).gt.tracr(k-1,ico,1)) then
         q=(tracr(k,ico,1)-thsnap(k1))/(tracr(k,ico,1)-tracr(k-1,ico,1))
         pvsnap(k1,ico)=q*potvor(k-1,ico)+(1.-q)*potvor(k,ico)
        end if
       end do
      end if
     end if
!<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>

#ifdef DEBUGPRINT
     if (vrbos)	then
      print '(i7,i3,2x,a,3f7.2,3f8.2)',perm(ico),k,			&
       'arg1/2/3,vor,dp,pv=',arg1,arg2,arg3,1.e5*worka(k,ico),		&
       .01*delp(k,ico),potvor(k,ico)
     end if
#endif
    end do		! k loop

! --- find PV surface 'PVparam'

#ifdef DEBUGPRINT
    if (vrbos) print '(i7,3x,a/(i13,2f12.3))',perm(ico),		&
      '(thsurf) theta      pot.vort',(k,tracr(k,ico,1),			&
      potvor(k,ico),k=1,nvl)
#endif

! --- find lowest and highest point in profile where pv = PVparam
! --- search only in layers above column minimum

    th_pvsrf(ico)=999.
    pr_pvsrf(ico)=999.
    us_pvsrf(ico)=999.
    vs_pvsrf(ico)=999.
    if (abs(deg_lat(ico)).gt.3.) then		! stay away from equator

! --- find minimum in column
     pvmin=1.e33
     k0=-1
     do k=2,nvl/2
      if (potvor(k,ico).lt.pvmin) then
       k0=k				! min PV value in column
       pvmin=potvor(k,ico)
      end if
     end do
     if (vrbos) print '(i7,a,f9.3,a,i3)',perm(ico),'  column PV min',	&
       pvmin,' at k=',k0

! --- find first crossover in column
     k1=-1
     do k=k0,nvl
      if (potvor(k-1,ico).lt.PVparam .and. potvor(k,ico).ge.PVparam) then
       k1=k
       exit
      end if
     end do
     if (k1.gt.0) then
      z=(potvor(k1-1,ico)-PVparam)/(potvor(k1-1,ico)-potvor(k1,ico))
      sumpos=(1.-z)*(potvor(k1  ,ico)-PVparam)
      z1=float(k1-1)+z
     end if

! --- find last crossover in column
     k2=-1
!    do k=nvl,2,-1
     do k=nvl*4/5,k0,-1		! bypass topmost layers where vort. may go crazy
      if (potvor(k,ico).gt.PVparam .and. potvor(k-1,ico).le.PVparam) then
       k2=k
       exit
      end if
     end do
     if (k2.gt.0) then
      z=(potvor(k2-1,ico)-PVparam)/(potvor(k2-1,ico)-potvor(k2,ico))
      sumneg=    z *(potvor(k2-1,ico)-PVparam)
      z2=float(k2-1)+z
     end if

     if (k1.lt.0 .or. k2.lt.0) then
      print '(i7,3x,a,f7.3,a,i3)',perm(ico),				&
       '(pvsurf) - no PV crossover point. use min PV=',pvmin,' at k=',k0
      th_pvsrf(ico)=tracr(k0,ico,1)
     else			! k1 > 0 .and. k2 > 0

! --- k1.ne.k2 indicates that the pv profile meanders around 'PVparam'. in that
! --- case, integrate both max(pv,0) and min(pv,0) over the interval (k1,k2).
! --- the rel. magnitude of the 2 integrals is used to construct a weighted
! --- average of k1,k2 which we take as the location of the pv=PVparam surface.

#ifdef DEBUGPRINT
      if (vrbos) print '(a,i8,3i4)','(pvsurf) i,k1,k2=',perm(ico),k1,k2
#endif
      if (k2.gt.k1) then	! multiple crossovers
       do k=k1+1,k2-1
        if (potvor(k,ico).ge.PVparam) then
         if (potvor(k-1,ico).ge.PVparam) then
          sumpos=sumpos+(potvor(k-1,ico)-PVparam)+(potvor(k,ico)-PVparam)
         else			! pv(k-1) < PVparam
          z=(potvor(k-1,ico)-PVparam)/(potvor(k-1,ico)-potvor(k,ico))
          sumpos=sumpos+(1.-z)*(potvor(k  ,ico)-PVparam)
          sumneg=sumneg+    z *(potvor(k-1,ico)-PVparam)
         end if
        else			! pv(k) < PVparam
         if (potvor(k-1,ico).le.PVparam) then
          sumneg=sumneg+(potvor(k-1,ico)-PVparam)+(potvor(k,ico)-PVparam)
         else			! pv(k-1) > PVparam
          z=(potvor(k-1,ico)-PVparam)/(potvor(k-1,ico)-potvor(k,ico))
          sumpos=sumpos+    z *(potvor(k-1,ico)-PVparam)
          sumneg=sumneg+(1.-z)*(potvor(k  ,ico)-PVparam)
         end if
        end if
       end do

       if (sumpos.ne.sumneg) then
        z=sumpos/(sumpos-sumneg)
        z=z*z1+(1.-z)*z2
       else
        z=.5*(z1+z2)		! arbitrary
       end if
#ifdef DEBUGPRINT
       if (vrbos) then
        print '(i7,3x,2(a,i3))',perm(ico),				&
          '(pvsurf) multiple crossovers between k=',k1,' and',k2
        print '(a,2es9.1,3f7.2)','sumpos,sumneg,z1,z2,final z=',	&
          sumpos,sumneg,z1,z2,z
       end if
#endif
       k=z+1.
       z=z-float(k-1)
      else			! k1 = k2, indicating single crossover
       k=k1
      end if
      th_pvsrf(ico)=(1.-z)*tracr(k-1,ico,1)+z*tracr(k,ico,1)
      pr_pvsrf(ico)=((1.-z)*exner(k-1)+z*exner(k))**(cp/rd)
      us_pvsrf(ico)=(1.-z)*us3d(k-1,ico)+z*us3d(k,ico)
      vs_pvsrf(ico)=(1.-z)*vs3d(k-1,ico)+z*vs3d(k,ico)
     end if			! k1,k2 > 0
    end if			! away from equator

#ifdef DEBUGPRINT
    if (vrbos) print '(i7,3x,a,i3,f6.3,4f7.2)',perm(ico),		&
      '(pvsurf) k,z,th/pr/us/vs_pvsrf=',k,z,th_pvsrf(ico),		&
      .01*pr_pvsrf(ico),us_pvsrf(ico),vs_pvsrf(ico)
#endif
   end do
!OMP END PARALLEL DO
!sms$ignore end

  write (string,'(f6.1,a2)') its2time(its),ArchvTimeUnit
#ifdef DEBUGPRINT
  call findmxmn1(th_pvsrf,nip,'th_pvsrf'//string)
#endif
  call stencl(th_pvsrf,1,1.,'th_pvsrf '//string)
  print *,'(pvsurf) ...PV diagnostics completed'
  return
  end subroutine pvsurf

!*********************************************************************
!     smoothr
!       Smoothing routine for data on icos grid
!	R.Bleck				July 2014
!*********************************************************************

! --- Method: Each icos point forms the center of a cluster of 6 (occasionally
! --- 5) triangles formed by triplets of icos points.  Averaging the 3 vertex
! --- values of each triangle yields an approximate areal average of 'field'
! --- in that triangle. This routine replaces the 'field' value at each grid
! --- point by the average of the 6 (or 5) surrounding triangular averages.

  subroutine smoothr(field,nvl,npass)
  use module_constants,only: prox,nprox
  use module_control  ,only: nip
  use findmaxmin2
  integer, intent(IN) :: nvl		! number of vertical levels
  integer, intent(IN) :: npass		! number of smoothing passes
!SMS$DISTRIBUTE (dh,2) BEGIN
  real,intent(INOUT) :: field(nvl,nip)	! field to be smoothed
  real    :: work(nvl,nip)
!SMS$DISTRIBUTE END
  real    :: sum(nvl)
  integer :: ipn,edg,iter
  character :: string*12
  do iter=1,npass
   work(:,:)=field(:,:)

#ifdef DEBUGPRINT
   do k=1,nvl
    write (string,'(a,i2)') 'smoothr k=',k
    call findmxmn2(work,nvl,nip,k,string)
   end do
#endif

!SMS$EXCHANGE(work)
!SMS$PARALLEL (dh,ipn) BEGIN
   do ipn=1,nip
    do k=1,nvl
     sum(k)=work(k,ipn)*.5*nprox(ipn)
     do edg=1,nprox(ipn)
      sum(k)=sum(k)+work(k,prox(edg,ipn))
     end do
     field(k,ipn)=sum(k)/(1.5*nprox(ipn))
    end do
   end do
!SMS$PARALLEL END

  end do		! iter
  return
  end subroutine smoothr
end module mdul_pvsurf
