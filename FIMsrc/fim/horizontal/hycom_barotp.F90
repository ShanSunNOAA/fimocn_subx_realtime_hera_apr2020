module hycom_barotp
use stencilprint
use stenedgprint
use findmaxmin2
contains

!*********************************************************************
!     barotp
!       Solves barotropic continuity and momentum equations
!	S.Sun, R.Bleck  	April 2010
!
!*********************************************************************

   subroutine barotp(nstep,leapn,utrop,vtrop,ptrop,btropfx,ubforc,vbforc)

use module_control  ,only: npp,nip
use module_constants,only: nprox,prox,proxs,rarea,sidevec_c,sidevec_e,	&
                           nedge,permedge,corio,deg_lat,cs,sn,sideln,perm
use hycom_constants ,only: wet,thref,dpth_edg,onem,batrop,slip
use hycom_control   ,only: bclin_frq
use fimnamelist     ,only: itest,diag_intvl,thkdff,veldff
use hycom_dffusn    ,only: dffusn_lev
use hycom_edgvar    ,only: uvedge
!use ersatzpipe

   implicit none
! External variables
   integer,intent (IN)    :: nstep			! model time step
   integer,intent (IN)    :: leapn			! leap frog time slot
!SMS$DISTRIBUTE (dh,1) BEGIN
   real   ,intent (INOUT) :: ubforc     (nip)	! eastw'd vel.forcing (m/s)
   real   ,intent (INOUT) :: vbforc     (nip)	! northw'd vel.forcing (m/s)
! Local variables
   real :: udissp(nip), vdissp(nip), pgfx(nip), pgfy(nip)
   real :: fld1d(nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
   real   ,intent (INOUT) :: utrop    (3,nip)	! eastw'd vel,m/sec
   real   ,intent (INOUT) :: vtrop    (3,nip)	! northw'd vel,m/sec
   real   ,intent (INOUT) :: ptrop    (3,nip)	! barotrop.prs
   real   ,intent (INOUT) :: btropfx(npp,nip)	! mass flux forcing (N/sec)
! Local variables
   real :: utrop_edg(npp,nip)	! utrop on edges, m/sec
   real :: vtrop_edg(npp,nip)	! vtrop on edges, m/sec
   real :: ptrop_edg(npp,nip)	! ptrop on edges, Pa
   real :: flx(npp,nip),uaux(1,nip),vaux(1,nip)
!SMS$DISTRIBUTE END
   integer :: i			! Index for icos point number
   integer :: edg,edgcount	! Index for icos edge number
   integer :: ipx		! neighbor across joint edge
   integer :: edgx		! joint edge index as seen by neighbor
   real    :: ungbor, vngbor, dfflxu, dfflxv, flxdv, rlstep, outgo,	&
              flxerr, valu1, valu2
   integer :: im1,ip1,lll,ml,nl,mn,nx
   logical :: vthenu,vrbos
   real,parameter :: tsmoo=0.25		! time smoothing weight
!  real,parameter :: tsmoo=0.125	! time smoothing weight
   real,parameter :: divby6=1./6.	! use in trapezoidal rule
   real,parameter :: divby18=1./18.	! use in simpson's rule

! --- explicit time integration of barotropic flow (forward-backward scheme)
! --- in order to combine forward-backward scheme with leapfrog treatment of
! --- coriolis term, v-eqn must be solved before u-eqn every other time step.
! --- bclin_frq must be even for final results to be returned in time slot leapn

   ml=leapn			! ml=slot for old time level
   nl=3				! nl=slot for new time level
   vthenu=.false.

   if (itest.gt.0 .and. mod(nstep,diag_intvl).eq.0) then
    fld1d(:)=ptrop(ml,:)
    call stenedg(fld1d,btropfx,1,'(ocn barotp) ptrop, btropfx')
    call stencl(ubforc,1,1000.,'(ocn barotp) ubforc (mm/s)')
    call stencl(vbforc,1,1000.,'(ocn barotp) vbforc (mm/s)')
   end if

!SMS$PARALLEL (dh,i) BEGIN
!SMS$HALO_COMP(<1,1>) BEGIN
   do i=1,nip
    flx   (:,i)=0.
    udissp  (i)=0.
    vdissp  (i)=0.
    pgfx    (i)=0.
    pgfy    (i)=0.
   end do
!SMS$HALO_COMP END

!sms$compare_var(btropfx,"barotp.F90 - btropfx1 ")
!sms$compare_var(ubforc, "barotp.F90 - ubforc1 ")
!sms$compare_var(vbforc, "barotp.F90 - vbforc1 ")
!sms$compare_var(utrop,  "barotp.F90 - utrop1 ")
!sms$compare_var(vtrop,  "barotp.F90 - vtrop1 ")

!SMS$PARALLEL END

   rlstep=1./float(bclin_frq)
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! --- test inherent stability of barotropic solver by turning off
! --- external forcing and running till the cows come home

!  if (nstep.eq.75) then
!   btropfx=0.
!   ubforc=0.
!   vbforc=0.
!   bclin_frq=1000
!   diag_intvl=1
!  end if
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

   do 80 lll=1,bclin_frq
!  print '(2(a,i7))','nstep =',nstep,'    starting barotropic time step',lll

! --- interpolate old -utrop,vtrop- to cell edges

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
   do i=1,nip
    uaux(1,i)=utrop(ml,i)
    vaux(1,i)=vtrop(ml,i)
   end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

   call uvedge(nstep,1,uaux,vaux,utrop_edg,vtrop_edg)

!sms$compare_var(utrop_edg, "barotp.F90 - utrop_edg2 ")
!sms$compare_var(vtrop_edg, "barotp.F90 - vtrop_edg2 ")

! --- -------------------
! --- continuity equation
! --- -------------------

!  print *,'solving -p- eqn... ml,nl    =',ml,nl

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,flxdv,edg,ipx,edgx,valu1,valu2,outgo)
   do i=1,nip			! horizontal loop
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    if (wet(i) > 0 ) then
! --- sum up edge fluxes to get tendency term
     flxdv=0.
     do edg=1,nprox(i)		! loop through edges
      ipx=prox(edg,i)
      if (wet(ipx) > 0 ) then
       edgx=proxs(edg,i)

! --- outgo=outging velocity component on edge times edge length
        valu1=sidevec_e(2,edg ,i)*utrop_edg(edg ,i)		&
            - sidevec_e(1,edg ,i)*vtrop_edg(edg ,i)		! m^2/s
        valu2=sidevec_e(2,edgx,ipx)*utrop_edg(edgx,ipx)		&
            - sidevec_e(1,edgx,ipx)*vtrop_edg(edgx,ipx)		! m^2/s
        outgo=.5*(valu1-valu2)		! avg.to prevent roundoff errors
       flx(edg,i)=outgo*dpth_edg(edg,i)*onem			! N/sec
       flxdv=flxdv + flx(edg,i)			! set btropfx to zero
!      flxdv=flxdv + flx(edg,i) + btropfx(edg,i)			! N/sec
      end if				! neighbor is ocean point
     end do				! loop through edges
     ptrop(nl,i)=(1.-tsmoo)*ptrop(ml,i) + tsmoo*ptrop(nl,i)	&
                  -(1.+tsmoo)*flxdv*rarea(i)*batrop
    end if				! ocean point
   end do				! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! --- check whether neighboring cells agree on how much mass to exchange

! !SMS$PARALLEL (dh,i) BEGIN
! !SMS$EXCHANGE (flx)
!    do i=1,nip
!     if (wet(i) > 0 ) then
!      do edg=1,nprox(i)
!       ipx=prox(edg,i)
!       edgx=proxs(edg,i)
!       flxerr=flx(edg,i)+flx(edgx,ipx)
!       if (abs(flxerr).gt.1.e-5*(abs(flx(edg,i))+abs(flx(edgx,ipx)))) then
!        valu1=sidevec_e(2,edg ,i  )*utrop_edg(edg ,i  )		&
!            - sidevec_e(1,edg ,i  )*vtrop_edg(edg ,i  )		! m^2/s
!        valu2=sidevec_e(2,edgx,ipx)*utrop_edg(edgx,ipx)		&
!            - sidevec_e(1,edgx,ipx)*vtrop_edg(edgx,ipx)		! m^2/s
!        print '(2(a,i7,i3))','flux mismatch at i,edg=',perm(i),edg,	&
!          ' communicating with',perm(ipx),edgx
!        print '(a)','   i  edg      vnorml            flx        btropfx'
!        print 102,perm(i),edg ,valu1,flx(edg ,i),btropfx(edg ,i)
!        print 102,perm(ipx),edgx,valu2,flx(edgx,ipx),btropfx(edgx,ipx)
!  102   format(i7,i3,3es15.6)
!       end if
!      end do
!     end if
!    end do
! !SMS$PARALLEL END

   call dffusn_lev(ptrop(1,nl),1,1,1,thkdff*batrop,mod(nstep,diag_intvl).eq.0)

!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE (ptrop)

! --- interpolate new -ptrop- to cell edges

!SMS$HALO_COMP(<1,1>) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,edgcount,edg,ipx,im1,ip1)
   do i=1,nip			! horizontal loop
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    ptrop_edg(:,i)=0.
    if (wet(i) > 0 ) then
     do edgcount=1,nedge(i)		! loop through edges
      edg=permedge(edgcount,i)
      ipx=prox(edg,i)          
      if (wet(ipx) == 0 ) ipx=i	! avoid stepping on land

      im1=mod(edg-2+nprox(i),nprox(i))+1
      ip1=mod(edg           ,nprox(i))+1
      im1=prox(im1,i)
      if (wet(im1) == 0 ) im1=i	! avoid stepping on land
      ip1=prox(ip1,i)
      if (wet(ip1) == 0 ) ip1=i	! avoid stepping on land
      ptrop_edg(edg,i)=(8.*(ptrop(nl,i  )+ptrop(nl,ipx))          	&
                          + ptrop(nl,im1)+ptrop(nl,ip1))*divby18
     end do
    end if
   end do				! horizontal loop
!$OMP END PARALLEL DO
!SMS$HALO_COMP END

!sms$compare_var(ptrop    , "barotp.F90 - ptrop3 ")
!sms$compare_var(ptrop_edg, "barotp.F90 - ptrop_edg3 ")

!SMS$PARALLEL END

   mn=ml
   if (vthenu) go to 91

! --- ---------------------
! --- u - momentum equation
! --- ---------------------

90 continue
!  print *,'solving -u- eqn... ml,nl,mn =',ml,nl,mn

!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE (utrop,vtrop)
!$OMP PARALLEL DO PRIVATE (vrbos,edg,dfflxu,ipx,ungbor)
   do i=1,nip	                ! horizontal loop
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    if (wet(i) > 0 ) then

     ! loop through edges and compute line integral of pressure gradient.
     !  Sidevec is a vectorial representation of the edge.

     pgfx(i)=0.
     do edg=1,nprox(i)		! loop through edges
      pgfx(i)=pgfx(i)-ptrop_edg(edg,i)*sidevec_c(2,edg,i)
     end do

! --- pressure force
     pgfx(i)=thref*pgfx(i)*rarea(i)

! --- lateral diffusive fluxes
     dfflxu=0.
     do edg=1,nprox(i)
      ipx=prox(edg,i)
      if (wet(ipx) > 0 ) then
! --- if either i or ipx is a pole point, must transform velocity components
       if (abs(deg_lat(i)).eq.90.) then                ! north or south pole
! --- express neighboring vector in coord.system used at pole
        ungbor= cs(1,edg,i)*utrop(ml,ipx)-sn(1,edg,i)*vtrop(ml,ipx)
       else if (abs(deg_lat(ipx)).eq.90.) then           ! north or south pole
! --- express pole velocity vector in coord.system used at -i-
        ungbor= cs(1,edg,i)*utrop(ml,ipx)+sn(1,edg,i)*vtrop(ml,ipx)
       else
        ungbor=utrop(ml,ipx)
       end if
       dfflxu=dfflxu+(utrop(ml,i)-ungbor)*sideln(edg,i)	! m^2/s
      else			! coastal edge. 
       dfflxu=dfflxu+utrop(ml,i)*(1.-slip)*sideln(edg,i)	! m^2/s
      end if
     end do

     udissp(i)=-veldff*dfflxu*rarea(i)			! m/s^2
     
     !  u/v tendcy is the sum of prs.grad., coriolis term, and ba'trop forcing

! --- coriolis term
     utrop(nl,i)=(1.-tsmoo)*utrop(ml,i) + tsmoo*utrop(nl,i)	&
                  +(1.+tsmoo)*(batrop*(pgfx(i)			&
                  + corio(i)*vtrop(mn,i) + udissp(i))		&
                  + ubforc(i)*rlstep)
    end if				! ocean point

    if (vrbos .and. vthenu) then
     if (lll.eq.1) write (*,100) nstep,itest
100  format ('b a r o t p    time step =',i8,'  i =',i8/	&
        '   uold   unew gradp corio frcng dissp',4x,		&
        '   vold   vnew gradp corio frcng dissp')
     write (*,101)						&
     utrop(ml,i),utrop(nl,i),pgfx(i)*batrop,			&
     corio(i)*vtrop(mn,i)*batrop,				&
     ubforc(i)*rlstep,udissp(i)*batrop,				&
     vtrop(ml,i),vtrop(nl,i),pgfy(i)*batrop,			&
     corio(i)*utrop(mn,i)*batrop,				&
     vbforc(i)*rlstep,vdissp(i)*batrop
101 format (3p,2f7.1,4f6.2,4x,2f7.1,4f6.2)
     call flush(6)
    end if

   end do		                ! horizontal loop
!$OMP END PARALLEL DO

!sms$compare_var(pgfx,  "barotp.F90 - pgfx4 ")
!sms$compare_var(udissp,"barotp.F90 - udissp4 ")
!sms$compare_var(utrop, "barotp.F90 - utrop4 ")

!SMS$PARALLEL END

   mn=nl
   if (vthenu) go to 92

! --- ---------------------
! --- v - momentum equation
! --- ---------------------

91 continue
!  print *,'solving -v- eqn... ml,nl,mn =',ml,nl,mn

!SMS$PARALLEL (dh,i) BEGIN
!SMS$EXCHANGE (utrop,vtrop)
!$OMP PARALLEL DO PRIVATE (vrbos,dfflxv,edg,ipx,vngbor)
   do i=1,nip	                ! horizontal loop
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    if (wet(i) > 0 ) then

     ! loop through edges and compute line integral of pressure gradient.
     !  Sidevec is a vectorial representation of the edge.

     pgfy(i)=0.
     do edg=1,nprox(i)		! loop through edges
      pgfy(i)=pgfy(i)+ptrop_edg(edg,i)*sidevec_c(1,edg,i)
     end do

! --- pressure force
     pgfy(i)=thref*pgfy(i)*rarea(i)

! --- lateral diffusive fluxes
     dfflxv=0.
     do edg=1,nprox(i)
      ipx=prox(edg,i)
      if (wet(ipx) > 0 ) then
! --- if either i or ipx is a pole point, must transform velocity components
       if (abs(deg_lat(i)).eq.90.) then                ! north or south pole
! --- express neighboring vector in coord.system used at pole
        vngbor= sn(1,edg,i)*utrop(ml,ipx)+cs(1,edg,i)*vtrop(ml,ipx)
       else if (abs(deg_lat(ipx)).eq.90.) then           ! north or south pole
! --- express pole velocity vector in coord.system used at -i-
        vngbor=-sn(1,edg,i)*utrop(ml,ipx)+cs(1,edg,i)*vtrop(ml,ipx)
       else
        vngbor=vtrop(ml,ipx)
       end if
       dfflxv=dfflxv+(vtrop(ml,i)-vngbor)*sideln(edg,i)	! m^2/s
      else			! coastal edge. 
       dfflxv=dfflxv+vtrop(ml,i)*(1.-slip)*sideln(edg,i)	! m^2/s
      end if
     end do

     vdissp(i)=-veldff*dfflxv*rarea(i)			! m/s^2
     
! --- coriolis term
     vtrop(nl,i)=(1.-tsmoo)*vtrop(ml,i) + tsmoo*vtrop(nl,i)	&
                  +(1.+tsmoo)*(batrop*(pgfy(i)			&
                  - corio(i)*utrop(mn,i) + vdissp(i))		&
                  + vbforc(i)*rlstep)
    end if				! ocean point

    if (vrbos .and. .not.vthenu) then
     if (lll.eq.1) write (*,100) nstep,itest
     write (*,101)						&
     utrop(ml,i),utrop(nl,i),pgfx(i)*batrop,			&
     corio(i)*vtrop(mn,i)*batrop,				&
     ubforc(i)*rlstep,udissp(i)*batrop,				&
     vtrop(ml,i),vtrop(nl,i),pgfy(i)*batrop,			&
     corio(i)*utrop(mn,i)*batrop,				&
     vbforc(i)*rlstep,vdissp(i)*batrop
     call flush(6)
    end if

   end do		                ! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

!sms$compare_var(pgfy,  "barotp.F90 - pgfy5 ")
!sms$compare_var(vdissp,"barotp.F90 - vdissp5 ")
!sms$compare_var(vtrop, "barotp.F90 - vtrop5 ")

   mn=nl
   if (vthenu) go to 90

! --- switch order in which -u,v- equations are solved
92 vthenu=.not.vthenu

! --- done with both u and v eqns. new time level becomes old time level
   nx=ml
   ml=nl
   nl=nx

80 continue

   call findmxmn2(utrop,3,nip,nl,'ptrop',wet)
   call findmxmn2(utrop,3,nip,nl,'utrop',wet)
   call findmxmn2(vtrop,3,nip,nl,'vtrop',wet)

   if (mod(nstep,diag_intvl).eq.0) then
    call stencl(ptrop,3,1./onem,'(ocn barotp) ptrop (m)')
    call stencl(utrop,3,1.e3,'(ocn barotp) utrop (mm/s)')
    call stencl(vtrop,3,1.e3,'(ocn barotp) vtrop (mm/s)')
    call stencl(pgfx ,1,1.e5,'(ocn barotp) pgfx x 10^5')
    call stencl(pgfy ,1,1.e5,'(ocn barotp) pgfy x 10^5')
   end if

!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!  if (bclin_frq.eq.1000) stop
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   return
   end subroutine barotp
end module hycom_barotp
