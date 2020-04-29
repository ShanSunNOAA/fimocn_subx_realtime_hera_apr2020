module hycom_fct3d
contains

  subroutine fct3d(nstep,fld3d,numfld,u,w,area,rarea,		&
                   thko,thkn,diagno)

! --- 3-d transport routine adapted from HYCOM

! --- fld3d - field(s) to be transported
! --- u,w   - mass fluxes (x time step) satisfying continuity equation
! ---         (w(k) > 0 means mass flows from layer k to layer k+1)
! --- area  - grid cell size
! --- rarea - inverse of area
! --- thko,thkn - layer thickness at beginning & end of time interval
  
  use fimnamelist  ,   only: kdm,itest,diag_intvl
  use hycom_control,   only: kgap
  use hycom_constants, only: wet,land_spval
  use hycom_variables, only: totmass
  use module_control,  only: npp,nip
  use module_constants,only: nprox,prox,proxs,grvity,perm
  use findmaxmin2
  use findmaxmin3
  use stencilprint

  implicit none
  integer,intent(IN) :: nstep		! model time step
  integer,intent(IN) :: numfld		! total number of fields in 'fld3d'
  logical,intent(IN) :: diagno		! activate diagnostic calculations
!SMS$DISTRIBUTE (dh,1) BEGIN
  real,   intent(IN) :: area(nip),rarea(nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
  real,intent(INOUT) :: fld3d(kdm,nip,numfld)
  real,   intent(IN) :: w(kdm,nip),thko(kdm,nip),thkn(kdm,nip)
! Local variables:
  real fld(kdm,nip),vertfx(kdm,nip),vertdv(kdm,nip),		&
       fmx(kdm,nip),fmn(kdm,nip),flp(kdm,nip),fln(kdm,nip),	&
       hordiv(kdm,nip),fldold(kdm,nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,3) BEGIN
  real,  intent(IN) :: u(kdm,npp,nip)
  real flx(kdm,npp,nip),uan(kdm,npp,nip)
!SMS$DISTRIBUTE END
  real*8 a(kdm),b(kdm),c(kdm),dx,fcdx,yl,yr,totin,totou,	&
         fldlo,q,clip,slab,dslab,thkchg,dxlft,dxmid,dxrgt,	&
         col1,col2,wtold,wtnew,bfore,after,amount,bounded,unbound
  integer i,edg,ix,edx,k,kp,type
  character string*24
  logical vrbos
  logical,parameter :: recovr=.true.		! enforce global conservation
! logical,parameter :: recovr=.false.		! enforce global conservation
  real   ,parameter :: athird=1./3.,epsil=1.e-11,onemu=1.e-6

! print *,'entering fct3d...'

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,thkchg)
  do i=1,nip
   if (wet(i) > 0 ) then
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0

    if (vrbos) then
! --- check mass conservation in test column
     write (*,'(i8,a/a)') perm(i),				&
      '  fct3d -- time-integrated continuity eqn diagnostics:',	&
      '     thknss_tndcy  horiz.flxdiv   vert.flxdiv      residuum'
     do k=1,kdm			! vertical loop
      thkchg=thkn(k,i)-thko(k,i)
      hordiv(k,i)=0.
      do edg=1,nprox(i)
       hordiv(k,i)=hordiv(k,i)+u(k,edg,i)
      end do
      hordiv(k,i)=hordiv(k,i)*rarea(i)
      if (k.eq.1) then
       write (*,103) k,thkchg,hordiv(k,i),w(k,i),		&
        thkchg+hordiv(k,i)+w(k,i)
      else
       write (*,103) k,thkchg,hordiv(k,i),w(k,i)-w(k-1,i),	&
        thkchg+hordiv(k,i)+w(k,i)-w(k-1,i)
      end if
     end do			! vertical loop
 103 format (i3,4es14.4)
    end if			! vrbos
   end if			! ocean point
  end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! --- optional: check mass conservation globally in select layers
  if (diagno) then

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
   do i=1,nip			! horiz. loop
   if (wet(i) > 0 ) then
     do k=1,kdm,kgap		! vertical loop
      hordiv(k,i)=0.
      do edg=1,nprox(i)
       hordiv(k,i)=hordiv(k,i)+u(k,edg,i)
      end do
      hordiv(k,i)=hordiv(k,i)*rarea(i)
      if (k.eq.1) then
        hordiv(k,i)=hordiv(k,i)					&
          +w(k,i)         +thkn(k,i)-thko(k,i)
      else
        hordiv(k,i)=hordiv(k,i)					&
          +w(k,i)-w(k-1,i)+thkn(k,i)-thko(k,i)
      end if
     end do			! vertical loop
    end if			! ocean point
   end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

   do k=1,kdm,kgap
    write (string,'(a,i2)') ' k=',k
    call findmxmn2(hordiv,kdm,nip,k,'o-fct3d hordiv'//string,wet)
    call findmxmn2(thko,kdm,nip,k,'o-fct3d thko'//string,wet)
    call findmxmn2(thkn,kdm,nip,k,'o-fct3d thkn'//string,wet)
   end do			! vertical loop
  end if			! diagno

!SMS$EXCHANGE(fld3d,thko,thkn)

  do 1 type=1,numfld		! loop over fields to be transported

  if (diagno) then
   do k=1,kdm,kgap
    write (string,'(a,i1,a,i2)') 'fct3d  IN fld',type,' k=',k
    call findmxmn3(fld3d,kdm,nip,numfld,k,type,trim(string),wet)
   end do
  end if

! --- get vertical flux by summing -fld- over upstream slab of thickness -w-

!SMS$PARALLEL (dh,i) BEGIN
!SMS$HALO_COMP(<1,1>) BEGIN
!$OMP PARALLEL DO
  do i=1,nip			! horiz. loop
   fld(:,i)=land_spval
   if (wet(i) > 0 ) then
    do k=1,kdm
     fld(k,i)=fld3d(k,i,type)
    end do
   end if
  end do			! horiz. loop
!$OMP END PARALLEL DO

!!$OMP PARALLEL DO
!  do i=1,nip			! horiz. loop
!   if (wet(i) > 0 ) then
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! --- optional:
! --- fill massless cells with data from layer(s) above or below
!   do k=kdm-1,1,-1
!    fld(k,i)=(fld(k,i)*thko(k,i)+fld(k+1,i)*onemu)	&
!            /(         thko(k,i)+           onemu)
!   end do
!   do k=2,kdm
!    fld(k,i)=(fld(k,i)*thko(k,i)+fld(k-1,i)*onemu)	&
!            /(         thko(k,i)+           onemu)
!   end do
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!   end if			! ocean point
!  end do			! horiz. loop
!!$OMP END PARALLEL DO

!SMS$HALO_COMP END

!$OMP PARALLEL DO PRIVATE (a,b,c,dxlft,dxmid,dxrgt,yl,yr,slab,amount,	&
!$OMP kp,dslab,dx,fcdx)
  do i=1,nip			! horiz. loop
   if (wet(i) > 0 ) then
! --- fit 0th, 1st, or 2nd deg. polynomial to tracer in each cell
    a(1  )=fld(1  ,i)
    b(1  )=0.
    c(1  )=0.
    a(kdm)=fld(kdm,i)
    b(kdm)=0.
    c(kdm)=0.
    do k=2,kdm-1		! vertical loop
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- piecewise constant method:
!    a(k)=fld(k,i)
!    b(k)=0.
!    c(k)=0.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- piecewise linear method:
! --- fit linear function a+bx to tracer in each cell (-.5 < x < +.5)
!    a(k)=fld(k,i)
!    b(k)=0.
!    if (fld(k,i).le.min(fld(k-1,i),fld(k+1,i)) .or.		&
!        fld(k,i).ge.max(fld(k-1,i),fld(k+1,i))) then
!      b(k)=0.
!    else if ((fld(k+1,i)-fld(k-1,i))*(fld(k-1,i)+fld(k+1,i)	&
!      -2.*fld(k,i)).gt.0.) then
!      b(k)=fld(k,i)-fld(k-1,i)
!    else
!      b(k)=fld(k+1,i)-fld(k,i)
!    end if
!    c(k)=0.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- piecewise parabolic method:
! --- construct parabola   a+bx+cx^2  whose integral over [-.5,+.5] equals
! --- fld(k) and which passes though points yl,yr at [-.5,+.5] resp.
     dxlft=max(epsil,thko(k-1,i))
     dxmid=max(epsil,thko(k  ,i))
     dxrgt=max(epsil,thko(k+1,i))
     yl=(dxlft*fld(k,i)+dxmid*fld(k-1,i))/(dxlft+dxmid)
     yr=(dxrgt*fld(k,i)+dxmid*fld(k+1,i))/(dxrgt+dxmid)

     a(k)=1.5*fld(k,i)-.25*(yl+yr)
     b(k)=yr-yl
     c(k)=6.*(.5*(yl+yr)-fld(k,i))
     if (abs(yr-yl) .lt. 6.*abs(.5*(yl+yr)-fld(k,i))) then
! --- apex of parabola o !urs inside interval [-.5,+.5], implying an over-
! --- or undershoot situation. change curve to prevent over/undershoots.
      if (abs(yr-yl) .gt. 2.*abs(.5*(yl+yr)-fld(k,i))) then
! --- put apex of parabola on edge of interval [-.5,+.5]
       if ((yr-yl)*(.5*(yl+yr)-fld(k,i)) .gt. 0.) then
! --- apex at x=-.5
        a(k)=.25*(3.*fld(k,i)+yl)
        c(k)=3.*(fld(k,i)-yl)
        b(k)=c(k)
       else
! --- apex at x=+.5
        a(k)=.25*(3.*fld(k,i)+yr)
        c(k)=3.*(fld(k,i)-yr)
        b(k)=-c(k)
       end if
      else			!  -1/6 < x < +1/6
! --- moving apex won't help. replace parabola by constant.
       a(k)=fld(k,i)
       b(k)=0.
       c(k)=0.
      end if
     end if
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    end do			! vertical loop

    do k=1,kdm-1		! vertical loop
     slab=onemu
     if (w(k,i).lt.0.) then			! interface moves down
      amount=slab*fld(k+1,i)
      kp=k
 24   kp=kp+1
      if (slab.ge.-w(k,i)) goto 23
      if (thko(kp,i).gt.0.) then
       dslab=min(slab+thko(kp,i),-w(k,i))	&
            -min(slab           ,-w(k,i))
       dx=dslab/thko(kp,i)
       fcdx=a(kp)				&
           +b(kp)*.5*(dx-1.)			& !  not needed in pcm
           +c(kp)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
       amount=amount+fcdx*dslab
       slab=slab+dslab
      end if
      if (kp.lt.kdm) go to 24
     else if (w(k,i).gt.0.) then		! interface moves up
      amount=slab*fld(k,i)
      kp=k+1
 25   kp=kp-1
      if (slab.ge.w(k,i)) goto 23
      if (thko(kp,i).gt.0.) then
       dslab=min(slab+thko(kp,i), w(k,i))	&
            -min(slab           , w(k,i))
       dx=dslab/thko(kp,i)
       fcdx=a(kp)				&
           +b(kp)*.5*(1.-dx)			& !  not needed in pcm
           +c(kp)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
       amount=amount+fcdx*dslab
       slab=slab+dslab
      end if
      if (kp.gt.2) go to 25
     end if
 23  vertfx(k,i)=w(k,i)*amount/slab
    end do			! vertical loop

    vertfx(kdm,i)=0.		!  don't allow flux thru bottom
    vertdv(1,i)=vertfx(1,i)
    do k=2,kdm
     vertdv(k,i)=vertfx(k,i)-vertfx(k-1,i)
    end do
   end if			! ocean point
  end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  if (recovr) then	! bitwise reproducibility not guaranteed if recovr=T
   bfore=0.
   wtold=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(col1,col2) REDUCTION(+:bfore,wtold)
   do i=1,nip			! horiz. loop
    if (wet(i) > 0 ) then
     col1=0.
     col2=0.
     do k=1,kdm
      col1=col1+         thko(k,i)*area(i)
      col2=col2+fld(k,i)*thko(k,i)*area(i)
     end do
     wtold=wtold+col1
     bfore=bfore+col2
    end if
   end do
!SMS$PARALLEL END
!SMS$REDUCE (bfore,wtold,SUM)
  end if			! recovr

! --- compute low-order & antidiffusive (high- minus low-order) fluxes
! --- find upper/lower limits (fmx,fmn) for transported field

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (ix,edx,q)
  do i=1,nip			! horiz. loop
   flx(:,:,i)=land_spval
   uan(:,:,i)=land_spval
   if (wet(i) > 0 ) then
    do k=1,kdm			! vertical loop
     fmx(k,i)=fld(k,i)
     fmn(k,i)=fld(k,i)

     do edg=1,nprox(i)
      flx(k,edg,i)=0.
      uan(k,edg,i)=0.
      ix=prox(edg,i)
      if (wet(ix) > 0 ) then
! +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+
!!     if (u(k,edg,i).ge.0.) then		! out-of cell > 0
!!      q=fld(k,i )
!!     else
!!      q=fld(k,ix)
!!     end if
!!     flx(k,edg,i)=u(k,edg,i)*q		! low order
! +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+
       edx=proxs(edg,i)		! index of joint edge as seen by neighbr
       flx(k,edg,i)=0.5*(		&	! low-order (out-of cell > 0)
           (u(k,edg,i )+abs(u(k,edg,i )))*fld(k,i )		&
         - (u(k,edx,ix)+abs(u(k,edx,ix)))*fld(k,ix))
! +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+
       q=.5*(fld(k,i)+fld(k,ix))		! high (2nd) order
       uan(k,edg,i)=q*u(k,edg,i)-flx(k,edg,i)	! antidiffusive
       fmx(k,i)=max(fmx(k,i),fld(k,ix))
       fmn(k,i)=min(fmn(k,i),fld(k,ix))
      end if			! neighbor wet/dry
     end do			! loop over edges
     if (k.lt.kdm) then
      if (w(k  ,i).lt.0.) then
       fmx(k,i)=max(fmx(k,i),vertfx(k,i)/w(k,i))
       fmn(k,i)=min(fmn(k,i),vertfx(k,i)/w(k,i))
      end if
     end if
     if (k.gt.1) then
      if (w(k-1,i).gt.0.) then
       fmx(k,i)=max(fmx(k,i),vertfx(k-1,i)/w(k-1,i))
       fmn(k,i)=min(fmn(k,i),vertfx(k-1,i)/w(k-1,i))
      end if
     end if
    end do			! vertical loop
   end if			! ocean point
  end do			! horiz. loop
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE (vrbos,unbound,bounded)
  do i=1,nip			! horiz. loop
   if (wet(i) > 0 ) then
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    do k=1,kdm			! vertical loop
     hordiv(k,i)=0.
     do edg=1,nprox(i)
      hordiv(k,i)=hordiv(k,i)+flx(k,edg,i)
     end do
     hordiv(k,i)=hordiv(k,i)*rarea(i)

     unbound=fld(k,i)*thko(k,i)-hordiv(k,i)-vertdv(k,i)
     bounded=max(0.,fmn(k,i)*thkn(k,i),min(unbound,fmx(k,i)*thkn(k,i)))
     fldold(k,i)=fld(k,i)

! --- step 1: low-order advection

     fld(k,i)=(fld(k,i)*onemu+bounded)/(onemu+thkn(k,i))
    end do			! vertical loop
   end if			! ocean point
  end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! --- at each grid point, determine the ratio of the largest permissible
! --- pos. (neg.) change in -fld- to the sum of all incoming (outgoing) fluxes

!SMS$EXCHANGE(fld)
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (totin,totou,ix)
  do i=1,nip			! horiz. loop
   if (wet(i) > 0 ) then
    do k=1,kdm			! vertical loop
     totin=0.
     totou=0.
     do edg=1,nprox(i)
      totin=totin-min(0.,uan(k,edg,i))
      totou=totou+max(0.,uan(k,edg,i))
      ix=prox(edg,i)
      if (wet(ix) > 0 ) then
! --- update upper/lower limits for transported field
       fmx(k,i)=max(fmx(k,i),fld(k,ix))
       fmn(k,i)=min(fmn(k,i),fld(k,ix))
      end if			! neighbor wet/dry
     end do			! loop over edges
     flp(k,i)=(fmx(k,i)-fld(k,i))*thkn(k,i)/(totin+epsil)*rarea(i)
     fln(k,i)=(fld(k,i)-fmn(k,i))*thkn(k,i)/(totou+epsil)*rarea(i)
    end do			! vertical loop
   end if			! ocean point
  end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! if (diagno) then
!  do k=1,kdm,kgap
!   write (string,'(a,i2.2,a,i2)') 'o-fct3d fld',type,' k=',k
!   call findmxmn2(fld,kdm,nip,k,string,wet)
!   write (string,'(a,i2.2,a,i2)') 'o-fct3d fmx',type,' k=',k
!   call findmxmn2(fmx,kdm,nip,k,string,wet)
!   write (string,'(a,i2.2,a,i2)') 'o-fct3d fmn',type,' k=',k
!   call findmxmn2(fmn,kdm,nip,k,string,wet)
!   write (string,'(a,i2.2,a,i2)') 'o-fct3d flp',type,' k=',k
!   call findmxmn2(flp,kdm,nip,k,string,wet)
!   write (string,'(a,i2.2,a,i2)') 'o-fct3d fln',type,' k=',k
!   call findmxmn2(fln,kdm,nip,k,string,wet)
!   write (string,'(a,i2.2,a,i2)') 'o-fct3d fmx',type,' k=',k
!   call stencl(fmx,kdm,1.,string)
!   write (string,'(a,i2.2,a,i2)') 'o-fct3d fmn',type,' k=',k
!   call stencl(fmn,kdm,1.,string)
!  end do
! end if

! --- limit antidiffusive fluxes

!SMS$EXCHANGE(flp,fln)

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,ix,clip,unbound,bounded,fldlo)
  do i=1,nip			! horiz. loop
   if (wet(i) > 0 ) then
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    do k=1,kdm			! vertical loop
     hordiv(k,i)=0.
     do edg=1,nprox(i)
      ix=prox(edg,i)
      clip=1.
      if (wet(ix) > 0 ) then
       if (uan(k,edg,i).ge.0.) then		! out-of cell > 0
        clip=min(1.,fln(k,i),flp(k,ix))
       else
        clip=min(1.,flp(k,i),fln(k,ix))
       end if
       flx(k,edg,i)=uan(k,edg,i)*clip
       hordiv(k,i)=hordiv(k,i)+flx(k,edg,i)
      end if			! neighbor is ocean point
     end do			! loop over edges
     hordiv(k,i)=hordiv(k,i)*rarea(i)		! antidiffusive divergence
     unbound=fld(k,i)*thkn(k,i)-hordiv(k,i)
     bounded=max(0.,fmn(k,i)*thkn(k,i),min(unbound,fmx(k,i)*thkn(k,i)))
     fldlo=fld(k,i)

! --- step 2: high-order advection

     fld(k,i)=(fld(k,i)*onemu+bounded)/(onemu+thkn(k,i))

     if (vrbos .and. mod(k,kgap).eq.1) then
      write (*,100) '(ocn fct3d) i,k=',perm(i),k,' fld',type,	&
        ' old/lo_ord/new:',fldold(k,i),fldlo,fld(k,i)
100   format (a,i7,i3,a,i2.2,a,3f10.2)
     end if

    end do			! vertical loop
   end if			! ocean point
  end do			! horiz. loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  if (recovr) then	! bitwise reproducibility not guaranteed if recovr=T
   after=0.
   wtnew=0.
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE(col1,col2) REDUCTION(+:after,wtnew)
   do i=1,nip			! horiz. loop
    if (wet(i) > 0 ) then
     col1=0.
     col2=0.
     do k=1,kdm
      col1=col1+         thkn(k,i)*area(i)
      col2=col2+fld(k,i)*thkn(k,i)*area(i)
     end do
     wtnew=wtnew+col1
     after=after+col2
    end if
   end do
!SMS$PARALLEL END
!SMS$REDUCE (after,wtnew,SUM)

!  bfore=bfore/(grvity*totmass)
!  after=after/(grvity*totmass)
   bfore=bfore/wtold
   after=after/wtnew
!  if (diagno)							&
     write (*,'(2a,i2,a,2f13.8)') '(ocn_fct3d) ',		&
     '  tracer',type,' concentr.  before,after=',bfore,after
   q=bfore-after
   write (*,'(a,es10.2,a,i2)') '(ocn_fct3d)  add',q,' to tracer',type

! --- modify transported fields to assure global conservation
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
    do i=1,nip
     if (wet(i) > 0 ) then
      do k=1,kdm
       fld(k,i)=fld(k,i)+q
      end do
     end if
    end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
  end if			! recovr

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
  do i=1,nip
   if (wet(i) > 0 ) then
    do k=1,kdm
     fld3d(k,i,type)=fld(k,i)
    end do
   end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  if (diagno) then
   do k=1,kdm,kgap
    write (string,'(a,i1,a,i2)') 'fct3d OUT fld',type,' k=',k
    call findmxmn3(fld3d,kdm,nip,numfld,k,type,trim(string),wet)
   end do
  end if

1 continue			! loop over tracers
! print *,'exiting fct3d...'

  return
  end subroutine fct3d
end module hycom_fct3d
