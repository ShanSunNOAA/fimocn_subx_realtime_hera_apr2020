module  hycom_transp3d
!*********************************************************************
!	3-dimensional tracer transport module designed for intermittent
!	(i.e., long time step) execution. there are 3 entries:
!
!	transp0 - initializes mass flux arrays and saves initial -dp-
!	transp1 - builds up time integral of horizontal mass fluxes
!	transp2 - performs the actual transport operation

!	R. Bleck	     March 2009
!*********************************************************************
  use findmaxmin1
  use findmaxmin2
  use stencilprint
  use edgmaxmin

contains
  subroutine transp0(nstep,leapn,cumuflx,dp,dpinit)
  use fimnamelist    ,only: kdm
  use hycom_constants,only: wet,onem,land_spval
  use module_control ,only: npp,nip
  implicit none

  integer,intent(IN)  :: nstep			! model time step
  integer,intent(IN)  :: leapn
!SMS$DISTRIBUTE (dh,2) BEGIN
  real,intent(IN)     :: dp     (kdm,    nip,2)	! layer thickness
  real,intent(INOUT)  :: dpinit (kdm,    nip)	! dp at start of time integr.
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,3) BEGIN
  real,intent(INOUT)  :: cumuflx(kdm,npp,nip)	! time-integrated mass flux
!SMS$DISTRIBUTE END
  integer    :: i,k

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
  do i=1,nip			! horizontal loop
   cumuflx(:,:,i)=land_spval
   if (wet(i) > 0 ) then
    do k=1,kdm
     cumuflx(k,:,i)=0.
     dpinit(k,i)=dp(k,i,leapn)
    end do
   end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  print *,'ocn tracer transport arrays initialized, time step',nstep
  return
  end subroutine transp0


  subroutine transp1(nstep,cumuflx,massflx)
  use fimnamelist     ,only: kdm
  use hycom_constants ,only: wet,batrop
  use module_control  ,only: npp,nip
  use module_constants,only: nprox,rarea
  implicit none
  integer,intent(IN)  :: nstep			! model time step

!SMS$DISTRIBUTE (dh,3) BEGIN
  real,intent(IN)     :: massflx(kdm,npp,nip)
  real,intent(INOUT)  :: cumuflx(kdm,npp,nip)
!SMS$DISTRIBUTE END
  integer    :: i		! Index for is point number
  integer    :: edg,k

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
  do i=1,nip			! horizontal loop
   if (wet(i) > 0 ) then
    do edg=1,nprox(i)
     do k=1,kdm
      cumuflx(k,edg,i)=cumuflx(k,edg,i)+massflx(k,edg,i)*2.*batrop
     end do
    end do
   end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  print *,'ocn mass fluxes added to time integral,  time step',nstep
  return
  end subroutine transp1


  subroutine transp2 (nstep,temp,saln,passv_tr,cumuflx,dpinit,dpfinl)
  use hycom_control,   only: numtr,bclin_frq,kgap
  use hycom_constants, only: wet,land_spval,batrop
  use module_control,  only: npp,nip
  use fimnamelist,     only: stencl_frst,kdm,itest,diag_intvl,temdff
  use module_constants,only: nprox,prox,area,rarea,perm
  use hycom_fct3d
  use hycom_dffusn,only: dffusn_lyr
  use hycom_diag,only: glob3d
  implicit none
#include <gptl.inc>

! External variables:
  integer,intent (IN)    :: nstep			! model time step
!SMS$DISTRIBUTE (dh,1) BEGIN
  real col_xpand(nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
  real   ,intent (INOUT) :: temp    (kdm,nip)		! temperature
  real   ,intent (INOUT) :: saln    (kdm,nip)		! salinity
  real   ,intent (INOUT) :: passv_tr(kdm,nip,numtr)	! passive tracer
  real   ,intent (IN)    :: dpinit  (kdm,nip)		! init'l lyr thknss
  real   ,intent (IN)    :: dpfinl  (kdm,nip)		! final lyr thknss
! Local variables:
  real :: vertfx (kdm,nip)
  real :: field  (kdm,nip)
  real :: all_trc(kdm,nip,2+numtr)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,3) BEGIN
  real   ,intent (IN)    :: cumuflx(kdm,npp,nip)	! time-integr. mass flx
!SMS$DISTRIBUTE END
  integer    :: k		! layer index
  integer    :: i		! Index for is point number
  integer    :: edg		! is edge number index
  integer    :: type		! tracer index
  integer    :: ret
  logical    :: vrbos
  character  :: text*24,string*24
  real       :: hordiv(kdm),vertdv(kdm)
! real*8     :: tstart,tstop,valmin,valmax
  real*8     :: totlin(2),totout(2),bias
  real,parameter :: toffset=10.	! used to keep temperature pos.definite

!sms$compare_var(temp, "transp3d.F90 - temp    ")
!sms$compare_var(saln, "transp3d.F90 - saln    ")
!sms$compare_var(cumuflx , "transp3d.F90 - cumuflx1 ")

! --- compute the various terms in the continuity equation integrated
! --- over time interval since last call to -transp0-
! --- the continuity eqn is split into horiz. and vert. terms as follows:
! ---        (dpfinl-dpinit) + hordiv + verdiv = 0

  write (text,'(a,i8)') 'ocn transp2 step',nstep

! call edgmxmn(cumuflx,kdm,'(ocn-transp2) cumuflx',wet)

! --- diagnose vertical mass fluxes

  vertfx(:,:)=land_spval
!SMS$EXCHANGE(cumuflx)
!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO PRIVATE (vrbos,hordiv,vertdv)
  do i=1,nip				! horizontal loop
   if (wet(i) > 0 ) then
    vrbos=i.eq.itest .and. mod(nstep,diag_intvl).eq.0
    hordiv   (:  )=0.
    col_xpand(  i)=0.
    vertfx   (:,i)=0.

    do edg=1,nprox(i)			! loop through edges
     do k=1,kdm				! loop through layers
      hordiv(k)=hordiv(k)+cumuflx(k,edg,i)
     end do
    end do

    do k=1,kdm				! loop through layers
     col_xpand(i)=col_xpand(i)+hordiv(k)
     vertdv(k)=(dpinit(k,i)-dpfinl(k,i))-hordiv(k)*rarea(i)
     if (k.eq.1) then
      vertfx(k,i)=              vertdv(k)
     else
      vertfx(k,i)=vertfx(k-1,i)+vertdv(k)
     end if

     if (vrbos .and. mod(k,kgap).eq.1) then
      write (*,'(i8,i3,a,3es12.4)') perm(i),k,				&
       ' transp2 hordiv,col_xpand,vertfx:',hordiv(k)*rarea(i),		&
       col_xpand(i)*rarea(i),vertfx(k,i)
      call flush(6)
     end if

    end do
   end if
  end do				! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! call findmxmn1(col_xpand,nip,'col_xpand',wet)
! do k=1,kdm,kgap
!  write (string,'(a,i3,a)') 'lyr',k,' vertfx'
!  call findmxmn2(vertfx,kdm,nip,k,string,wet)
! end do
! print *
! call flush(6)

! --- having determined the vertical flux term in the time-integrated
! --- continuity eqn, we can now perform the actual tracer transport

! --- combine active and passive tracers in single array

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
  do i=1,nip
   all_trc(:,i,:)=land_spval
   if (wet(i) > 0 ) then
    all_trc(:,i,1)=temp(:,i)
    all_trc(:,i,2)=saln(:,i)
    do type=1,numtr
     all_trc(:,i,2+type)=passv_tr(:,i,type)
    end do
   end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  if (mod(nstep,diag_intvl).eq.0) then
   call stencl(temp,kdm,1.,text//'  temperature IN')
   call stencl(saln,kdm,1.,text//'  salinity IN')
   do type=1,numtr
     write (string,'(a,i2,a)') '  trcr',type,' IN'
    call stencl(passv_tr(:,:,type),kdm,1.,text//trim(string))
   end do
  end if

! do type=1,2+numtr
!  totlin(:)=glob3d(all_trc(:,:,type),dpinit)
! end do

!SMS$PARALLEL (dh,i) BEGIN
!JR Could thread but not enough work?
  do i=1,nip
   if (wet(i) > 0 ) all_trc(:,i,1)=all_trc(:,i,1)+toffset	! trcr 1 = temp
  end do
!SMS$PARALLEL END

! call StartTimer(tstart)

  ret = gptlstart('hycom_transp3d')
  call fct3d(nstep,all_trc,2+numtr,cumuflx,vertfx,area,rarea,		&
             dpinit,dpfinl,mod(nstep,diag_intvl).eq.0)
  ret = gptlstop ('hycom_transp3d')

! tstop=0
!  call IncrementTimer(tstart,tstop)
! valmin=tstop*1.e6
! valmax=valmin
!!SMS$REDUCE(valmin,MIN)
!!SMS$REDUCE(valmax,MAX)
!  print '(a,i2,a,2(i9,a))','time spent in subr fct3d (',2+numtr,	&
!   ' tracers) ',nint(valmin),' -',nint(valmax),' usec'

!SMS$PARALLEL (dh,i) BEGIN
  do i=1,nip
   if (wet(i) > 0 ) all_trc(:,i,1)=all_trc(:,i,1)-toffset	! trcr 1 = temp
  end do
!SMS$PARALLEL END

! --- recover active and passive tracers from combined array

!SMS$PARALLEL (dh,i) BEGIN
!$OMP PARALLEL DO
  do i=1,nip
   if (wet(i) > 0 ) then
    temp(:,i)=all_trc(:,i,1)
    saln(:,i)=all_trc(:,i,2)
    do type=1,numtr
     passv_tr(:,i,type)=all_trc(:,i,2+type)
    end do
   end if
  end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END

! do type=1,2+numtr
!  totout(type)=glob3d(all_trc(:,:,type),dpfinl)
! end do                       ! loop over tracers

! --- transport done. now diffuse  a c t i v e  tracer laterally

  print *,'now diffusing temp,saln...'

  if (temdff.gt.0.) then
!JR Ensure vrbos is set consistently
   vrbos= itest > 0 .and. mod(nstep,diag_intvl).eq.0
   call dffusn_lyr(temp,dpfinl,kdm,1,kdm,temdff*batrop*bclin_frq,vrbos)
   call dffusn_lyr(saln,dpfinl,kdm,1,kdm,temdff*batrop*bclin_frq,vrbos)
  end if

  print *,'...back from dffusn_lyr'

  if (mod(nstep,diag_intvl).eq.0) then
   call stencl(temp,kdm,1.,text//'  temperature OUT')
   call stencl(saln,kdm,1.,text//'  salinity OUT')
   do type=1,numtr
     write (string,'(a,i2,a)') '  trcr',type,' OUT'
    call stencl(passv_tr(:,:,type),kdm,1.,text//trim(string))
   end do
  end if

! print 100,'(transp2) INITL tracer mean',totlin
! print 100,'(transp2) INTMD tracer mean',totout
100 format (a,6f13.8)

! field(:,:)=land_spval
! do type=1,2+numtr
!  bias=totlin(type)-totout(type)
!  print '(a,es10.2,a,i2.2)','(transp2) adding',bias,' to tracer',type
!!SMS$PARALLEL (dh,i) BEGIN
!  do i=1,nip
!   if (wet(i) > 0 ) then
!    do k=1,kdm
!     all_trc(k,i,type)=all_trc(k,i,type)+bias
!     field(k,i)=activ_tr(k,i,type)
!    end do
!   end if
!  end do
!!SMS$PARALLEL END
!  totout(type)=glob3d(field,dpfinl)
! end do		! loop over tracers
! print 100,'(transp2) FINAL tracer mean',totout

  print *,'ocn tracers transported,  time step',nstep
  return
  end subroutine transp2

end module hycom_transp3d
