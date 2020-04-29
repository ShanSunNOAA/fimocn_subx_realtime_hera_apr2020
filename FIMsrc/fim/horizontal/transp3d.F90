module  module_transp3d
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

  implicit none

#include <gptl.inc>

contains

  subroutine transp0(its,cumufx,dp3d,dpinit)
  use module_control ,only: npp,nip
  use fimnamelist    ,only: nvl

  integer,intent(IN)  :: its		! model time step
!sms$distribute (dh,2) begin
  real,intent(IN)     :: dp3d  (nvl,    nip)	! layer thickness
  real,intent(OUT)    :: dpinit(nvl,    nip)	! dp at start of time integr.
!sms$distribute end

!sms$distribute (dh,3) begin
  real,intent(OUT)    :: cumufx(nvl,npp,nip)	! time-integrated mass flux
!sms$distribute end

  cumufx(:,:,:)=0.
  dpinit(:,:)=dp3d(:,:)

  print *,'tracer transport arrays initialized, time step',its
  return
  end subroutine transp0


  subroutine transp1(its,nf,of,vof,adbash1,adbash2,adbash3,	&
                     cumufx,massfx)
  use module_control ,only: npp,nip
  use fimnamelist    ,only: nvl
  use module_constants,only: nprox,rarea

  implicit none
  integer,intent(IN)  :: its		! model time step
  integer,intent(IN)  :: nf,of,vof	! time slots: new,old,very old
  real   ,intent(IN)  :: adbash1,adbash2,adbash3

!sms$distribute (dh,3) begin
  real,intent(IN)     :: massfx(nvl,npp,nip,3)
  real,intent(INOUT)  :: cumufx(nvl,npp,nip)
!sms$distribute end
  integer    :: ico		! Index for icos point number
  integer    :: edg,k
  logical    :: vrbos

!SMS$PARALLEL (dh,ico) BEGIN
  do ico=1,nip			! horizontal loop
   do edg=1,nprox(ico)
    do k=1,nvl
     cumufx(k,edg,ico)=cumufx(k,edg,ico)			&
       +adbash1*massfx(k,edg,ico, nf)				&
       +adbash2*massfx(k,edg,ico, of)				&
       +adbash3*massfx(k,edg,ico,vof)
    end do
   end do
  end do
!SMS$PARALLEL END

! print *,'mass fluxes added to time integral,  time step',its
  return
  end subroutine transp1


  subroutine transp2 (its,             &
                      tracr, cumufx,   &
                      dpinit, dpfinl,  &
                      TimingBarriers )

  use module_control  ,only: nvlp1,npp,nip,dt,ntra,ntrb
  use fimnamelist     ,only: nvl,PrintIpnDiag
  use module_constants,only: nprox,prox,area,rarea,perm
  use module_fct3d

  implicit none
! External variables:
  integer,intent (IN)    :: its			! model time step
! integer,intent (IN)    :: ntr			! number of tracer fields
  logical,intent (IN)    :: TimingBarriers	! measure task skew when .true.
!sms$distribute (dh,1) begin
  real :: col_xpand(nip)
!sms$distribute end
!sms$distribute (dh,2) begin
  real   ,intent (INOUT) :: tracr (nvl,nip,ntra+ntrb) ! tracer
  real   ,intent (IN)    :: dpinit(nvl,nip)           ! init'l lyr thknss
  real   ,intent (IN)    :: dpfinl(nvl,nip)           ! final lyr thknss
  real :: vertfx(nvl,nip)
  real :: field(nvl,nip)
!sms$distribute end
!sms$distribute (dh,3) begin
  real   ,intent (IN)    :: cumufx(nvl,npp,nip)       ! time-integr. mass flx
!sms$distribute end
  integer    :: k		! layer index
  integer    :: ico		! Index for icos point number
  integer    :: edg		! icos edge number index
  integer    :: type		! tracer index
  logical    :: vrbos
  character  :: string*32
  real       :: hordiv(nvl),vertdv(nvl)
  real*8  :: tfct3d, valmin, valmax
  integer :: ret

!sms$compare_var(tracr  , "transp3d.F90 - tracr1   ")
!sms$compare_var(cumufx , "transp3d.F90 - cumufx1  ")

  ret = gptlstart ('transp2')

! --- compute the various terms in the continuity equation integrated
! --- over time interval since last call to -transp0-
! --- the continuity eqn is split into horiz. and vert. terms as follows:
! ---        (dpfinl-dpinit) + hordiv + verdiv = 0

  if (TimingBarriers) then
    ret = gptlstart ('transp2_barrier')
!SMS$BARRIER
    ret = gptlstop ('transp2_barrier')
  end if
  ret = gptlstart ('transp2_exchange')
!SMS$EXCHANGE(cumufx,dpinit,dpfinl)
  ret = gptlstop ('transp2_exchange')

!SMS$PARALLEL (dh,ico) BEGIN
!SMS$HALO_COMP(<1,1>) BEGIN
  do ico=1,nip				! horizontal loop
   vrbos=ico.eq.PrintIpnDiag
   hordiv   (:    )=0.
   col_xpand(  ico)=0.
   vertfx   (:,ico)=0.

   do edg=1,nprox(ico)			! loop through edges
    do k=1,nvl				! loop through layers
     hordiv(k)=hordiv(k)+cumufx(k,edg,ico)
    end do
   end do

   do k=1,nvl				! loop through layers
    col_xpand(ico)=col_xpand(ico)+hordiv(k)
    vertdv(k)=(dpinit(k,ico)-dpfinl(k,ico))-hordiv(k)*rarea(ico)
    if (k.eq.1) then
     vertfx(k,ico)=                vertdv(k)
    else
     vertfx(k,ico)=vertfx(k-1,ico)+vertdv(k)
    end if

    if (vrbos) then
    write (*,'(i7,i3,a,3es12.4)') perm(ico),k,			&
    ' transp2 hordiv,col_xpand,vertfx:',hordiv(k)*rarea(ico),	&
      col_xpand(ico)*rarea(ico),vertfx(k,ico)
    call flush(6)
    end if

   end do
  end do				! horizontal loop
!SMS$HALO_COMP END
!SMS$PARALLEL END

#ifdef DEBUGPRINT
  call findmxmn1(col_xpand,nip,'col_xpand')
  do k=1,nvl,7
   write (string,'(a,i3,a)') 'lyr',k,' vertfx'
   call findmxmn2(vertfx,nvl,nip,k,string)
  end do
  print *
  call flush(6)
#endif

! --- having determined the vertical flux term in the time-integrated
! --- continuity eqn, we can now perform the actual tracer transport

  if (ntrb.gt.0) then

   do type=1,ntrb	! loop over class B tracers
!SMS$PARALLEL (dh,ico) BEGIN
    do ico=1,nip
     do k=1,nvl
      field(k,ico)=tracr(k,ico,ntra+type)
     end do
    end do
!SMS$PARALLEL END

    if (PrintIpnDiag > 0) then
      write (string,'(a,i6,a,i2,a)') '(transp2) step',its,'  trcr',type,' in'
      call stencl(field,nvl,1.,trim(string))
    end if

   end do		! loop over tracers

   ret = gptlstart ('fct3d')

   call fct3d(tracr,ntra+ntrb,ntra+1,ntra+ntrb,cumufx,vertfx,		&
     area,rarea,dpinit,dpfinl,.false.)

   ret = gptlstop ('fct3d')
   ret = gptlget_wallclock ('fct3d', 0, tfct3d)  ! The "0" is thread number
   valmin = tfct3d*1.e3
   valmax = valmin
!SMS$REDUCE(valmin,MIN)
!SMS$REDUCE(valmax,MAX)
   print '(a,i2,a,2(i7,a))','time spent in subr fct3d (',ntrb,		&
    ' tracers)',nint(valmin),' -',nint(valmax),' msec'

   do type=1,ntrb	! loop over class B tracers
!SMS$PARALLEL (dh,ico) BEGIN
    do ico=1,nip
     do k=1,nvl
      field(k,ico)=tracr(k,ico,ntra+type)
     end do
    end do
!SMS$PARALLEL END

    if (PrintIpnDiag > 0) then
      write (string,'(a,i6,a,i2,a)') '(transp2) step',its,'  trcr',type,' out'
      call stencl(field,nvl,1.,trim(string))
    end if
   end do		! loop over tracers

  end if		!  ntrb > 0

!sms$compare_var(tracr  , "transp3d.F90 - tracr2   ")
!sms$compare_var(cumufx , "transp3d.F90 - cumufx2  ")

  ret = gptlstop ('transp2')

  return
  end subroutine transp2
end module module_transp3d
