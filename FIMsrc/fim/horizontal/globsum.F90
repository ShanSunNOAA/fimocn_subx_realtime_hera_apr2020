module module_globsum
  use module_control,   only: nvlp1, nip, dt, ntra, ntrb
  use fimnamelist,      only: nvl
  use module_constants, only: grvity, area, cp, rd
  use global_bounds,    only: ims, ime, ips, ipe

  implicit none

  real    :: qmstrold   ! save previous value
  real    :: qmstrcold  ! save previous value
  real    :: qmstrnold  ! save previous value
  logical :: qdtr_set = .false.

!JR These things were moved from dyn_run to here for restart capability.

  real :: qmass
  real :: qmsqv
  real :: qmsqw
  real :: qmso3
  real :: qmste
  real :: qmstr  = 0.
  real :: qmstrn = 0.
  real :: qmstrc = 0.
  real :: qdtr
  real :: qdtrn
  real :: qdtrc

contains

!*********************************************************************
!     globsum
!       Calculate global sums of various useful quantities
!       Alexander E. MacDonald  11/14/2005
!       J. Lee                  01/04/2006
!       J. Rosinski             10/03/2011
!         Modified for restart capability
!*********************************************************************

  subroutine globsum (its, dp3d, tr, rn2d, rc2d, pr3d, ex3d, qf2d, qtrcr)
    integer, intent(in) :: its
    real, intent(in) :: dp3d(nvl,ims:ime)
    real, intent(in) :: tr(nvl,ims:ime,ntra+ntrb)
    real, intent(in) :: rn2d(ims:ime)         ! from PHY via CPL
    real, intent(in) :: rc2d(ims:ime)         ! from PHY via CPL
    real, intent(in) :: pr3d(nvlp1,ims:ime)
    real, intent(in) :: ex3d(nvlp1,ims:ime)
    real, intent(in) :: qf2d(ims:ime)         ! from PHY via CPL
    real, intent(out) :: qtrcr(ntra+ntrb)

! Local variables
    integer :: ipn  ! Index for icos point number
    integer :: k    ! Index vertical level
    integer :: t    ! Index for tracer type
    real    :: den
    real*8  :: summ,smqv,smqw,smo3,sumr,sume,sumrc,sumrn,sumtr(ntra+ntrb)


!  Initialize all sums as zero
    summ  = 0.
    smqv  = 0.
    smqw  = 0.
    smo3  = 0.
    sumr  = 0.
    sume  = 0.
    sumrc = 0.
    sumrn = 0.
    sumtr(:) = 0.

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (k,t,den) &
!$OMP             REDUCTION(+: summ,smqv,smqw,smo3,sumr,sumrc,sumrn,sume,sumtr) &
!$OMP             SCHEDULE (static)
    do ipn=ips,ipe
      do k=1,nvl
        summ = summ + area(ipn)*dp3d(k,ipn)                  ! integrated mass
        smqv = smqv + area(ipn)*dp3d(k,ipn)*tr(k,ipn,2)      ! integrated water vapor
        smqw = smqw + area(ipn)*dp3d(k,ipn)*tr(k,ipn,3)      ! integrated condensate
        smo3 = smo3 + area(ipn)*dp3d(k,ipn)*tr(k,ipn,4)      ! integrated ozone
      end do

      sumr  = sumr + area(ipn)*rn2d(ipn)                     ! integrated precipitation
      sumrc = sumrc + area(ipn)*rc2d(ipn)                    ! integrated sub-gridscale (conv) precip
      sumrn = sumrn + area(ipn)*(rn2d(ipn)-rc2d(ipn))        ! integrated resolved precip
      den   = pr3d(1,ipn)/((rd/cp)*ex3d(1,ipn)*tr(1,ipn,1))
      sume  = sume + area(ipn)*den*qf2d(ipn)*dt/(1.25*2.5e6) ! integrated evaporation

      do t=1,ntra+ntrb
        do k=1,nvl
          sumtr(t) = sumtr(t) + area(ipn)*dp3d(k,ipn)*tr(k,ipn,t) ! tracer
        end do
      end do
    end do
!$OMP END PARALLEL DO
!sms$ignore end

    summ = summ/grvity
    smqv = smqv/grvity
    smqw = smqw/grvity
    smo3 = smo3/grvity

!sms$reduce(summ,smqv,smqw,smo3,sumr,sume,sumrc,sumrn,SUM)

    qmass = summ 
    qmsqv = smqv
    qmsqw = smqw
    qmso3 = smo3
    qmste = sume
    qtrcr(:) = sumtr(:)
! save previous values
    qmstrold  = qmstr
    qmstrcold = qmstrc
    qmstrnold = qmstrn
! store new values
    qmstr  = sumr
    qmstrc = sumrc
    qmstrn = sumrn

    if (qdtr_set) then
      qdtr  = qmstr  - qmstrold
      qdtrn = qmstrn - qmstrnold
      qdtrc = qmstrc - qmstrcold
    else
      qdtr  = 0.
      qdtrn = 0.
      qdtrc = 0.
      qdtr_set = .true.
    end if

    return
  end subroutine globsum
end module module_globsum
