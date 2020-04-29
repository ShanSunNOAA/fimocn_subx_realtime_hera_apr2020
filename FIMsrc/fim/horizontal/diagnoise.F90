module module_diagnoise
  use findmaxmin1
  use module_control,      only: nip
  use fimnamelist,         only: nvl, TimingBarriers
  use module_constants,    only: nedge, prox, nprox, deg_lat, deg_lon
  use module_outqv_mn,     only: outqv_mn
  use module_outqv_mn_lat, only: outqv_mn_lat
  use global_bounds,       only: ims, ime, ips, ipe

  implicit none

#include <gptl.inc>

contains

!*********************************************************************
!       diagnoise
!	A routine to diagnose external gravity wave noise in FIM
!       (expressed as rms of 2nd time derivative of surface pressure)
!	R. Bleck                July 2006
!*********************************************************************

  subroutine diagnoise (its, ptdcy)
! Arguments
    integer,intent (IN) :: its		    ! model time step
!sms$distribute (dh,1) begin
    real,   intent (IN) :: ptdcy(nip,2) ! d(surf.prs.)/dt at 2 consecutive time levels
!sms$distribute end

! Local variables
    real    :: work(ims:ime)
    real    :: deriv2		! rms of 2nd derivative of surface pressure
    real    :: derivma		! mean abs of 2nd derivative of surface pressure
    real    :: laplac_mean      ! mean absolute value of Laplacian of surface presure tendency
    real    :: scale = 5.	! arbitrary scale for plotting noise parameter
    integer :: ivl		! layer index
    integer :: ipn		! index for icosahedral grid
    real    :: valmin,valmax
    real    :: sum1
    integer :: ipx              ! neighbor across joint edge
    integer :: edg              ! joint edge index as seen by neighbor
    integer :: ret              ! return code from gptl routines

    if (its < 2) return

    ret = gptlstart ('diagnoise')

    deriv2 = 0.
    derivma = 0.

!sms$ignore begin
!$OMP PARALLEL DO REDUCTION(+: deriv2,derivma) SCHEDULE (static)
    do ipn=ips,ipe
      deriv2 = deriv2 + (ptdcy(ipn,1)-ptdcy(ipn,2))**2
      derivma = derivma + abs (ptdcy(ipn,1)-ptdcy(ipn,2))
    end do
!$OMP END PARALLEL DO
!sms$ignore end

!sms$reduce(deriv2,derivma,SUM)

! --- plot time series of noise index in stdout
! --- (type 'grep =+= stdout' to display the time series)
    deriv2 = sqrt(deriv2/nip)
    derivma = derivma/nip
    call linout(scale * deriv2,'x',its)

#ifdef DEBUGPRINT
    call findmxmn1 (ptdcy(:,1), nip, 'ptdcy(:,1)')
    call findmxmn1 (ptdcy(:,2), nip, 'ptdcy(:,2)')
#endif

    laplac_mean = 0.

    if (TimingBarriers) then
      ret = gptlstart ('diagnoise_barrier')
!SMS$BARRIER
      ret = gptlstop ('diagnoise_barrier')
    end if
    ret = gptlstart ('diagnoise_exchange')
!SMS$EXCHANGE(ptdcy)
    ret = gptlstop ('diagnoise_exchange')

    ! calculate Laplacian of sfc pressure tendency
!sms$ignore begin
!$OMP PARALLEL DO REDUCTION(+: laplac_mean) PRIVATE (sum1,edg,ipx)
    do ipn=ips,ipe
      sum1 = 0.
      do edg=1,nedge(ipn)
        ipx = prox(edg,ipn)
        sum1 = sum1 + ptdcy(ipx,1)
      end do
      work(ipn) = ptdcy(ipn,1) - sum1/nedge(ipn)
      laplac_mean = laplac_mean + abs(work(ipn))
    end do
!$OMP END PARALLEL DO
!sms$ignore end

!sms$reduce(laplac_mean,SUM)

    laplac_mean = (laplac_mean/nip)
    write (6,'(a,i10,3f8.3)')'rms-d(psfc)^2/dt^2 & deriv-mean abs & laplac-mean abs =',	&
         its,deriv2,derivma,laplac_mean
#ifdef DEBUGPRINT
    call findmxmn1(work,nip,'laplac-ptdcy')
#endif

    call outqv_mn (work,     1,                   1.,      'Laplac-ptdcy')
    call outqv_mn_lat (work, 1, deg_lat, deg_lon, 30., 1., 'Laplac-ptdcy')

    ret = gptlstop ('diagnoise')
 
    return
  end subroutine diagnoise

  subroutine linout(value,char,labl)
!
! --- print single characters in a manner mimicking a curve plot in x,y space
! --- abscissa: down the page; ordinate: across the page
!
    real,intent(IN)        :: value
    character*1,intent(IN) :: char
    integer,intent(IN)     :: labl

    integer,parameter      :: length=72
    character*1 line(length)
    integer l,n
 
! --- replace n-th element of array 'line' by character 'char', where
! --- n = 'value' modulo 'length'
! --- initialize 'line' by blanks before adding 'char'
! --- output 'line' after adding 'char'
! --- labl       -- abscissa value (integer), added to output line
!
    line(:) = ' '
    if (value.gt.0.) then
      n = int(mod(value + float(length-1),float(length))) + 1
      line(n) = char
    end if
    write (*,'(''=+='',i6,80a1)') labl,(line(l),l=1,length)
    return
  end subroutine linout
end module module_diagnoise
