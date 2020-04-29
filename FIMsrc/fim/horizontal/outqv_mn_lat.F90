module module_outqv_mn_lat
  use global_bounds,  only: ims, ime, ips, ipe

  implicit none

#include <gptl.inc>

contains
!
!$$$   SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    OUTQV_MN    PRINT MEAN VALUE FOR EACH LAYER FROM ICOS-ARRAY
!   PRGMMR:  JIN, ADAPTED FROM OUTQV ORIGINALLY BY S.BENJAMIN  DATE: 07-06-20
!
! ABSTRACT:  PRINT MEAN LAYER VALUE OF 2-D ARRAY
!
! PROGRAM HISTORY LOG:
!    2007/06        J. Lee          adapted codes from outqv.F90
!    2013/05        J. Rosinski     vectorized/optimized (10X speedup)
!
! USAGE:   call outqv_mn_lat (qva, nlev, lat, lon, xlim_lat, factor, string)
!
!   INPUT ARGUMENT LIST:
!     QVA      - REAL     2-D ICOSAHEDRAL TRACER ARRAY
!     NLEV     - INTEGER  NO. OF POINTS IN VERTICAL DIRECTION
!
! OUTPUT ARGUMENT LIST: none
!
! REMARKS: NONE
!

  subroutine outqv_mn_lat (qva, nlev, lat, lon, xlim_lat, factor, string)
    integer,intent(in) :: nlev              ! number of vertical levels in qva
    real, intent(in)   :: qva(nlev,ims:ime)
    real, intent(in)   :: lat(ims:ime)
    real, intent(in)   :: lon(ims:ime)
    real, intent(in)   :: xlim_lat
    real, intent(in)   :: factor
    character(len=*), intent(in) :: string  ! string to be printed

    real*8  :: sum(nlev), sum1(nlev)
    real    :: qvamn(nlev), qvamn1(nlev)
    integer :: ipn, k
    integer :: isum(nlev), isum1(nlev)

    integer :: ret       ! return code from gptl timing routines

    ret = gptlstart ('outqv_mn_lat')

    sum(:)   = 0.d0
    sum1(:)  = 0.d0
    isum(:)  = 0
    isum1(:) = 0
!JR The following OMP directive works, but on CPU is usually somewhat slower than not threading
!JR May want to enable for Xeon Phi
!!!$OMP PARALLEL DO PRIVATE(K) REDUCTION(+:SUM,SUM1,ISUM,ISUM1)
!sms$ignore begin
    do ipn=ips,ipe
!DIR$ VECTOR
      do k=1,nlev
        if (abs(lat(ipn)) > xlim_lat) then
          sum(k) = sum(k) + qva(k,ipn)
          isum(k) = isum(k) + 1 
        else 
          sum1(k) = sum1(k) + qva(k,ipn)
          isum1(k) = isum1(k) + 1
        end if
      end do
    end do
!sms$ignore end
!SMS$REDUCE(ISUM,ISUM1,SUM,SUM1,SUM)
    do k=1,nlev
      qvamn(k)  = sum(k)/isum(k)
      qvamn1(k) = sum1(k)/isum1(k)
    end do

    write(6,*) string
    write(6,118)xlim_lat,factor
118 format ('  Latitude-limit=',f6.1,'  Scaling-factor=',G10.2)
    do k=1,nlev
      write (6,120)k,xlim_lat,isum(k),qvamn(k)*factor,isum1(k),qvamn1(k)*factor
120   format ('  K=',i3,'  lat GT/LE ',f6.2, 2(' npts=',i8,'   MeanVal=',f12.4))
    end do

    ret = gptlstop ('outqv_mn_lat')

    return
  end subroutine outqv_mn_lat
end module module_outqv_mn_lat
