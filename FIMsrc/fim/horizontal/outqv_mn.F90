module module_outqv_mn
  use module_control, only: nip
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
! USAGE:   call outqv_mn (qva, nlev, factor, string)
!
!   INPUT ARGUMENT LIST:
!     QVA      - REAL     2-D ICOSAHEDRAL TRACER ARRAY
!     NLEV     - INTEGER  NO. OF POINTS IN VERTICAL DIRECTION
!     factor   - real
!     string   - character string to be printed
!
! OUTPUT ARGUMENT LIST: none
!
! REMARKS: Due to the SMS$REDUCE(SUM), this routine does NOT produce bitwise-identical answers
!          across arbitrary task counts.
!

  subroutine outqv_mn (qva, nlev, factor, string)
    integer, intent(in) :: nlev             ! number of vertical levels in qva
    real, intent(in)    :: qva(nlev,ims:ime)
    real, intent(in)    :: factor
    character(len=*), intent(in) :: string  ! string to be printed

    real*8  :: sum(nlev)
    real    :: qvamn(nlev)
    integer :: ipn,k
    integer :: ret       ! return code from gptl timing routines

    ret = gptlstart ('outqv_mn')

    sum(:) = 0.d0
!JR The following OMP directive works, but on CPU is usually somewhat slower than not threading
!JR May want to enable for Xeon Phi
!!!$OMP PARALLEL DO PRIVATE(K) REDUCTION(+:SUM)
!sms$ignore begin
    do ipn=ips,ipe
      do k=1,nlev
        sum(k) = sum(k) + qva(k,ipn)
      end do
    end do
!sms$ignore end

    do k=1,nlev
      qvamn(k) = sum(k)/nip
    end do
!SMS$REDUCE(QVAMN,SUM)
    write(6,'(2a)') '(outqv_mn) ',string
    do k=1,nlev
      write (6,120)k,qvamn(k)*factor
120   format ('  K=',i3,'   MeanVal=',f12.4 )
    end do

    ret = gptlstop ('outqv_mn')

    return
  end subroutine outqv_mn
end module module_outqv_mn
