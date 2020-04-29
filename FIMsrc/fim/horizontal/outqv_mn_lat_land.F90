module module_outqv_mn_lat_land
  use module_control, only: nip

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
!
! USAGE:   call outqv_mn_lat_land (qva, nlev, lat, lon, xlim_lat, slmsk2d, iland, factor, string)
!
!   INPUT ARGUMENT LIST:
!     QVA      - REAL     2-D ICOSAHEDRAL TRACER ARRAY
!     NLEV     - INTEGER  NO. OF POINTS IN VERTICAL DIRECTION
!
! OUTPUT ARGUMENT LIST: none
!
! REMARKS: TODO: Rewrite this routine similarly to what was done for
!          outqv.F90. Vectorization and other optimizations should
!          result in around 10X speedup. This routine is currently
!          not called often enough to warrant rewrite though

subroutine outqv_mn_lat_land (qva, nlev, lat, lon, xlim_lat, slmsk2d, iland, factor, string)
    integer, intent(in) :: nlev             ! number of vertical levels in qva
    character(len=*), intent(in) :: string  ! string to be printed
!sms$distribute (dh,1) begin
    real, intent(in) :: lat(nip)
    real, intent(in) :: lon(nip)
    real, intent(in) :: slmsk2d(nip)
!sms$distribute end
!sms$distribute (dh,2) begin
    real, intent(in) :: qva(nlev,nip)
!sms$distribute end
    real, intent(in) :: factor
    real, intent(in) :: xlim_lat 
    integer, intent(in) :: iland

    REAL*8 sum, sum1
    REAL qvamn, qvamn1
    INTEGER IPN,IVL, isum, isum1
    integer :: ret       ! return code from gptl timing routines

    ret = gptlstart ('outqv_mn_lat_land')
    write (6,*) string
    write (6,118) xlim_lat, factor
118 format ('  Latitude-limit=',f6.1, '  Scaling-factor=',G10.2 )

!SMS$PARALLEL (dh,ipn) BEGIN
    DO IVL=1,NLEV
      sum=0.d0
      sum1=0.d0
      isum = 0
      isum1 = 0
      DO IPN=1,NIP
        if (abs(lat(ipn)).gt.xlim_lat.and.iland.eq.int(slmsk2d(ipn)+0.5)) then
          sum = sum + QVA(ivl,ipn)
          isum = isum + 1 
        else if (abs(lat(ipn)).le.xlim_lat.and.iland.eq.int(slmsk2d(ipn)+0.5)) then
          sum1 = sum1 + QVA(ivl,ipn)
          isum1= isum1+ 1
        end if
      ENDDO
!SMS$reduce(isum,isum1,sum,sum1,SUM)
      qvamn=sum  /float(isum)
      qvamn1=sum1/float(isum1)
!
      write (6,120)IVL,xlim_lat,isum,qvamn*factor,isum1,qvamn1*factor
120   format ('  K=',i3, '  lat GT/LE ',f6.2, 2(' npts=',i8, '   MeanVal=',f12.4))
    ENDDO
!SMS$PARALLEL END
    ret = gptlstop ('outqv_mn_lat_land')
    RETURN
  end subroutine outqv_mn_lat_land
end module module_outqv_mn_lat_land
