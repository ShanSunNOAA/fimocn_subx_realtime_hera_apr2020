module module_outqv_mn_theta_abs
  use module_control, only: nip

  implicit none

#include <gptl.inc>

contains
!
!$$$   SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    OUTQV_MN_theta_abs    PRINT MEAN abs VALUE FOR EACH LAYER FROM ICOS-ARRAY
!                                      for isentropic or non-isentropic 3-d points. 
!   PRGMMR:  JIN, ADAPTED FROM OUTQV ORIGINALLY BY S.BENJAMIN  DATE: 07-06-20
!
! ABSTRACT:  PRINT MEAN LAYER VALUE OF 2-D ARRAY
!
! PROGRAM HISTORY LOG:
!    2007/06        J. Lee          adapted codes from outqv.F90
!
! USAGE:   call outqv_mn (qva, nlev, lat, lon, xlim_theta, factor, string)
!
!   INPUT ARGUMENT LIST:
!     QVA      - REAL     2-D ICOSAHEDRAL TRACER ARRAY
!     nlev     - INTEGER  NO. OF POINTS IN VERTICAL DIRECTION
!
! OUTPUT ARGUMENT LIST: none
!
! REMARKS: TODO: Rewrite this routine similarly to what was done for
!          outqv.F90. Vectorization and other optimizations should
!          result in around 10X speedup. This routine is currently
!          not called often enough to warrant rewrite though

  subroutine outqv_mn_theta_abs (qva, tr,ntra,thetac, nlev, nlev2, lat, lon, xlim_theta, factor, string)
    integer,intent(in) :: nlev,nlev2        ! number of vertical levels
    integer,intent(in) :: ntra              ! number of A-type tracers
    character(len=*), intent(in) :: string  ! string to be printed
!sms$distribute (dh,1) begin
    real, intent(in) :: lat(nip),lon(nip)
!sms$distribute end
!sms$distribute (dh,2) begin
    real, intent(in) :: qva(nlev,nip)
    real, intent(in) :: tr(nlev2,nip,ntra)
!sms$distribute end
    real, intent(in) :: factor
    real, intent(in) :: xlim_theta
    real, intent(in) :: thetac(nlev2)

    REAL*8 xsum, sum1
    REAL qvamn, qvamn1
    INTEGER IPN,IVL, isum, isum1
    integer :: ret       ! return code from gptl timing routines

    ret = gptlstart ('outqv_mn_theta_abs')
    write(6,*) string
    write (6,118)xlim_theta,factor
118 format ('  DTheta-limit=',f6.1, '  Scaling-factor=',G10.2)

!   do ivl=1,nlev2
!   write (6,*)' theta/thetac',tr(ivl,1,1),thetac(ivl)
!   end do

!SMS$PARALLEL (dh,ipn) BEGIN
    DO IVL=2,NLEV2
      xsum=0.d0
      sum1=0.d0
      isum = 0
      isum1 = 0
      DO IPN=1,NIP
        if (max( abs(tr(ivl,ipn,1)-thetac(ivl)), &
          abs(tr(ivl-1,ipn,1)-thetac(ivl-1))).lt.xlim_theta) then

          xsum = xsum + abs(QVA(ivl,ipn) )
          isum = isum + 1 
        else 
          sum1 = sum1 + abs( QVA(ivl,ipn) )
          isum1= isum1+ 1
        end if
      ENDDO
!SMS$reduce(isum,isum1,xsum,sum1,SUM)
      qvamn=xsum  /max(1.,float(isum))
      qvamn1=sum1/max(1.,float(isum1))
!
      write (6,120)IVL,xlim_theta,isum,qvamn*factor,isum1,qvamn1*factor
120   format ('  K=',i3,'  theta LE/GT ',f6.2, 2(' npts=',i8,'   MeanVal=',f12.4 ) )
    ENDDO
!SMS$PARALLEL END
    ret = gptlstop ('outqv_mn_theta_abs')

    RETURN
  end subroutine outqv_mn_theta_abs
end module module_outqv_mn_theta_abs
