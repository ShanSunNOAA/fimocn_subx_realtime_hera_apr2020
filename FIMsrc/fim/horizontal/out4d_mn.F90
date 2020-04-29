module module_out4d_mn
  use module_control, only: nip

  implicit none

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
! USAGE:   CALL OUTQV_MN(QVA,NVL,NIP,ITS)
!
!   INPUT ARGUMENT LIST:
!     QVA      - REAL     2-D ICOSAHEDRAL TRACER ARRAY
!     NVL      - INTEGER  NO. OF POINTS IN VERTICAL DIRECTION
!
! OUTPUT ARGUMENT LIST: none
!
! REMARKS: NONE
!

  SUBROUTINE OUT4D_MN (QVA, NLEV, NTR, factor, n, string)
    INTEGER,intent(IN) ::  NLEV,NTR,n
    character(len=*), intent(in) :: string
!sms$distribute (dh,2) begin
    real, intent(in) :: qva(nlev,nip,ntr)
!sms$distribute end
    real, intent(in) :: factor
    REAL*8 sum
    REAL qvamn
    INTEGER IPN,IVL

    write(6,*)string
!SMS$PARALLEL (dh,ipn) BEGIN
    DO IVL=1,NLEV
      sum=0.d0
      DO IPN=1,NIP
        sum = sum + QVA(ivl,ipn,n)
      ENDDO
      qvamn=sum/float(nip)
!SMS$reduce(qvamn,SUM)
!
      write (6,120)IVL,qvamn*factor
120   format ('  K=',i3,'   MeanVal=',f12.4 )
    ENDDO

!SMS$PARALLEL END
    RETURN
  end subroutine out4d_mn
end module module_out4d_mn
