subroutine GetIpnGlobal (ipn, ipnGlobal, DiagPrint)

! For an input ipn this routine returns the global ipn (ipnGLobal) and DiagPrint.
! DiagPrint=T means that the input ipn matches PrintIpnDiag from the namelist.

  use fimnamelist,only: PrintIpnDiag

  implicit none

  integer,intent(in)  :: ipn
  integer,intent(out) :: ipnGlobal
  logical,intent(out) :: DiagPrint

  ipnGlobal = ipn
  DiagPrint = ipnGlobal==PrintIpnDiag

end subroutine GetIpnGlobal
