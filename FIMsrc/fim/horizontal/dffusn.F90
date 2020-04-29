module dffusn
  use module_control  ,only: nip
  use fimnamelist     ,only: nvl, PrintIpnDiag
  use module_constants,only: nprox, prox, rarea, sideln, perm
  use global_bounds,   only: ims, ime, ips, ipe

  implicit none

  private

#ifdef USE_INTERFACE_BLOCKS
  public :: dffusn_lev, dffusn_lyr

  interface dffusn_lev
     module procedure dffusn_lev, dffusn8_lev
  end interface
     
  interface dffusn_lyr
     module procedure dffusn_lyr, dffusn8_lyr
  end interface
#else
  public :: dffusn_lev, dffusn8_lev, dffusn_lyr, dffusn8_lyr
#endif
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nothing below needs to be processed by SMS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!*********************************************************************
!     dffusn_lev
!	Diffuse level variable (no thickness weighting)
!	S. Sun                           September 2009
!*********************************************************************

  subroutine dffusn_lev (fld, dfflen, kdim, k1, k2)

! Arguments
    integer,intent (IN)    :: kdim		! vert.dim. of fld and dfflen
    integer,intent (IN)    :: k1,k2		! operate on levels k1 ... k2
    real   ,intent (IN)    :: dfflen(kdim)	! diffusion length scale (m)
    real   ,intent (INOUT) :: fld(kdim,ims:ime)     ! field(s) to be diffused

! Local variables
    real    :: flxdv(kdim,ims:ime) ! line integral of dffus.flux across the edge

#include "dffusn_lev.dec"

!sms$ignore begin

#include "dffusn_lev.exe"

    return

!sms$ignore end

  end subroutine dffusn_lev

  subroutine dffusn8_lev (fld, dfflen, kdim, k1, k2)

! Arguments
    integer,intent (IN)    :: kdim		! vert.dim. of fld and dfflen
    integer,intent (IN)    :: k1,k2		! operate on levels k1 ... k2
    real   ,intent (IN)    :: dfflen(kdim)	! diffusion length scale (m)
    real*8 ,intent (INOUT) :: fld(kdim,ims:ime)     ! field(s) to be diffused

! Local variables
    real*8  :: flxdv(kdim,ims:ime) ! line integral of dffus.flux across the edge

#include "dffusn_lev.dec"

!sms$ignore begin  

#include "dffusn_lev.exe"

    return

!sms$ignore end

  end subroutine dffusn8_lev

!*********************************************************************
!     dffusn_lyr
!	Diffuse layer variable (thickness-weighted for conservation)
!	S. Sun                           September 2009
!*********************************************************************

  subroutine dffusn_lyr (fld, delp, dfflen)
! Arguments
    real   ,intent (IN)    :: dfflen              ! diffusion length scale (m)
    real   ,intent (INOUT) :: fld(nvl,ims:ime)    ! field to be diffused
    real   ,intent (IN)    :: delp(nvl,ims:ime)   ! lyr thknss, Pa

! Local variables
    real :: flxdv(nvl,ims:ime)	! line integral of dffus.flux across the edge

#include "dffusn_lyr.dec"

!sms$ignore begin

#include "dffusn_lyr.exe"

    return

!sms$ignore end

  end subroutine dffusn_lyr

  subroutine dffusn8_lyr (fld, delp, dfflen)
! Arguments
    real   ,intent (IN)    :: dfflen              ! diffusion length scale (m)
    real*8 ,intent (INOUT) :: fld(nvl,ims:ime)    ! field to be diffused
    real   ,intent (IN)    :: delp(nvl,ims:ime)   ! lyr thknss, Pa

! Local variables
    real*8 :: flxdv(nvl,ims:ime)	! line integral of dffus.flux across the edge

#include "dffusn_lyr.dec"

!sms$ignore begin

#include "dffusn_lyr.exe"

    return

!sms$ignore end

  end subroutine dffusn8_lyr

end module dffusn
