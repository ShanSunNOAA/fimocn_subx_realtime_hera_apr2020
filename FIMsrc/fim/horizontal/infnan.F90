module infnan

  implicit none

  private
  public :: inf, nan, negint, isinf, isnan

! TODO: add big endian ifdef for e.g. IBM
! TODO: Should NaN be signaling or non-signaling?

  integer, parameter :: negint = -999       ! Bad integer value
#ifdef NAG
  integer            :: iinf,inan
  real               :: inf,nan

  equivalence (iinf,inf)
  equivalence (inan,nan)

  data iinf/z'7f800000'/    ! Infinity
  data inan/z'ffc00000'/    ! NaN
#else
  real,    parameter :: inf = z'7f800000'   ! Infinity
  real,    parameter :: nan = z'ffc00000'   ! NaN
#endif

contains

! These functions are experimental! If they work they can provide a nice means to test
! if something has become an inf or a nan

  logical function isnan (x)
    real :: x

    isnan = (x < 0. .eqv. x >= 0.)
    return

  end function isnan

  logical function isinf (x)
    real :: x

    isinf = abs(x) > huge (x)
    return

  end function isinf

end module infnan
