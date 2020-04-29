module module_outqv
  use global_bounds, only: ims, ime, ips, ipe
  use module_constants,only: perm

  implicit none

#include <gptl.inc>

contains
!
!$$$   SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    OUTQV       PRINT MAX VALUE OF ARRAY
!   PRGMMR:  BENJAMIN, STAN ORG: ERL/PROFS      DATE: 93-01-18
!
! ABSTRACT:  PRINT MAX VALUE AND IPN,IVL OF 2-D ARRAY
!
! PROGRAM HISTORY LOG:
!    88/05/31       S. BENJAMIN     ORIGINAL VERSION
!    2006/02        J. Lee          convert from F77 to F90
!    2006/02        J. Middlecoff   converted to icos notation and parallelized
!    2013/05        J. Rosinski     vectorized/optimized (10X speedup)
!
! USAGE:   call outqv (qva, nlev, lat, lon, string)
!
!   INPUT ARGUMENT LIST:
!     QVA      - REAL     2-D ARRAY
!     nlev     - INTEGER  NO. OF POINTS IN VERTICAL DIRECTION
!     lat      - REAL     1-D array
!     lon      - REAL     1-D array
!     string   - character string to be printed
!
! OUTPUT ARGUMENT LIST: none
!
! REMARKS: NONE
!

  subroutine outqv (qva, nlev, lat, lon, factor, string)
    integer, intent(in) :: nlev             ! number of vertical levels in qva
    real, intent(in) :: qva(nlev,ims:ime)
    real, intent(in) :: lat(ims:ime)
    real, intent(in) :: lon(ims:ime)
    real, intent(in) :: factor
    character(len=*), intent(in) :: string  ! string to be printed
! Local workspace
    real :: qvamax(nlev)      ! max values per level
    real :: my_qvamax(nlev)   ! local max values per level
    real :: latimax(nlev)     ! latitude of max values
    real :: lonimax(nlev)     ! longitude of max values
    integer :: imax(nlev)     ! index of max values
    integer :: glbmax(nlev)   ! permutated index of max values
    integer :: imaxk(1)       ! Fortran "maxloc" intrinsic returns array

    real :: qvamin(nlev)      ! min values per level
    real :: my_qvamin(nlev)   ! local min values per level
    real :: latimin(nlev)     ! latitude of min values
    real :: lonimin(nlev)     ! longitude of min values
    integer :: imin(nlev)     ! index of min values
    integer :: glbmin(nlev)   ! permutated index of min values
    integer :: imink(1)       ! Fortran "minloc" intrinsic returns array

    integer :: ipn, k         ! indices
    integer :: ret            ! return code from gptl timing routines

    ret = gptlstart ('outqv')
    qvamax(:) = -1.e30
    qvamin(:) = +1.e30
    glbmax(:) = -9999
    glbmin(:) = -9999

!sms$ignore begin
    do ipn=ips,ipe
!DIR$ VECTOR
      do k=1,nlev
        if (qva(k,ipn) >= qvamax(k)) then
          qvamax(k) = qva(k,ipn)
          imax(k) = ipn
        end if
        if (qva(k,ipn) <= qvamin(k)) then
          qvamin(k) = qva(k,ipn)
          imin(k) = ipn
        end if
      end do
    end do
!sms$ignore end
!JR Save off local max/min before SMS reduction to discover who owns the values
!JR and what the ipn index is
    my_qvamax(:) = qvamax(:)
    my_qvamin(:) = qvamin(:)
!SMS$REDUCE (QVAMAX,MAX)
!SMS$REDUCE (QVAMIN,MIN)

!JR Set indices to zero if my rank is not the owner. That way the owner
!JR of the global max/min can discover who he is.
!DIR$ VECTOR
    do k=1,nlev
      if (my_qvamax(k) /= qvamax(k)) then
        imax(k) = 0
      end if
      if (my_qvamin(k) /= qvamin(k)) then
        imin(k) = 0
      end if
    end do
!JR Reduction of indices on "MAX" guarantees the owner will have non-zero values
!SMS$REDUCE(IMAX,IMIN,MAX)

! exit if data is totally hosed
!JR Fortran standard weirdness requires output from maxloc be an array
    imaxk(:) = maxloc (qvamax)
    if (qvamax(imaxk(1)) >= 1.e30) then
      print *,'error outqv:  qvamax, k= ',qvamax(imaxk(1)), k
      stop
    end if

    imink(:) = minloc (qvamin)
    if (qvamin(imink(1)) <= -1.e30) then
      print *,'error outqv:  qvamin, k= ',qvamin(imink(1)), k
      stop
    end if
!JR Discover the lat/lon points of the owner of max/min values
    latimax(:) = -9999.
    lonimax(:) = -9999.
    latimin(:) = -9999.
    lonimin(:) = -9999.
!DIR$ VECTOR
    do k=1,nlev
      if (imax(k) >= ips .and. imax(k) <= ipe) then
        latimax(k) = lat(imax(k))
        lonimax(k) = lon(imax(k))
        glbmax(k) = perm(imax(k))
      end if
      if (imin(k) >= ips .and. imin(k) <= ipe) then
        latimin(k) = lat(imin(k))
        lonimin(k) = lon(imin(k))
        glbmin(k) = perm(imin(k))
      end if
      if (lonimax(k) > 180.) lonimax(k) = lonimax(k) - 360.
      if (lonimin(k) > 180.) lonimin(k) = lonimin(k) - 360.
    end do
!JR Since latimax, etc. were initialized to a negative number, reduction on "MAX"
!JR will set the value based on the owner of the actual max/min
!SMS$REDUCE(LATIMAX,LONIMAX,LATIMIN,LONIMIN,glbmax,glbmin,MAX)
    do k=1,nlev
      qvamax(k) = qvamax(k)*factor
      qvamin(k) = qvamin(k)*factor
    end do
    
    write(6,'(2a)') '(outqv) ',string
    print 100
100 format ('    IVL    IPN    LAT     LON    MAX VALUE       IPN    LAT     LON    MIN VALUE')
    do k=1,nlev
      print 101, k, glbmax(k), latimax(k), lonimax(k), qvamax(k), glbmin(k), latimin(k), lonimin(k), qvamin(k)
101   format (i6, i8, 2f8.1, 1pe12.3, i10, 2f8.1, 1pe12.3)
    end do
    ret = gptlstop ('outqv')

    return
  end subroutine outqv
end module module_outqv
