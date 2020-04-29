!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read 32-bit real array and distribute to other MPI tasks.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readarr32_r (arr, dim1siz, unitno)
  use module_control, only: nip

  implicit none

  integer, intent(in) :: dim1siz  ! size of 1st dimension of arr
  integer, intent(in) :: unitno   ! unit number to read from
  integer             :: ierr

!SMS$DISTRIBUTE(dh,2) BEGIN
  real*4, intent(inout) :: arr(dim1siz,nip)
!SMS$DISTRIBUTE END

  read (unitno,iostat=ierr) arr
  if (ierr.ne.0) then
    write(*,'(a,i0,a)') 'readarr32_r: error reading from unit ',unitno,' Stopping'
    call flush(6)
    stop
  endif

end subroutine readarr32_r
