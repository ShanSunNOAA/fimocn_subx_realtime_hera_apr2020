!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read 64-bit real array and distribute to other MPI tasks.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readarr64t_r (arr, dim1siz, unitno)
  use module_control, only: nip

  implicit none

  integer, intent(in) :: dim1siz  ! size of 1st dimension of arr
  integer, intent(in) :: unitno   ! unit number to read from
  integer             :: ierr

!SMS$DISTRIBUTE(dh,2) BEGIN
  real*8, intent(inout) :: arr(dim1siz,nip)
!SMS$DISTRIBUTE END

  read (unitno,iostat=ierr) arr
  if (ierr.ne.0) then
    write(*,'(a,i0,a)') 'readarr64t_r: error reading from unit ',unitno,' Stopping'
    call flush(6)
    stop
  endif

end subroutine readarr64t_r
