!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read 64-bit real array and distribute to other MPI tasks.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readarr64_r (arr, dim2siz, unitno)
  use module_control, only: nip

  implicit none

  integer, intent(in) :: dim2siz  ! size of 2nd dimension of arr
  integer, intent(in) :: unitno   ! unit number to read from
  integer             :: ierr

!SMS$DISTRIBUTE(dh,1) BEGIN
  real*8, intent(inout) :: arr(nip,dim2siz)
!SMS$DISTRIBUTE END

  read (unitno,iostat=ierr) arr
  if (ierr.ne.0) then
    write(*,'(a,i0,a)') 'readarr64_r: error reading from unit ',unitno,' Stopping'
    call flush(6)
    stop
  endif

end subroutine readarr64_r
