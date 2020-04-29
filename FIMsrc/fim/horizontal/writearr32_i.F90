!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write 32-bit integer array after gathering from other MPI tasks.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine writearr32_i (arr, dim1siz, unitno)
  use module_control, only: nip
  
  implicit none

  integer, intent(in) :: dim1siz  ! size of 1st dimension of arr
  integer, intent(in) :: unitno   ! unit number to write to
  integer             :: ierr

!SMS$DISTRIBUTE(dh,2) BEGIN
  integer, intent(in) :: arr(dim1siz,nip)
!SMS$DISTRIBUTE END

  write (unitno,iostat=ierr) arr
  if (ierr.ne.0) then
    write(*,'(a,i0,a)') 'writearr32_i: error writing to unit ',unitno,' Stopping'
    call flush(6)
    stop
  endif

end subroutine writearr32_i
