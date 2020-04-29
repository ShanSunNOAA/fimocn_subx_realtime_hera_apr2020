!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write 32-bit real array after gathering from other MPI tasks.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine writearr32_r (arr, dim1siz, unitno)
  use module_control, only: nip
  
  implicit none

  integer, intent(in) :: dim1siz  ! size of 1st dimension of arr
  integer, intent(in) :: unitno   ! unit number to write to
  integer             :: ierr

!SMS$DISTRIBUTE(dh,2) BEGIN
  real*4, intent(in) :: arr(dim1siz,nip)
!SMS$DISTRIBUTE END

  write (unitno,iostat=ierr) arr
  if (ierr.ne.0) then
    write(*,'(a,i0,a)') 'writearr32_r: error writing to unit ',unitno,' Stopping'
    call flush(6)
    stop
  endif

end subroutine writearr32_r
