!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write 64-bit real array after gathering from other MPI tasks.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine writearr64_r (arr, dim2siz, unitno)
  use module_control, only: nip
  
  implicit none

  integer, intent(in) :: dim2siz  ! size of 2nd dimension of arr
  integer, intent(in) :: unitno   ! unit number to write to
  integer             :: ierr

!SMS$DISTRIBUTE(dh,1) BEGIN
  real*8, intent(in) :: arr(nip,dim2siz)
!SMS$DISTRIBUTE END

  write (unitno,iostat=ierr) arr
  if (ierr.ne.0) then
    write(*,'(a,i0,a)') 'writearr64_r: error writing to unit ',unitno,' Stopping'
    call flush(6)
    stop
  endif

end subroutine writearr64_r
