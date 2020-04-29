subroutine GetLastError(tag)
!This routine returns the last error on the GPU.
!This routine does nothing on the CPU.
!See horizontal/cuda/GetLastError.cu for the real routine.
! Author:  Jacques Middlecoff
! Date:    August 2011

implicit none
integer,intent(IN) :: tag
return
end
