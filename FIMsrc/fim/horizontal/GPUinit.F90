! Routine to initialize the GPU
! Author:  Jacques Middlecoff
! Date:  September 2010 
! For Fortran this routine does nothing except return error=0
! See horizontal/cuda/GPUinit.cu for the real routine.

subroutine GPUinit(npes,me,error)
integer, intent(IN ) :: npes
integer, intent(IN ) :: me
integer, intent(OUT) :: error

error = 0

return
end subroutine GPUinit
