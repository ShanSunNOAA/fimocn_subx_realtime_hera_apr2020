! Routine to initialize the GPU
! Author:  Jacques Middlecoff
! Date:    August 2011
! For Fortran this routine does nothing except return status=0
! See horizontal/cuda/SMS_GetPointer.cu for the real routine.

subroutine SMS_GetPointer(varname, var, status)
INTEGER,INTENT(IN) :: VarName !The name of the variable to exchange(really character but F2C doesn't support character)
INTEGER,INTENT(IN) :: var     !The variable to exchange
INTEGER,INTENT(OUT):: status  !Non-zero means error

status = 0

return
end subroutine SMS_GetPointer
