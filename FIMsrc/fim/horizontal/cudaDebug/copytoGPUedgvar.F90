! Routine to copy variables from the CPU to the GPU before running edgvar on the GPU
! This is a no-op routine for the CPU
! Author: Jacques Middlecoff
! Date:   April 2012 

subroutine copytoGPUedgvar(nvl,ims,ime,ntra,nvlp1,u_vel,v_vel,delp,pres,montg,tracr)

implicit none
integer,intent(in) ::       nvl  ,ims,ime,ntra,nvlp1
real   ,intent(in) :: u_vel(nvl  ,ims:ime)      ! west wind on s level
real   ,intent(in) :: v_vel(nvl  ,ims:ime)      ! south wind on s level
real   ,intent(in) :: delp (nvl  ,ims:ime)      ! layer thickness
real   ,intent(in) :: pres (nvlp1,ims:ime)      ! pressure on interfaces
real   ,intent(in) :: montg(nvl  ,ims:ime)      ! montgomery potential
real   ,intent(in) :: tracr(nvl  ,ims:ime,ntra) ! mass field tracers

!ACC$DATA(<u_vel,v_vel,delp,pres,montg,tracr:in,global>)

return
end subroutine copytoGPUedgvar
