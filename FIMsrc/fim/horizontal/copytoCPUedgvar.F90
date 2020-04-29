! Routine to copy variables from the GPU to the CPU after running edgvar on the CPU
! This is a no-op routine for the CPU
! Author: Jacques Middlecoff
! Date:   April 2012 

subroutine copytoCPUedgvar(nvl,npp,ims,ime,ntra,u_edg,v_edg,dp_edg,lp_edg,bnll_edg,trc_edg)

implicit none
integer,intent(in) ::          nvl,npp,ims,ime,ntra
real   ,intent(in) :: u_edg   (nvl,npp,ims:ime)	     ! west wind on edges
real   ,intent(in) :: v_edg   (nvl,npp,ims:ime)	     ! south wind on edges
real   ,intent(in) :: dp_edg  (nvl,npp,ims:ime)      ! layer thickness on edges
real   ,intent(in) :: lp_edg  (nvl,npp,ims:ime)      ! midlayer pressure on edges
real   ,intent(in) :: trc_edg (nvl,npp,ims:ime,ntra) ! tracers on edges
real   ,intent(in) :: bnll_edg(nvl,npp,ims:ime)	     ! bernoulli fct on edges

!ACC$DATA(<u_edg,v_edg,dp_edg,lp_edg,bnll_edg,trc_edg:out,extern>)

return
end subroutine copytoCPUedgvar
