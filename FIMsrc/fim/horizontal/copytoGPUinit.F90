! Routine to copy constants from the CPU to the GPU once and for all.
! This is a no-op routine for the CPU
! Author: Jacques Middlecoff
! Date:   April 2012 

subroutine copytoGPUinit(ims,ime,ips,ipe,ihe,nvl,nvlp1,npp,ntra,PrintIpnDiag,&
                                nprox,proxs,prox,cs,sn,nedge,permedge,perm )

implicit none
integer,intent(in) :: ims, ime, ips, ipe, ihe        ! bounds
integer,intent(in) :: nvl, nvlp1, npp, ntra, PrintIpnDiag
integer,intent(in) :: nprox(ims:ime), proxs(npp,ims:ime), prox(npp,ims:ime)
real   ,intent(in) :: cs(4,npp,ims:ime),sn(4,npp,ims:ime)
integer,intent(in) :: nedge(ims:ime), permedge(npp,ims:ime), perm(ims:ime)

!ACC$DATA(<ims,ime,ips,ipe,ihe,nvl,nvlp1,npp,ntra,PrintIpnDiag,   &
!ACC$>     nprox,proxs,prox,cs,sn,nedge,permedge,perm:in,global>)

return
end subroutine copytoGPUinit
