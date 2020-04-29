!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines in this module are computational only. They contain no SMS directives
! except for IGNORE. Thus they are good candidates for migration to MIC or GPU.
!
! massadv_compute: contains all the computational routines for mass transport  
! calculations.
!
! Ning Wang, Dec 2015 - Original version
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module massadv_compute
  implicit none

  private
  public :: massflux_newdepth, update_pres

#include <gptl.inc>

contains

   subroutine massflux_newdepth(dp_edg,massfx,dp_tndcy,delp,psurf) 
    use module_control,   only: npp, nabl
    use fimnamelist,      only: nvl
    use module_variables, only: nf, of, vof, adbash1, adbash2, adbash3
    use module_constants, only: nprox, prox, rarea, nedge, permedge
    use global_bounds,    only: ims, ime, ips, ipe, ihe
    use Fluxlimiter,      only : vnorm

    real, intent(in)    :: dp_edg(nvl,npp,ims:ime)      ! fluid depth at edges
    real, intent(inout) :: massfx(nvl,npp,ims:ime,3)    ! N/sec
    real, intent(inout) :: dp_tndcy  (nvl,ims:ime,nabl) ! Pa/sec
    real, intent(inout) :: delp(nvl,ims:ime)
    real, intent(out)   :: psurf(ims:ime)               ! surface pressure

! Local variables
    integer :: ipn
    integer :: edg
    integer :: edgcount            ! count of icos edges
    integer :: k
    integer :: ret
    
    ret = gptlstart ('massflux_newdepth')
!$OMP PARALLEL DO PRIVATE (k,edgcount,edg) SCHEDULE (runtime)
    do ipn=ips,ihe
      psurf(ipn) = 0.
      do k=nvl,1,-1      ! top-down for psurf
        psurf(ipn) = psurf(ipn) + delp(k,ipn)
      end do

      do edgcount=1,nedge(ipn)
        edg = permedge(edgcount,ipn)
        do k=1,nvl
          massfx(k,edg,ipn,nf) = vnorm(k,edg,ipn)*dp_edg(k,edg,ipn)
        end do
      end do
    end do
!$OMP END PARALLEL DO

!!sms$ignore begin
!    do ipn=ips,ihe
!if (ipn > 0 .and. ipn < 10)then
!do k=nvl,1,-1      ! top-down for psurf
!print*, 'delp(k,ipn)', delp(k,ipn)
!enddo
!print*, 'ipn,  psurf(ipn)', ipn,  psurf(ipn)
!endif
!    end do 
!!sms$ignore begin

!$OMP PARALLEL DO PRIVATE (k,edg) SCHEDULE (runtime)
    do ipn=ips,ipe
      do k=1,nvl
        dp_tndcy(k,ipn,nf) = 0.
      end do
      do edg=1,nprox(ipn)
        do k=1,nvl            
          dp_tndcy(k,ipn,nf) = dp_tndcy(k,ipn,nf) + massfx(k,edg,ipn,nf)
        end do
      end do
      do k=1,nvl
        dp_tndcy(k,ipn,nf) = -dp_tndcy(k,ipn,nf)*rarea(ipn)
      end do

      do k=1,nvl
        ! advance delp to new time step
        delp(k,ipn) = delp(k,ipn) + adbash1*dp_tndcy(k,ipn,nf) + &
                                    adbash2*dp_tndcy(k,ipn,of) + &
                                    adbash3*dp_tndcy(k,ipn,vof)
        delp(k,ipn) = max(delp(k,ipn), 0.0)
      end do
    end do
!$OMP END PARALLEL DO
    ret = gptlstop ('massflux_newdepth')
  end subroutine massflux_newdepth

  subroutine update_pres (pres, delp, exner, dp_tndcy, omega, &
                      lp_edg, u_velo, v_velo)
    use module_control,   only: nvlp1, npp, nabl, dt
    use fimnamelist,      only: nvl, PrintIpnDiag
    use module_constants, only: rarea, nedge, permedge, sidevec_c,	&
                                p1000, cp, rd, perm
    use module_variables, only: nf, of, vof, adbash1, adbash2, adbash3
    use global_bounds,    only: ims, ime, ips, ipe, ihe

! Arguments
    real, intent(inout) :: pres(nvlp1,ims:ime)	       ! interface pressure, Pa
    real, intent(in)    :: delp(nvl,ims:ime)           ! layer thickness, Pa
!JR exner needs to be inout: why is only 1-nvl set here, not exner(nvlp1)???
    real, intent(inout) :: exner(nvlp1,ims:ime)	       ! exner function.
    real, intent(in)    :: dp_tndcy(nvl,ims:ime,nabl)  ! Pa/sec
    real, intent(out)   :: omega(nvl,ims:ime)	       ! n/sec
    real, intent(in)    :: lp_edg(nvl,npp,ims:ime)     ! mid-layer pressure on edges, Pa
    real, intent(in)    :: u_velo(nvl,ims:ime)         ! west wind, m/sec
    real, intent(in)    :: v_velo(nvl,ims:ime)         ! west wind, m/sec
      
! Local variables
    integer :: ipn
    integer :: k
    integer :: edgcount   ! count of icos edges
    integer :: edg        ! Index for icos edge number

    real    :: coltend    ! column pressure tendency
    real    :: lyrtend    ! layer pressure tendency
    real    :: dpdx,dpdy  ! pressure gradient in x,y direction
    real    :: old
    integer :: ret

    ret = gptlstart ('update_pres')
!$OMP PARALLEL DO PRIVATE (coltend,k,lyrtend,dpdx,dpdy,edgcount,edg,old) SCHEDULE (runtime)
    do ipn=ips,ipe
      coltend = 0.		! integrated mass flux convergence
      do k=nvl,1,-1		! loop through layers (top-down for p,omega)
        ! update pressure and Exner fct by vertically summing up thickness values
        pres(k,ipn) = pres(k+1,ipn) + delp(k,ipn)
        exner(k,ipn) = cp*(pres(k,ipn)/p1000)**(rd/cp)

        ! evaluate omega = dp/dt as
        ! (partial_p/partial_t) + (v_dot_grad_p) + (s_dot partial_p/ partial_s)

        lyrtend = (adbash1*dp_tndcy(k,ipn,nf) +			&
                   adbash2*dp_tndcy(k,ipn,of) +			&
                   adbash3*dp_tndcy(k,ipn,vof)) / dt
        omega(k,ipn) = coltend + 0.5*lyrtend          ! evaluate at mid depth
        coltend = coltend + lyrtend                   ! flux conv. intgral

        ! pressure gradient
        dpdx = 0.
        dpdy = 0.
        do edgcount=1,nedge(ipn)	! loop through edges
          edg = permedge(edgcount,ipn)
          dpdx = dpdx + lp_edg(k,edg,ipn)*sidevec_c(2,edg,ipn)
          dpdy = dpdy - lp_edg(k,edg,ipn)*sidevec_c(1,edg,ipn)
        end do

        old = omega(k,ipn)
        omega(k,ipn) = omega(k,ipn) +			&
          (u_velo(k,ipn)*dpdx + v_velo(k,ipn)*dpdy)*rarea(ipn)
#ifdef DEBUGPRINT
        if (ipn == PrintIpnDiag .and. mod(k,7) == 1) then
          print '(a,i8,i4,a,3f9.3)','ipn,k  = ',perm(ipn),k,		&
              '  omega terms 1+3,term 2,total:',old,(u_velo(k,ipn)*dpdx	&
              + v_velo(k,ipn)*dpdy)*rarea(ipn),omega(k,ipn)
        end if
#endif
      end do
    end do
!$OMP END PARALLEL DO
    ret = gptlstop ('update_pres')
  end subroutine update_pres
end module massadv_compute
