!------------------------------------------------------------------------------------------------
! This module contains the implementation of the flux limiter algorithm
! described in Thuburn 1996 and Miura 2013.
!  
! Dec. 2015, N. Wang, original version. 
! 
! May 2016, N. Wang, modified the original algorithm to work better with
! Adams-bashforth III time integration scheme.
!
!------------------------------------------------------------------------------------------------
#define M1
module Fluxlimiter
  use module_control,   only: npp, nabl, dt
  use fimnamelist,      only: nvl, TimingBarriers
  use module_variables, only: nf, of, vof, q_min_out, q_max_out, adbash1, adbash2, adbash3
  use global_bounds,    only: ims, ime, ips, ipe, ihe
  use module_constants, only: prox,proxs,nedge,nprox,sidevec_e,sideln,permedge,rarea

  implicit none

  real,public,allocatable,save :: vnorm(:,:,:) ! velocity normal to edges
  real,public,allocatable,save :: c_num(:,:,:) ! Courant number for the flow normal to edges
  real,public,allocatable,save :: c_num_t(:,:,:) ! Courant number for the flow normal to edges
  
  private
  logical,save :: mem_allocated = .false.

  public :: init_mass_limiter, init_tracer_limiter, limiter
#include <gptl.inc>

contains
  subroutine init_mass_limiter(u_edg, v_edg, dp_cent, dp_edge)
    real, intent(in)    :: u_edg (nvl,npp,ims:ime)   ! u on edges, m/sec
    real, intent(in)    :: v_edg (nvl,npp,ims:ime)   ! v on edges, m/sec
    real, intent(in)    :: dp_cent(nvl,ims:ime)        ! depth value
    real, intent(in)    :: dp_edge(nvl,npp,ims:ime)  ! depth edge value
    
    if (.not. mem_allocated ) then
      allocate(vnorm(nvl,npp,ims:ime))
      allocate(c_num(nvl,npp,ims:ime))
      allocate(c_num_t(nvl,npp,ims:ime))
      mem_allocated  = .true.
    endif

    call edgeNormVel(u_edg, v_edg)

  end subroutine init_mass_limiter
     
  subroutine init_tracer_limiter(u_edg, v_edg, dp_cent, dp_edge)
    real, intent(in)    :: u_edg (nvl,npp,ims:ime)   ! u on edges, m/sec
    real, intent(in)    :: v_edg (nvl,npp,ims:ime)   ! v on edges, m/sec
    real, intent(in)    :: dp_cent(nvl,ims:ime)        ! depth value
    real, intent(in)    :: dp_edge(nvl,npp,ims:ime)  ! depth edge value
    
    if (.not. mem_allocated) then
      allocate(vnorm(nvl,npp,ims:ime))
      allocate(c_num(nvl,npp,ims:ime))
      allocate(c_num_t(nvl,npp,ims:ime))
      mem_allocated  = .true.
      call edgeNormVel(u_edg, v_edg)
    endif

    call c_num_tilde(dp_cent,dp_edge)

  end subroutine init_tracer_limiter
     

  subroutine edgeNormVel(u_edg, v_edg)
    real, intent(in)    :: u_edg (nvl,npp,ims:ime)   ! u on edges, m/sec
    real, intent(in)    :: v_edg (nvl,npp,ims:ime)   ! v on edges, m/sec

    integer :: ipn        ! index for cell
    integer :: edg	  ! index for icos edge number
    integer :: edgcount   ! count of icos edges
    integer :: k          ! index for vertical level
    integer :: ret

    ret = gptlstart ('edgeNormVel')
!$OMP PARALLEL DO PRIVATE (k,edgcount,edg) SCHEDULE (runtime)
    do ipn=ips,ihe
      do edgcount=1,nedge(ipn)	! loop through edges
        edg = permedge(edgcount,ipn)
        do k=1,nvl
          vnorm(k,edg,ipn) = sidevec_e(2,edg,ipn)*u_edg(k,edg,ipn) - &
                             sidevec_e(1,edg,ipn)*v_edg(k,edg,ipn)  
          c_num(k,edg,ipn) = abs(vnorm(k,edg,ipn))*dt*rarea(ipn)
        end do
      end do
    end do
!$OMP END PARALLEL DO
    ret = gptlstop ('edgeNormVel')
  end subroutine edgeNormVel
  
  subroutine c_num_tilde(dp_cent,dp_edge)
    real, intent(in)    :: dp_cent(nvl,ims:ime)        ! depth value
    real, intent(in)    :: dp_edge(nvl,npp,ims:ime)  ! depth edge value

    integer :: ipn        ! index for cell
    integer :: edg	  ! index for icos edge number
    integer :: edgcount   ! count of icos edges
    integer :: k          ! index for vertical level
    integer :: ret

    ret = gptlstart ('c_num_tilde')
!$OMP PARALLEL DO PRIVATE (k,edgcount,edg) SCHEDULE (runtime)
    do ipn=ips,ihe
      do edgcount=1,nedge(ipn)	! loop through edges
        edg = permedge(edgcount,ipn)
        do k=1,nvl
          c_num_t(k,edg,ipn) = c_num(k,edg,ipn) * dp_edge(k,edg,ipn) / &
                               dp_cent(k,ipn) 
        end do
      end do
    end do
!$OMP END PARALLEL DO
    ret = gptlstop ('c_num_tilde')
  end subroutine c_num_tilde


  subroutine limiter(gridVal, edgeVal, c_num_inp, dp_tndcy)
! Arguments
    real, intent(in)    :: gridVal(nvl,    ims:ime)   ! grid value 
    real, intent(inout) :: edgeVal(nvl,npp,ims:ime)   ! edge value to be limited
    real, intent(inout) :: c_num_inp(nvl,npp,ims:ime) ! either c_num or c_num_t 
    real, intent(in)    :: dp_tndcy(nvl,ims:ime,nabl) ! Pa/sec

! Local variables
    real :: q_min, q_max, q_star(nvl,npp,ims:ime)
    real :: c_num_in_sum(nvl), c_num_out_sum(nvl)
    real :: cq_min_in_sum(nvl), cq_max_in_sum(nvl)
    real :: q_min_in_min(nvl), q_max_in_max(nvl)
    real :: q_min_next, q_max_next
    real :: ev, ev2, qs

    integer :: ipn
    integer :: ipx	  ! neighbor across joint edge
    integer :: edg	  ! Index for icos edge number
    integer :: edgcount   ! count of icos edges
    integer :: edgm1,edgp1,nbm1,nbp1
    integer :: k
    integer :: ret
    real, parameter :: max_out_bound = 100000.0

    ret = gptlstart ('limiter')
! q_star can be computed in the halo, thus no need to exchange it
!$OMP PARALLEL DO PRIVATE (edgcount, edg, edgm1, nbm1, edgp1, nbp1, ipx, k, q_min, q_max, &
!$OMP                      q_min_in_min, q_max_in_max, cq_min_in_sum, cq_max_in_sum,      &
!$OMP                      c_num_in_sum, c_num_out_sum, q_min_next, q_max_next)
    do ipn=ips,ihe
      do edgcount=1,nedge(ipn)
        edg = permedge(edgcount,ipn)
        edgm1 = mod(edg-2+nprox(ipn),nprox(ipn))+1
        nbm1 = prox(edgm1,ipn)
        edgp1 = mod(edg,nprox(ipn))+1
        nbp1 = prox(edgp1,ipn)
        ipx = prox(edg,ipn)          ! neighbor across shared edge
        do k=1,nvl
          q_min = min(gridVal(k,ipn), gridVal(k,ipx), gridVal(k,nbm1), gridVal(k,nbp1)) 
          q_max = max(gridVal(k,ipn), gridVal(k,ipx), gridVal(k,nbm1), gridVal(k,nbp1)) 
          q_star(k,edg,ipn) = max(q_min, min(q_max, edgeVal(k,edg,ipn)))
        end do
      end do

! q_min_out and q_max_out computed at the end of the loop will need to be exchanged
      if (ipn <= ipe) then
        q_min_in_min(:)  = huge(q_min_in_min)
        q_max_in_max(:)  = 0.0 
        cq_min_in_sum(:) = 0.0
        cq_max_in_sum(:) = 0.0
        c_num_in_sum(:)  = 0.0
        c_num_out_sum(:) = 0.0
        do edgcount=1,nedge(ipn)
          edg = permedge(edgcount,ipn)
          edgm1 = mod(edg-2+nprox(ipn),nprox(ipn))+1
          nbm1 = prox(edgm1,ipn)
          edgp1 = mod(edg,nprox(ipn))+1
          nbp1 = prox(edgp1,ipn)
          ipx = prox(edg,ipn)          ! neighbor across shared edge
          do k=1,nvl
            if (vnorm(k,edg,ipn) > 0) then
#ifdef M1
              q_min = min(q_star(k,edg,ipn),gridVal(k,ipn))
              q_max = max(q_star(k,edg,ipn),gridVal(k,ipn))
#else
              q_min = min(q_star(k,edg,ipn),gridVal(k,ipn),gridVal(k,nbm1), gridVal(k,nbp1))
              q_max = max(q_star(k,edg,ipn),gridVal(k,ipn),gridVal(k,nbm1), gridVal(k,nbp1))
#endif
              c_num_out_sum(k) = c_num_out_sum(k)+c_num_inp(k,edg,ipn)
            else
#ifdef M1
              q_min = min(q_star(k,edg,ipn),gridVal(k,ipx))
              q_max = max(q_star(k,edg,ipn),gridVal(k,ipx))
#else
              q_min = min(q_star(k,edg,ipn),gridVal(k,ipx),gridVal(k,nbm1), gridVal(k,nbp1))
              q_max = max(q_star(k,edg,ipn),gridVal(k,ipx),gridVal(k,nbm1), gridVal(k,nbp1))
#endif
              c_num_in_sum(k) = c_num_in_sum(k)+ c_num_inp(k,edg,ipn)

              q_min_in_min(k) = min(q_min_in_min(k),q_min)
              q_max_in_max(k) = max(q_max_in_max(k),q_max)
              cq_min_in_sum(k) = cq_min_in_sum(k) + c_num_inp(k,edg,ipn)*q_min
              cq_max_in_sum(k) = cq_max_in_sum(k) + c_num_inp(k,edg,ipn)*q_max
            endif
          end do
        end do
! outflow bounds -------
! Commented out code implements Thuburn 1996 and Miura 2013. Actual code reflects modifications for
! Adams-Bashforth 3 time integration scheme
        do k=1,nvl
          q_min_next = min(gridVal(k,ipn),q_min_in_min(k)) 
          q_max_next = max(gridVal(k,ipn),q_max_in_max(k)) 
          if (c_num_out_sum(k) /= 0.0) then
!          q_min_out(k,ipn) = (gridVal(k,ipn) + cq_max_in_sum(k) - &
!            q_max_next * (1.0 + c_num_in_sum(k) - c_num_out_sum(k)) ) / &
!            c_num_out_sum(k)
            q_min_out(k,ipn) = (gridVal(k,ipn) + adbash1/dt*cq_max_in_sum(k) - &
                 q_max_next * (1.0 + c_num_in_sum(k) - c_num_out_sum(k)) + &
                 adbash2*dp_tndcy(k,ipn,of) + adbash3*dp_tndcy(k,ipn,vof) ) / &
                 c_num_out_sum(k) / adbash1 * dt
!          q_max_out(k,ipn) = (gridVal(k,ipn) + cq_min_in_sum(k) - &
!            q_min_next * (1.0 + c_num_in_sum(k) - c_num_out_sum(k))) / &
!            c_num_out_sum(k)
            q_max_out(k,ipn) = (gridVal(k,ipn) + adbash1/dt*cq_min_in_sum(k) - &
                 q_min_next * (1.0 + c_num_in_sum(k) - c_num_out_sum(k)) + &
                 adbash2*dp_tndcy(k,ipn,of) + adbash3*dp_tndcy(k,ipn,vof) ) / &
                 c_num_out_sum(k) / adbash1 * dt
          else
            q_min_out(k,ipn) = -max_out_bound
            q_max_out(k,ipn) = max_out_bound
          endif
        end do
      end if
    end do
!$OMP END PARALLEL DO

    if (TimingBarriers) then
      ret = gptlstart ('limiter_barrier')
!SMS$BARRIER
      ret = gptlstop ('limiter_barrier')
    end if

    ret = gptlstart ('limiter_exchange')
!SMS$EXCHANGE(q_min_out,q_max_out)
    ret = gptlstop ('limiter_exchange')

!$OMP PARALLEL DO PRIVATE (k,edgcount,edg,ipx) SCHEDULE (runtime)
    do ipn=ips,ihe
      do edgcount=1,nedge(ipn)
        edg = permedge(edgcount,ipn)
        ipx = prox(edg,ipn)
        do k=1,nvl
          if (vnorm(k,edg,ipn) > 0) then
            edgeVal(k,edg,ipn) = max(q_min_out(k,ipn),min(q_max_out(k,ipn),q_star(k,edg,ipn))) 
          else 
            edgeVal(k,edg,ipn) = max(q_min_out(k,ipx),min(q_max_out(k,ipx),q_star(k,edg,ipn))) 
          endif
       end do
      end do 
    end do
!$OMP END PARALLEL DO
    ret = gptlstop ('limiter')
  end subroutine limiter
end module FluxLimiter
