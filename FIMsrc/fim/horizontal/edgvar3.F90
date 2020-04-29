module module_edgvar3

  use stenedgprint
  use global_bounds,   only: ims, ime, ips, ipe, ihe
  use module_control  ,only: nvlp1, npp, ntra
  use fimnamelist     ,only: nvl, PrintIpnDiag, TimingBarriers
  use module_constants,only: cs, sn, nedge, permedge, perm, prox, proxs, &
                             nprox, hfcoef
!  use globalutils
  implicit none
  
#include <gptl.inc>
  integer :: ret

contains

!*******************************************************************************
! edgvar3
! Compute Voronoi cell edge values for high-order scheme mass and tracer flux.   
!
! N. Wang                 Jan 2014 Implementation of a higher oder flux
!                                  computation, according to Miura's 2012 paper.
!*******************************************************************************

  subroutine edgvar3 (its, u_vel, v_vel, delp, pres, exner,             &
                     geop, montg, tracr, u_edg, v_edg, dp_edg,          &
                     lp_edg, geop_edg, mont_edg, uvsq_edg, trc_edg)
! Arguments
    integer, intent(in) :: its                          ! model time step
    real, intent(in) :: u_vel (nvl,ims:ime)             ! west wind on s level
    real, intent(in) :: v_vel (nvl,ims:ime)             ! south wind on s level
    real, intent(in) :: delp  (nvl,ims:ime)             ! layer thickness
    real, intent(in) :: pres  (nvlp1,ims:ime)           ! pres on interfaces
    real, intent(in) :: geop  (nvlp1,ims:ime)           ! geopotential
    real, intent(in) :: exner (nvlp1,ims:ime)           ! exner fct on interfc
    real, intent(in) :: montg (nvl,ims:ime)             ! montgomery potential
    real, intent(in) :: tracr (nvl,ims:ime,ntra)        ! mass field tracers

!SMS$DISTRIBUTE(dh,3) BEGIN
    real, intent(inout)   :: u_edg   (nvl,npp,ims:ime)  ! west wind on edges
    real, intent(inout)   :: v_edg   (nvl,npp,ims:ime)  ! south wind on edges
    real, intent(inout)   :: dp_edg  (nvl,npp,ims:ime)  ! layer thkns on edges
    real, intent(inout)   :: lp_edg  (nvl,npp,ims:ime)  ! midlayer prs on edges
    real, intent(inout)   :: mont_edg(nvl,npp,ims:ime)  ! montg.pot. on edges
    real, intent(inout)   :: geop_edg(nvl,npp,ims:ime)  ! geopot. on edges
    real, intent(inout)   :: trc_edg (nvl,npp,ims:ime,ntra) ! tracers on edges
!SMS$DISTRIBUTE END

    real, intent(inout)   :: uvsq_edg(nvl,npp,ims:ime)  ! vel.-squared on edges

! Local workspace
    real    :: work_edg1
    real    :: work_edg2
    real    :: work_edg_tr
    real    :: var_nb(nvl,7),edge_val(nvl,6)
    integer :: nb, edg_os
    integer :: edg                   ! icos edge index
    integer :: edgcount              ! icos edge count
    integer :: ipn                   ! icos index
    integer :: k,kx,ka,kb            ! layer Index
    integer :: nt                    ! tracer index: 1=theta, 2=qv, ....
    logical :: vrbos

    ret = gptlstart ('edgvar3_loop1')
! compute the edge values for interior (non-halo) grid points. 
!$OMP PARALLEL DO PRIVATE (vrbos,k,edg,nb,var_nb,edge_val,nt)
    do ipn=ips,ipe
      vrbos = ipn == PrintIpnDiag .and. ipn == perm(ipn)

! interpolate layer thickness to edges,
      do edg = 1,nprox(ipn)
        nb = prox(edg,ipn)
        do k=1,nvl
          var_nb(k,edg) = delp(k,nb) 
        end do
      end do
      do k=1,nvl
        var_nb(k,nprox(ipn)+1) = delp(k,ipn)
      end do
      CALL edgvalEval (var_nb, dp_edg(:,:,ipn), nprox(ipn), ipn)
! and tracer variables ...
      do nt=2,ntra
        do edg = 1,nprox(ipn)
          nb = prox(edg,ipn)
          do k=1,nvl
            var_nb(k,edg) = tracr(k,nb,nt) 
          end do
        end do
        do k=1,nvl
          var_nb(k,nprox(ipn)+1) = tracr(k,ipn,nt)
        end do
        CALL edgvalEval (var_nb, trc_edg(:,:,ipn,nt), nprox(ipn), ipn)
      end do
    enddo
    ret = gptlstop ('edgvar3_loop1')

    if (TimingBarriers) then
      ret = gptlstart ('edgvar3_barrier')
!SMS$BARRIER
      ret = gptlstop ('edgvar3_barrier')
    end if

!exchange the computed edge values 
    ret = gptlstart ('edgvar3_exchange')
!SMS$EXCHANGE(dp_edg,trc_edg(:,:,:,2:ntra))
    ret = gptlstop ('edgvar3_exchange')

!assign the edge values for the interior (non-halo) cells 
    ret = gptlstart ('edgvar3_loop2')
!$OMP PARALLEL DO PRIVATE (edgcount,edg,nb,edg_os,k,work_edg1,nt,work_edg_tr)
    do ipn=ips,ihe
      do edgcount=1,nedge(ipn)
        edg = permedge(edgcount,ipn)
        nb = prox(edg,ipn)
        edg_os = proxs(edg,ipn)
! The following "if" test enables us to only store averaged values once
        if (ipn < nb) then
          do k=1,nvl
            work_edg1 = 0.5*(dp_edg(k,edg,ipn)+dp_edg(k,edg_os,nb))
            dp_edg(k,edg,ipn) = work_edg1
            dp_edg(k,edg_os,nb) = work_edg1
          end do
          do nt = 2, ntra
            do k=1,nvl
              work_edg_tr = 0.5*(trc_edg(k,edg,ipn,nt)+trc_edg(k,edg_os,nb,nt))
              trc_edg(k,edg,ipn,nt) = work_edg_tr
              trc_edg(k,edg_os,nb,nt) = work_edg_tr
            enddo
          end do
        end if
      end do
    end do
    ret = gptlstop ('edgvar3_loop2')
  end subroutine edgvar3

  subroutine edgvalEval(b, edgval, n_nb, ipn)
    implicit none
    real, intent(in) :: b(:,:)
    real, intent(out) :: edgval(:,:)
    integer, intent(in) :: n_nb, ipn 

    integer i, j, k

! Ensure pentagons have zero stored in 6th position
    do i=n_nb+1,6
      do k=1,nvl
        edgval(k,i) = 0.
      end do
    end do

    do i=1,n_nb
!DIR$ VECTOR ALWAYS
      do k=1,nvl
        edgval(k,i) = b(k,n_nb+1)
      end do
      do j=1,n_nb  
!DIR$ VECTOR ALWAYS
        do k=1,nvl
          edgval(k,i) = edgval(k,i) + (b(k,j)-b(k,n_nb+1))*hfcoef(j,i,ipn)
        end do
      enddo
    end do
  end subroutine edgvalEval
end module module_edgvar3
