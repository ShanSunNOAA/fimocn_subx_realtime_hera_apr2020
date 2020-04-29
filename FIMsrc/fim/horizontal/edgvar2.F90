module module_edgvar2

  use stenedgprint
  use global_bounds,   only: ims, ime, ips, ipe, ihe
  use module_control  ,only: nvlp1, npp, ntra
  use fimnamelist     ,only: nvl, PrintIpnDiag, janjic, eqwgt
  use module_constants,only: cs, sn, nedge, permedge, perm, prox, proxs, &
                             nprox, idxtrl, idxldg, wgttrl, wgtldg, &
                             sidevec_c, rarea, corner_xy
  use module_variables,only: grad
!  use globalutils
  implicit none
  
contains

!*********************************************************************
! edgvar2
! Step 2 of the edge variable interpolation
! N. Wang		  Jan 2014  Implementation of piece-wise 
!                                   linear interpolation between 
!                                   hexagonal cells.
!*********************************************************************

  subroutine edgvar2 (its, u_vel, v_vel, delp, pres, exner,              &
                     geop, montg, tracr, u_edg, v_edg, dp_edg,          &
                     lp_edg, geop_edg, mont_edg, uvsq_edg, trc_edg)
! Arguments
    integer, intent(in) :: its				! model time step
    real, intent(in) :: u_vel (nvl,ims:ime)		! west wind on s level
    real, intent(in) :: v_vel (nvl,ims:ime)		! south wind on s level
    real, intent(in) :: delp  (nvl,ims:ime)		! layer thickness
    real, intent(in) :: pres  (nvlp1,ims:ime)		! pres on interfaces
    real, intent(in) :: geop  (nvlp1,ims:ime)		! geopotential
    real, intent(in) :: exner (nvlp1,ims:ime)		! exner fct on interfc
    real, intent(in) :: montg (nvl,ims:ime)		! montgomery potential
    real, intent(in) :: tracr (nvl,ims:ime,ntra)	! mass field tracers
    real, intent(inout)   :: u_edg   (nvl,npp,ims:ime)	! west wind on edges
    real, intent(inout)   :: v_edg   (nvl,npp,ims:ime)	! south wind on edges
    real, intent(inout)   :: dp_edg  (nvl,npp,ims:ime)	! layer thkns on edges
    real, intent(inout)   :: lp_edg  (nvl,npp,ims:ime)	! midlayer prs on edges
    real, intent(inout)   :: trc_edg (nvl,npp,ims:ime,ntra) ! tracers on edges
    real, intent(inout)   :: geop_edg(nvl,npp,ims:ime)	! geopot. on edges
    real, intent(inout)   :: mont_edg(nvl,npp,ims:ime)	! montg.pot. on edges
    real, intent(inout)   :: uvsq_edg(nvl,npp,ims:ime)	! vel.-squared on edges

! Local workspace
    real    :: exmid(nvlp1)
    real    :: grad_conr(nvl,npp,ims:ime)	! gradient for cell corners
    real    :: uv_spd(nvl,ims:ime), uv_spd_edg(nvl,npp,ims:ime)
    real    :: uvm1,uvm2,r_m
    integer :: edg                   ! icos edge index
    integer :: ipx                   ! neighbor across edge with index 'edg'
    integer :: em1,ep1               ! edges edg-1 and edg+1 (cyclic)
    integer :: im1,ip1               ! neighbors across edg-1 and edg+1
    integer :: ipn                   ! icos index
    integer :: k,kx,ka,kb 	     ! layer Index
    integer :: nt                    ! tracer index: 1=theta, 2=qv, ....
    integer :: edgcount              ! count of icos edges
!    These are u and v at neighboring icos points, NOT on edge
    real    :: u1trl,u2trl,u3trl,v1trl,v2trl,v3trl
    real    :: u1ldg,u2ldg,u3ldg,v1ldg,v2ldg,v3ldg
    real    :: var,exbot,extop,thbot,thtop,thav
    integer :: conr1,conr2,nbaconr1,nbaconr2,nbm1conr,nbp1conr
    integer :: nbaidx,nbaedg,nbm1idx,nbm1edg,nbp1idx,nbp1edg,edgm1,edgp1
    real    :: uv1conr1,uv2conr1,uv3conr1,uv1conr2,uv2conr2,uv3conr2
    real    :: var1conr1,var2conr1,var3conr1,var1conr2,var2conr2,var3conr2 
    logical :: vrbos
    character(len=32) :: string
    real, parameter :: divby18=1./18.
    real, parameter :: divby6=1./6.
    real, parameter :: divby4=1./4.
    real, parameter :: divby3=1./3.
    real, parameter :: eps = 0.00001

! --- variable u, v  
  CALL uv_speed_edge(u_edg,v_edg,uv_spd_edg,nvl)  
  CALL gradient(uv_spd_edg,grad_conr,nvl)
  CALL uv_speed(u_vel,v_vel,uv_spd,nvl)  
!print*, 'min max of grad_conr', minval(grad_conr), maxval(grad_conr)

!Halo comp means horizontal loop indices are ips,ihe
!$OMP PARALLEL DO PRIVATE (vrbos,edgcount,edg,em1,ep1,ipx,im1,ip1,k,	   &
!$OMP                      conr1,conr2,nbaidx,nbaedg,nbaconr1,nbaconr2,	   &
!$OMP                      edgm1,nbm1idx,nbm1edg,nbm1conr,                 &
!$OMP                      edgp1,nbp1idx,nbp1edg,nbp1conr,                 &
!$OMP                      uv1conr1,uv2conr1,uv3conr1,uv1conr2,uv2conr2,uv3conr2,&
!$OMP                      uvm1,uvm2,r_m) SCHEDULE (runtime)
  do ipn=ips,ihe
    vrbos = ipn == PrintIpnDiag .and. ipn == perm(ipn)
    do edgcount=1,nedge(ipn)
!     For each edge, compute all relevant indices
      edg = permedge(edgcount,ipn)
      conr1 = edg                !first corner index for the edge
      conr2 = mod(edg,nprox(ipn))+1 !second corner index for the edge 
 
      nbaidx = prox(edg,ipn)     !cell index for neighbor across (nba)
      nbaedg = proxs(edg,ipn)    !edge index from the perspective of nba
      nbaconr1 = mod(nbaedg,nprox(nbaidx))+1 ! first corner index
      nbaconr2 = nbaedg          !second corner index 

      edgm1 = mod(edg-2+nprox(ipn),nprox(ipn))+1
      nbm1idx = prox(edgm1,ipn)  !cell index for edge-1
      nbm1edg = proxs(edgm1,ipn) !edge index from the perspective of nbm1
      nbm1conr = nbm1edg         !nbm1 corner index 

      edgp1 = mod(edg,nprox(ipn))+1
      nbp1idx = prox(edgp1,ipn) !cell index for edge+1
      nbp1edg = proxs(edgp1,ipn)  !edge index from the perspective of nbp
      nbp1conr = mod(nbp1edg,nprox(nbp1idx))+1 !nbp1 corner index

!     Interpolate the wind speed to cell's corners with gradients along the
!     center-corner lines, from each of 4 icos cells.
      do k=1,nvl
        uv1conr1 = uv_spd(k,ipn) + grad_conr(k,conr1,ipn)
        uv2conr1 = uv_spd(k,nbaidx) + grad_conr(k,nbaconr1,nbaidx)
        uv3conr1 = uv_spd(k,nbm1idx) + grad_conr(k,nbm1conr,nbm1idx)

        uv1conr2 = uv_spd(k,ipn) + grad_conr(k,conr2,ipn)
        uv2conr2 = uv_spd(k,nbaidx) + grad_conr(k,nbaconr2,nbaidx)
        uv3conr2 = uv_spd(k,nbp1idx) + grad_conr(k,nbp1conr,nbp1idx)

! --- trapezoidal:
        uvm1 = sqrt(u_edg(k,edg,ipn)**2 + v_edg(k,edg,ipn)**2)
        uvm2 = (uv1conr1+uv2conr1+uv3conr1+uv1conr2+uv2conr2+uv3conr2)*divby6
        if (uvm1 < eps) then
          r_m = 1.0
        else
          r_m = uvm2 / uvm1 
        endif
        u_edg(k,edg,ipn) = u_edg(k,edg,ipn) * r_m 
        v_edg(k,edg,ipn) = v_edg(k,edg,ipn) * r_m 


!       interpolate kinetic energy to edges
        uvsq_edg(k,edg,ipn) = .5*(u_edg(k,edg,ipn)**2		&
                                + v_edg(k,edg,ipn)**2)

      enddo ! end loop for vertical levels
    enddo  ! end loop for edges
  enddo !end loop for grid points
!$OMP END PARALLEL DO

#ifdef MORE_VARS
! more variables -- lp, mont  
  CALL gradient(lp_edg,grad_conr_v1,nvl)
  CALL gradient(mont_edg,grad_conr_v2,nvl)
!$OMP PARALLEL DO PRIVATE (vrbos,edgcount,edg,k,                           &
!$OMP                      conr1,conr2,nbaidx,nbaedg,nbaconr1,nbaconr2,	   &
!$OMP                      edgm1,nbm1idx,nbm1edg,nbm1conr,                 &
!$OMP                      edgp1,nbp1idx,nbp1edg,nbp1conr,                 &
!$OMP                      var1conr1,var2conr1,var3conr1,                  &
!$OMP                      var1conr2,var2conr2,var3conr2) SCHEDULE (runtime)
  do ipn=ips,ihe
    vrbos = ipn == PrintIpnDiag .and. ipn == perm(ipn)
    do edgcount=1,nedge(ipn)
!     for each edge, compute all relevant indices
      edg = permedge(edgcount,ipn)
      conr1 = edg                !first corner index for the edge
      conr2 = mod(edg,nprox(ipn))+1 !second corner index for the edge 
 
      nbaidx = prox(edg,ipn)     !cell index for neighbor across (nba)
      nbaedg = proxs(edg,ipn)    !edge index from the perspective of nba
      nbaconr1 = mod(nbaedg,nprox(nbaidx))+1 ! first corner index
      nbaconr2 = nbaedg          !second corner index 

      edgm1 = mod(edg-2+nprox(ipn),nprox(ipn))+1
      nbm1idx = prox(edgm1,ipn)  !cell index for edge-1
      nbm1edg = proxs(edgm1,ipn) !edge index from the perspective of nbm1
      nbm1conr = nbm1edg         !nbm1 corner index 
     
      edgp1 = mod(edg,nprox(ipn))+1
      nbp1idx = prox(edgp1,ipn) !cell index for edge+1
      nbp1edg = proxs(edgp1,ipn)  !edge index from the perspective of nbp
      nbp1conr = mod(nbp1edg,nprox(nbp1idx))+1 !nbp1 corner index

      do k = 1, nvl
!       interpolate mid-layer pressure to edges
        var1conr1 = .5*(pres(k,ipn)+pres(k+1,ipn))                       &
                  + grad_conr_v1(k,conr1,ipn)
        var2conr1 = .5*(pres(k,nbaidx)+pres(k+1,nbaidx))                &
                  + grad_conr_v1(k,nbaconr1,nbaidx)
        var3conr1 = .5*(pres(k,nbm1idx)+pres(k+1,nbm1idx))              &
                  + grad_conr_v1(k,nbm1conr,nbm1idx)

        var1conr2 = .5*(pres(k,ipn)+pres(k+1,ipn))                      &
                  + grad_conr_v1(k,conr2,ipn)
        var2conr2 = .5*(pres(k,nbaidx)+pres(k+1,nbaidx))                &
                  + grad_conr_v1(k,nbaconr2,nbaidx)
        var3conr2 = .5*(pres(k,nbp1idx)+pres(k+1,nbp1idx))              &
                  + grad_conr_v1(k,nbp1conr,nbp1idx)
 
        lp_edg(k,edg,ipn) = (var1conr1+var2conr1+var3conr1+              &
                             var1conr2+var2conr2+var3conr2)*divby6
 
      !   interpolate montgomery potential to edges
        var1conr1 = montg(k,ipn) + grad_conr_v2(k,conr1,ipn)
        var2conr1 = montg(k,nbaidx) + grad_conr_v2(k,nbaconr1,nbaidx)
        var3conr1 = montg(k,nbm1idx) + grad_conr_v2(k,nbm1conr,nbm1idx)

        var1conr2 = montg(k,ipn) + grad_conr_v2(k,conr2,ipn)
        var2conr2 = montg(k,nbaidx) + grad_conr_v2(k,nbaconr2,nbaidx)
        var3conr2 = montg(k,nbp1idx) + grad_conr_v2(k,nbp1conr,nbp1idx)

        mont_edg(k,edg,ipn) = (var1conr1+var2conr1+var3conr1+           &
                               var1conr2+var2conr2+var3conr2)*divby6

      end do			! vertical loop
    end do			! loop over edges
  end do			! horizontal loop
!$OMP END PARALLEL DO
#endif
  return

end subroutine edgvar2

! Subroutine to compute cell gradient and their projection on center-corner lines
! from first guess at edges
subroutine gradient(edge_v, grad_conr,nvl)
  integer nvl
  real, intent(in) :: edge_v(nvl,npp,ims:ime)
  real grad_conr(nvl,npp,ims:ime)
  integer ipn, edgcount,edg, k

! compute gradients of the interior (non-halo) cells, then do a halo-exchange.  
!$OMP PARALLEL DO PRIVATE (edgcount,edg,k) SCHEDULE (runtime)  
  do ipn=ips,ipe
    grad(1:nvl,1:2,ipn) = 0.0
    do edgcount=1,nedge(ipn)
      edg = permedge(edgcount,ipn)
      do k=1,nvl
        grad(k,1,ipn) = grad(k,1,ipn) + edge_v(k,edg,ipn)*sidevec_c(2,edg,ipn)
        grad(k,2,ipn) = grad(k,2,ipn) - edge_v(k,edg,ipn)*sidevec_c(1,edg,ipn)
      enddo
    enddo
    grad(1:nvl,1:2,ipn) = grad(1:nvl,1:2,ipn) * rarea(ipn)
  enddo
!$OMP END PARALLEL DO

!SMS$EXCHANGE(grad)

! compute the projections of cell gradient onto center-corner lines for each cell, including halo-cells
! (corner_xy are constant, therefore no need to do halo-exchange every time step.
!$OMP PARALLEL DO PRIVATE (edg,k) SCHEDULE (runtime)  
  do ipn=ips,ihe
    do edg = 1,nprox(ipn)
      do k=1,nvl
        grad_conr(k,edg,ipn) = 0.5*(grad(k,1,ipn)*corner_xy(edg,1,ipn) &
                             + grad(k,2,ipn)*corner_xy(edg,2,ipn)) 
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

end subroutine gradient

! Subroutine to convert u,v at edge from the projection at the middle of edge
! to the projection at the cell center 
subroutine edge2center(u_edge,v_edge,nvl)
  integer nvl
  real, intent(inout) ::  u_edge(nvl,npp,ims:ime),v_edge(nvl,npp,ims:ime)
  integer ipn, edgcount,edg, k
  real u, v

! calculation of new u,v 
!$OMP PARALLEL DO PRIVATE (edgcount,edg,k,u,v) SCHEDULE (runtime)  
  do ipn=ips,ihe
    do edgcount=1,nedge(ipn)
      edg = permedge(edgcount,ipn)
      do k=1,nvl
        u = cs(1,edg,ipn)*u_edge(k,edg,ipn) - sn(1,edg,ipn)*v_edge(k,edg,ipn)
        v = sn(1,edg,ipn)*u_edge(k,edg,ipn) + cs(1,edg,ipn)*v_edge(k,edg,ipn)
        u_edge(k,edg,ipn) = u
        v_edge(k,edg,ipn) = v
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

end subroutine edge2center

subroutine uv_speed_edge(u_edge, v_edge, uv_spd_edge, nvl)  
  integer, intent(in) :: nvl
  real, intent(in) :: u_edge(nvl,npp,ims:ime),v_edge(nvl,npp,ims:ime)
  real, intent(out) :: uv_spd_edge(nvl,npp,ims:ime)

  integer ipn,edgcount,edg,k

!$OMP PARALLEL DO PRIVATE (edgcount,edg,k) SCHEDULE (runtime)  
  do ipn=ips,ihe
    do edgcount=1,nedge(ipn)
      edg = permedge(edgcount,ipn)
      do k=1,nvl
        uv_spd_edge(k,edg,ipn) = sqrt(u_edge(k,edg,ipn)**2 + &
                                      v_edge(k,edg,ipn)**2) 
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

end subroutine uv_speed_edge

subroutine uv_speed(u_velo, v_velo, uv_spd, nvl)  
  integer, intent(in) :: nvl
  real, intent(in) :: u_velo(nvl,ims:ime), v_velo(nvl,ims:ime)
  real, intent(out) :: uv_spd(nvl,ims:ime)

  integer ipn, k

!$OMP PARALLEL DO PRIVATE (k) SCHEDULE (runtime)  
  do ipn=ips,ihe
    do k=1,nvl
      uv_spd(k,ipn) = sqrt(u_velo(k,ipn)**2 + &
                           v_velo(k,ipn)**2) 
    enddo
  enddo
!$OMP END PARALLEL DO

end subroutine uv_speed


end module module_edgvar2
