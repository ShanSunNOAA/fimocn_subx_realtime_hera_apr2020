!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines in this module are computational only. They contain no SMS directives
! except for IGNORE. Thus they are good candidates for migration to MIC or GPU.
!
! cnuity_compute: contains all the computational routines for continuity equation 
! calculations.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cnuity_compute
  implicit none

  private
  public :: cnuity1, cnuity2, cnuity3, cnuity4, cnuity5

contains

  subroutine cnuity1 (antifx, massfx, u_edg, v_edg, dp_edg, delp)
    use module_control,   only: npp
    use fimnamelist,      only: nvl
    use module_variables, only: nf
    use module_constants, only: prox, proxs, nedge, sidevec_e, permedge
    use global_bounds,    only: ims, ime, ips, ipe, ihe
! Arguments
    real, intent(out)   :: antifx(nvl,npp,ims:ime)   ! N/sec
    real, intent(inout) :: massfx(nvl,npp,ims:ime,3) ! N/sec
    real, intent(in)    :: u_edg (nvl,npp,ims:ime)   ! u on edges, m/sec
    real, intent(in)    :: v_edg (nvl,npp,ims:ime)   ! v on edges, m/sec
    real, intent(in)    :: dp_edg(nvl,npp,ims:ime)   ! dp on edges, Pa
    real, intent(in)    :: delp  (nvl,    ims:ime)   ! layer thickness, Pa

! Local variables
    real    :: vnorm(nvl,npp,ims:ime)    ! outward-directed velocity on edge x edge length
    real    :: flxhi
    integer :: ipn
    integer :: ipx	  ! neighbor across joint edge
    integer :: edg	  ! Index for icos edge number
    integer :: edgcount   ! count of icos edges
    integer :: edx
    integer :: k

! Halo comp means horizontal loop indices are ips,ihe
!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (k,edgcount,edg) SCHEDULE (runtime)
    do ipn=ips,ihe
      do edgcount=1,nedge(ipn)	! loop through edges
        edg = permedge(edgcount,ipn)
        do k=1,nvl
          vnorm(k,edg,ipn) = sidevec_e(2,edg,ipn)*u_edg(k,edg,ipn) - &
                             sidevec_e(1,edg,ipn)*v_edg(k,edg,ipn)  
        end do
      end do
    end do
!$OMP END PARALLEL DO
!sms$ignore end

! Halo comp means horizontal loop indices are ips,ihe
!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (edgcount,edg,ipx,edx,k,flxhi) SCHEDULE (runtime)
    do ipn=ips,ihe
      ! Initialize these local and INTENT(OUT) arrays so COMPARE_VAR does not get 
      ! confused by uninitialized edg = 6 edges of pentagonal grid cells. NOTE
      ! that this only *partially* initializes massfx, which must therefore be
      ! declared intent inout, not out.
      do k=1,nvl
        antifx(k,npp,ipn)    = 0.
        massfx(k,npp,ipn,nf) = 0.
      end do

      do edgcount=1,nedge(ipn)
        edg = permedge(edgcount,ipn)
        ipx = prox(edg,ipn)          ! neighbor across shared edge
        edx = proxs(edg,ipn)	     ! neighbors index of shared edge
        do k=1,nvl
! high-order mass flux (2nd order centered, out-of cell > 0)
          flxhi = vnorm(k,edg,ipn)*dp_edg(k,edg,ipn)

! low-order mass flux (donor-cell, out-of cell > 0)
          massfx(k,edg,ipn,nf) =				&
               0.5*((vnorm(k,edg,ipn) +				&
                 abs(vnorm(k,edg,ipn)))*delp(k,ipn)		&
                  - (vnorm(k,edx,ipx) +				&
                 abs(vnorm(k,edx,ipx)))*delp(k,ipx))

! anti-diffusive flux = difference between high order flux (flxhi)
! and low-order flux (massfx), from Zalesek, p336, Eqn (3):
          antifx(k,edg,ipn) = flxhi - massfx(k,edg,ipn,nf)   	! N/sec
        end do
      end do
    end do
!$OMP END PARALLEL DO
!sms$ignore end

    return

  end subroutine cnuity1

  subroutine cnuity2 (dplo_tndcy, massfx, delp, delp_lo)
    use module_control,   only: npp, nabl
    use fimnamelist,      only: nvl
    use module_variables, only: nf, of, vof, adbash1, adbash2, adbash3
    use module_constants, only: nprox, rarea
    use global_bounds,    only: ims, ime, ips, ipe
! Arguments
    real, intent(inout) :: dplo_tndcy(nvl,ims:ime,nabl) ! Pa/sec (inout not out due to halo)
    real, intent(in)    :: massfx(nvl,npp,ims:ime,3)    ! N/sec
    real, intent(in)    :: delp(nvl,ims:ime)            ! layer thickness, Pa
    real, intent(inout) :: delp_lo(nvl,ims:ime)         ! Pa (inout not out due to halo)

! Local variables
    integer :: ipn  ! horizontal index
    integer :: edg  ! edge index
    integer :: k    ! layer index

! sum up edge fluxes to get low-order tendency term 
!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (edg,k) SCHEDULE (runtime)
    do ipn=ips,ipe
      do k=1,nvl
        dplo_tndcy(k,ipn,nf) = 0.
      end do

      do edg=1,nprox(ipn)
        do k=1,nvl            
          dplo_tndcy(k,ipn,nf) = dplo_tndcy(k,ipn,nf) + massfx(k,edg,ipn,nf)
        end do
      end do

!.............................................................
!  Sec. 2.  Calculate new low-order dp using full Adams Bashforth
!.............................................................

      do k=1,nvl
        ! divide tendency by cell area to convert to  Pa/sec
        dplo_tndcy(k,ipn,nf) = -dplo_tndcy(k,ipn,nf)*rarea(ipn)
        
        ! get new value for the low-order delp field using Adams Bashforth
        ! 3 time levels, the one just calculated (nf), and the two prev ones,
        ! marked of (old field) and vof (very old field):
        
        delp_lo(k,ipn) = delp(k,ipn) + adbash1*dplo_tndcy(k,ipn, nf) + &
                                       adbash2*dplo_tndcy(k,ipn, of) + &
                                       adbash3*dplo_tndcy(k,ipn,vof)
      end do
    end do
!$OMP END PARALLEL DO
!sms$ignore end

    return

  end subroutine cnuity2

  subroutine cnuity3 (p_plus, p_mnus, antifx, delp_lo, delp, dp_tndcy, r_plus, r_mnus)
    use module_control,   only: npp, nabl
    use fimnamelist,      only: nvl
    use module_constants, only: nprox, prox, area
    use module_variables, only: nf, of, vof, adbash1, adbash2, adbash3
    use global_bounds,    only: ims, ime, ips, ipe

! Arguments
    real, intent(out) :: p_plus  (nvl,ims:ime)       ! Zalesak's p_plus, N/sec
    real, intent(out) :: p_mnus  (nvl,ims:ime)       ! Zalesak's p_minus, N/sec
    real, intent(in)  :: antifx  (nvl,npp,ims:ime)   ! N/sec
    real, intent(in)  :: delp_lo (nvl,ims:ime)       ! Pa
    real, intent(in)  :: delp    (nvl,ims:ime)       ! layer thickness, Pa
    real, intent(in)  :: dp_tndcy(nvl,ims:ime,nabl)  ! Pa/sec
    real, intent(out) :: r_plus  (nvl,ims:ime)	     ! Zalesak's r_plus, dimensionless
    real, intent(out) :: r_mnus  (nvl,ims:ime)	     ! Zalesak's r_minus, dimensionless

! Local variables
    integer :: ipn
    integer :: k
    integer :: edg

    real, parameter :: thshld = 1.e-11
    real :: q_plus,q_mnus               ! Zalesak's q_plus,q_mnus (N/sec)
    real :: dpmax,dpmin	                ! Zalesak's phi_max,phi_min

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (edg,k,dpmax,dpmin,q_plus,q_mnus) SCHEDULE (runtime)
    do ipn=ips,ipe
      do k=1,nvl
        p_plus(k,ipn) = 0. 
        p_mnus(k,ipn) = 0.
      end do

      do edg=1,nprox(ipn)
        do k=1,nvl
          if (antifx(k,edg,ipn) <= 0.) then		      ! flux into cell
            p_plus(k,ipn) = p_plus(k,ipn) - antifx(k,edg,ipn)
          else						      ! flux out of cell
            p_mnus(k,ipn) = p_mnus(k,ipn) + antifx(k,edg,ipn)    
          end if
        end do
      end do

!   At this stage we have calculated the low-order delp_lo, and the 
!   unclipped antidiffusive flux for the entire icos global grid.

!............................................................
!   Sec 3.  Monotonicity Limit on Fluxes
!............................................................

!   Determine the amount of antidiffusive fluxes that can be
!   added to the low-order solution without violating monotonicity.

      do k=1,nvl

        ! According to Zal (17), we must limit according to max of
        ! dp from any gazebo direction, in either current or prev time
        ! step.  Likewise for the minimum.
        ! For the pentagons prox(6,ipn) is set to prox(5,ipn) in init.F90.
     
        dpmax = max (delp_lo(k,ipn), delp(k,ipn),		&
             delp_lo(k,prox(1,ipn)), delp_lo(k,prox(2,ipn)),	&
             delp_lo(k,prox(3,ipn)), delp_lo(k,prox(4,ipn)),	&
             delp_lo(k,prox(5,ipn)), delp_lo(k,prox(6,ipn)),	&
             delp(k,prox(1,ipn)), delp(k,prox(2,ipn)),	&
             delp(k,prox(3,ipn)), delp(k,prox(4,ipn)),	&
             delp(k,prox(5,ipn)), delp(k,prox(6,ipn))  )

        dpmin = min (delp_lo(k,ipn), delp(k,ipn),		&
             delp_lo(k,prox(1,ipn)), delp_lo(k,prox(2,ipn)),	&
             delp_lo(k,prox(3,ipn)), delp_lo(k,prox(4,ipn)),	&
             delp_lo(k,prox(5,ipn)), delp_lo(k,prox(6,ipn)),  &
             delp(k,prox(1,ipn)), delp(k,prox(2,ipn)),	&
             delp(k,prox(3,ipn)), delp(k,prox(4,ipn)),	&
             delp(k,prox(5,ipn)), delp(k,prox(6,ipn))  )

        dpmax = max(0.,dpmax)	! cannot allow negative dpmax
        dpmin = max(0.,dpmin)	! cannot allow negative dpmin	

        ! q_plus/q_mnus are the upper/lower limits on antidiffusive dp tendencies
        ! q_plus is Zalesak (8):
        q_plus = (dpmax - delp_lo(k,ipn)		& ! N/sec
               -(adbash2*min(0.,dp_tndcy(k,ipn, of))	&
               + adbash3*max(0.,dp_tndcy(k,ipn,vof))))	&
               /adbash1*area(ipn)

        ! q_mnus is Zalesak (11):
        q_mnus = (delp_lo(k,ipn) - dpmin		& ! N/sec
               +(adbash2*max(0.,dp_tndcy(k,ipn, of))	&
               + adbash3*min(0.,dp_tndcy(k,ipn,vof))))	&
               /adbash1*area(ipn)

        !  Having p_plus(k,ipn) and q_plus, we can calc r_plus, Zal (9):

        !  reduce fluxes to stay within limits posed by q_plus,q_mnus.
        !  r_plus,r_mnus are dimensionless
          
        r_plus(k,ipn) = max(0.,min(1.,q_plus/max(p_plus(k,ipn),thshld)))
        r_mnus(k,ipn) = max(0.,min(1.,q_mnus/max(p_mnus(k,ipn),thshld)))
      end do			! loop through layers
    end do			! horizontal loop
!$OMP END PARALLEL DO
!sms$ignore end

    return

  end subroutine cnuity3

  subroutine cnuity4 (psurf, antifx, r_plus, r_mnus, massfx, &
                      anti_tndcy, dplo_tndcy, dp_tndcy, delp)
    use module_control,   only: npp, nabl
    use fimnamelist,      only: nvl
    use module_variables, only: nf, of, vof, adbash1, adbash2, adbash3
    use module_constants, only: nprox, prox, rarea, nedge, permedge
    use global_bounds,    only: ims, ime, ips, ipe, ihe

    real, intent(out)   :: psurf(ims:ime)               ! surface pressure
    real, intent(inout) :: antifx(nvl,npp,ims:ime)      ! N/sec
    real, intent(in)    :: r_plus(nvl,ims:ime)          ! Zalesak's r_plus, dimensionless
    real, intent(in)    :: r_mnus(nvl,ims:ime)          ! Zalesak's r_minus, dimensionless
    real, intent(inout) :: massfx(nvl,npp,ims:ime,3)    ! N/sec
    real, intent(out)   :: anti_tndcy(nvl,ims:ime)      ! Pa/sec
    real, intent(in)    :: dplo_tndcy(nvl,ims:ime,nabl) ! Pa/sec
    real, intent(inout) :: dp_tndcy  (nvl,ims:ime,nabl) ! Pa/sec
    real, intent(inout) :: delp(nvl,ims:ime)

! Local variables
    integer :: ipn
    integer :: edg
    integer :: edgcount            ! count of icos edges
    integer :: k

!ss    real    :: recovr(npp,ims:ime) ! fluxes discarded in clipping process
    real    :: clip
    
! Halo comp means horizontal loop indices are ips,ihe
!sms$ignore begin
!DIR$ vector always
!$OMP PARALLEL DO PRIVATE (k,edgcount,edg,clip) SCHEDULE (runtime)
    do ipn=ips,ihe
      psurf(ipn) = 0.
      do k=nvl,1,-1      ! top-down for psurf
        psurf(ipn) = psurf(ipn) + delp(k,ipn)
      end do

      do edgcount=1,nedge(ipn)
        edg = permedge(edgcount,ipn)
!ss        recovr(edg,ipn) = 0.
        do k=1,nvl
          if (antifx(k,edg,ipn) >= 0.) then              ! outgoing
            clip = min (r_mnus(k,ipn), r_plus(k,prox(edg,ipn)))
          else						! incoming
            clip = min (r_plus(k,ipn), r_mnus(k,prox(edg,ipn)))
          end if

          ! limit antidiffusive fluxes
          ! set aside vertical integral of discarded fluxes for later recovery
!ss          recovr(edg,ipn) = recovr(edg,ipn) + antifx(k,edg,ipn)*(1.-clip)
          antifx(k,edg,ipn) = antifx(k,edg,ipn)*clip

          ! add clipped antidiff to low-order flux to obtain total mass flux
          massfx(k,edg,ipn,nf) = massfx(k,edg,ipn,nf) + antifx(k,edg,ipn)
        end do
      end do
    end do
!$OMP END PARALLEL DO
!sms$ignore end

!ss   ! include fluxes discarded during clipping process as barotropic
!ss   !  corrections to antidiffusive fluxes
!ss
!ss   do ipn=ips,ihe
!ss    do edgcount = 1,nedge(ipn)	! loop through edges
!ss     edg = permedge(edgcount,ipn)
!ss     if (recovr(edg,ipn).ge.0.) then		! outgoing
!ss      do k = 1,nvl
!ss       antifx(k,edg,ipn) = antifx(k,edg,ipn)			&
!ss         +recovr(edg,ipn)*delp   (k,ipn )/psurf(ipn )
!ss!        +recovr(edg,ipn)*delp_lo(k,ipn )/psurf(ipn )
!ss      end do
!ss     else					! incoming
!ss      ipx = prox(edg,ipn)
!ss      do k = 1,nvl
!ss       antifx(k,edg,ipn) = antifx(k,edg,ipn)			&
!ss         +recovr(edg,ipn)*delp   (k,ipx)/psurf(ipx)
!ss!        +recovr(edg,ipn)*delp_lo(k,ipx)/psurf(ipx)
!ss      end do
!ss     end if			! recovr > or < 0
!ss    end do			! loop through edges
!ss   end do			! horizontal loop

!JR What does "vector always" do?
!sms$ignore begin
!DIR$ vector always
!$OMP PARALLEL DO PRIVATE (k,edg) SCHEDULE (runtime)
    do ipn=ips,ipe
      do k=1,nvl
        anti_tndcy(k,ipn) = 0.

        ! sum up edge fluxes to get antidiffusive tendency term
        do edg=1,nprox(ipn)
          anti_tndcy(k,ipn) = anti_tndcy(k,ipn) + antifx(k,edg,ipn)
        end do

        ! divide antidiff tendency by cell area to convert to  Pa/sec
        anti_tndcy(k,ipn) = -anti_tndcy(k,ipn)*rarea(ipn)

        ! combine low-order with clipped antidiffusive tendency
        dp_tndcy(k,ipn,nf) = dplo_tndcy(k,ipn,nf) + anti_tndcy(k,ipn)

        ! advance delp to new time step
        delp(k,ipn) = delp(k,ipn) + adbash1*dp_tndcy(k,ipn,nf) + &
                                    adbash2*dp_tndcy(k,ipn,of) + &
                                    adbash3*dp_tndcy(k,ipn,vof)
      end do
    end do
!$OMP END PARALLEL DO
!sms$ignore end

    return

  end subroutine cnuity4

  subroutine cnuity5 (pres, delp, exner, dp_tndcy, omega, &
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

!sms$ignore begin
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
!sms$ignore end

    return

  end subroutine cnuity5

end module cnuity_compute
