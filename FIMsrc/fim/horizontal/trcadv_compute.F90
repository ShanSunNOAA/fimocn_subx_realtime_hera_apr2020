!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines in this module are computational only. They contain no SMS directives
! except for IGNORE. Thus they are good candidates for migration to MIC or GPU.
!
! trcadv_compute: contains all the computational routines for tracer advection
! calculations.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module trcadv_compute
  implicit none

  private
  real, parameter :: thshld = 1.e-11

  public :: trcadv1, trcadv2, trcadv3

contains

  subroutine trcadv1 (t, trclo_tdcy, trc_tdcy, massfx, trcr_edg,	&
                      tracr, trcdp, antiflx, trmin, trmax, trcr_lo,	&
                      flx_lo, delp)
    use module_control,   only: npp, nabl, ntra
    use fimnamelist,      only: nvl
    use module_variables, only: nf, of, vof, adbash1, adbash2, adbash3
    use module_constants, only: nprox, prox, proxs, rarea, nedge, permedge
    use global_bounds,    only: ims, ime, ips, ipe, ihe

! Arguments
    integer, intent(in) :: t                                  ! tracer index
    real, intent(inout) :: trclo_tdcy(nvl,ims:ime,nabl,ntra)  ! low order tdcy
    real, intent(inout) :: trc_tdcy  (nvl,ims:ime,nabl,ntra)  ! final trcr tdcy
    real, intent(in)    :: massfx    (nvl,npp,ims:ime,3)      ! mass flux across edges
    real, intent(in)    :: trcr_edg  (nvl,npp,ims:ime,ntra)   ! tracer on edges
    real, intent(in)    :: tracr     (nvl,ims:ime,ntra)       ! tracer
    real, intent(in)    :: trcdp     (nvl,ims:ime,ntra)       ! tracr*dp
    real, intent(inout) :: antiflx   (nvl,npp,ims:ime)        ! antidiffusive tracer flux
    real, intent(out)   :: trmax(nvl,ims:ime)		! regional tracer max for flux clipping
    real, intent(out)   :: trmin(nvl,ims:ime)		! regional tracer min for flux clipping
    real, intent(out)   :: trcr_lo(nvl,ims:ime)		! tracer after low-order transport
    real, intent(out)   :: flx_lo(nvl,npp,ims:ime)	! tracer flux, low order
    real, intent(in)    :: delp(nvl,ims:ime)

! Local variables
    integer :: ipn
    integer :: k
    integer :: edgcount		! count of icos edges
    integer :: edg		! icos edge number index
    integer :: ipx		! neighbor across joint edge
    integer :: edx		! neighbor's index of joint edge
    real :: flx_hi		! tracer flux, high order
    real :: trcdp_lo		! tracer*dp after low-order transport

! Halo comp so horizontal loop indices are ips,ihe
!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (edgcount,edg,ipx,edx,k,flx_hi,trcdp_lo) SCHEDULE (runtime)
    do ipn=ips,ihe
      do k=1,nvl
        trclo_tdcy(k,ipn,nf,t) = 0.
        trc_tdcy  (k,ipn,nf,t) = 0.

        ! Initialize these local arrays so COMPARE_VAR does not get 
        ! confused by uninitialized edg=6 edges of pentagonal grid cells.  This 
        ! code could be safely omitted if COMPARE_VAR were not used.  
        antiflx(k,npp,ipn) = 0.
      end do

      do edgcount=1,nedge(ipn)
        edg = permedge(edgcount,ipn)
        ipx = prox(edg,ipn)		! neighbor across edge
        edx = proxs(edg,ipn)		! index of joint edge as seen by nghbor
        do k=1,nvl
! --- high-order tracer flux (2nd order centered, out of cell > 0):
          flx_hi = massfx(k,edg,ipn,nf)*trcr_edg(k,edg,ipn,t)	! trcdim x N/s

! --- low-order tracer flux (donor-cell, out of cell > 0):
          flx_lo(k,edg,ipn) =					&
            0.5*((massfx(k,edg,ipn,nf) +			&
              abs(massfx(k,edg,ipn,nf)))*tracr(k,ipn,t)		&
               - (massfx(k,edx,ipx,nf) +			&
              abs(massfx(k,edx,ipx,nf)))*tracr(k,ipx,t))

          ! get anti-diffusive flux as the difference between high order
          ! (flx_hi) and low order flux (flx_lo), from Zalesak, p336, Eqn (3):
          antiflx(k,edg,ipn) = flx_hi - flx_lo(k,edg,ipn)   ! trcdim x N/sec
        end do
      end do
! End of halo comp
      if (ipn <= ipe) then
        do edg=1,nprox(ipn)
          do k=1,nvl
            ! sum up edge fluxes to get low-order tracer tendency
            trclo_tdcy(k,ipn,nf,t) = trclo_tdcy(k,ipn,nf,t) + flx_lo(k,edg,ipn)
          end do
        end do
      
      !.....................................................................
      !  Sec. 2.  Calculate new low order trcr*dp using full Adams Bashforth
      !.....................................................................

        do k=1,nvl
          !  divide tendency by cell area to convert to (trcdim x Pa/sec)
          trclo_tdcy(k,ipn,nf,t) = -trclo_tdcy(k,ipn,nf,t)*rarea(ipn)
          
          ! get new value for the low order tracer*dp field using Adams Bashforth
          ! 3 time levels, the one just calculated (nf), and the two prev ones,
          ! marked 'of' (old field) and 'vof' (very old field):
          
          trcdp_lo = trcdp(k,ipn,t)			&
               +adbash1*trclo_tdcy(k,ipn, nf,t)		&
               +adbash2*trclo_tdcy(k,ipn, of,t)		&
               +adbash3*trclo_tdcy(k,ipn,vof,t)

          ! initiate setting upper/lower bounds for flux clipping
          trmax(k,ipn) = max(tracr(k,ipn,t),			&
               tracr(k,prox(1,ipn),t),tracr(k,prox(2,ipn),t),	&
               tracr(k,prox(3,ipn),t),tracr(k,prox(4,ipn),t),	&
               tracr(k,prox(5,ipn),t),tracr(k,prox(6,ipn),t)  )

          trmin(k,ipn) = min(tracr(k,ipn,t),			&
               tracr(k,prox(1,ipn),t),tracr(k,prox(2,ipn),t),	&
               tracr(k,prox(3,ipn),t),tracr(k,prox(4,ipn),t),	&
               tracr(k,prox(5,ipn),t),tracr(k,prox(6,ipn),t)  )

          trmax(k,ipn) = max(0.,trmax(k,ipn))	! cannot allow negative trmax
          trmin(k,ipn) = max(0.,trmin(k,ipn))	! cannot allow negative trmin	

          ! now divide trcdp_lo by delp to get new low order tracer field 
          ! (low-ord solution is diffusive and hence will not exceed range of
          ! neighboring values; thus, they can be used as bounds when dividing.)
          trcr_lo(k,ipn) = max(trmin(k,ipn),			&
                           min(trmax(k,ipn),			&
                           trcdp_lo/max(thshld,delp(k,ipn))))
        end do			! vertical loop
      end if
    end do			! horiz. loop
!$OMP END PARALLEL DO
!sms$ignore end

    return

  end subroutine trcadv1

  subroutine trcadv2 (t, trc_tdcy, delp, antiflx, trmin, &
                      trmax, trcr_lo, r_plus, r_mnus)
    use module_variables, only: nf, of, vof, adbash1, adbash2, adbash3
    use module_control,   only: npp, nabl, ntra
    use fimnamelist,      only: nvl
    use module_constants, only: nprox, prox, area
    use global_bounds,    only: ims, ime, ips, ipe

! Arguments
    integer, intent(in) :: t                                ! tracer index
    real, intent(in)    :: trc_tdcy(nvl,ims:ime,nabl,ntra)  ! forcing
    real, intent(in)    :: delp    (nvl,ims:ime)            ! layer thickness
    real, intent(in)    :: antiflx(nvl,npp,ims:ime)	! antidiff tracer flux
    real, intent(inout) :: trmax(nvl,ims:ime)	! regional trcr max for flux clipping
    real, intent(inout) :: trmin(nvl,ims:ime)	! regional trcr min for flux clipping
    real, intent(in)    :: trcr_lo(nvl,ims:ime)	! tracer after low-ord transport
    real, intent(out)   :: r_plus(nvl,ims:ime)	! Zalesak's r_plus
    real, intent(out)   :: r_mnus(nvl,ims:ime)	! Zalesak's r_minus 

! Local workspace
    integer :: ipn
    integer :: k
    integer :: edg
    real :: q_plus	  ! Zalesak's q_plus
    real :: q_mnus	  ! Zalesak's q_minus
    real :: s_plus(nvl)   ! Zalesak's p_plus
    real :: s_mnus(nvl)   ! Zalesak's p_minus

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (k,edg,s_plus,s_mnus,q_plus,q_mnus) SCHEDULE (runtime)
    do ipn=ips,ipe
      do k=1,nvl
        s_plus(k) = 0.
        s_mnus(k) = 0.
      end do
      do edg=1,nprox(ipn)
        do k=1,nvl
          s_plus(k) = s_plus(k) - min(0., antiflx(k,edg,ipn))
          s_mnus(k) = s_mnus(k) + max(0., antiflx(k,edg,ipn))
        end do
      end do
      
      !............................................................
      !  Sec 3.  Monotonicity Limit on Fluxes
      !............................................................
      
      !  Determine the amount of antidiffusive fluxes that can be
      !  added to the low order solution without violating monotonicity.

      do k=1,nvl
        ! According to Zal (17), we must limit according to max of
        ! tracer from any gazebo direction, in either current or prev time
        ! step.  Likewise for the minimum.
        ! For the pentagons prox(6,ipn) is set to prox(5,ipn) in init.F90.
        
        ! relax upper/lower bounds by incorporating new low-order field
        trmax(k,ipn) = max(trmax(k,ipn),trcr_lo(k,ipn),		&
             trcr_lo(k,prox(1,ipn)),trcr_lo(k,prox(2,ipn)),	&
             trcr_lo(k,prox(3,ipn)),trcr_lo(k,prox(4,ipn)),	&
             trcr_lo(k,prox(5,ipn)),trcr_lo(k,prox(6,ipn))  )
          
        trmin(k,ipn) = min(trmin(k,ipn),trcr_lo(k,ipn),		&
             trcr_lo(k,prox(1,ipn)),trcr_lo(k,prox(2,ipn)),	&
             trcr_lo(k,prox(3,ipn)),trcr_lo(k,prox(4,ipn)),	&
             trcr_lo(k,prox(5,ipn)),trcr_lo(k,prox(6,ipn))  )
          
        trmax(k,ipn) = max(0.,trmax(k,ipn))	! cannot allow negative trmax
        trmin(k,ipn) = max(0.,trmin(k,ipn))	! cannot allow negative trmin	
          
        ! q_plus/q_mnus are the upper/lower limits on antidiff trcr tendencies,
        ! dimensioned trdim x N/sec

        ! q_plus is Zalesak (8):
        q_plus = ((trmax(k,ipn) - trcr_lo(k,ipn))*delp(k,ipn)	&
             -(adbash2*min(0.,trc_tdcy(k,ipn, of,t))		&
             + adbash3*max(0.,trc_tdcy(k,ipn,vof,t))))		&
             /adbash1*(area(ipn))
          
        ! q_mnus is Zalesak (11):
        q_mnus = ((trcr_lo(k,ipn) - trmin(k,ipn))*delp(k,ipn)	&
             +(adbash2*max(0.,trc_tdcy(k,ipn, of,t))		&
             + adbash3*min(0.,trc_tdcy(k,ipn,vof,t))))		&
             /adbash1*(area(ipn))
          
        ! Having s_plus(k,ipn) and q_plus, we can calc r_plus, Zal (9):

        !  reduce fluxes to stay within limits posed by q_plus/q_mnus
        !  r_plus/r_mnus are dimensionless

        r_plus(k,ipn) = max (0., min(1., q_plus / max (s_plus(k), thshld)))
        r_mnus(k,ipn) = max (0., min(1., q_mnus / max (s_mnus(k), thshld)))
      end do
    end do
!$OMP END PARALLEL DO
!sms$ignore end

    return

  end subroutine trcadv2

  subroutine trcadv3 (t, trc_tdcy, trclo_tdcy, trcdp, tracr,	&
                      delp, antiflx, r_mnus, r_plus, trmin,trmax)
    use module_control,   only: npp, nabl, ntra
    use fimnamelist,      only: nvl
    use module_variables, only: nf, of, vof, adbash1, adbash2, adbash3
    use module_constants, only: nprox, prox, rarea
    use global_bounds,    only: ims, ime, ips, ipe
! Arguments
    integer, intent(in) :: t                                  ! thread index
    real, intent(inout) :: trc_tdcy  (nvl,ims:ime,nabl,ntra)  ! forcing
    real, intent(in)    :: trclo_tdcy(nvl,ims:ime,nabl,ntra)  ! low order tdcy 
    real, intent(inout) :: trcdp     (nvl,ims:ime,ntra)       ! tracr*dp
    real, intent(inout) :: tracr     (nvl,ims:ime,ntra)       ! tracer
    real, intent(in)    :: delp      (nvl,ims:ime)            ! layer thickness
    real, intent(inout) :: antiflx(nvl,npp,ims:ime) ! antidiffusive tracer flux
    real, intent(in) :: r_plus(nvl,ims:ime)         ! Zalesak's r_plus
    real, intent(in) :: r_mnus(nvl,ims:ime)         ! Zalesak's r_minus 
    real, intent(in) :: trmax(nvl,ims:ime)	    ! regional tracer max for flux clipping
    real, intent(in) :: trmin(nvl,ims:ime)	    ! regional tracer min for flux clipping

! Local workspace
    integer :: ipn
    integer :: k
    integer :: edg
    real :: anti_tdcy(nvl)          ! antidiff trcr tendency

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (edg,k,anti_tdcy) SCHEDULE (runtime)
    do ipn=ips,ipe
      do k=1,nvl
        anti_tdcy(k) = 0.
      end do
      
          ! clip antidiffusive fluxes
      do edg=1,nprox(ipn)
        do k=1,nvl
          if (antiflx(k,edg,ipn) >= 0.) then	! outgoing
            antiflx(k,edg,ipn) = antiflx(k,edg,ipn)*		&
              min(r_mnus(k,ipn),r_plus(k,prox(edg,ipn)))
          else					! incoming
            antiflx(k,edg,ipn) = antiflx(k,edg,ipn)*		&
              min(r_plus(k,ipn),r_mnus(k,prox(edg,ipn)))
          end if
          ! sum up edge fluxes to get antidiff tracer tendency
          anti_tdcy(k) = anti_tdcy(k) + antiflx(k,edg,ipn)
        end do
      end do

      do k=1,nvl
        ! divide tendency by cell area to convert to  trcdim x Pa/sec
        anti_tdcy(k) = -anti_tdcy(k)*rarea(ipn)

        ! combine low order with clipped antidiff tendency
        trc_tdcy(k,ipn,nf,t) = trclo_tdcy(k,ipn,nf,t) + anti_tdcy(k)
        
        ! advance (tracer x thickness) field in time
        trcdp(k,ipn,t) = trcdp(k,ipn,t)		&
             + adbash1*trc_tdcy(k,ipn, nf,t)	&
             + adbash2*trc_tdcy(k,ipn, of,t)	&
             + adbash3*trc_tdcy(k,ipn,vof,t)
          
        ! finally divide new trcr*dp by new delp to get new tracer field 
        tracr(k,ipn,t) = max(trmin(k,ipn),			&
                         min(trmax(k,ipn),			&
                         trcdp(k,ipn,t)/max(thshld,delp(k,ipn))))
      end do
    end do                  ! ipn
!$OMP END PARALLEL DO
!sms$ignore end

    return

  end subroutine trcadv3

end module trcadv_compute
