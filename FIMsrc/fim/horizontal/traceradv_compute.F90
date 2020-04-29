!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines in this module are computational only. They contain no SMS directives
! except for IGNORE. Thus they are good candidates for migration to MIC or GPU.
!
! trcadv_compute: contains all the computational routines for tracer advection
! calculations.
!
! Ning Wang, Jan 2016 - original version
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module traceradv_compute
  implicit none

  private
  real, parameter :: thshld = 1.e-11

  public :: tracer_tendency, tracer_update

contains

  subroutine tracer_tendency (t, trc_tdcy, tracer_edge, massfx)
    use module_variables, only: nf
    use module_control,   only: npp, nabl, ntra
    use fimnamelist,      only: nvl
    use module_constants, only: nprox, prox, rarea
    use global_bounds,    only: ims, ime, ips, ipe

! Arguments
    integer, intent(in) :: t                                ! tracer index
    real, intent(out)   :: trc_tdcy(nvl,ims:ime,nabl,ntra)  ! forcing
    real, intent(in)    :: tracer_edge(nvl,npp,ims:ime,ntra)  ! tracer on edge
    real,intent (in)    :: massfx (nvl,npp,ims:ime,3)	! mass flux across edges (N/sec)

! Local workspace
    integer :: ipn
    integer :: k
    integer :: edg

!$OMP PARALLEL DO PRIVATE (k,edg) SCHEDULE (runtime)
    do ipn=ips,ipe
      do k=1,nvl
        trc_tdcy (k,ipn,nf,t) = 0.
      end do

      do edg=1,nprox(ipn)
        do k=1,nvl
          ! sum up edge mass fluxes and tracer values to get tracer tendency
          trc_tdcy(k,ipn,nf,t) = trc_tdcy(k,ipn,nf,t) + massfx(k,edg,ipn,nf) * &
                                 tracer_edge(k,edg,ipn,t)
        end do
      end do

      do k=1,nvl
        trc_tdcy (k,ipn,nf,t) = -trc_tdcy(k,ipn,nf,t) * rarea(ipn)
      end do
    end do
!$OMP END PARALLEL DO

  end subroutine tracer_tendency

  subroutine tracer_update (t, trc_tdcy, trcdp, tracr,	&
                      delp, trmin, trmax)
    use module_control,   only: npp, nabl, ntra
    use fimnamelist,      only: nvl
    use module_variables, only: nf, of, vof, adbash1, adbash2, adbash3
    use module_constants, only: nprox, prox, rarea
    use global_bounds,    only: ims, ime, ips, ipe
! Arguments
    integer, intent(in) :: t                                  ! thread index
    real, intent(inout) :: trc_tdcy  (nvl,ims:ime,nabl,ntra)  ! forcing
    real, intent(inout) :: trcdp     (nvl,ims:ime,ntra)       ! tracr*dp
    real, intent(inout) :: tracr     (nvl,ims:ime,ntra)       ! tracer
    real, intent(in)    :: delp      (nvl,ims:ime)            ! layer thickness
    real, intent(inout) :: trmax(nvl,ims:ime)	    ! regional tracer max for flux clipping
    real, intent(inout) :: trmin(nvl,ims:ime)	    ! regional tracer min for flux clipping

! Local workspace
    integer :: ipn
    integer :: k
    integer :: edg
    real :: anti_tdcy(nvl)          ! antidiff trcr tendency

!$OMP PARALLEL DO PRIVATE (k) SCHEDULE (runtime)
    do ipn=ips,ipe
      do k=1,nvl
          trmax(k,ipn) = max(tracr(k,ipn,t),                    &
               tracr(k,prox(1,ipn),t),tracr(k,prox(2,ipn),t),   &
               tracr(k,prox(3,ipn),t),tracr(k,prox(4,ipn),t),   &
               tracr(k,prox(5,ipn),t),tracr(k,prox(6,ipn),t)  )

          trmin(k,ipn) = min(tracr(k,ipn,t),                    &
               tracr(k,prox(1,ipn),t),tracr(k,prox(2,ipn),t),   &
               tracr(k,prox(3,ipn),t),tracr(k,prox(4,ipn),t),   &
               tracr(k,prox(5,ipn),t),tracr(k,prox(6,ipn),t)  )
      end do
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE (k) SCHEDULE (runtime)
    do ipn=ips,ipe
      do k=1,nvl
        ! advance (tracer x thickness) field in time
        trcdp(k,ipn,t) = trcdp(k,ipn,t)		&
             + adbash1*trc_tdcy(k,ipn, nf,t)	&
             + adbash2*trc_tdcy(k,ipn, of,t)	&
             + adbash3*trc_tdcy(k,ipn,vof,t)
          
        ! finally divide new trcr*dp by new delp to get new tracer field 
        tracr(k,ipn,t) = max(trmin(k,ipn),			&
                         min(trmax(k,ipn),			&
                         trcdp(k,ipn,t)/max(thshld,delp(k,ipn))))
!        tracr(k,ipn,t) = trcdp(k,ipn,t)/max(thshld,delp(k,ipn))
      end do
    end do                  ! ipn
!$OMP END PARALLEL DO

    return

  end subroutine tracer_update

end module traceradv_compute
