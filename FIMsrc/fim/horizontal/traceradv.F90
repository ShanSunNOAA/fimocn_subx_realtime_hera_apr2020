module traceradv_driver
  implicit none

contains
!*********************************************************************
!     traceradv
!
!     A new tracer transport scheme based on flux-limit scheme of
!     Thuburn 1995, 1996 and Miura 2012.
!    
!     The implementation of this module, its interface and code
!     structure, is partially adapted from trcadv.F90.
!
!     N. Wang, Jan 2016. 
!
!*********************************************************************

  subroutine traceradv (its, u_edg, v_edg, trcr_edg, tracr, trcdp, trc_tdcy, &
                        massfx, delp, dp_edg, TimingBarriers )
    use stencilprint,      only: stencl
    use stenedgprint,      only: stenedg
    use module_control,    only: npp, nip, nabl, ntra
    use fimnamelist,       only: nvl, PrintIpnDiag
    use traceradv_compute, only: tracer_tendency, tracer_update
    use fluxlimiter,       only: init_tracer_limiter, limiter, c_num_t
    use findmaxmin2
    use findmaxmin3

!..............................................................
!	variable declarations
!..............................................................

! Arguments
    integer, intent(in) :: its                             ! model time step
    logical, intent(in) :: TimingBarriers                  ! measure task skew when .true.
!sms$distribute (dh,2) begin
    real, intent(inout) :: tracr     (nvl,nip,ntra)       ! tracer
    real, intent(inout) :: trcdp     (nvl,nip,ntra)       ! tracr*dp
    real, intent(inout) :: trc_tdcy  (nvl,nip,nabl,ntra)  ! forcing
    real, intent(in)    :: delp      (nvl,nip)            ! layer thickness
    real, intent(in)    :: dp_edg    (nvl,npp,nip)            ! layer thickness
    real, intent(in)    :: u_edg     (nvl,npp,nip)            ! layer thickness
    real, intent(in)    :: v_edg     (nvl,npp,nip)            ! layer thickness
    real :: trmax(nvl,nip)		! regional tracer max for flux clipping
    real :: trmin(nvl,nip)		! regional tracer min for flux clipping
    real :: work (nvl,nip)
!sms$distribute end
!sms$distribute (dh,3) begin
    real, intent(inout) :: trcr_edg  (nvl,npp,nip,ntra)   ! tracer on edges
    real, intent(in)    :: massfx    (nvl,npp,nip,3)      ! mass flux across edges
!sms$distribute end

    character(len=48)  :: string

    integer :: k	! layer index
    integer :: t	! tracer index (1=theta; 2=specif.hum., ...)
    integer :: ret      ! return code from gptl routines

    real    :: del_dp	! delta_p used for upstream integral of flux

#include <gptl.inc>

!.............................................................
!   Calculate and enforec the flux limits for all tracers 
!.............................................................

!sms$compare_var(tracr  , "trcadv.F90 - tracr1   ")
!sms$compare_var(massfx , "trcadv.F90 - massfx1   ")
!sms$compare_var(trc_tdcy, "trcadv.F90 - trc_tdcy1  ")
!sms$compare_var(trcr_edg, "trcadv.F90 - trcr_edg1  ")

#ifdef DEBUGPRINT
    do k=1,nvl,7
     write (string,'(a,i3,a)') 'k',k,' trcadv:theta'
     call findmxmn3(tracr,nvl,nip,ntra,k,1,string)
    end do
    print *
#endif

! initialize mass flux limiter 
    call init_tracer_limiter(u_edg, v_edg, delp, dp_edg)

    do t=2,ntra                   ! loop through tracers

      if (PrintIpnDiag > 0 .and. t.eq.1) then
        write (string,'(a,i5,a,i2)') '(trcadv) step',its,'  tracer',t
        work(:,:)=tracr(:,:,t)
        call stencl(work,nvl,1.,trim(string)//' (old)')
      end if

      ret = gptlstart ('trcadv_init_limiter' )

! clip tracer fluxes if necessary
      call limiter(tracr(:,:,t),trcr_edg(:,:,:,t),c_num_t,trc_tdcy(:,:,:,t))
      ret = gptlstop ('trcadv_init_limiter' )

      if (TimingBarriers) then
        ret = gptlstart ('trcadv_barrier')
!SMS$BARRIER
        ret = gptlstop ('trcadv_barrier')
      end if

! compute tracer tendiency 
      call tracer_tendency(t, trcr_edg, trc_tdcy, massfx)

! computes high-order tracer at a new time step
      call tracer_update (t, trc_tdcy, trcdp, tracr,	&
                      delp, trmin, trmax)

    end do                	! loop through tracers

    return
  end subroutine traceradv
end module traceradv_driver
