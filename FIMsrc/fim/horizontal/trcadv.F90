module trcadv_driver
  implicit none

contains
!*********************************************************************
!     trcadv
! 	trcadv = flux corrected transport for mass field tracers
!	J. Lee		Author				September 2005
!	A. E. MacDonald Documentor			November  2005
!	R. Bleck	major rewrite			April     2008
!       R. Bleck        fixed bug in high-order flux	April     2011
!
! 	This routine is based on Zalesak, JOURNAL OF COMPUTATIONAL
!       PHYSICS, 31, 335-362, 1979.  Dale Durran provides an
!       excellent discussion of flux corrected transport in his book
!     NUMERICAL METHODS FOR WAVE EQUATIONS IN GEOPHYSICAL FLUID DYNAMICS.
!*********************************************************************

  subroutine trcadv (its, trcr_edg, tracr, trcdp, trc_tdcy, &
                     trclo_tdcy, massfx, delp, TimingBarriers )
    use stencilprint,   only: stencl
    use stenedgprint,   only: stenedg
    use module_control, only: npp, nip, nabl, ntra
    use fimnamelist,    only: nvl, PrintIpnDiag
    use trcadv_compute, only: trcadv1, trcadv2, trcadv3
    use findmaxmin2
    use findmaxmin3

!..............................................................
!	Sec. 0  Dimension and Type
!..............................................................

! Arguments
    integer, intent(in) :: its                             ! model time step
    logical, intent(in) :: TimingBarriers                  ! measure task skew when .true.
!sms$distribute (dh,2) begin
    real, intent(inout) :: tracr     (nvl,nip,ntra)       ! tracer
    real, intent(inout) :: trcdp     (nvl,nip,ntra)       ! tracr*dp
    real, intent(inout) :: trc_tdcy  (nvl,nip,nabl,ntra)  ! forcing
    real, intent(inout) :: trclo_tdcy(nvl,nip,nabl,ntra)  ! low order forcing
    real, intent(in)    :: delp      (nvl,nip)            ! layer thickness
    real :: r_plus(nvl,nip)		! Zalesak's r_plus
    real :: r_mnus(nvl,nip)		! Zalesak's r_minus 
    real :: trcr_lo(nvl,nip)		! tracer after low-order transport
    real :: trmax(nvl,nip)		! regional tracer max for flux clipping
    real :: trmin(nvl,nip)		! regional tracer min for flux clipping
    real :: work (nvl,nip)
!sms$distribute end
!sms$distribute (dh,3) begin
    real, intent(in)    :: trcr_edg  (nvl,npp,nip,ntra)   ! tracer on edges
    real, intent(in)    :: massfx    (nvl,npp,nip,3)      ! mass flux across edges
    real :: flx_lo (nvl,npp,nip)	! low-order tracer flux
    real :: antiflx(nvl,npp,nip)	! antidiffusive tracer flux
!sms$distribute end

    logical, parameter :: low_ord = .false. ! if true, skip antidiffusive part
    character(len=48)  :: string

    integer :: k	! layer index
    integer :: t	! tracer index (1=theta; 2=specif.hum., ...)
    integer :: ret      ! return code from gptl routines

    real    :: del_dp	! delta_p used for upstream integral of flux

#include <gptl.inc>

!.............................................................
!  Sec. 1 Calculate Low and Antidiffusive Flux
!.............................................................

! Calculates the FCT low and high order fluxes, and defines
! the antidiffusive flux as the difference between the high and
! low order fluxes.  The low order flux is computed based on 
! the assumption of piecewise continuity, with a constant value 
! in each cell.  The higher order uses the "gazebo" with a
! sloped line used for the flux integral.

! Avoid exchange via HALO_COMP in cnuity.F90 and edgvar.F90
!!!SMS$EXCHANGE(massfx,trcr_edg)
!!SMS$EXCHANGE(tracr) exchanged in edgvar

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

    do t=1,ntra                   ! loop through tracers

      if (PrintIpnDiag > 0 .and. t.eq.1) then
        write (string,'(a,i5,a,i2)') '(trcadv) step',its,'  tracer',t
        work(:,:)=tracr(:,:,t)
        call stencl(work,nvl,1.,trim(string)//' (old)')
      end if

      ret = gptlstart ('trcadv' )
      ret = gptlstart ('trcadv1')
! --- trcadv1 computes low-ord & antidff fluxes, low-ord tracer at new time step
      call trcadv1 (t, trclo_tdcy, trc_tdcy, massfx, trcr_edg,	&
                    tracr, trcdp, antiflx, trmin, trmax,	&
                    trcr_lo, flx_lo, delp)
      ret = gptlstop ('trcadv1')
      ret = gptlstop ('trcadv' )

    if (PrintIpnDiag > 0 .and. t.eq.1) then
      call stenedg(work,flx_lo,nvl,trim(string)//' (old), low-ord flx')
      call stencl(trcr_lo,nvl,1.,trim(string)//' (low-ord)')
    end if

!JR Dont thread this because of small amount of work and low_ord is false
      if (low_ord) then
        tracr(:,:,t) = trcr_lo(:,:)
      else				! evaluate and apply antidiffusive fluxes

        ! Dale Durran's book indicates that condition Zal (14) should be
        ! satisfied, although Zalesak says it's cosmetic.  We believe Dale:
        ! Also Calculate s_plus Zalesak (7) and s_minus Zalesak (10).
        ! (aggregate of incoming and outgoing fluxes, trcdim x N/sec)

        if (TimingBarriers) then
          ret = gptlstart('trcadv_barrier')
!SMS$BARRIER
          ret = gptlstop ('trcadv_barrier')
        end if

        ret = gptlstart ('trcadv_exchange')
!SMS$EXCHANGE(trcr_lo)
        ret = gptlstop  ('trcadv_exchange')

        ret = gptlstart ('trcadv' )
        ret = gptlstart ('trcadv2')
! --- trcadv2 computes flux limiters
        call trcadv2 (t, trc_tdcy, delp, antiflx, trmin, &
                      trmax, trcr_lo, r_plus, r_mnus)
        ret = gptlstop  ('trcadv2')
        ret = gptlstop  ('trcadv' )

#ifdef DEBUGPRINT
   do k=1,nvl,7
    write (string,'(a3,i3)') ' k=',k
    call findmxmn2(r_plus,nvl,nip,k,'(trcadv) r_plus'//string(1:6))
    call findmxmn2(r_mnus,nvl,nip,k,'(trcadv) r_mnus'//string(1:6))
   end do
   if (PrintIpnDiag > 0) then
     call stencl(r_plus,nvl,1000.,'(trcadv) r_plus x 1000')
     call stencl(r_mnus,nvl,1000.,'(trcadv) r_mnus x 1000')
   end if
#endif

        ! Now we have a flux limiter that will assure monotonicity.
        ! Next, we do the clipping.

        !.......................................................
        !  Sec. 4. Flux Clipping
        !.......................................................
        
        ! As explained by Durran, once you have the r_plus and
        ! r_mnus over the whole grid, you can assure that the clipping
        ! can be done so that it does not cause a problem in the center
        ! cell, NOR IN ANY OF THE NEIGHBORING CELLS.  The clipping
        ! is from Zalesak (13):

        if (TimingBarriers) then
          ret = gptlstart ('trcadv_barrier')
!SMS$BARRIER
          ret = gptlstop ('trcadv_barrier')
        end if

        ret = gptlstart ('trcadv_exchange')
!SMS$EXCHANGE(r_plus,r_mnus)
        ret = gptlstop ('trcadv_exchange')

!sms$compare_var(r_plus   , "trcadv.F90 - r_plus4 ")
!sms$compare_var(r_mnus   , "trcadv.F90 - r_mnus4 ")
!sms$compare_var(antiflx  , "trcadv.F90 - antiflx4 ")
!sms$compare_var(trcr_lo  , "trcadv.F90 - trcr_lo3 ")

      if (PrintIpnDiag > 0 .and. t.eq.1) then
        call stenedg(trcr_lo,antiflx,nvl,			&
                  trim(string)//' (low-ord), unclipped anti.flx')
      end if

        ret = gptlstart ('trcadv' )
        ret = gptlstart ('trcadv3')
! --- trcadv3 clips antidiff fluxes, computes high-ord tracer at new time step
        call trcadv3 (t, trc_tdcy, trclo_tdcy, trcdp, tracr,	&
                      delp, antiflx, r_mnus, r_plus, trmin,	&
                      trmax)
        ret = gptlstop  ('trcadv3')
        ret = gptlstop  ('trcadv' )
      end if			! low_ord = true or false

      if (PrintIpnDiag > 0 .and. t.eq.1) then
        work(:,:)=tracr(:,:,t)
        call stenedg(work,antiflx,nvl,				&
                  trim(string)//' (new), clipped antidiff flx')
        call stencl(work,nvl,1.,trim(string)//' (new)')
      end if

    end do                	! loop through tracers
!sms$compare_var(antiflx   , "trcadv.F90 - antiflx5 ")
!sms$compare_var(trcdp     , "trcadv.F90 - trcdp5 ")
!sms$compare_var(trclo_tdcy, "trcadv.F90 - trclo_tdcy5  ")

    return
  end subroutine trcadv
end module trcadv_driver
