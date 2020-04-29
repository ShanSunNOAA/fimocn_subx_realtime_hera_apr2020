module cnuity_driver
  implicit none

  private
  public :: cnuity

contains
!*********************************************************************
!     cnuity
! 	based on fct = flux corrected transport
!	J. Lee		Author				September, 2005
!	A. E. MacDonald Documentor			November,  2005
!	R. Bleck	major rewrite			April      2008
!	R. Bleck	revised omega diagnostics	August     2009
!	R. Bleck	discarded-flux recovery		November   2009
!
! 	This routine is based on Zalesak, JOURNAL OF COMPUTATIONAL
!       PHYSICS, 31, 335-362, 1979.  Dale Durran provides an
!       excellent discussion of flux corrected transport in his book
!     NUMERICAL METHODS FOR WAVE EQUATIONS IN GEOPHYSICAL FLUID DYNAMICS.
!*********************************************************************

  subroutine cnuity (its, u_velo, v_velo, u_edg, v_edg, &
                     dp_edg, lp_edg, delp, pres, exner, &
                     dp_tndcy, dplo_tndcy, massfx, omega, TimingBarriers)

    use module_control,   only: nvlp1, npp, nip, nabl
    use fimnamelist,      only: nvl, PrintIpnDiag
    use module_variables, only: nf
    use module_constants, only: sidevec_e
    use findmaxmin2
    use stencilprint
    use stenedgprint
    use cnuity_compute, only: cnuity1, cnuity2, cnuity3, cnuity4, cnuity5
!..............................................................
!	Sec. 0  Dimension and Type
!..............................................................

! Arguments
    integer,intent (IN)    :: its              ! model time step
    logical,intent (IN)    :: TimingBarriers   ! measure task skew when .true.
!sms$distribute (dh,1) begin
    real    :: psurf(nip)		   ! surface pressure
!sms$distribute end
!sms$distribute (dh,2) begin
    real   ,intent (IN)    :: u_velo    (nvl,nip)	! west wind (m/s)
    real   ,intent (IN)    :: v_velo    (nvl,nip)	! south wind (m/s)
    real   ,intent (INOUT) :: delp      (nvl,nip)	! layer thickness (Pa)
    real   ,intent (INOUT) :: pres      (nvlp1,nip)	! pressure on interfaces (Pa)
    real   ,intent (INOUT) :: exner     (nvlp1,nip)	! Exner function
    real   ,intent (INOUT) :: dp_tndcy  (nvl,nip,nabl)	! Forcing for dp (Pa/sec)
    real   ,intent (INOUT) :: dplo_tndcy(nvl,nip,nabl)	! Forcing for low order dp (Pa/sec)
    real   ,intent (OUT)   :: omega     (nvl,nip)	! dp/dt (N/sec)
    real    :: p_plus(nvl,nip)	           ! Zalesak's p_plus, N/sec
    real    :: p_mnus(nvl,nip)	           ! Zalesak's p_minus, N/sec
    real    :: r_plus(nvl,nip)	           ! Zalesak's r_plus, dimensionless
    real    :: r_mnus(nvl,nip)	           ! Zalesak's r_minus, dimensionless
    real    :: delp_lo(nvl,nip)		   ! Pa
    real    :: anti_tndcy(nvl,nip)	   ! Pa/sec
!sms$distribute end
!sms$distribute (dh,3) begin
    real   ,intent (IN)    :: u_edg     (nvl,npp,nip)	! u on edges (m/s)
    real   ,intent (IN)    :: v_edg     (nvl,npp,nip)	! v on edges (m/s)
    real   ,intent (IN)    :: dp_edg    (nvl,npp,nip)	! delta-p on edges (Pa)
    real   ,intent (IN)    :: lp_edg    (nvl,npp,nip)	! mid-layer pressure on edges (Pa)
    real   ,intent (INOUT) :: massfx    (nvl,npp,nip,3)	! mass flux across edges (N/sec)
    real    :: antifx(nvl,npp,nip)	   ! N/sec
    real    :: work(nvl,npp,nip)
!sms$distribute end

    integer :: k	 		   ! layer index
    integer :: ipn		           ! Index for icos cell number
    character :: string*32
    logical,parameter :: low_ord = .false. ! if true, skip antidiffusive part
    integer :: ret

#include <gptl.inc>

!.............................................................
!  Sec. 1 Calculate Anti Diffusive Flux, Low order forcing
!.............................................................

! Calculates the FCT low and high order fluxes, and defines
! the antidiffusive flux as the difference between the high and
! low-order fluxes.  The low order flux is computed based on 
! the assumption of piecewise continuity, with a constant value 
! in each cell.  The higher order uses the "gazebo" with a
! sloped line used for the flux integral.

#ifdef DEBUGPRINT
   do k=1,nvl,7
    write (string,'(a,i2)') 'cnuity: old dp k=',k
    call findmxmn2(delp,nvl,nip,k,string)
   end do
#endif

!sms$compare_var(sidevec_e, "cnuity.F90 - sidevec_e6 ")
!sms$compare_var(u_edg    , "cnuity.F90 - u_edg6 ")
!sms$compare_var(v_edg    , "cnuity.F90 - v_edg6 ")
!sms$compare_var(dp_edg   , "cnuity.F90 - dp_edg6 ")
!sms$compare_var(delp     , "cnuity.F90 - delp6 ")

    ret = gptlstart ('cnuity')

    ret = gptlstart ('cnuity1')
! --- cnuity1 computes low-order and antidiffusive fluxes
    call cnuity1 (antifx, massfx, u_edg, v_edg, dp_edg, delp)
    ret = gptlstop ('cnuity1')
!sms$compare_var(antifx, "cnuity.F90 - antifx7 ")

    if (PrintIpnDiag > 0) then
      write (string,'(a,i6,2x)') 'step',its
!JR Could merge cnuity1 and cnuity2 if this stencl call could be eliminated
      call stencl(delp,nvl,.01,string(1:12)//'(cnuity) old dp')
      work(:,:,:)=massfx(:,:,:,nf)
      call stenedg(delp,work,nvl,				&
                   string(1:12)//'(cnuity) old dp, low-ord flx')
    end if

    ret = gptlstart ('cnuity2')
! --- cnuity2 computes low-order dp at new time step
    call cnuity2 (dplo_tndcy, massfx, delp, delp_lo)
    ret = gptlstop ('cnuity2')

!JR Dont thread this because of small amount of work and low_ord is false
    if (low_ord) then		! use low-order mass fluxes only
!sms$parallel (dh,ipn) begin
      do ipn=1,nip
        do k=1,nvl
          delp(k,ipn) = delp_lo(k,ipn)
          dp_tndcy(k,ipn,nf) = dplo_tndcy(k,ipn,nf)
        end do
      end do
!sms$parallel end
    else               ! evaluate and apply antidiffusive fluxes

! Dale Durrans book indicates that condition Zal (14) should be
! satisfied, although Zalesak says its cosmetic.  We believe Dale:
! Also calculate p_plus Zalesak (7) and p_minus Zalesak (10)
! (aggregate of incoming and outgoing fluxes, N/sec) 

      ret = gptlstop ('cnuity')

      if (TimingBarriers) then
        ret = gptlstart ('cnuity_barrier')
!SMS$BARRIER
        ret = gptlstop ('cnuity_barrier')
      end if

      ret = gptlstart ('cnuity_exchange')
!SMS$EXCHANGE(delp_lo)
!!!SMS$EXCHANGE(delp) exchanged in edgvar
      ret = gptlstop ('cnuity_exchange')

!sms$compare_var(antifx, "cnuity.F90 - antifx8 ")
!sms$compare_var(delp_lo,"cnuity.F90 - delp_lo8 ")

      ret = gptlstart ('cnuity')
      ret = gptlstart ('cnuity3')

! --- cnuity3 computes limiters for antidiffusive fluxes
      call cnuity3 (p_plus, p_mnus, antifx, delp_lo, delp, dp_tndcy, r_plus, r_mnus)
      ret = gptlstop ('cnuity3')
!
!   Now we have chosen flux limiters that will assure monotonicity.
!   Next, we do the clipping.

!.......................................................
!   Sec. 4. Flux Clipping
!.......................................................

!   As explained by Durran, once you have the r_plus and
!   r_mnus over the whole grid, you can assure that the clipping
!   can be done so that it does not cause a problem in the center
!   cell, NOR IN ANY OF THE NEIGHBORING CELLS.  The clipping
!   is from Zalesek (13):
    
      ret = gptlstop ('cnuity')

      if (PrintIpnDiag > 0) then
        write (string,'(a,i6,2x)') 'step',its                        
        call stenedg (delp_lo,antifx,nvl,				&
                      string(1:12)//'(cnuity) low-ord dp, unclipped antidiff flx')
      end if

      if (TimingBarriers) then
        ret = gptlstart ('cnuity_barrier')
!SMS$BARRIER
        ret = gptlstop ('cnuity_barrier')
      end if

      ret = gptlstart ('cnuity_exchange')
!SMS$EXCHANGE(r_plus,r_mnus)
      ret = gptlstop ('cnuity_exchange')

!sms$compare_var(antifx    , "cnuity.F90 - antifx9 ")
      ret = gptlstart ('cnuity')
      ret = gptlstart ('cnuity4')
! --- cnuity4 clips antidiff fluxes and computes new high-order dp
      call cnuity4 (psurf, antifx, r_plus, r_mnus, massfx, &
                    anti_tndcy, dplo_tndcy, dp_tndcy, delp)
      ret = gptlstop ('cnuity4')
    end if			! low_ord = true or false

#ifdef DEBUGPRINT
   anti_tndcy(:,:)=dp_tndcy(:,:,nf)
   do k=1,nvl,7
    write (string,'(a,i2)') 'dp_tndcy (nf), k=',k
    call findmxmn2(anti_tndcy,nvl,nip,k,string)
    write (string,'(a,i2)') 'cnuity: new dp k=',k
    call findmxmn2(delp,nvl,nip,k,string)
   end do
#endif

    if (PrintIpnDiag > 0) then
      write (string,'(a,i6,2x)') 'step',its                        
      call stenedg (delp_lo,antifx,nvl,				&
                  string(1:12)//'(cnuity) low-ord dp, clipped antidiff flx')
      call stencl (delp,nvl,.01,string(1:12)//'(cnuity) new dp')
    end if

    ret = gptlstart ('cnuity5')
! --- cnuity5 
    call cnuity5 (pres, delp, exner, dp_tndcy, omega, &
                  lp_edg, u_velo, v_velo)
    ret = gptlstop ('cnuity5')

#ifdef DEBUGPRINT
  do k = 1,nvl,7
   write (string,'(a,i3,a)') 'k',k,' cnuity:omega'
   call findmxmn2(omega,nvl,nip,k,string)
  end do
  print *
#endif

    if (.not.low_ord) then
!sms$compare_var(r_plus    , "cnuity.F90 - r_plus10 ")
!sms$compare_var(r_mnus    , "cnuity.F90 - r_mnus10 ")
!sms$compare_var(antifx    , "cnuity.F90 - antifx10 ")
    end if
!sms$compare_var(massfx    , "cnuity.F90 - massfx10")
!sms$compare_var(dplo_tndcy, "cnuity.F90 - dplo_tndcy10")
!sms$compare_var(delp_lo   , "cnuity.F90 - delp_lo10")

    ret = gptlstop ('cnuity')

    return
  end subroutine cnuity
end module cnuity_driver
