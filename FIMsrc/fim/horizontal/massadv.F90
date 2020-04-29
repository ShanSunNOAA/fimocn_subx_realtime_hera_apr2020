module massadv_driver
  implicit none

  private
  public :: mass_adv

contains
!*********************************************************************
!     massadv
!
!     A new mass transport scheme based on a flux-limitied advection
!     scheme of Thuburn 1995, 1996 and Miura 2012.
!    
!     The implementation of this module, its interface and code
!     structure, is partially adapted from cnuity.F90.
!
!     N. Wang, Dec 2015. 
!
!*********************************************************************

  subroutine mass_adv (its, u_velo, v_velo, u_edg, v_edg, &
                     dp_edg, lp_edg, delp, pres, exner, &
                     dp_tndcy, massfx, omega, TimingBarriers)

    use module_control,   only: nvlp1, npp, nip, nabl
    use fimnamelist,      only: nvl, PrintIpnDiag
    use module_variables, only: nf
    use module_constants, only: sidevec_e
    use findmaxmin2
    use stencilprint
    use stenedgprint
    use fluxlimiter, only: init_mass_limiter, limiter, c_num
    use massadv_compute, only: update_pres, massflux_newdepth
 
    use global_bounds,    only: ims, ime, ips, ipe, ihe
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
    real   ,intent (OUT)   :: omega     (nvl,nip)	! dp/dt (N/sec)
!sms$distribute end
!sms$distribute (dh,3) begin
    real   ,intent (IN)    :: u_edg     (nvl,npp,nip)	! u on edges (m/s)
    real   ,intent (IN)    :: v_edg     (nvl,npp,nip)	! v on edges (m/s)
    real   ,intent (INOUT) :: dp_edg    (nvl,npp,nip)	! delta-p on edges (Pa)
    real   ,intent (IN)    :: lp_edg    (nvl,npp,nip)	! mid-layer pressure on edges (Pa)
    real   ,intent (INOUT) :: massfx    (nvl,npp,nip,3)	! mass flux across edges (N/sec)
    real                   :: limter    (nvl,npp,nip)
!sms$distribute end

    integer :: k	 		   ! layer index
    integer :: ipn		           ! Index for icos cell number
    character :: string*32
    logical,parameter :: low_ord = .false. ! if true, skip antidiffusive part
    integer :: ret


#include <gptl.inc>

    ret = gptlstart ('mass_adv')

#ifdef DEBUGPRINT
   do k=1,nvl,7
    write (string,'(a,i2)') 'mass_adv: old dp k=',k
    call findmxmn2(delp,nvl,nip,k,string)
   end do
#endif

!sms$compare_var(u_velo    , "enter massadv.F90 - u_velo ")
!sms$compare_var(v_velo    , "enter massadv.F90 - v_velo ")
!sms$compare_var(u_edg    , "enter massadv.F90 - u_edg ")
!sms$compare_var(v_edg    , "enter massadv.F90 - v_edg ")
!sms$compare_var(lp_edg    , "enter massadv.F90 - lp_edg ")
!sms$compare_var(delp    , "enter massadv.F90 - delp ")
!sms$compare_var(pres    , "enter massadv.F90 - pres ")
!sms$compare_var(exner   , "enter massadv.F90 - exner ")
!sms$compare_var(dp_tndcy   , "enter massadv.F90 - dp_tndcy ")
!sms$compare_var(omega   , "enter massadv.F90 - omega ")
!sms$compare_var(massfx   , "enter massadv.F90 - massfx ")
!sms$compare_var(dp_edg   , "enter massadv.F90 - dp_edg ")


! initialize mass flux limiter 
   call init_mass_limiter(u_edg, v_edg, delp, dp_edg)

! compute and apply limits to mass flux 
   call limiter(delp, dp_edg, c_num, dp_tndcy)

! compute the mass flux and new depth field
   call massflux_newdepth(dp_edg,massfx,dp_tndcy,delp,psurf) 

#ifdef DEBUGPRINT
   do k=1,nvl,7
    write (string,'(a,i2)') 'mass_adv: new dp k=',k
    call findmxmn2(delp,nvl,nip,k,string)
   end do
#endif

! update pressure, exner function, and omega function
    call update_pres(pres, delp, exner, dp_tndcy, omega, &
                  lp_edg, u_velo, v_velo)

#ifdef DEBUGPRINT
  do k = 1,nvl,7
   write (string,'(a,i3,a)') 'k',k,' mass_adv:omega'
   call findmxmn2(omega,nvl,nip,k,string)
  end do
  print *
#endif

#ifdef DEBUGPRINT
  do k = 1,nvl,7
   write (string,'(a,i3,a)') 'k',k,' mass_adv:exner'
   call findmxmn2(exner,nvl,nip,k,string)
  end do
  print *
#endif

#ifdef DEBUGPRINT
  do k = 1,nvl,7
   write (string,'(a,i3,a)') 'k',k,' mass_adv:pres'
   call findmxmn2(pres,nvl,nip,k,string)
  end do
  print *
#endif

    ret = gptlstop ('mass_adv')

!!$OMP PARALLEL DO PRIVATE (k) SCHEDULE (runtime)
!    do ipn = ips, ipe
!      if (ipn < 100) then
!        do k = 1, nvl
!          print*, dp_edg(k, 1, ipn)
!        enddo
!      endif
!    enddo  

!sms$compare_var(delp    , "exit massadv.F90 - delp ")
!sms$compare_var(pres    , "exit massadv.F90 - pres ")
!sms$compare_var(exner   , "exit massadv.F90 - exner ")
!sms$compare_var(dp_tndcy   , "exit massadv.F90 - dp_tndcy ")
!sms$compare_var(omega   , "exit massadv.F90 - omega ")
!sms$compare_var(massfx   , "exit massadv.F90 - massfx ")
!sms$compare_var(dp_edg   , "exit massadv.F90 - dp_edg ")

    return
  end subroutine mass_adv

end module massadv_driver
