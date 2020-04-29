module momtum_driver
  implicit none

contains

!*********************************************************************
!     momtum
!	Solves momentum equations
!	Alexander E. MacDonald			12/22/2004
!	J. Lee					September  2005
!	R. Bleck   major rewrite		April      2008
!	R. Bleck   removed omega diagnostics	August     2009
!	R. Bleck   biharmonic dissipation	April      2012
!*********************************************************************

  subroutine momtum (its, u_velo, v_velo, exner, relvor, 	 	&
                     u_edg, v_edg, trcr_edg, geop_edg, mont_edg,	&
                     uvsq_edg, u_tndcy, v_tndcy, dp3d)

    use findmaxmin2
    use module_constants, only: nprox, prox, rarea, sideln, sidevec_c,	&
                                sidevec_e, corio, cs, sn, dfflen_max,	&
                                biharw, grndspd
    use momtum_compute,   only: del4prep, dissip, momtum1
    use stencilprint,     only: stencl
    use module_variables, only: nf, of, vof, adbash1, adbash2, adbash3
    use module_control,   only: nvlp1, npp, nip, nabl, dt, ntra
    use fimnamelist,      only: nvl, PrintIpnDiag, TimingBarriers
    use global_bounds,    only: ims, ime, ips, ipe, ihe

!  External type and dimension:
    integer,intent (IN)  :: its			! model time step
!sms$distribute (dh,2) begin
    real   ,intent (INOUT) :: u_velo  (nvl,nip)
    real   ,intent (INOUT) :: v_velo  (nvl,nip)
    real   ,intent (IN)    :: exner   (nvlp1,nip)
    real   ,intent (INOUT) :: u_tndcy (nvl,nip,nabl)
    real   ,intent (INOUT) :: v_tndcy (nvl,nip,nabl)
    real   ,intent (OUT)   :: relvor  (nvl,nip)
    real   ,intent (IN)    :: dp3d    (nvl,nip)     ! layer thickness
    real :: worku(nvl,nip)
    real :: workv(nvl,nip)
!sms$distribute end
!sms$distribute (dh,3) begin
    real   ,intent (IN)    :: u_edg   (nvl,npp,nip)
    real   ,intent (IN)    :: v_edg   (nvl,npp,nip)
    real   ,intent (IN)    :: trcr_edg(nvl,npp,nip,ntra)
    real   ,intent (IN)    :: geop_edg(nvl,npp,nip)
    real   ,intent (IN)    :: mont_edg(nvl,npp,nip)
    real   ,intent (IN)    :: uvsq_edg(nvl,npp,nip)
!sms$distribute end

    integer   :: k		! layer index
    integer   :: ipn		! icos point index
    character :: string*24
    integer   :: ret

#include <gptl.inc>

#ifdef DEBUGPRINT
    do k=1,nvl,7
     write (string,'(a,i3,a)') 'k=',k,' u_velo'
     call findmxmn2(u_velo,nvl,nip,k,string)
  
     write (string,'(a,i3,a)') 'k=',k,' v_velo'
     call findmxmn2(v_velo,nvl,nip,k,string)
    end do
    print *
#endif

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! --- momentum dissipation (enhanced in the upper stratosphere to dampen
! --- gravity waves before they get reflected at model top)

! --- 'diffusion length' dfflen = (diffusivity) * (time step) / (mesh size)
! ---                           = (diffusion velocity) x (time step)
! --- dfflen_max = max value of dfflen

#ifdef DEBUGPRINT
    if (PrintIpnDiag > 0) then
      call stencl(u_velo,nvl,1.,'(atm dissip) -u- input')
      call stencl(v_velo,nvl,1.,'(atm dissip) -v- input')
    end if
#endif

    if (dfflen_max > 0.) then
      if (TimingBarriers) then
        ret = gptlstart ('momtum_barrier')
!SMS$BARRIER
        ret = gptlstop ('momtum_barrier')
      end if
      ret = gptlstart ('momtum_exchange')
!SMS$EXCHANGE (u_velo,v_velo,dp3d)
      ret = gptlstop ('momtum_exchange')

! Halo comp means horizontal loop indices are ips,ihe
!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (k)
      do ipn=ips,ihe
        do k=1,nvl
          worku(k,ipn) = u_velo(k,ipn)
          workv(k,ipn) = v_velo(k,ipn)
! --- use the following line to achieve angular momentum conservation
          worku(k,ipn) = worku(k,ipn)+grndspd(ipn)
        end do
      end do
!$OMP END PARALLEL DO
!sms$ignore end

! --- del4prep subtracts average of neighboring values from -u,v-
      if (maxval(biharw).gt.0.) then
        call del4prep(worku,workv)

        if (TimingBarriers) then
          ret = gptlstart ('momtum_barrier')
!SMS$BARRIER
          ret = gptlstop ('momtum_barrier')
        end if
        ret = gptlstart ('momtum_exchange')
!SMS$EXCHANGE (worku,workv)
        ret = gptlstop ('momtum_exchange')
      end if

! --- del^2 based dissipation of worku,workv (which can be either -u,v- or
! --- a quantity proportional to their negative laplacian, depending on
! --- whether del4prep has been called)

      call dissip (worku, workv, dp3d)

!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (k)
      do ipn=ips,ipe
        do k=1,nvl
          u_velo(k,ipn) = u_velo(k,ipn) + worku(k,ipn)
          v_velo(k,ipn) = v_velo(k,ipn) + workv(k,ipn)
        end do
      end do
!$OMP END PARALLEL DO
!sms$ignore end
    end if			! dfflen_max > 0
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

! --- evaluate remaining terms in momentum equation

!sms$compare_var(u_tndcy , "momtum.F90 - u_tndcy5 ")
!sms$compare_var(v_tndcy , "momtum.F90 - v_tndcy5 ")
!sms$compare_var(exner   , "momtum.F90 - exner5 ")
!sms$compare_var(mont_edg, "momtum.F90 - mont_edg5 ")

    call momtum1 (its, u_velo, v_velo, exner, relvor,		&
                  u_edg, v_edg, trcr_edg, geop_edg, mont_edg,	&
                  uvsq_edg, u_tndcy, v_tndcy, dp3d)

!sms$compare_var(u_tndcy, "momtum.F90 - u_tndcy6 ")
!sms$compare_var(v_tndcy, "momtum.F90 - v_tndcy6 ")
!sms$compare_var(u_velo , "momtum.F90 - u_velo6 ")
!sms$compare_var(v_velo , "momtum.F90 - v_velo6 ")

    return
  end subroutine momtum
end module momtum_driver
