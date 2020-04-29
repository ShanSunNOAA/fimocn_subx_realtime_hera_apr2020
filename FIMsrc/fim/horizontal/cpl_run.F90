module module_fim_cpl_run

  implicit none
  save

!TODO:  move to a new cpl_internal_state module
!sms$distribute(dh,2) begin
real*8,allocatable :: u_tdcy_phy  (:,:)         ! physics forcing of u
real*8,allocatable :: v_tdcy_phy  (:,:)         ! physics forcing of v
real*8,allocatable :: trc_tdcy_phy(:,:,:)       ! physics forcing of tracers
!sms$distribute end
!sms$distribute(dh,1) begin
real*8,allocatable :: gfs_u_old(:,:)
real*8,allocatable :: gfs_v_old(:,:)
real*8,allocatable :: gfs_t_old(:,:)
real*8,allocatable :: gfs_q_old(:,:)
real*8,allocatable :: gfs_oz_old(:,:)
real*8,allocatable :: gfs_cld_old(:,:)
!sms$distribute end

#include <gptl.inc>

contains

!*********************************************************************
!       "Run" method for the FIM DYN-PHY coupler component.  
!       Argument its is time step count.  
!       Argument dyn_to_phy controls the direction of coupling:  
!         dyn_to_phy == .true.    Couple from DYN to PHY
!         dyn_to_phy == .false.   Couple from PHY to DYN
!       T. Henderson            February, 2009  - Moved code here from physics()
!       R. Bleck                July, 2010      - fixed layer pressure formula
!*********************************************************************

  subroutine cpl_run(its, dyn_to_phy)

    USE gfs_physics_internal_state_mod, only: gfs_physics_internal_state, gis_phy
    use module_control  ,only: CallPhysics,itsDFI
    use fimnamelist     ,only: nts,itsStart,digifilt
    use module_variables,only: us3d,vs3d,pr3d,tr3d,ws3d
    use module_sfc_variables

! Arguments
    integer, intent(in) :: its
    logical, intent(in) :: dyn_to_phy

! Local variables
    integer :: laststep
    logical :: callphy  ! .true. iff physics is called during this timestep
    integer :: ret

    ret = gptlstart ('cpl_run')

    !...........................................................
    ! Couple components unless this is the last (nts+1) 
    ! iteration (in which DYN just finishes).  
    ! This complexity is required for the NCEP ESMF approach 
    ! in which single-phase DYN and PHY components alternate 
    ! execution during each time step.  
    !

    if (digifilt) then
      laststep=itsDFI+nts
    else
      laststep=itsStart+nts
    end if
    if (its < laststep) then

      !TODO:  Eliminate duplication by encapsulating this logic
      callphy = mod(its,CallPhysics)==0.or.its==1.or.(its==itsDFI.and.digifilt)

!sms$compare_var(st3d   , 'begin cpl_run')
!sms$compare_var(sm3d   , 'begin cpl_run')
!sms$compare_var(rn2d   , 'begin cpl_run')
!sms$compare_var(rc2d   , 'begin cpl_run')
!sms$compare_var(ts2d   , 'begin cpl_run')
!sms$compare_var(us2d   , 'begin cpl_run')
!sms$compare_var(hf2d   , 'begin cpl_run')
!sms$compare_var(rsds   , 'begin cpl_run')
!sms$compare_var(slmsk2d, 'begin cpl_run')
!sms$compare_var(strs2d,  'begin cpl_run')

      if (dyn_to_phy) then
        ! Subroutine cpl_dyn_to_phy() converts FIM values to GFS 
        ! values.  
        ! All arrays passed via the ESMF coupler are passed as 
        ! arguments, allowing this subroutine to be called from
        ! the ESMF coupler too.  
        call cpl_dyn_to_phy(its, callphy,                         &
! IN args
             pr3d, us3d, vs3d, ws3d, tr3d,                        &
! OUT args
             gis_phy%ps, gis_phy%dp, gis_phy%p,                   &
             gis_phy%u , gis_phy%v , gis_phy%dpdt,                &
             gis_phy%q , gis_phy%oz, gis_phy%cld,                 &
             gis_phy%t )
      else
        ! Subroutine cpl_phy_to_dyn() converts GFS values to FIM 
        ! values.  
        ! All arrays passed via the ESMF coupler are passed as 
        ! arguments, allowing this subroutine to be called from
        ! the ESMF coupler too.  

        call cpl_phy_to_dyn(its, callphy, gis_phy%deltim,		&
! IN args
             gis_phy%p, gis_phy%u  , gis_phy%v,				&
             gis_phy%q, gis_phy%cld, gis_phy%oz, gis_phy%t,		&
             ! these GFS PHY fields are passed to FIM DYN for 
             ! output and diagnostics only.  
             gis_phy%flx_fld%GESHEM,    gis_phy%flx_fld%RAINC,		&
             gis_phy%sfc_fld%TSEA,      gis_phy%sfc_fld%UUSTAR,		&
             gis_phy%flx_fld%HFLX,      gis_phy%flx_fld%EVAP,		&
             gis_phy%sfc_fld%SHELEG,    gis_phy%sfc_fld%CANOPY,		&
             gis_phy%sfc_fld%HICE,      gis_phy%sfc_fld%FICE,		&
             gis_phy%sfc_fld%STC,       gis_phy%sfc_fld%SMC,		&
             gis_phy%flx_fld%SFCDSW,    gis_phy%flx_fld%SFCDLW,		&
             gis_phy%flx_fld%SFCUSW,    gis_phy%flx_fld%SFCULW,		&   !hli
             gis_phy%flx_fld%TOPDSW,    gis_phy%flx_fld%TOPUSW,		&   !hli
             gis_phy%flx_fld%TOPULW,                     		&   !hli
             gis_phy%sfc_fld%T2M,       gis_phy%sfc_fld%Q2M,		&
             gis_phy%flx_fld%u10m,      gis_phy%flx_fld%v10m,		&
             gis_phy%sfc_fld%SLMSK,     gis_phy%HPRIME,			&
             gis_phy%FLUXR,             gis_phy%SFALB,			&
             gis_phy%sfc_fld%STRESS,    gis_phy%flx_fld%runoff,		&
! OUT args
             us3d, vs3d, tr3d, rn2d, rc2d, ts2d, us2d, hf2d, qf2d,	&
             sheleg2d, canopy2d, hice2d, fice2d, st3d, sm3d,		&
             rsds, rlds, rsus, rlus, rsdt, rsut, rlut,                  &
             t2m2d, q2m2d, u10m, v10m,					&
	     slmsk2d, hprm2d, alb2d, strs2d, rain, runoff)
    endif

!sms$compare_var(st3d   , 'end cpl_run')
!sms$compare_var(sm3d   , 'end cpl_run')
!sms$compare_var(rn2d   , 'end cpl_run')
!sms$compare_var(rc2d   , 'end cpl_run')
!sms$compare_var(ts2d   , 'end cpl_run')
!sms$compare_var(us2d   , 'end cpl_run')
!sms$compare_var(hf2d   , 'end cpl_run')
!sms$compare_var(rsds   , 'end cpl_run')
!sms$compare_var(slmsk2d, 'end cpl_run')

    endif

    ret = gptlstop ('cpl_run')

    return
  end subroutine cpl_run

! Couple from DYN->PHY.  
  subroutine cpl_dyn_to_phy(its, callphy, &
               pr3d,us3d,vs3d,ws3d,tr3d,  &
               gfs_ps, gfs_dp, gfs_p,     &
               gfs_u , gfs_v , gfs_dpdt,  &
               gfs_q , gfs_oz, gfs_cld,   &
               gfs_t )

! TODO:  Pass elements of gfs_physics_internal_state_mod via argument 
! TODO:  list and remove use of gfs_physics_internal_state_mod.  
!TBH:  NOTE removal of "only" clause.  This was forced upon us by *bug* in the 
!TBH:  ifort 11.1 compiler on njet.  Restore the "only" clause when the broken 
!TBH:  compiler is fixed!  
!USE gfs_physics_internal_state_mod, only: gfs_physics_internal_state, gis_phy
    USE gfs_physics_internal_state_mod
    use module_constants,only: cp, rd, p1000, lat, lon, qvmin
    use module_control  ,only: nip, CallRadiation
    use fimnamelist     ,only: nvl, itsStart
    use module_sfc_variables, only: zorl2d,srflag2d
    USE MACHINE         ,only: kind_evod

    integer, intent(in) :: its
    logical, intent(in) :: callphy
!sms$distribute (dh,2) begin
    real, intent(in) :: pr3d(:,:)
    real, intent(in) :: us3d(:,:)
    real, intent(in) :: vs3d(:,:)
    real, intent(in) :: ws3d(:,:)
    real, intent(in) :: tr3d(:,:,:)
!sms$distribute end
!sms$distribute (dh,1) begin
    real(kind=kind_evod) , intent(out) :: gfs_ps(:)
    real(kind=kind_evod) , intent(out) :: gfs_dp(:,:)
    real(kind=kind_evod) , intent(out) :: gfs_p(:,:)
    real(kind=kind_evod) , intent(out) :: gfs_u(:,:)
    real(kind=kind_evod) , intent(out) :: gfs_v(:,:)
    real(kind=kind_evod) , intent(out) :: gfs_dpdt(:,:)
    real(kind=kind_evod) , intent(out) :: gfs_q(:,:)
    real(kind=kind_evod) , intent(out) :: gfs_oz(:,:)
    real(kind=kind_evod) , intent(out) :: gfs_cld(:,:)
    real(kind=kind_evod) , intent(out) :: gfs_t(:,:)
!sms$distribute end

!  local variables
    integer :: ipn
    integer :: k
    real    :: rocp1
    real    :: rocpr
!sms$distribute (dh,2) BEGIN
    real :: tr_2(nvl,nip)
    real :: theta_nv(nvl,nip)
!sms$distribute end

    if (callphy) then
      rocp1 = rd/cp + 1.
      rocpr = cp/rd

!SMS$PARALLEL (dh,ipn) BEGIN
!$OMP PARALLEL DO PRIVATE (k) SCHEDULE (runtime)
      do ipn=1,nip
        do k=1,nvl
          tr_2(k,ipn) = max(qvmin, tr3d(k,ipn,2))
        end do

!NOTE:  ZORL is held constant in FIM-GFS coupling.  Bao confirms that this 
!NOTE:  is what we want.  
        gis_phy%sfc_fld%ZORL(ipn,1) = zorl2d(ipn)
        !NOTE:  This logic replicates Bao's original logic in which 
        !NOTE:  initial values of TRPCP was read from a file in phy_init() but 
        !NOTE:  overwritten from GESHEM after the first call to phy_run().  
        !NOTE:  Bao has checked the original logic and verified that it 
        !NOTE:  behaved as he intended.  
        if (its > 1) then
          gis_phy%sfc_fld%TPRCP(ipn,1) = max(0.0d0, gis_phy%flx_fld%GESHEM(ipn,1))
        end if
!NOTE:  SRFLAG is held constant in FIM-GFS coupling.  Bao says "OK".  
!       gis_phy%sfc_fld%SRFLAG(ipn,1) = srflag2d(ipn)
        ! Bao confirms that overwrite with TG3 is intentional here.  
        gis_phy%flx_fld%TMPMIN(ipn,1) = gis_phy%sfc_fld%TG3(ipn,1)
        gis_phy%flx_fld%TMPMAX(ipn,1) = gis_phy%sfc_fld%TG3(ipn,1)

        gfs_ps(ipn)       =  0.001*pr3d(1,ipn)
        do k=1,nvl
          gfs_dp(ipn,k)   = (0.001*pr3d(k,ipn))-(0.001*pr3d(k+1,ipn))
!       gfs_p(ipn,k)    = 0.5*(pr3d(k,ipn)+pr3d(k+1,ipn))
! get energetically consistent mid-lyr prs from partial[p^(kap+1)]/partial[p]
          gfs_p(ipn,k)    = ((pr3d(k,ipn)**rocp1-pr3d(k+1,ipn)**rocp1)/	&
                            ((pr3d(k,ipn)       -pr3d(k+1,ipn)       )*rocp1))**rocpr
          gfs_u(ipn,k)    = us3d(k,ipn)
          gfs_v(ipn,k)    = vs3d(k,ipn)
          gfs_dpdt(ipn,k) = 0.001*ws3d(k,ipn)
          gfs_q(ipn,k)    = tr_2(k,ipn)
          gfs_oz(ipn,k)   = tr3d(k,ipn,4)
          gfs_cld(ipn,k ) = tr3d(k,ipn,3)
          theta_nv(k,ipn) = tr3d(k,ipn,1) / (1. + 0.6078*max(qvmin,tr_2(k,ipn)))
          gfs_t(ipn,k)    = theta_nv(k,ipn)*(gfs_p(ipn,k)/p1000)**(rd/cp)
          ! save values prior to physics call for use in computing 
          ! tendencies after physics is called
          gfs_u_old(ipn,k) = gfs_u(ipn,k)
          gfs_v_old(ipn,k) = gfs_v(ipn,k)
          gfs_t_old(ipn,k) = gfs_t(ipn,k)
          gfs_q_old(ipn,k) = gfs_q(ipn,k)
          gfs_oz_old(ipn,k) = gfs_oz(ipn,k)
          gfs_cld_old(ipn,k) = gfs_cld(ipn,k)
          ! Save DYN fields to avoid INOUT INTENT for these arrays in 
          ! cpl_phy_to_dyn() which breaks ESMF data model.  See rant below.  
        end do
      end do
!$OMP END PARALLEL DO
!SMS$PARALLEL END
    end if

    return
  end subroutine cpl_dyn_to_phy

! Couple from PHY->DYN.  
  subroutine cpl_phy_to_dyn(its, callphy, dtp,		&
! IN args
               gfs_p, gfs_u  , gfs_v,			&
               gfs_q, gfs_cld, gfs_oz, gfs_t,		&
               ! these GFS PHY fields are passed to FIM DYN for 
               ! output and diagnostics only.  
               gfs_geshem, gfs_rainc,			&
               gfs_tsea,   gfs_uustar,			&
               gfs_hflx,   gfs_evap,			&
               gfs_sheleg, gfs_canopy,			&
               gfs_hice,   gfs_fice,			&
               gfs_stc,    gfs_smc,			&
               gfs_sfcdsw, gfs_sfcdlw,			&
               gfs_sfcusw, gfs_sfculw,               &  !hli 09/2014
               gfs_topdsw, gfs_topusw,               &  !hli 09/2014
               gfs_topulw,                           &  !hli 09/2014
               gfs_t2m,    gfs_q2m,			&
               gfs_u10m,   gfs_v10m,			&
               gfs_slmsk,  gfs_hprime,			&
               gfs_fluxr,  gfs_alb,			&
               gfs_stress, gfs_runoff,			&
! OUT args
             us3d,  vs3d,    tr3d,			&
               ! these GFS PHY fields are passed to FIM DYN for 
               ! output and diagnostics only.  
             rn2d, rc2d, ts2d, us2d, hf2d, qf2d,        &
             sheleg2d, canopy2d, hice2d, fice2d,        &
             st3d, sm3d, rsds, rlds, rsus, rlus,        &
             rsdt, rsut, rlut,                          &
             t2m2d, q2m2d, u10m, v10m, slmsk2d,         &
	     hprm2d, alb2d, strs2d, rain, runoff )

!tgs
    USE dffusn, only: dffusn8_lyr
    use module_constants,only: cp, rd, p1000, qvmin, qwmin
    USE gfs_physics_internal_state_mod, only: gis_phy
    use module_control  ,only: nip, ntra, ntrb, dt
    use fimnamelist     ,only: nts, nvl, addtend, smoothtend, TimingBarriers
    USE MACHINE         ,only: kind_evod,kind_phys,kind_rad
! <begin rant>
! dp3d should be an intent(in) argument but it is contrary to ESMF architecture 
! to make it available to be read from the DYN import state (because import 
! state should be write-only).  So we must get it from the DYN module directly. 
! This illustrates a fundamental flaw in the ESMF component design.  
! <end rant>
    use module_sfc_variables,only: cld_hi,cld_md,cld_lo,cld_tt,cld_bl
    use module_variables,only: dp3d
    
    integer, intent(in) :: its
    logical, intent(in) :: callphy
    real(kind=kind_phys), intent(in) :: dtp
!sms$distribute (dh,1) begin
    real(kind=kind_evod), intent(in) :: gfs_p(:,:)
    real(kind=kind_evod), intent(in) :: gfs_u(:,:)
    real(kind=kind_evod), intent(in) :: gfs_v(:,:)
    real(kind=kind_evod), intent(in) :: gfs_q(:,:)
    real(kind=kind_evod), intent(in) :: gfs_cld(:,:)
    real(kind=kind_evod), intent(in) :: gfs_oz(:,:)
    real(kind=kind_evod), intent(in) :: gfs_t(:,:)
    real(kind=kind_phys), intent(in) :: gfs_geshem(:,:)
    real(kind=kind_phys), intent(in) :: gfs_rainc(:,:)
    real(kind=kind_phys), intent(in) :: gfs_tsea(:,:)
    real(kind=kind_phys), intent(in) :: gfs_uustar(:,:)
    real(kind=kind_phys), intent(in) :: gfs_hflx(:,:)
    real(kind=kind_phys), intent(in) :: gfs_evap(:,:)
    real(kind=kind_phys), intent(in) :: gfs_sheleg(:,:)
    real(kind=kind_phys), intent(in) :: gfs_canopy(:,:)
    real(kind=kind_phys), intent(in) :: gfs_hice  (:,:)
    real(kind=kind_phys), intent(in) :: gfs_fice  (:,:)
    real(kind=kind_phys), intent(in) :: gfs_sfcdsw(:,:)
    real(kind=kind_phys), intent(in) :: gfs_sfcdlw(:,:)
    real(kind=kind_phys), intent(in) :: gfs_sfcusw(:,:) !hli 09/2014
    real(kind=kind_phys), intent(in) :: gfs_sfculw(:,:) !hli 09/2014
    real(kind=kind_phys), intent(in) :: gfs_topdsw(:,:) !hli 09/2014
    real(kind=kind_phys), intent(in) :: gfs_topusw(:,:) !hli 09/2014
    real(kind=kind_phys), intent(in) :: gfs_topulw(:,:) !hli 09/2014
    real(kind=kind_phys), intent(in) :: gfs_t2m   (:,:)
    real(kind=kind_phys), intent(in) :: gfs_q2m   (:,:)
    real(kind=kind_phys), intent(in) :: gfs_u10m  (:,:)
    real(kind=kind_phys), intent(in) :: gfs_v10m  (:,:)
    real(kind=kind_phys), intent(in) :: gfs_slmsk (:,:)
    real(kind=kind_phys), intent(in) :: gfs_alb   (:,:)
    real(kind=kind_phys), intent(in) :: gfs_stress(:,:)
    real(kind=kind_phys), intent(in) :: gfs_runoff(:,:)
!sms$distribute end

!sms$distribute (dh,2) begin
    real(kind=kind_phys), intent(in) :: gfs_stc   (:,:,:)
    real(kind=kind_phys), intent(in) :: gfs_smc   (:,:,:)
    real(kind=kind_rad),  intent(in) :: gfs_hprime(:,:,:)
    real(kind=kind_rad),  intent(in) :: gfs_fluxr (:,:,:)
    real, intent(out) :: us3d(:,:)
    real, intent(out) :: vs3d(:,:)
  ! TBH:  tr3d must be inout instead of out because tr3d(:,:,4) is 
  ! TBH:  never set.  
    real, intent(inout) :: tr3d(:,:,:)
!sms$distribute end

!sms$distribute (dh,1) begin
    real, intent(inout) :: rn2d(:)
    real, intent(inout) :: rc2d(:)
    real, intent(out)   :: ts2d(:)
    real, intent(out)   :: us2d(:)
    real, intent(out)   :: hf2d(:)
    real, intent(out)   :: qf2d(:)
    real, intent(out)   :: sheleg2d(:)
    real, intent(out)   :: canopy2d(:)
    real, intent(out)   :: hice2d(:)
    real, intent(out)   :: fice2d(:)
    real, intent(out)   :: rsds(:)
    real, intent(out)   :: rlds(:)
    real, intent(out)   :: rsus(:)
    real, intent(out)   :: rlus(:)
    real, intent(out)   :: rsdt(:)
    real, intent(out)   :: rsut(:)
    real, intent(out)   :: rlut(:)
    real, intent(out)   :: t2m2d(:)
    real, intent(out)   :: q2m2d(:)
    real, intent(out)   :: u10m(:)
    real, intent(out)   :: v10m(:)
    real, intent(out)   :: slmsk2d(:)
    real, intent(out)   :: alb2d(:)
    real, intent(out)   :: strs2d(:)
    real, intent(out)   :: rain(:)
    real, intent(out)   :: runoff(:)
!sms$distribute end

!sms$distribute (dh,2) begin
    real, intent(out)   :: st3d(:,:)
    real, intent(out)   :: sm3d(:,:)
    real, intent(out)   :: hprm2d(:,:)
!sms$distribute end

!  local variables
    integer :: ipn
    integer :: k
    integer :: t   ! tracer index
    integer :: ret ! gptl return code
    real (kind=kind_phys) :: rn2dten
    real (kind=kind_phys) :: rc2dten

    if (callphy) then
!SMS$PARALLEL (dh,ipn) BEGIN

!----------------------------------------------------------------------
!   Move all output into correct FIM arrays
!----------------------------------------------------------------------

      IF (addtend) THEN
!        write(6,*) 'Addtend, its, dt',addtend, its, dt
        do ipn=1,nip
          do k=1,nvl
            u_tdcy_phy(k,ipn) = (gfs_u(ipn,k)-gfs_u_old(ipn,k))/dtp
            v_tdcy_phy(k,ipn) = (gfs_v(ipn,k)-gfs_v_old(ipn,k))/dtp
            trc_tdcy_phy(k,ipn,1) = (gfs_t(ipn,k)-gfs_t_old(ipn,k))/dtp
!tgs    trc_tdcy_phy(k,ipn,1) = (gfs_t(ipn,k)-gfs_t(ipn,k))/dtp &
!       *(p1000/gfs_p(ipn,k))**(rd/cp)
            trc_tdcy_phy(k,ipn,2) = (gfs_q(ipn,k)-gfs_q_old(ipn,k))/dtp
            trc_tdcy_phy(k,ipn,3) = (gfs_cld(ipn,k)-gfs_cld_old(ipn,k))/dtp
            trc_tdcy_phy(k,ipn,4) = (gfs_oz(ipn,k)-gfs_oz_old(ipn,k))/dtp
          enddo
        enddo
      ENDIF

      IF(addtend .and. smoothtend) then
        write(6,*) ' Smooth tendencies'
        if (TimingBarriers) then
          ret = gptlstart ('cpl_run_barrier')
!SMS$BARRIER
          ret = gptlstop ('cpl_run_barrier')
        end if
        ret = gptlstart ('cpl_run_exchange')
!SMS$EXCHANGE(dp3d, u_tdcy_phy, v_tdcy_phy, trc_tdcy_phy(:,:,1:4))
        ret = gptlstop ('cpl_run_exchange')

        ret = gptlstart ('dffusn8_lyr')
        call dffusn8_lyr (u_tdcy_phy, dp3d, dt*20.)
        call dffusn8_lyr (v_tdcy_phy, dp3d, dt*20.)
        do t=1,4
          call dffusn8_lyr (trc_tdcy_phy(:,:,t), dp3d, dt*20.)
        enddo
        ret = gptlstop ('dffusn8_lyr')
      ENDIF

!      write(6,*) 'Physics to dynamics, addtend, dt',addtend, dt
!$OMP PARALLEL DO PRIVATE (k,rn2dten,rc2dten) SCHEDULE (runtime)
      do ipn=1,nip
        IF(addtend) then
        ! TBH:  It should be possible to remove this and instead call 
        ! TBH:  addphytend() every time step.  When this was tried, 
        ! TBH:  prognostic variables diverged gradually.  
          do k=1,nvl
            us3d(k,ipn) = gfs_u_old(ipn,k) + u_tdcy_phy(k,ipn)*dt
            vs3d(k,ipn) = gfs_v_old(ipn,k) + v_tdcy_phy(k,ipn)*dt
          ! Replace values for prognostic variables
! Tanya:  why use a different equation here than addphytend() uses?  
            tr3d(k,ipn,1) = (gfs_t_old(ipn,k) + trc_tdcy_phy(k,ipn,1)*dt) &
                 *(p1000/gfs_p(ipn,k))**(rd/cp)*(1.+0.6078*max(REAL(qvmin,kind_evod), &
                 (gfs_q_old(ipn,k)+ trc_tdcy_phy(k,ipn,2)*dt)))
            tr3d(k,ipn,2) = max(REAL(qvmin,kind_evod),(gfs_q_old(ipn,k)+ trc_tdcy_phy(k,ipn,2)*dt)) ! qv
            tr3d(k,ipn,3) = max(REAL(qwmin,kind_evod),(gfs_cld_old(ipn,k)+trc_tdcy_phy(k,ipn,3)*dt))
            tr3d(k,ipn,4) = max(REAL(qvmin,kind_evod),(gfs_oz_old(ipn,k) + trc_tdcy_phy(k,ipn,4)*dt))
          enddo
        ELSE
          do k=1,nvl
            us3d(k,ipn) = gfs_u(ipn,k)
            vs3d(k,ipn) = gfs_v(ipn,k)
            ! Replace values for prognostic variables
            tr3d(k,ipn,1) = gfs_t(ipn,k)*(p1000/gfs_p(ipn,k))**(rd/cp)* &
                            (1. + 0.6078*max(REAL(qvmin,kind_evod),gfs_q(ipn,k)))
            tr3d(k,ipn,2) = max(REAL(qvmin,kind_evod),gfs_q(ipn,k))
            tr3d(k,ipn,3) = max(REAL(qwmin,kind_evod),gfs_cld(ipn,k))
            tr3d(k,ipn,4) = max(REAL(qvmin,kind_evod),gfs_oz(ipn,k))
          end do
        END IF
!TBH:  Remaining fields are needed only for FIM diagnostics.  
        rn2dten = max(0.0_kind_phys, gfs_GESHEM(ipn,1))
        rc2dten = max(0.0_kind_phys, gfs_RAINC(ipn,1))
!TBH:  I assume here that NEMS allows the modified state to be either 
!TBH:  intent(out) *or* intent(inout) ...  
        rn2d(ipn)	= rn2d(ipn) + rn2dten*1000.
        rc2d(ipn)	= rc2d(ipn) + rc2dten*1000.
        us2d(ipn)	=  gfs_UUSTAR(ipn,1)
        hf2d(ipn)	= -gfs_HFLX  (ipn,1)	! sensible heatflux, switch to positive down
        qf2d(ipn)	= -gfs_EVAP  (ipn,1)	!   latent heatflux, switch to positive down
        rain(ipn)       =  rn2dten		! m/per time step
        runoff(ipn)     =  gfs_runoff(ipn,1)	! m/per time step
        sheleg2d(ipn)	=  gfs_SHELEG(ipn,1)
        canopy2d(ipn)	=  gfs_CANOPY(ipn,1)
        st3d(:,ipn)	=  gfs_STC (:,ipn,1)
        sm3d(:,ipn)	=  gfs_SMC (:,ipn,1)
        t2m2d(ipn)	=  gfs_T2M   (ipn,1)
        q2m2d(ipn)	=  gfs_Q2M   (ipn,1)
        u10m(ipn)	=  gfs_u10m(ipn,1)
        v10m(ipn)	=  gfs_v10m(ipn,1)
        slmsk2d(ipn)	=  gfs_SLMSK (ipn,1)
        hprm2d(:,ipn)	=  gfs_HPRIME(:,ipn,1)
        strs2d(ipn)	=  gfs_STRESS(ipn,1)	! sfc wind stress
        alb2d(ipn)	=  gfs_alb   (ipn,1)	! sfc albedo
        ts2d(ipn)	=  gfs_TSEA  (ipn,1)
        hice2d(ipn)	=  gfs_HICE  (ipn,1)
        fice2d(ipn)	=  gfs_FICE  (ipn,1)
        slmsk2d(ipn)	=  gfs_SLMSK (ipn,1)

! 2 identical ways to get rsds/rlds
        rsds(ipn)       =  gfs_SFCDSW(ipn,1)    ! sfc downward SW, pos down
        rlds(ipn)       =  gfs_SFCDLW(ipn,1)    ! sfc downward LW, pos down
        rsdt(ipn)       =  gfs_TOPDSW(ipn,1)    ! top downward SW, pos down
        rsus(ipn)       = -gfs_SFCUSW(ipn,1)    ! sfc   upward SW, switch to positive down
        rlus(ipn)       = -gfs_SFCULW(ipn,1)    ! sfc   upward LW, switch to positive down
        rsut(ipn)       = -gfs_TOPUSW(ipn,1)    ! top   upward SW, switch to positive down
        rlut(ipn)       = -gfs_TOPULW(ipn,1)    ! top   upward LW, switch to positive down
        cld_hi(ipn)     =  gfs_FLUXR( 5,ipn,1)
        cld_md(ipn)     =  gfs_FLUXR( 6,ipn,1)
        cld_lo(ipn)     =  gfs_FLUXR( 7,ipn,1)
        cld_tt(ipn)     =  gfs_FLUXR(17,ipn,1)
        cld_bl(ipn)     =  gfs_FLUXR(18,ipn,1)
      end do !ipn
!$OMP END PARALLEL DO
!SMS$PARALLEL END
    else   ! callphy is .false.
        IF(addtend) then
      ! Add physics tendencies for time steps between the physics calls
       write(6,*) 'Call Addphytend, addtend, its',addtend, its
          call addphytend(u_tdcy_phy,v_tdcy_phy,trc_tdcy_phy, &
                      us3d,vs3d,tr3d)
        ENDIF
    end if  ! callphy

    return
  end subroutine cpl_phy_to_dyn
end module module_fim_cpl_run
