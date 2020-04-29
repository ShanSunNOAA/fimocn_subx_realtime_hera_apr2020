module module_fim_dyn_run
  use module_constants
  use module_control  ,only: nvlp1,nip,ntra,ntrb,nabl,TestDiagProgVars,    &
                             TestDiagNoise,PrintDiagProgVars,dt,ArchvStep, &
                             ArchvStep0,itsDFI,digifilt,npp
  use fimnamelist     ,only: nts,nvl,itsStart,PrintIpnDiag,TimingBarriers, &
                             dtratio,ArchvTimeUnit,gribout,hiOrderFluxComp,&
                             hiOrderUVComp,flux_conserv_schm

  use module_variables,only: u_tdcy,v_tdcy,dp_tdcy,dpl_tdcy,               &
                             massfx,cumufx,dpinit,nf,of,vof,               &
                             adbash1,adbash2,adbash3,                      &
                             us3d,vs3d,dp3d,pr3d,ph3d,ex3d,mp3d,sdot,      &
                             tk3d,u_edg,v_edg,dp_edg,diaga,diagb,          &
                             lp_edg,geop_edg,mont_edg,uvsq_edg,            &
                             tr3d,trc_edg,trdp,trl_tdcy,trc_tdcy,          &
                             ws3d,sdot,rh3d,pw2d,pq2d,relvor,psrf,ptdcy,   &
                             worka,curr_write_time,wt_int_rev,k_int_rev,   &
                             flxavg

  use module_initial_chem_namelists, only: chem_opt,cu_physics,mp_physics
  USE gfs_physics_internal_state_mod, only:gis_phy
  use module_wrf_variables,only: phys3dwrf,phys2dwrf,exch
  use module_sfc_variables   ,only: aod2d,rn2d,rc2d,ts2d,us2d,st3d,sn2d,	&
                                    rsds,rlds,rsus,rlus,rsdt,rsut,rlut,qf2d,    &
                                    hf2d,sm3d,t2m2d,q2m2d,canopy2d,fice2d,	&
                                    hice2d,sheleg2d,slmsk2d,rn2d0,rc2d0,        &
                                    rg2d0,u10m,v10m
  use module_abstart         ,only: abstart
  use module_chem_output     ,only: chem_output
  use module_wrf_output      ,only: wrf_output
  use module_output          ,only: output
  use module_edgvar1         ,only: edgvar1
  use module_edgvar2         ,only: edgvar2
  use module_edgvar3         ,only: edgvar3
  use cnuity_driver          ,only: cnuity
  use massadv_driver         ,only: mass_adv
  use trcadv_driver          ,only: trcadv
  use traceradv_driver       ,only: traceradv
  use module_transp3d        ,only: transp0,transp1,transp2
  use momtum_driver          ,only: momtum
  use module_hystat          ,only: hystat
  use module_hybgen          ,only: hybgen 
  use module_diagnoise       ,only: diagnoise
  use module_globsum         ,only: globsum, qmass, qmsqv, qmsqw, qmso3, qmste, qmstr, qmstrn, qmstrc, &
                                    qdtr, qdtrn, qdtrc
  use module_profout         ,only: profout
  use module_outDiags        ,only: outDiags
  use module_chem_variables  ,only: p10,pm25,trfall,tr1_tavg
  use findmaxmin1
  use findmaxmin3
  use module_core_setup      ,only: iam_compute_root,use_write_tasks
  use global_bounds          ,only: ims, ime, ips, ipe
  use mdul_transects         ,only: trnsec1

  implicit none

#include <gptl.inc>

contains

  subroutine dyn_run(its)
!*********************************************************************
!	"Run" method for fim global model dynamics
!	Alexander E. MacDonald  12/24/2005
!	J. LEE                  12/28/2005
!*********************************************************************


! Arguments
    integer, intent(in) :: its

! Local variables
    integer :: itsm1        ! its - 1
    integer :: k,ipn,t,laststep
!SMS$DISTRIBUTE (dh,1) BEGIN
    real work(nip)
!SMS$DISTRIBUTE END
! diagnostic sums
    real    :: qtrcr(ntra+ntrb)
    logical, save :: qdtr_print=.false.
    real :: time
    integer :: ret  ! return code from GPTL routines
    real, external :: its2time
    logical::post_init_file_called=.false.

    ret = gptlstart ('Dynamics')
    itsm1 = its - 1

  !...........................................................
  ! Finish up previous dynamics time step unless this is the 
  ! first time step.  
  ! This complexity is required for the NCEP ESMF approach 
  ! in which single-phase DYN and PHY components alternate 
  ! execution during each time step.  
  !
!TBH:  Restore if statement label once Mark fixes PPP
!TBH skip_first_time1: if (its > itsStart ) then
    if (its > 1 ) then

!sms$compare_var(st3d   , "dyn_run begin second half of iteration ")
!sms$compare_var(sm3d   , "dyn_run begin second half of iteration ")
!sms$compare_var(rn2d   , "dyn_run begin second half of iteration ")
!sms$compare_var(rc2d   , "dyn_run begin second half of iteration ")
!sms$compare_var(ts2d   , "dyn_run begin second half of iteration ")
!sms$compare_var(us2d   , "dyn_run begin second half of iteration ")
!sms$compare_var(hf2d   , "dyn_run begin second half of iteration ")
!sms$compare_var(rsds   , "dyn_run begin second half of iteration ")

      !...........................................................
      ! Hybrid sigma-theta grid maintenance
      !
      if (TimingBarriers) then
        ret = gptlstart ('hybgen_barrier')
!SMS$BARRIER
        ret = gptlstop ('hybgen_barrier')
      end if

      call hybgen (itsm1,		&
           thetac,			& ! target pot.temp.
           us3d,vs3d,tr3d,		& ! zonal & merid wind, mass field tracers
           sdot,ex3d,dp3d,pr3d, TimingBarriers)	  ! interface displ, exner, layer thknss, prs
      !...........................................................
      ! Hydrostatic calculations

!sms$compare_var(ex3d, "dyn_run before hystat - ex3d5 ")
!sms$compare_var(ph3d, "dyn_run before hystat - ph3d5 ")

      ret = gptlstart ('hystat')
      call hystat (itsm1, ph3d,	ex3d, mp3d, dp3d, &
                   tr3d, trdp, psrf,ptdcy)
      ret = gptlstop ('hystat')

      !...........................................................
      ! Evaluate noise parameter
      !
      if( mod(itsm1,TestDiagNoise)==0 .or. itsm1==itsStart .or. &
         (itsm1==itsDFI.and.digifilt)) then
        call diagnoise (itsm1,ptdcy) ! sfc.pres.tdcy at 2 consec. time levels
      ! Diagnose wind speed near model top
!$OMP PARALLEL DO schedule (static)
        do ipn=ips,ipe
          work(ipn)=sqrt(us3d(nvl-2,ipn)**2+vs3d(nvl-2,ipn)**2)
        end do
!$OMP END PARALLEL DO
        call findmxmn1(work,nip,'wind speed in layer nvl-2')
      end if
      if (mod (itsm1,TestDiagProgVars) == 0 .or. itsm1 == itsStart .or. its == itsStart+nts .or. &
           ((itsm1==itsDFI.or.itsm1==itsDFI+nts).and.digifilt)) then
        call globsum (itsm1, dp3d, tr3d, rn2d, rc2d, pr3d, ex3d, qf2d, qtrcr)
      end if

      if (mod(itsm1,TestDiagProgVars) == 0 .or. ((itsm1 == itsStart .or. &
           ((itsm1==itsDFI.or.itsm1==itsDFI+nts).and.digifilt)) .and. &
           PrintDiagProgVars /= -2)) then
        call outDiags (itsm1, nvl, nvlp1, ntra+ntrb, pr3d, ph3d, tr3d, rn2d, rc2d, &
                       sdot, dp3d, us3d, vs3d, rh3d, rlds, rsds, hf2d, qf2d)
      end if
!sms$compare_var(mp3d, "dyn_run.F90 - mp3d10 ")

!sms$compare_var(st3d   , "dyn_run end iteration ")
!sms$compare_var(sm3d   , "dyn_run end iteration ")
!sms$compare_var(rn2d   , "dyn_run end iteration ")
!sms$compare_var(rc2d   , "dyn_run end iteration ")
!sms$compare_var(ts2d   , "dyn_run end iteration ")
!sms$compare_var(us2d   , "dyn_run end iteration ")
!sms$compare_var(hf2d   , "dyn_run end iteration ")
!sms$compare_var(rsds   , "dyn_run end iteration ")

    endif		! its > 1

!TBH:  Restore if statement label once Mark fixes PPP
!TBH endif skip_first_time1

    ret = gptlstop ('Dynamics')
    ret = gptlstart ('Output')
! Always call output routine, even for its-1 == 0
!............................................................

! If gribout is enabled and I am the compute root, prepare the grib output file
! for writing. The unsupported case where gribout is enabled and more than one
! write task is specified has already been handled in dyn_init. Note that it is
! assumed here that all routines calling icosio_out are doing so on the same
! time step. If any routine needs to write history to disk on a different
! schedule, the first "if" conditional below will need to change.

    if (mod(itsm1,ArchvStep) == 0) then
      if (gribout.and..not.use_write_tasks.and.iam_compute_root()) then
        call post_init_file(nint(its2time(itsm1)),ret)
        if (ret /= 0) then
          write(6,*) 'output: bad return from post_init_file: stopping'
          call flush(6)
          stop
        end if
        post_init_file_called = .true.
      end if
    end if

! Generate the output field of main variables:
 if (itsm1-ArchvStep0.ge.0) 			& ! skip time-lag prior to official start
      call output (itsm1, nts,                  &
             us3d, vs3d, dp3d, sdot,            & ! west & south wind, layer thickness, vert.velocity
             pr3d, ex3d, mp3d,                  & ! pressure, exner fct, montg.pot.,
             tr3d, rh3d, relvor, ws3d,          & ! field tracers, rel. humidity, vorticity, omega
             chem_opt,diaga, diagb,             & ! diagnostic arrays
             ph3d, tk3d, rn2d, sn2d,rc2d,pw2d,pq2d, &
             ts2d, us2d, rsds, rlds, rsus, rlus,&
             rsdt, rsut, rlut, qf2d, hf2d,      &
             st3d, sm3d, t2m2d, q2m2d, canopy2d,&
             fice2d, hice2d, sheleg2d, slmsk2d, &
             u10m, v10m, rn2d0, rc2d0, rg2d0,   &
             TimingBarriers, curr_write_time,   &
             wt_int_rev,k_int_rev)

! if done like this it should be separated into wrfphys and chem
    if (chem_opt > 0) then ! output chem variables
      call chem_output(itsm1,nts,aod2d,exch,p10,pm25,pr3d,tk3d,tr3d,trfall,phys2dwrf,tr1_tavg)
    end if
! if done like this it should be separated into wrfphys and chem
    if (mp_physics > 0 .or. cu_physics > 0) then ! output wrfphys variables
      call wrf_output(itsm1,nts,pr3d,tk3d,tr3d,phys2dwrf)
    end if

! See the comment accompanying the post_init_file() call.

    if (post_init_file_called) then
      call post_finalize_file(ret)
      if (ret.ne.0) then
        write(6,*) 'output: bad return from post_finalize_file: stopping'
        stop
      endif
      post_init_file_called = .false.
    end if

    ret = gptlstop ('Output')
    ret = gptlstart ('Dynamics')

!TBH:  Restore if statement label once Mark fixes PPP
!TBH skip_first_time2: if (its > itsStart ) then
    if (its > 1) then
      ret = gptlstart ('profout')
      call profout(                 &
           itsm1,PrintIpnDiag,          & ! index time step
           us3d,vs3d,dp3d,              & ! west & south wind, layer thickness
           pr3d,tr3d(:,:,1),mp3d,tk3d,  & ! pressure, pot.temp., montg.pot, temp (k)
           tr3d(:,:,2),rh3d,            & ! specific and relative humidity
           ph3d,                        & ! geopot.
           ts2d,us2d,hf2d,qf2d  )         ! skin temp., ustar, snsbl heatflx, vapor flx
      ret = gptlstop ('profout')

      if (mod (itsm1,TestDiagProgVars) == 0 .or. itsm1 == 1 .or. itsm1 == itsStart+nts	&
	.or. ((itsm1==itsDFI.or.itsm1==itsDFI+nts).and.digifilt)) then
        time = its2time(itsm1)
        write (6,*)'precip - total/nonconv/conv=',qdtr,qdtrn,qdtrc
        write (6,80)'Global 3D mass          =',qmass,' at ',time,ArchvTimeUnit,', time step=',itsm1
        write (6,80)'Global 3D water vapor   =',qmsqv,' at ',time,ArchvTimeUnit,', time step=',itsm1
        write (6,80)'Global 3D ozone         =',qmso3,' at ',time,ArchvTimeUnit,', time step=',itsm1
        write (6,80)'Global 3D cloud water   =',qmsqw,' at ',time,ArchvTimeUnit,', time step=',itsm1
        write (6,80)'Global integ acc precip =',qmstr,' at ',time,ArchvTimeUnit,', time step=',itsm1

        if (qdtr_print) then
          write (6,100) PrintDiagProgVars,qdtr,time,ArchvTimeUnit,itsm1
100       format ('Global precip, last',i3,'h =',es14.7,' at ',f5.1,1x,a,', time step=',i8)
        else
          qdtr_print = .true.
        end if

        write (6,80)'Global integ evaporation=',qmste,' at ',time,ArchvTimeUnit,', time step=',itsm1
        do t=1,ntrb
          write (6,'(a,i2.2,a,es15.7,a,i12)') 'Tracer B',t,        &
               ' global amount =',qtrcr(ntra+t),'  at time step=',itsm1
        end do
      end if
80    format (a,es14.7,a,f5.1,1x,2a,i8)

!TODO:  move this down into a subroutine...  
    end if
!TBH:  Restore if statement label once Mark fixes PPP
!TBH endif skip_first_time2

    !...........................................................
    ! Begin current dynamics time step unless this is the last 
    ! (itsStart+nts+1) iteration.  
    ! This complexity is required for the NCEP ESMF approach 
    ! in which single-phase DYN and PHY components alternate 
    ! execution during each time step.  
    !
    !TBH:  Restore if statement label once Mark fixes PPP
    !TBH skip_last_iteration: if (its <= nts ) then

    if (digifilt) then
      laststep = itsDFI + nts
    else
      laststep = itsStart + nts
    end if

    if (its < laststep) then
!     if (mp_physics /= 0 .or. cu_physics /= 0) then
        !
        ! store this time level t,qv
        !
!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (k) schedule (static)
        do ipn=ips,ipe
          do k=1,nvl
            phys3dwrf(k,ipn,3) = tr3d(k,ipn,2)
            phys3dwrf(k,ipn,7) = tr3d(k,ipn,1)*ex3d(k,ipn)/(1.+.6078*tr3d(k,ipn,2))/cp
!     if(ipn.eq.5000)then
!        write(6,*)'dyn0',k,phys3dwrf(k,ipn,7),ex3d(k,ipn)
!     endif
          end do
        end do
!$OMP END PARALLEL DO
!sms$ignore end
!     end if

!sms$compare_var(st3d   , "dyn_run begin iteration ")
!sms$compare_var(sm3d   , "dyn_run begin iteration ")
!sms$compare_var(rn2d   , "dyn_run begin iteration ")
!sms$compare_var(rc2d   , "dyn_run begin iteration ")
!sms$compare_var(ts2d   , "dyn_run begin iteration ")
!sms$compare_var(us2d   , "dyn_run begin iteration ")
!sms$compare_var(hf2d   , "dyn_run begin iteration ")
!sms$compare_var(rsds   , "dyn_run begin iteration ")

      !.............................................................
      ! Start the Adams Bashforth 3rd order time diff:
      ret = gptlstart ('abstart')
      call abstart (its,itsStart, &
           nf,of,vof,			& ! time slots for Adams Bashforth
           adbash1,adbash2,adbash3	) ! Adams Bashforth time dif. weights
      ret = gptlstop ('abstart')

      if (TimingBarriers) then
        ret = gptlstart ('dyn_run_barrier')
!SMS$BARRIER
        ret = gptlstop ('dyn_run_barrier')
      end if

      ret = gptlstart ('dyn_run_exchange')
!SMS$EXCHANGE(us3d,vs3d,dp3d,tr3d(:,:,1:ntra),mp3d,ph3d,pr3d,ex3d)
      ret = gptlstop ('dyn_run_exchange')

      !.........................................................
      ! interpolate variables to cell edges
      ret = gptlstart ('edgvar1')
      call edgvar1 (its, us3d, vs3d, dp3d, pr3d, ex3d,           &
                   ph3d, mp3d, tr3d, u_edg, v_edg, dp_edg,      &
                   lp_edg, geop_edg, mont_edg, uvsq_edg,	&
                   trc_edg)
      ret = gptlstop ('edgvar1')

      ret = gptlstart ('edgvar2')
      if (hiOrderUVComp) then
        call edgvar2 (its, us3d, vs3d, dp3d, pr3d, ex3d,           &
                   ph3d, mp3d, tr3d, u_edg, v_edg, dp_edg,      &
                   lp_edg, geop_edg, mont_edg, uvsq_edg,	&
                   trc_edg)
      endif
      ret = gptlstop ('edgvar2')

      ret = gptlstart ('edgvar3')
      if (hiOrderFluxComp) then 
        call edgvar3 (its, us3d, vs3d, dp3d, pr3d, ex3d,           &
                   ph3d, mp3d, tr3d, u_edg, v_edg, dp_edg,      &
                   lp_edg, geop_edg, mont_edg, uvsq_edg,	&
                   trc_edg)
      end if
      ret = gptlstop ('edgvar3')

!SMS$COMPARE_VAR(u_edg,"in dyn_run")
!SMS$COMPARE_VAR(v_edg,"in dyn_run")
!SMS$COMPARE_VAR(dp_edg,"in dyn_run")

      !........................................................
      ! Solve continuity equation
      if (flux_conserv_schm == 1) then
        call cnuity (its, us3d, vs3d, u_edg, v_edg,    &
                     dp_edg, lp_edg, dp3d, pr3d, ex3d, &
                     dp_tdcy, dpl_tdcy, massfx, ws3d, TimingBarriers)
      else 
        call mass_adv (its, us3d, vs3d, u_edg, v_edg,    &
                       dp_edg, lp_edg, dp3d, pr3d, ex3d, &
                       dp_tdcy, massfx, ws3d, TimingBarriers)
      endif
  !........................................................
  ! Build up time integral of mass fluxes
      if (ntrb > 0) then
        call transp1(its,             &
             nf,of,vof,                  & ! time slots for Adams Bashforth
             adbash1,adbash2,adbash3,    & ! Adams Bashforth time dif. weights
             cumufx,massfx)              ! time-integrated and instant. massfx

        !........................................................
        ! Solve transport equation for class B tracers
        if (mod(its,dtratio) == 0) then
          call transp2 (its,            &
               tr3d, cumufx,   & ! tracer, time-integrated mass flux
               dpinit, dp3d,   & ! initial & final thickness
               TimingBarriers )

          ! re-initialize time integrals
          call transp0(its,cumufx,dp3d,dpinit)
        end if
      end if

  !........................................................
  ! Solve tracer transport equation for class A tracers
      call trcadv (its, trc_edg, tr3d, trdp, trc_tdcy, &
                   trl_tdcy, massfx, dp3d, TimingBarriers)

!      call traceradv (its, u_edg, v_edg,trc_edg, tr3d, trdp, trc_tdcy, &
!                        massfx, dp3d, dp_edg, TimingBarriers )
      !.........................................................   
      ! Solve momentum equation

!sms$compare_var(u_tdcy  , "before momtum - u_tndcy5 ")
!sms$compare_var(v_tdcy  , "before momtum - v_tndcy5 ")
!sms$compare_var(ex3d    , "before momtum - exner5 ")
!sms$compare_var(mont_edg, "before momtum - mont_edg5 ")

      ret = gptlstart ('momtum')
      call momtum (its, us3d, vs3d, ex3d, relvor,		&
                   u_edg, v_edg, trc_edg, geop_edg, mont_edg,	&
                   uvsq_edg, u_tdcy, v_tdcy, dp3d)
      ret = gptlstop  ('momtum')

!sms$compare_var(u_tdcy  , "after momtum - u_tdcy6 ")
!sms$compare_var(v_tdcy  , "after momtum - v_tdcy6 ")
!sms$compare_var(us3d    , "after momtum - us3d ")
!sms$compare_var(vs3d    , "after momtum - vs3d ")
!sms$compare_var(trdp    , "dyn_run - trdp5 ")
!sms$compare_var(trc_tdcy, "dyn_run - trc_tdcy5 ")


      !.........................................................
      ! Save theta at end of dyn_run but before physics_run
      !JR This needs to be threaded
      worka(:,:)=tr3d(:,:,1)

!sms$compare_var(st3d   , "dyn_run end first half of iteration ")
!sms$compare_var(sm3d   , "dyn_run end first half of iteration ")
!sms$compare_var(rn2d   , "dyn_run end first half of iteration ")
!sms$compare_var(rc2d   , "dyn_run end first half of iteration ")
!sms$compare_var(ts2d   , "dyn_run end first half of iteration ")
!sms$compare_var(us2d   , "dyn_run end first half of iteration ")
!sms$compare_var(hf2d   , "dyn_run end first half of iteration ")
!sms$compare_var(rsds   , "dyn_run end first half of iteration ")

    end if
!TBH:  Restore if statement label once Mark fixes PPP
!TBH endif skip_last_iteration

!TODO:  move this down into a subroutine...  
!   if (mp_physics /=0 .or. cu_physics /= 0) then
      !
      ! get dynamic tendencies for th and qv
!sms$ignore begin
!$OMP PARALLEL DO PRIVATE (k) schedule (static)
      do ipn=ips,ipe
        do k=1,nvl
!     if(ipn.eq.2)then
!        write(6,*)'dyn',dt,tr3d(k,ipn,1)* (ex3d(k,ipn)/(1.+.6078*tr3d(k,ipn,2))/cp),phys3dwrf(k,ipn,7)
!        write(6,*)'dyn2',k,tr3d(k,ipn,2),phys3dwrf(k,ipn,3),ex3d(k,ipn)
!     endif
          phys3dwrf(k,ipn,3) = (tr3d(k,ipn,2) - phys3dwrf(k,ipn,3))/dt
          phys3dwrf(k,ipn,7) = (tr3d(k,ipn,1)*(ex3d(k,ipn)/(1. + .6078*tr3d(k,ipn,2)))/cp - &
               phys3dwrf(k,ipn,7))/dt
!        if(ipn.eq.2)write(6,*)'dyn2',ipn,k,dt,phys3dwrf(k,ipn,7)
        end do
      end do
!$OMP END PARALLEL DO
!sms$ignore end
!   endif

!sms$compare_var(phys3dwrf   , "dyn_run end first half of iteration ")

    ret = gptlstop ('Dynamics')

    return
  end subroutine dyn_run
end module module_fim_dyn_run
