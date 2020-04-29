module module_FIM_INTEGRATE
  use esmf_mod
  use module_err_msg

  use module_fim_cpl_run           ,only: cpl_run
  use module_fim_dyn_run           ,only: dyn_run
  use module_fim_phy_run           ,only: phy_run

  ! TODO:  move this to internal state
  use module_core_setup     ,only: use_write_tasks
  ! TODO:  not sure if this needs to move
  use icosio,only:icosio_out_stop,icosio_stop

  implicit none
  private

  public :: fim_integrate

CONTAINS

! only FIM compute tasks execute this routine
  subroutine fim_integrate ( clock_fim, &
                             rc_integrate) 

    type(esmf_clock),    intent(inout) :: clock_fim
    integer,             intent(  out) :: rc_integrate
!
! Local variables
!
    integer :: rc
    integer(esmf_kind_i8) :: ntimestep_esmf
    integer :: its,itsm1
    type(esmf_timeinterval) :: timestep
    type(esmf_time) :: stoptime, newstoptime

     ! Run the clock one more time step (i.e. stop after its=nts+1), then 
     ! back up one step to mimic run.F90.  
     !  * set stoptime = stoptime+dt
     !  * run "integrate" loop
     !  * set ESMF_MODE_REVERSE
     !  * advance backwards one time step
     !  * set ESMF_MODE_FORWARD
     !  * reset stoptime to its original value
     !NOTE:  This hackery works around the fact that the original 
     !NOTE:  FIM run.F90 executes one extra time step in which the 
     !NOTE:  dynamics component finishes its final computations.  This 
     !NOTE:  was required by early versions of NEMS which did not 
     !NOTE:  allow multiple run phases.  See run.F90 for a very 
     !NOTE:  detailed discussion of this issue.  
     !NOTE:  This complexity could be avoided if we allowed a 2-phase 
     !NOTE:  run method for the DYN component -- and run.F90 would 
     !NOTE:  also be simplified.  However, interoperability with other 
     !NOTE:  components would be more difficult due to potential 
     !NOTE:  mismatches in numbers of phases.  

     call esmf_clockget (clock=clock_fim, stoptime=stoptime, rc=rc)
     call err_msg (rc,'esmf_clockget(stoptime)', rc_integrate)
     call esmf_clockget (clock=clock_fim, timestep=timestep, rc=rc)
     call err_msg (rc,'esmf_clockget(timestep)', rc_integrate)
     newstoptime = stoptime + timestep
     call esmf_clockset(clock=clock_fim, stoptime=newstoptime, rc=rc)
     call err_msg (rc,'esmf_clockset(newstoptime)', rc_integrate)

     integrate: do while (.not. esmf_clockisstoptime (clock_fim, rc=rc))
      call err_msg (rc,'esmf_clockisstoptime', rc_integrate)
!
!-----------------------------------------------------------------------
!***  Advance the model one timestep
!*** TODO:  Restore wrf physics option, chem option, nophysics, and write tasks as in run.F90
!-----------------------------------------------------------------------
      call esmf_clockget (clock=clock_fim,advancecount=ntimestep_esmf,rc=rc)
      call err_msg (rc, 'get time step from clock', rc_integrate)
      itsm1 = ntimestep_esmf
      its = itsm1 + 1

!TBH:  DRY to remove duplication with run.F90
      call dyn_run (its)                     ! Dynamics run method.  
      call cpl_run (its, dyn_to_phy=.true.)  ! Coupler run method: dyn->phy
      call phy_run (its)                   ! Physics run method.
      call cpl_run (its, dyn_to_phy=.false.) ! Coupler run method: phy->dyn
!
!-----------------------------------------------------------------------
!***  Flush buffered output to write tasks and/or clear list of files
!***  written to during this output frame.      
!-----------------------------------------------------------------------
!
      call icosio_out_stop(itsm1) ! signal end of output event
!
!-----------------------------------------------------------------------
!***  Advance clock to next time step.
!-----------------------------------------------------------------------
!
      call esmf_clockadvance (clock=clock_fim, rc=rc)
      call err_msg (rc, 'advance clock', rc_integrate)

     end do integrate    ! time step loop

     call icosio_stop(itsm1) ! shut down write task(s)

     ! reset clock to state expected by caller upon return
     call esmf_clockset(clock=clock_fim, direction=ESMF_MODE_REVERSE, rc=rc)
     call err_msg (rc,'esmf_clockset(ESMF_MODE_REVERSE)', rc_integrate)
     call esmf_clockadvance(clock=clock_fim, rc=rc)
     call err_msg (rc,'esmf_clockadvance(one step backwards)', rc_integrate)
     call esmf_clockset(clock=clock_fim, direction=ESMF_MODE_FORWARD, rc=rc)
     call err_msg (rc,'esmf_clockset(ESMF_MODE_FORWARD)', rc_integrate)
     call esmf_clockset(clock=clock_fim, stoptime=stoptime, rc=rc)
     call err_msg (rc,'esmf_clockset(restore original stoptime)', rc_integrate)

  end subroutine fim_integrate
end module module_FIM_INTEGRATE
