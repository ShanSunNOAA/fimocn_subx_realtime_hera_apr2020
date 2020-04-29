module module_fim_init

implicit none

contains

!*********************************************************************
!       Loads the initial variables and constants to start fim
!       Alexander E. MacDonald  11/27/04
!       J. Lee                  September, 2005
!*********************************************************************

subroutine init(call_core_setup_in)
  use module_core_setup   ,only: core_setup_fim, iam_fim_task, iam_write_task
  use module_fim_chem_init,only: chem_init
  use module_fim_cpl_init ,only: cpl_init
  use module_fim_dyn_init ,only: dyn_init
  use module_fim_phy_init ,only: phy_init
  use module_wrf_phy_init ,only: wrf_phy_init
  use fimnamelist         ,only: readrestart
  use restart             ,only: read_restart
  use units               ,only: initunit
  use hycom_init          ,only: ocn_init

  implicit none
! Optional arguments
  logical, intent(in), optional :: call_core_setup_in

! Local variables
  integer :: big_endian_units(5) = (/11,21,28,30,82/)
  logical :: call_core_setup           ! call core_setup_fim?

!JR: Sit in a spin-wait loop so a debugger can attach, halt the process,
!JR  reset the variable spinwait, and then continue.
!JR: Placed here because MPI is now active (if enabled), but the model is
!JR  still in the startup phase.

#ifdef ATTACH_DEBUGGER
  integer :: spinwait = 1

  do while (spinwait == 1)
  end do
#endif

! Initialize units handler with list of big-endian units
  call initunit (big_endian_units(:))

! When MPI is used, set up communicators for compute tasks and 
! optional write tasks.  Mirrors NEMS approach.  
! NOTE:  Executable SMS directives must not be placed before 
! NOTE:  this call!  
! NOTE:  This includes writes or prints without !SMS$ignore because they 
! NOTE:  cause SMS to generate code.

! If call_core_setup is .true., call core_setup_fim from here for FIM.
! If call_core_setup is .false., core_setup_fim is called by NEMS.

  call_core_setup = .true.
  if (present(call_core_setup_in)) call_core_setup = call_core_setup_in


  if (call_core_setup) then
    call core_setup_fim ()
  end if

  if (iam_fim_task .or. iam_write_task) then
    call dyn_init ()
  end if

  if (iam_fim_task) then
    call phy_init ()
    call chem_init ()
    call wrf_phy_init ()
    call cpl_init ()
    call ocn_init ()		! check atmonly/coupled/ocnonly & do ocean initialization
    if (readrestart) then
      call read_restart ()      ! read dynamics, physics & ocean data from restart file
    end if
  end if
  
  print *,'... exiting init'
  return
end subroutine init

end module module_fim_init
