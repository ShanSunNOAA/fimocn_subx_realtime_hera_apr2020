!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module restart contains these public entities:
!   write_restart: writes a restart file
!   read_restart: reads a restart file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module restart

  USE fimnamelist,    only: itsstart,atmonly,coupled
  use units,          only: getunit, returnunit
  use hycom_variables,          only: nstep_ocn

  implicit none

#include <gptl.inc>

  private
  public :: write_restart, read_restart

contains

! write_restart: open the output restart file, then call write_restart_dyn,
! write_restart_cpl, and write_restart_phy

  subroutine write_restart (itsm1)

    use restart_cpl,  only: write_restart_cpl
    use restart_dyn,  only: write_restart_dyn
    use restart_phy,  only: write_restart_phy
    use hycom_restart,only: write_restart_ocn

    integer, intent(in) :: itsm1     ! current time step to write to restart file

    character(len=16) :: rfilename   ! Restart file name
    integer           :: unitno      ! Fortran unit number to write to
    integer           :: ret         ! function return
    character(len=8)  :: funcval     ! for Lahey

    character(len=8), external :: its2string ! Converts an integer to a string

    write (6,*) 'Entered write_restart itsm1=', itsm1
    ret = gptlstart ('write_restart')

!JR The following causes Lahey to fail! Looks like a compiler bug, so code around it:
!    write (rfilename, "('restart_',a8)") its2string (itsm1)
    funcval = its2string (itsm1)
    write (rfilename, "(a8,a8)") 'restart_', funcval

    unitno = getunit ()
    if (unitno < 0) then
      write(6,*)'write_restart: Bad return from getunit'
      stop
    end if

    open (unitno, file=rfilename, form='unformatted', action='write', err=70)
    write (unitno, err=90) itsm1

    ret = gptlprint_memusage ('before write_restart')
    call write_restart_dyn (unitno)
    call write_restart_cpl (unitno)
    call write_restart_phy (unitno)
    if (.not. atmonly) call write_restart_ocn (unitno)
    ret = gptlprint_memusage ('after write_restart')

    close (unitno)

    write(6,*)'write_restart: opening file "rpointer" on unit ', unitno
    open (unit=unitno, file='rpointer', form='formatted', action='write', err=70)
    write (unitno, '(a16)', err=90) rfilename
    write(6,*)'write_restart: wrote text file rpointer containing:', rfilename
    close (unitno)

    call returnunit (unitno)

    ret = gptlstop ('write_restart')
    return

70  write(6,*)'write_restart: Error opening unit ', unitno, '. Stopping'
    stop

90  write(6,*)'write_restart: Error writing restart file. Stopping'
    stop
  end subroutine write_restart

! read_restart: open the input restart file, then call read_restart_dyn,
! read_restart_cpl, and read_restart_phy
  subroutine read_restart ()

    use restart_cpl,  only: read_restart_cpl
    use restart_dyn,  only: read_restart_dyn
    use restart_phy,  only: read_restart_phy
    use hycom_restart,only: read_restart_ocn

    integer :: itsm1       ! time step read from restart file
    integer :: unitno      ! unit number for restart file
    integer :: ret         ! function return
    character(len=16) :: rfilename   ! Restart file name

    ret = gptlstart ('read_restart')

    unitno = getunit ()
    if (unitno < 0) then
      write(6,*)'read_restart: Bad return from getunit'
      stop
    end if

    write(6,*)'read_restart: trying to open file "rpointer" on unit ', unitno

    open (unitno, file='rpointer', form='formatted', action='read', err=70)
    write(6,*)'read_restart: successfully opened restart file rpointer on unit ', unitno
    read (unitno, '(a16)', err=90) rfilename
    close (unitno)

    open (unitno, file=rfilename, form='unformatted', action='read', err=70)
    write (6,*)'read_restart: opened restart file:', rfilename
    read (unitno, err=90) itsm1

    write (6,*) 'read_restart: starting up from itsm1=', itsm1

    ret = gptlprint_memusage ('before read_restart')
    call read_restart_dyn (unitno)
    call read_restart_cpl (unitno)
    call read_restart_phy (unitno)
    if (.not. atmonly) call read_restart_ocn (unitno)
    ret = gptlprint_memusage ('after read_restart')

    close (unitno)
    call returnunit (unitno)
    write(6,*)'read_restart: returned unit ', unitno, ' to the list of available units'

    itsstart = itsm1 + 1     ! Set starting iteration counter for the restart run
    nstep_ocn = itsm1        ! set starting iteration counter for the ocean
    ret = gptlstop ('read_restart')
    return

70  write(6,*)'read_restart: Error opening rpointer. Stopping'
    stop

90  write(6,*)'read_restart: Error reading restart file. Stopping'
    stop
  end subroutine read_restart

end module restart
