!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module units contains these public entities
!   initunit: initialize the units module, passing an optionl list of (probably big-endian)
!             special unit numbers that the user must ask for.
!   getunit: obtain an available fortran unit number
!   returnunit: return a fortran unit number to the list of available units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module units

  implicit none

  private
  public :: initunit   ! unit number initialization (must be called before getunit/returnunit)
  public :: getunit    ! obtain a Fortran unit number for use
  public :: returnunit ! return a Fortran unit number to the pot

  integer, parameter :: maxunits = 99 ! Not all Fortran compilers allow more than 2 digits for units
  integer :: i                        ! index for "isinuse" array which follows
  logical, save :: isinuse(maxunits) = (/(.false.,i=1,maxunits)/) ! internal state of unit assignments
! The following units will be given only to callers that specifically ask for them, e.g. for situations
! where invoking scripts need to specify special handling via unit numbers for byte swapping
  integer, allocatable :: mustask(:)    ! list of unit numbers user must ask for (generally big-endian)
  integer :: num_mustask = -1           ! size of array mustask
  logical, save :: initcalled = .false. ! must call init routine before any other in the module

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initunit: initialize unit handler
!   Arguments: user_mustask: Optional: If present, list of Fortran unit numbers which are not
!              to be given out unless the user specifically asks for them. Generally the list
!              is unit numbers which are mentioned in calling scriptery, e.g. big-endian units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initunit (user_mustask)
    integer, intent(in), optional :: user_mustask(:) ! List of special (big-endian?) unit numbers

!sms$ignore begin
    if (initcalled) then
      write(6,*) 'initunit: initunit was already called--not doing anything'
      return
    end if

    if (present(user_mustask)) then
      num_mustask = size (user_mustask)
      allocate (mustask (num_mustask))
      mustask(:) = user_mustask(:)
    end if

    initcalled = .true.
    return
!sms$ignore end

  end subroutine initunit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! getunit: provide the caller a unit number they can use for I/O
!   Arguments:
!     unitno: Optional: If present, see if the requested unit number is available. If not available,
!             return an error code.
!   Return value: unit number to use, or -1 if a unit number cannot be provided
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function getunit (unitno)
    integer, intent(in), optional :: unitno  ! requested unit number (if present)

    integer :: i   ! index over unit numbers

!sms$ignore begin
    getunit = -1   ! initialize to bad return value

    if (.not. initcalled) then
      write(6,*) 'getunit: initunit() has not yet been called. Invoking now with no list of must-ask unit numbers.'
      call initunit ()
    end if

! If optional argument "unitno" is present, give the requestor that unit if it is available.

    if (present (unitno)) then
      if (unitno > maxunits .or. unitno < 1 .or. unitno == 5 .or. unitno == 6) then
        write(6,*) 'getunit: Unit ', unitno, ' is not valid'
        return
      end if
      
      if (isinuse (unitno)) then
        write(6,*) 'getunit: Unit ', unitno, ' is already in use'
        return
      end if
        
      isinuse (unitno) = .true.
      getunit = unitno
      return
    end if

! Do not allocate units 5 (stdin), 6 (stdout), or any of the special units normally reserved
! for byte-swapping.

    do i=1,maxunits
      if (.not. isinuse (i) .and. .not. isspecial(i)) then
        isinuse(i) = .true.
        getunit = i
        return
      end if
    end do

    write(6,*) 'getunit: No more Fortran unit numbers available!'
    return
!sms$ignore end

  end function getunit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! isspecial: private function returns whether or not unit requested is "special", i.e. 
!            whether it is 5, 6, or is a member of a user-specified list of units which 
!            must be specifically requested.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function isspecial (unitno)
    integer, intent(in) :: unitno

    integer :: i

!sms$ignore begin
    isspecial = .false.
    if (unitno == 5 .or. unitno == 6) then
      isspecial = .true.
      return
    end if

    do i=1,num_mustask
      if (unitno == mustask(i)) then
        isspecial = .true.
        return
      end if
    end do

    return
!sms$ignore end

  end function isspecial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! returnunit: return unitno to the pile of available units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine returnunit (unitno)
    integer, intent(in) :: unitno

!sms$ignore begin
    if (.not. initcalled) then
      write(6,*) 'returnunit: initunit() must be called before returnunit'
      return
    end if

    if (unitno > maxunits .or. unitno < 1) then
      write(6,*) 'returnunit: Unit ', unitno, ' is not valid'
      return
    end if

    if (.not. isinuse(unitno)) then
      write(6,*) 'returnunit: WARNING--unit ', unitno, ' is not in use'
      return
    end if

    isinuse(unitno) = .false.
    return
!sms$ignore end

  end subroutine returnunit

end module units
