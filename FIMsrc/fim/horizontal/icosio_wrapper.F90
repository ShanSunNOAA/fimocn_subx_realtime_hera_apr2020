module icosio_wrapper
  use post,           only: nvars, post_init_readnl_called
  use icosio,         only: icosio_out
  use module_control, only: nip, filename_len
  use fimnamelist,    only: fimout, gribout, var_list, only_write_var_list_entries
  use module_header,  only: header

  implicit none
  
  private
  public :: maybe_write

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! maybe_write
!   Wrapper code determines whether the data need to be written to disk 
!   before actually performing the icosio_out call. This can avoid significant
!   MPI overhead associated with gathering the data for output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine maybe_write (its, name, field, nlev, scalefactor, accum_start, twodfile)
    integer, intent(in) :: its           ! iteration
    character(len=*), intent(in) :: name ! 4-character field name
    integer, intent(in) :: nlev          ! number of vertical levels
!sms$distribute (dh,2) begin
    real, intent(in) :: field(nlev,nip)  ! data to be written
!sms$distribute end
    real, intent(in), optional :: scalefactor    ! some GRIB fields need a scale factor applied
    integer, intent(in), optional :: accum_start ! some GRIB fields need an accumulation count
    logical, intent(in), optional :: twodfile    ! put on 2d file
  
    logical :: twodfile_local          ! whether to write to file prefix "2D__"
    character(len=filename_len) :: fn  ! file name to write to

    integer :: time                    ! its converted to ArchvTimeUnit units
! These 2 things are here to avoid permutation hell in deciding whether to pass
! optional arguments down to icosio_out
    real :: scalefactor_local
    integer :: accum_start_local

    real, external :: its2time
    character(len=filename_len), external :: filename

    if (.not.fimout.and..not.gribout) return

    scalefactor_local = 1.
    accum_start_local = -1
    twodfile_local = .false.
    
    if (present (scalefactor)) then
      scalefactor_local = scalefactor
    end if

    if (present (accum_start)) then
      accum_start_local = accum_start
    end if

    if (present (twodfile)) then
      twodfile_local = twodfile
    end if

    time = nint(its2time (its))

    if (isonlist (name) .or. (fimout .and. .not. only_write_var_list_entries)) then
      if (twodfile_local) then
        fn = filename ('2D__', its)
      else
        fn = filename (name, its)
      end if

      write(6,*)'maybe_write: calling icosio_out for field ', name
      call icosio_out (its, time, name, field, fn, header (name, nlev, its), &
                       scalefactor=scalefactor_local, accum_start=accum_start_local)
    else
      write(6,*)'maybe_write: skipping ', name, ' because:'
      write(6,*)'.not. (isonlist .or. (fimout .and. .not. only_write_var_list_entries))'
    end if
    
    return

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! isonlist: 
!   determine whether the input field name is on the list of GRIB variables
!   to be written.
!   NOTE: This is a stupid linear search. If performance is a problem,
!   it should be rewritten to do a binary search through a sorted list!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical function isonlist (name)
      character(len=*), intent(in) :: name
      
      integer :: i

      if (.not. post_init_readnl_called) then
        write(6,*) 'isonlist: post_init_readnl_called is false. Must be true before var_list can be checked'
        stop
      end if

      isonlist = .false.
      do i=1,nvars
        if (name == var_list(i)) then
          isonlist = .true.
        end if
      end do
      
      return
    end function isonlist
  end subroutine maybe_write
end module icosio_wrapper
