program speinfo
  use sigio_module

  implicit none

  integer, parameter :: lusig = 82   ! Big endian!!!
  integer :: imax, jmax, nvp         ! gaussian grid dimensions
  integer :: maxwv                   ! max number of waves
  integer :: i, j, k                 ! gaussian grid indices
  integer :: iret                    ! Return code from sigio routines
  type(sigio_head) :: head           ! NCEP spectral data header
  type(sigio_data) :: data           ! NCEP spectral data
  character(len=120) :: sanlfile     ! spectral coefficient file

  do while (.true.)
    write(6,*) 'Enter NCEP IC file: Remember: unit 82 must be set up for BIG ENDIAN'
    read (5,*) sanlfile
    write(6,*)'Reading spectal data from file=',trim(sanlfile)
    call sigio_srohdc (lusig, sanlfile, head, data, iret)
    if (iret /= 0) then
      write(6,*)'Error from sigio_srohdc'
      continue
    end if

    nvp  = head%levs
    imax = head%lonb
    jmax = head%latb
    maxwv = head%jcap
    write (*,'(a)'   ) 'Input spectral grid is as follows:'
    write (*,'(a,i0)') 'nvp (vertical levels)   = ',nvp
    write (*,'(a,i0)') 'imax (num longitudes)   = ',imax
    write (*,'(a,i0)') 'jmax (num latitudes)    = ',jmax
    write (*,'(a,i0)') 'maxwv (max wavenumbers) = ',maxwv
    call sigio_axdata (data, iret)     ! deallocate array
  end do
end program speinfo
