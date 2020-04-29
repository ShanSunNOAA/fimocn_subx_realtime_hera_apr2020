module module_header

  use module_control,only:nip
  use fimnamelist   ,only:ArchvTimeUnit,glvl,yyyymmddhhmm

  implicit none
  public header,headero
  private

! Parameters

  integer,parameter::header_cols=80
  integer,parameter::header_rows=10
  integer,parameter::header_len=header_cols*header_rows

! Module variables

  character(len=header_cols)::line
  character(len=header_len)::h
  integer::pos

contains

  function header(varname,levels,its)
    implicit none
    character*(*),intent(in)::varname
    character(len=header_len)::header
    integer,intent(in)::its,levels
    real::its2time
    pos=1
    h=' '
    write (line,'(a,a,a,a)') 'FIM ',varname,' Forecast initial time YYYYMMDDHHMM: ',yyyymmddhhmm
    call append
    write (line,'(a,i0,a,i0,a,i0,a,i0,a,a)') 'Level ',levels,', GLVL= ',glvl,', Step ',its,', ',nint(its2time(its)),' ',ArchvTimeUnit
    call append
    write (line,'(a,i0,a,i0)') 'dim1=',levels,', nip=',nip
    call append
    write (line,'(i0)') 4
    call append
    write (line,'(i0)') 5
    call append
    write (line,'(i0)') 6
    call append
    write (line,'(i0)') 7
    call append
    write (line,'(i0)') 8
    call append
    write (line,'(i0)') 9
    call append
    write (line,'(i0)') 10
    call append
    header=h
  end function header

  function headero(varname,levels,its)
    implicit none
    character*(*),intent(in)::varname
    character(len=header_len)::headero
    integer,intent(in)::its,levels
    real::its2time
    pos=1
    h=' '
    write (line,'(a,a,a,a)') 'OCN ',varname,' Forecast initial time YYYYMMDDHHMM: ',yyyymmddhhmm
    call append
    write (line,'(a,i0,a,i0,a,i0,a,i0,a,a)') 'Level ',levels,', GLVL= ',glvl,', Step ',its,', ',nint(its2time(its)),' ',ArchvTimeUnit
    call append
    write (line,'(a,i0,a,i0)') 'dim1=',levels,', nip=',nip
    call append
    write (line,'(i0)') 4
    call append
    write (line,'(i0)') 5
    call append
    write (line,'(i0)') 6
    call append
    write (line,'(i0)') 7
    call append
    write (line,'(i0)') 8
    call append
    write (line,'(i0)') 9
    call append
    write (line,'(i0)') 10
    call append
    headero=h
  end function headero

  subroutine append
    implicit none
    integer::i,j
    if (pos.ge.header_len) then
      write (*,'(a)') 'ERROR: Attempt to write past end of header.'
      stop
    endif
    do i=1,len(line)
      j=pos+i-1
      h(j:j)=line(i:i)
    enddo
    pos=(int((pos+header_cols)/header_cols)*header_cols)+1
  end subroutine append

end module module_header
