character(len=*) function filename(tag,its)
  use fimnamelist,only:ArchvTimeUnit
  implicit none
  character(len=*),intent(in)::tag
  integer,intent(in)::its
  character(len=6)::timestr
  real::its2time
!sms$ignore begin
  write (timestr,'(i6.6)') nint(its2time(its))
!sms$ignore end
  filename='fim_out_'//trim(tag)//timestr//ArchvTimeUnit
end function filename

character(len=*) function flexflnm(tag,its)
! --- same as subr.filename but without the hardwired 'fim_out_' part.
! --- file name returned by flexflnm starts with string 'tag'.
  use fimnamelist,only:ArchvTimeUnit
  implicit none
  character(len=*),intent(in)::tag
  integer,intent(in)::its
  character(len=6)::timestr
  real::its2time
!sms$ignore begin
  write (timestr,'(i6.6)') nint(its2time(its))
!sms$ignore end
  flexflnm=trim(tag)//timestr//ArchvTimeUnit
end function flexflnm
