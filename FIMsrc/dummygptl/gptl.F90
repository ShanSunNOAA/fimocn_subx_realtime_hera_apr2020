integer function gptlinitialize ()
  implicit none

  gptlinitialize = 0
  return
end function gptlinitialize

integer function gptlstart (name)
  implicit none

  character(len=*), intent(in) :: name

  gptlstart = 0
  return
end function gptlstart

integer function gptlstart_handle (name, handle)
  implicit none

  character(len=*), intent(in) :: name
  integer(8), intent(in) :: handle

  gptlstart_handle = 0
  return
end function gptlstart_handle

integer function gptlstop (name)
  implicit none

  character(len=*), intent(in) :: name

  gptlstop = 0
  return
end function gptlstop

integer function gptlstop_handle (name, handle)
  implicit none

  character(len=*), intent(in) :: name
  integer(8), intent(in) :: handle

  gptlstop_handle = 0
  return
end function gptlstop_handle

integer function gptlget_nregions (t, nregions)
  implicit none
  
  integer, intent(in) :: t
  integer, intent(out) :: nregions

  nregions = 0
  gptlget_nregions = 0
  return
end function gptlget_nregions

integer function gptlget_regionname (t, n, name)
  implicit none

  integer, intent(in) :: t
  integer, intent(in) :: n
  character(len=*), intent(out) :: name

  name = ' '
  gptlget_regionname = 0
  return
end function gptlget_regionname

integer function gptlget_wallclock (name, t, value)
  implicit none

  character(len=*), intent(in) :: name
  integer, intent(in) :: t
  real(8), intent(out) :: value

  value = 0.
  gptlget_wallclock = 0
  return
end function gptlget_wallclock

integer function gptlpr (mype)
  implicit none

  integer, intent(in) :: mype

  gptlpr = 0
  return
end function gptlpr

integer function gptlpr_file (file)
  implicit none

  character(len=*), intent(in) :: file

  gptlpr_file = 0
  return
end function gptlpr_file

integer function gptlprint_memusage (string)
  implicit none

  character(len=*), intent(in) :: string
  gptlprint_memusage = 0
  return
end function gptlprint_memusage

integer function gptlget_memusage (size, rss, share, text, datastack)
  implicit none

  integer, intent(out) :: size, rss, share, text, datastack
  size = 0
  rss = 0
  share = 0
  text = 0
  datastack = 0
  gptlget_memusage = 0
  return
end function gptlget_memusage

integer function gptlstamp (wall, usr, sys)
  implicit none

  real*8, intent(out) :: wall, usr, sys
  wall = 0.
  usr = 0.
  sys = 0.
  gptlstamp = 0
  return
end function gptlstamp

integer function gptlpr_summary (mpicomm)
  implicit none

  integer, intent(in) :: mpicomm
  gptlpr_summary = 0
  return
end function gptlpr_summary
