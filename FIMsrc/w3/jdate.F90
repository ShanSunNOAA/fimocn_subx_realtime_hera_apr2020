program jdate
#ifdef NAG
USE F90_UNIX_ENV, ONLY : getarg
#endif
implicit none
CHARACTER(len=12) :: yyyymmddhhmm
CHARACTER(len=9) :: JulianDate 

call getarg(1,yyyymmddhhmm)
call GetJdate(yyyymmddhhmm,JulianDate)
#ifdef NAG
write(*, '(A9)', ADVANCE = "NO") JulianDate
#else
print'(A9,$)',JulianDate
#endif
stop
end program jdate
