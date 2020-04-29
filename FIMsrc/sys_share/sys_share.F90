! Substitutes for non-standard system intrinsic subroutines.  
!JR Why is an empty module with a (also possibly empty) subroutine appended 
!JR to the end needed? Mac complained about empty file, so added a stub.

module module_sys_share
  implicit none
end module module_sys_share

! Substitute for non-standard flush() intrinsic subroutine.  
!
! If your system does not support flush(), #define NOFLUSH 
! and build this proxy.  

#if (defined NAG || defined GFORTRAN)
real function timef()
  implicit none
  timef = 0.
end function timef
#endif

#ifdef NAG
subroutine flush(lun)
  implicit none
  integer,intent(in)::lun
  flush(lun)
end subroutine flush

subroutine abort
  implicit none
  stop 999
end subroutine abort

! NOTE: The following erf() function is needed only for pre-v6.0 NAG compilers,
!       which do not supply an ERF intrinsic. It should be removed when a newer
!       NAG compiler is available.

real function erf(x)

! See https://en.wikipedia.org/wiki/Error_function#Numerical_approximation

  implicit none

  real,intent(in)::x
  real::t,tau

  t=1/(1+0.5*abs(x))
  tau=t*exp(         &
    -x**2            &
    -1.26551223      &
    +1.00002368*t    &
    +0.37409196*t**2 &
    +0.09678418*t**3 &
    -0.18628806*t**4 &
    +0.27886807*t**5 &
    -1.13520398*t**6 &
    +1.48851587*t**7 &
    -0.82215223*t**8 &
    +0.17087277*t**9 &
    )
  if (x<0) then
    erf=tau-1
  else
    erf=1-tau
  end if

end function erf

#endif /* NAG compiler */

#ifdef IBM
subroutine flush(lun)
  implicit none
  integer,intent(in)::lun
  call flush_(lun) ! for IBM, generalize later if needed
end subroutine flush
#endif

subroutine stub_to_satisfy_linkers_that_cant_handle_empty_doto_files
end subroutine stub_to_satisfy_linkers_that_cant_handle_empty_doto_files
