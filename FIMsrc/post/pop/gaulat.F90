!      PROGRAM glat
!      IMPLICIT NONE
       
!      INTEGER, PARAMETER :: num_lats = 576
!      INTEGER, PARAMETER :: num_lats = 10001
!      REAL*8 gl(num_lats)
!      REAL gl4(num_lats)

!      CALL gauss_lats(gl, num_lats)
!      gl4 = gl 
!      PRINT*, gl4 

!      END PROGRAM glat
!==============================================================================
!  Subroutine gaussian_lats(gl, n_l)
!
!  N. Wang, Aug. 9, 2011.
!
!  Created partially from the code snippets in GFS code, with the following
!  changes:
!  
!  1. The stop condition for the Newton iteration has been changed to a test
!     of the maximum value of the Legendre polynomials at the current 
!     approximations of the roots.
!  2. The initial cos values of the Gaussian latitutes are computed with a 
!     simpler formula. No pre-defined array of constants are necessary.  
!  3. Some cosmtic changes are made and comments are added to the code snippets 
!     to make it more readable.
!     
!==============================================================================  
      SUBROUTINE gaussian_lats(gl,n_l)
      IMPLICIT NONE 
      INTEGER n_l
      REAL*8 gl(n_l)

      REAL cgl4(n_l)
      REAL*8 cgl(n_l)
      REAL*8 pk(n_l/2), pk_1(n_l/2), pk_2(n_l/2)
      REAL*8 dd, r, oneovdpdx, pkmax, eps, z

!      REAL*8, PARAMETER:: pi = 4.0 * atan(1.0), r2d = 180.0 / pi
!      REAL*8, PARAMETER:: c=(1.-(2./pi)**2)*0.25
      REAL*8 pi, r2d, c

      INTEGER j,jh,jhe, jmax, n

      pi = 4.0 * atan(1.0)
      r2d = 180.0 / pi
      c=(1.-(2./pi)**2)*0.25

      eps = 1.0E-09

      jmax = n_l 
      jh=jmax/2
      jhe=(jmax+1)/2
      r=1./jmax

      cgl = 0.0

      cgl(1) = cos(2.405*r)
      cgl(2) = cos(5.520*r)
      cgl(3) = cos(8.654*r)
      z = 8.654
      DO j=4,jh
        z = z + pi
        cgl(j)=cos(z*r)
      ENDDO
      pkmax=1.
      ! Use Newton method to find the roots iteratively 
      DO WHILE(pkmax > eps)
        pkmax=0.
        DO j=1,jh
          pk_1(j)=1.
          pk(j)=cgl(j)
        ENDDO
        ! Evaluate the Legendre ploynomials recursively
        DO n=2,jmax
          DO j=1,jh
            pk_2(j)=pk_1(j)
            pk_1(j)=pk(j)
            pk(j)=((2*n-1)*cgl(j)*pk_1(j)-(n-1)*pk_2(j))/n
          ENDDO
        ENDDO
        ! Derivative dPn/dx = [Pn-1(x) - xPn(x)] * n / (1 - x^2) 
        DO j=1,jh
          oneovdpdx = (1.-cgl(j)**2)/(jmax*(pk_1(j)-cgl(j)*pk(j)))
          dd=pk(j)*oneovdpdx
          cgl(j)=cgl(j)-dd
          pkmax=max(pkmax,abs(pk(j)))
        ENDDO
        PRINT*, pkmax, eps
      ENDDO

      ! To match GFS Gaussian grid values, use single precision cgl
      cgl4 = cgl
      DO j = 1, jh
        gl(j) = acos(cgl4(j)) * r2d - 90.0
        gl(jmax+1-j) = -gl(j) 
      END DO
      IF(jhe > jh) THEN
        gl(jhe)=0.0
      ENDIF

      END SUBROUTINE gaussian_lats
