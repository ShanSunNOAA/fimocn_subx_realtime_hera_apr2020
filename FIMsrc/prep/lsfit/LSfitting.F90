!-------------------------------------------------------------------------------
!  This file contains a module that implement the least square fitting, 
!  using polynomial terms of bi-quadratic function as basis functions.
!
!  The least square problem is solved through QR decomposition.
!  
!  For computational efficiency, QR decomposition is precomputed and
!  saved for multiple right hand side values.
!
!
!  Ning Wang, Feb.. 2009
!  
!-------------------------------------------------------------------------------
!#define ORIG
MODULE LSfitting
   IMPLICIT NONE

   INTEGER diag

   INTEGER nBasisFuns, nObs

   TYPE Site2d 
     REAL :: x
     REAL :: y
   END TYPE Site2d

!SMS$DISTRIBUTE(dh,2) BEGIN
   REAL, ALLOCATABLE :: Amtxs(:,:)
!SMS$DISTRIBUTE END

   PUBLIC :: initLSmod, closeLSmod, LSQRdecomp, LSsolveEvalShift,&
             pseudoInverse, solve,  solveiThLS, check_solution, Site2d

   PUBLIC :: Amtxs

   PRIVATE  
        
CONTAINS

   SUBROUTINE initLSmod(maxNobs, maxNbfs, npts)

     INTEGER maxNobs, maxNbfs, npts 

     ALLOCATE(Amtxs((maxNobs+1) * maxNbfs, npts))
     Amtxs = 0.0

   END SUBROUTINE initLSmod

   SUBROUTINE closeLSmod()

     DEALLOCATE(Amtxs)

   END SUBROUTINE closeLSmod

   SUBROUTINE LSQRdecomp(xy, num_obs, num_basisfuns, ith)
     INTEGER, INTENT(IN):: num_obs, num_basisfuns, ith
     REAL, INTENT(IN) :: xy(2,num_obs)

     INTEGER :: j
     TYPE(Site2d) :: st2d(num_obs)

     DO j=1,num_obs
       st2d(j)%x = xy(1,j)
       st2d(j)%y = xy(2,j)
     END DO

     CALL setnOBSnBF(num_obs,num_basisfuns)
     CALL constructiThObsMat(st2d,ith)
     CALL preSolveiThLS(ith)

   END SUBROUTINE LSQRdecomp

   SUBROUTINE pseudoInverse(xy, num_obs, num_basisfuns, inv)
     INTEGER, INTENT(IN):: num_obs, num_basisfuns
     REAL, INTENT(IN) :: xy(2,num_obs)
     REAL, INTENT(OUT) :: inv(num_basisfuns, num_obs)

     INTEGER :: i, j, k
     TYPE(Site2d) :: st2d(num_obs)

     REAL :: Amtx(num_obs,num_basisfuns)
     REAL :: sing(num_basisfuns), cond
     INTEGER :: rank, lwork, info

     REAL :: u(num_obs,num_obs), vt(num_basisfuns,num_basisfuns)
     REAL :: vsigma(num_basisfuns, num_obs)
     INTEGER, PARAMETER :: sz_wkarray = 200 ! optimal work array size 
     REAL :: work(sz_wkarray)


     DO j=1,num_obs
       st2d(j)%x = xy(1,j)
       st2d(j)%y = xy(2,j)
     END DO

! construct matrix A
     CALL setnOBSnBF(num_obs,num_basisfuns)
     CALL constructObsMat(st2d,Amtx)

     lwork = sz_wkarray     ! size of the work array
     CALL sgesvd ('A', 'A', num_obs, num_basisfuns, Amtx, num_obs, &
                  sing, u, num_obs, vt, num_basisfuns, work, lwork, info)

#ifdef DIAG
     PRINT*, 'info=', info, 'lwork=',work(1)
     PRINT*, 'sing=',sing
     PRINT*, 'u:'
     PRINT*, u
     PRINT*, 'vt:'
     PRINT*, vt
#endif

! construct pseudo inverse matrix V*Sigma^{-T}*U^T
     vsigma = 0.0
     DO i = 1, num_basisfuns
       DO j = 1, num_basisfuns
         vsigma(i,j) = vt(j, i) / sing(j) 
       ENDDO
     ENDDO

     DO i = 1, num_basisfuns
       DO j = 1, num_obs
         inv(i,j) = 0.0
         DO k = 1, num_obs 
           inv(i,j) = inv(i,j) + vsigma(i, k) * u(j, k) 
         ENDDO
       ENDDO
     ENDDO
     
   END SUBROUTINE pseudoInverse

   SUBROUTINE solve(inv, b, x)
     REAL, INTENT(IN) :: inv(nBasisFuns, nobs)
     REAL, INTENT(IN) :: b(nobs)
     REAL, INTENT(OUT) :: x(nBasisFuns)
    
     INTEGER :: i, j

     DO i = 1, nBasisFuns
        x(i) = 0.0
        DO j = 1, nobs
          x(i) = x(i) + inv(i,j)*b(j)
        ENDDO
      ENDDO
    END SUBROUTINE solve


   SUBROUTINE LSsolveEvalShift(obsval, edgval, n_nb, ith, pxy)
     INTEGER, INTENT(IN):: ith, n_nb
     REAL, INTENT(IN) :: obsval(n_nb+1)
     REAL, INTENT(IN) :: pxy(5,n_nb)
     REAL, INTENT(OUT) :: edgval(n_nb)

     INTEGER :: j
     REAL :: x, y, b(6)
    
     nObs = n_nb 
     nBasisFuns = 5
    
     DO j = 1, nObs
       b(j) = obsval(j)-obsval(n_nb+1)
!       b(j) = obsval(j)
     END DO

! solve the least square fitting problem
     CALL solveiThLS(b, ith)

! evaluate the values for all edges
     DO j = 1, n_nb
#ifndef ORIG
       edgval(j) = obsval(n_nb+1)+b(1)*pxy(1,j)+b(2)*pxy(2,j)+ &
         b(3)*pxy(3,j)+b(4)*pxy(4,j)+b(5)*pxy(5,j) 
#else
       edgval(j) = obsval(n_nb+1)+b(1)*pxy(3,j)+b(2)*pxy(5,j)+ &
         b(3)*pxy(4,j)+b(4)*pxy(1,j)+b(5)*pxy(2,j) 
#endif
     END DO

   END SUBROUTINE LSsolveEvalShift

! private subroutines for internal use 
   SUBROUTINE setnOBSnBF(nOb, nbf)
     INTEGER nOb, nbf, ith
 
     nObs = nOb
     nBasisFuns = nbf

   END SUBROUTINE setnOBSnBF
     
   SUBROUTINE constructObsMat(st2d, Amtx) 
     TYPE(Site2d) st2d(nObs) 
     REAL :: Amtx(:,:)

     INTEGER i, j, sizeA

     sizeA = nObs * nBasisFuns
     DO i = 1, nObs
       DO j = 1, nBasisFuns
         Amtx(i,j) = Quadratic(st2d(i), j);
       ENDDO
     ENDDO

   END SUBROUTINE constructObsMat

   SUBROUTINE constructiThObsMat(st2d, ith) 
     TYPE(Site2d) st2d(nObs) 
     INTEGER ith

     INTEGER i, j, sizeA

     sizeA = nObs * nBasisFuns
     DO i = 1, nObs
       DO j = 1, nBasisFuns
         Amtxs(i + (j-1) * nObs, ith) = Quadratic(st2d(i), j);
       ENDDO
     ENDDO

   END SUBROUTINE constructiThObsMat

   REAL FUNCTION Quadratic(site, bi)
      TYPE(Site2d) site
      INTEGER bi

      SELECT CASE (bi)
          CASE (1)   
#ifndef ORIG
              Quadratic = site%x 
#else
              Quadratic = site%x * site%x
#endif
   
          CASE (2)   
#ifndef ORIG
              Quadratic = site%y
#else
              Quadratic = site%y * site%y
#endif

          CASE (3)   
#ifndef ORIG
              Quadratic = site%x * site%x
#else
              Quadratic = site%x * site%y
#endif

          CASE (4)
#ifndef ORIG
              Quadratic = site%x * site%y
#else
              Quadratic = site%x 
#endif
   
          CASE (5)
#ifndef ORIG
              Quadratic = site%y * site%y
#else
              Quadratic = site%y 
#endif
   
          CASE (6)
              Quadratic = 1

          CASE default  
              PRINT*,  'Wrong basis function index!'

      END SELECT

   END FUNCTION Quadratic

   SUBROUTINE preSolveiThLS(ith)  
     INTEGER ith
     INTEGER i, m, n, lda, lwk, info, sizeA
     REAL, ALLOCATABLE :: wk(:) 
     REAL cond

     m = nObs;
     n = nBasisFuns;
     lda = m;
     lwk = m * n
     sizeA = m * n
     ALLOCATE(wk(lwk)) 

     CALL sgeqr2(m, n, Amtxs(1,ith), lda, wk(1), wk(n+1), info)
     DO i = 1, n
       Amtxs(sizeA + i, ith) = wk(i)
     ENDDO
     DEALLOCATE(wk)

   END SUBROUTINE preSolveiThLS

   SUBROUTINE solveiThLS(b, ith) 
     REAL, INTENT(INOUT) :: b(:)
     INTEGER, INTENT(IN) :: ith
     INTEGER m, n, lda, ldb, lwk, info
     INTEGER sizeA 

     INTEGER i, j
     REAL, ALLOCATABLE :: wk(:)

     lwk = nObs * nBasisFuns
     ALLOCATE(wk(lwk)) 

     m = nObs;
     n = nBasisFuns;
     lda = m;
     ldb = m;
     lwk =  n * m;

     sizeA = nObs * nBasisFuns
     CALL mySormqr(m, n, Amtxs(1, ith), lda, &
              Amtxs(sizeA+1, ith), b, ldb)
     CALL myStrtrs(n, Amtxs(1, ith), lda, b)  
     IF (diag == 2) THEN
       PRINT*, 'Solved a system of ',  m,  ' x ',  n
       PRINT*, 'info = ', info
       DO i = 1, nBasisFuns
         PRINT*, 'b(i)= ', b(i) 
       ENDDO
     ENDIF
     DEALLOCATE(wk)

   END SUBROUTINE solveiThLS 
 
   SUBROUTINE mySormqr(m, n, A, lda, tau, b, ldb)
     INTEGER m, n, lda, ldb
     INTEGER i, j, mi
     REAL A(lda, *), b(*), tau(*)
     REAL aii, wk

     DO i = 1, n
       mi = m - i + 1
       aii = A(i, i)
       A(i, i) = 1.0
       wk = 0.0
       DO j = 1, mi
         wk = wk + A(i+j-1, i) * b(i+j-1)
       ENDDO
       b(i:m) = b(i:m) -  tau(i) * wk * A(i:m, i) 
       A(i, i) = aii
     ENDDO
   END SUBROUTINE mySormqr

   SUBROUTINE myStrtrs(n, A, lda, b)
     INTEGER n, lda
     REAL A(lda, *), b(*)
     INTEGER i, j
     
     DO j = n, 1, -1
       IF (b(j) /= 0.0) THEN
         b(j) = b(j) / A(j,j) 
         DO i = 1, j-1
           b(i) = b(i) - b(j)*A(i, j)
         ENDDO
       ENDIF
     ENDDO 
   END SUBROUTINE myStrtrs

   SUBROUTINE mulHb(m, v, tau, b, wk) 
     INTEGER m, i
     REAL v(*), b(*), wk(*) 
     REAL tau

     wk(1) = 0.0
     DO i = 1, m
       wk(1) = wk(1) + v(i) * b(i)
     ENDDO

     b(1:m) = b(1:m) - tau * wk(1) * v(1:m) 

   END SUBROUTINE mulHb        

   SUBROUTINE check_solution(st2d,b,b2)
      TYPE(Site2d) st2d(nObs)
      REAL b(nBasisFuns),b2(nObs)

      REAL s, c(nObs)
      INTEGER i, j

      PRINT*, 'Number of basis functions:',nBasisFuns
      PRINT*, 'Nunmber of observations:',nObs
      s = 0.0
      DO i = 1, nObs
        c(i) = 0;
        DO j = 1, nBasisFuns
          c(i) = c(i) +  Quadratic(st2d(i), j) * b(j) 
        ENDDO
        PRINT*, 'c(i)= ', c(i), ', b(i)= ', b2(i)
        s = s + (c(i) - b2(i)) * (c(i) - b2(i));
      ENDDO
      PRINT*, 'L2 Residual = ',  sqrt(s)
   END SUBROUTINE check_solution

END MODULE
