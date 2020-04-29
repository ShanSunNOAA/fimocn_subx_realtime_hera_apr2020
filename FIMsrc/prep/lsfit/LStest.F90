#define NEW_S1
      PROGRAM tstLS

      USE LSfitting
      IMPLICIT NONE

      INTEGER i, j, nob, nbf, nedge
      REAL start_time, end_time
      REAL x(7), y(7), f(7), f2(7), xy(2,7), b(7)
      TYPE(Site2d) st2d(7)
      REAL inv(5,6)

      REAL a(5)
      
      INTEGER, PARAMETER :: npts = 1
      INTEGER, PARAMETER :: n_itr = 1

      DATA x /0.0, -0.866, -0.866, 0.0, 0.866, 0.866, 0.0/

      DATA y /1.0, 0.5, -0.5, -1.0,  -0.5, 0.5, 0.0/

      DATA f /27.0, 26.5, 26.3, 27.2, 26.2, 26.4, 26.0/

      nedge = 6
      f2 = f

      DO i = 1, 7
        xy(1,i) = x(i)
        xy(2,i) = y(i)
        st2d(i)%x = x(i)
        st2d(i)%y = y(i)
      ENDDO

#ifdef NEW_S1
      nob = 6
      nbf = 5
#else
      nob = 7
      nbf = 6
#endif

      CALL initLSmod(nob, nbf, npts)
     
      CALL cpu_time(start_time)
      DO j = 1, npts
        CALL pseudoInverse(xy, nob, nbf, inv)
      ENDDO
      CALL cpu_time(end_time)
      PRINT*, "Time taken to compute pseduo inverse (single prec):" ,  &
                end_time - start_time, " seconds"

#ifdef DIAG
      PRINT*, 'inverse matrix'
      PRINT*, inv
#endif

      CALL cpu_time(start_time)
      DO j = 1,  n_itr
        b(1:6) = f(1:6)-f(7) 
        CALL solve(inv, b, a)
      ENDDO
      CALL cpu_time(end_time)

      PRINT*, 'Solution from solve() - through pseduo inverse:'
      PRINT*, a

      CALL cpu_time(start_time)
      DO j = 1, npts
        CALL LSQRdecomp(xy, nob, nbf, j)
      ENDDO
      CALL cpu_time(end_time)
      PRINT*, "Time taken to do a QR deomposition (single prec):" ,  &
                end_time - start_time, " seconds"

      CALL cpu_time(start_time)
      DO j = 1,  n_itr
#ifdef NEW_S1
        b(1:nob) = f2(1:nob) - f2(nob+1)
#else
        b = f2
#endif
        CALL solveiThLS(b, 1)
      ENDDO
      CALL cpu_time(end_time)

      PRINT*, 'Solution from solveiThLS() - through QR decomposition:'
      PRINT*, b

#ifdef NEW_S1
      f(1:nob) = a(1:nob)
#else
      f = b
#endif


#ifdef NEW_S1
      CALL check_solution(st2d, f,f2-f2(nob+1))
#else
      CALL check_solution(st2d, f,f2)
#endif

      CALL cpu_time(end_time)

      PRINT*, 'Time taken to solve',n_itr,'systems (single prec):' , &
                 end_time - start_time, " seconds"
      
      CALL closeLSmod()

      END PROGRAM tstLS

