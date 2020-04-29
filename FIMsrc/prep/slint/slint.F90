!------------------------------------------------------------
! This file contains the routines needed to perform linear
! interpolation on sphere.
!
!
! Ning Wang, Jan 2007, init version
! This file contains the routines needed to perform linear
! interpolation on sphere.
!
!
! Ning Wang, Jan 2007, initial version
! Ning Wang, Jan 2011, Added some subroutines and comments.
!
! General purpose subroutines:
!  (1) slint_init(grid1, n1, grid2, n2)
!      grid1, grid2: array of lat-lons that specifies source
!                    and target grid;
!       n1, n2, gripoint numners of source and target grids.
!  (2) slint_init_fn(grid_file1, n1, grid_file2, n2, nn)
!    grid_file1: file name for the source grid specification;
!    grid_file2: file name for the target grid specification;
!    n1, n2: grid point numbers of source and target grids; 
!  slint_init (...) initialize the associated data structures
!  and computes the 
!        
!  (3) bilinear_interp (src_data, tgt_data)
!    src_data: an array of n1 elements that contains the data
!              at source grid points.
!    tgt_data: an array of n2 elements that contains the data
!              at target grid points.
!  bilinear_interp (...)  interpolates the src_data bilinearly 
!  to tgt_dat. 
!
!  (4) nn_interp (src_data, tgt_data)
!  Same as bilinear_interp, except it assigns the nearest 
!  neighbor's value in the src_data to the tgt_data.  
!
! Special and legacy subroutines:
!    Following the similar naming convention as those used 
!    in general subroutines. Details see in-line comments.    
!
! Ning Wang, July 2011, important revision to the package.
!
!-------------------------------------------------------------
MODULE slint

    TYPE GRID
      INTEGER :: type
      INTEGER :: ngp, mx, my
      REAL, ALLOCATABLE :: latlon(:,:)
      REAL, ALLOCATABLE :: coeffs(:,:)
      INTEGER, ALLOCATABLE :: mask(:)
      REAL, ALLOCATABLE :: data(:)
      INTEGER, ALLOCATABLE :: nn(:,:)
       
    END TYPE GRID

    TYPE(GRID) src_grid, tgt_grid
    LOGICAL use_mask
!    integer, external :: gptlinitialize, gptlstart, gptlstop, gptlpr

CONTAINS

! General purpose subroutines, 
SUBROUTINE slint_init(grid1, n1, grid2, n2)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n1, n2
    REAL,INTENT(IN) :: grid1(n1, 2), grid2(n2, 2)

    CALL init_intern_array(grid1, n1, grid2, n2)

END SUBROUTINE slint_init

SUBROUTINE slint_init_ocn(grid1, n1, grid2, n2, icos_ocn_mask, nns_r)
    USE kd, ONLY: set_nns_radius 
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: n1, n2
    REAL,INTENT(IN) :: grid1(n1, 2), grid2(n2, 2)
    INTEGER, INTENT(IN):: icos_ocn_mask(*)
    REAL,INTENT(IN) :: nns_r 
   
    CALL set_nns_radius(nns_r)
    CALL init_intern_array_w_mask(grid1, n1, grid2, n2, icos_ocn_mask)

END SUBROUTINE slint_init_ocn

SUBROUTINE slint_init_fn(grid_file1, n1, grid_file2, n2)
    IMPLICIT NONE
    INTEGER :: n1, n2
    CHARACTER *(*),intent(in) :: grid_file1, grid_file2

    CALL init_intern_fn(grid_file1, n1, grid_file2, n2, 0)

END SUBROUTINE slint_init_fn

SUBROUTINE slint_init_save(fn, grid1, n1, grid2, n2)
   
    IMPLICIT NONE

    CHARACTER *(*), INTENT(IN) :: fn
    INTEGER :: n1, n2
    REAL,INTENT(IN) :: grid1(n1, 2), grid2(n2, 2)
    integer :: ret


!    ret = gptlinitialize ()
!    ret = gptlstart ('total')
!    ret = gptlstart ('init_intern_array')
    CALL init_intern_array(grid1, n1, grid2, n2)
!    ret = gptlstop ('init_intern_array')
    
!    ret = gptlstart ('write_array')
    OPEN (10,file=fn,form='unformatted')
    WRITE(10) n1, n2
    WRITE(10) tgt_grid%nn, tgt_grid%coeffs
    WRITE(10) src_grid%latlon, tgt_grid%latlon
!    ret = gptlstop ('write_array')
   
    CLOSE(10)
!    ret = gptlstop ('total')
!    ret = gptlpr (0)

END SUBROUTINE slint_init_save

SUBROUTINE slint_init_read(fn)

    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: fn

    INTEGER :: n1, n2

    
    OPEN (10,file=fn,status='old', form='unformatted')
    READ(10) n1, n2

    src_grid%type = 1
    src_grid%ngp = n1
    tgt_grid%type = 1
    tgt_grid%ngp = n2

    IF (ALLOCATED(tgt_grid%nn)) THEN
      DEALLOCATE(tgt_grid%nn)
    ENDIF
    IF (ALLOCATED(tgt_grid%coeffs)) THEN
      DEALLOCATE(tgt_grid%coeffs)
    ENDIF
    IF (ALLOCATED(src_grid%latlon)) THEN
      DEALLOCATE(src_grid%latlon)
    ENDIF
    IF (ALLOCATED(tgt_grid%latlon)) THEN
      DEALLOCATE(tgt_grid%latlon)
    ENDIF
    ALLOCATE(tgt_grid%nn(3, n2))
    ALLOCATE(tgt_grid%coeffs(3, n2))
    ALLOCATE(src_grid%latlon(2, n1))
    ALLOCATE(tgt_grid%latlon(2, n2))

    READ(10) tgt_grid%nn, tgt_grid%coeffs
    READ(10) src_grid%latlon, tgt_grid%latlon

    CLOSE(10)

END SUBROUTINE slint_init_read

SUBROUTINE bilinear_interp (src_data, tgt_data)
    IMPLICIT NONE

    REAL, INTENT(IN):: src_data(*)
    REAL, INTENT(OUT):: tgt_data(*)

    REAL v(3)
    REAL*8 c(3)
    INTEGER i

    DO i = 1, tgt_grid%ngp
      c = tgt_grid%coeffs(1:3,i)
      v(1) = src_data(tgt_grid%nn(1, i))
      v(2) = src_data(tgt_grid%nn(2, i))
      v(3) = src_data(tgt_grid%nn(3, i))
      tgt_data(i) = c(1) * v(1) + c(2) * v(2) + c(3) * v(3)
    END DO

END SUBROUTINE bilinear_interp

SUBROUTINE bilinear_interp_ocn (src_data, tgt_data, icos_ocn_mask)
    IMPLICIT NONE

    REAL, INTENT(IN):: src_data(*)
    REAL, INTENT(OUT):: tgt_data(*)
    INTEGER, INTENT(IN):: icos_ocn_mask(*)

    REAL v(3)
    REAL*8 c(3)
    INTEGER i

    DO i = 1, tgt_grid%ngp
      IF (icos_ocn_mask(i)==1) THEN
          c = tgt_grid%coeffs(1:3,i)
          IF (c(1) /= -1.0) THEN 
            v(1) = src_data(tgt_grid%nn(1, i))
            v(2) = src_data(tgt_grid%nn(2, i))
            v(3) = src_data(tgt_grid%nn(3, i))
            tgt_data(i) = c(1) * v(1) + c(2) * v(2) + c(3) * v(3)
          ELSE
            PRINT*, 'Cannot find any source grid point !!'
            PRINT*, 'Target grid point index:', i
            PRINT*, 'Target grid point location (lat,lon):', &
                     tgt_grid%latlon(1, i), tgt_grid%latlon(2, i) 
            STOP
          ENDIF
      ELSE
        tgt_data(i) = -1.e33
      ENDIF 
    END DO

    return
END SUBROUTINE bilinear_interp_ocn

SUBROUTINE nn_interp (src_data, tgt_data)
    IMPLICIT NONE

    REAL, INTENT(IN):: src_data(*)
    REAL, INTENT(OUT):: tgt_data(*)

    REAL v(3)
    REAL*8 c(3)
    INTEGER i

    DO i = 1, tgt_grid%ngp
      v(1) = src_data(tgt_grid%nn(1, i))
      tgt_data(i) = v(1)
    END DO

END SUBROUTINE nn_interp

SUBROUTINE nn_interp_ocn (src_data, tgt_data, icos_ocn_mask)
    IMPLICIT NONE

    REAL, INTENT(IN):: src_data(*)
    REAL, INTENT(OUT):: tgt_data(*)
    INTEGER, INTENT(IN):: icos_ocn_mask(*)

    REAL v(3)
    REAL*8 c(3)
    INTEGER i

    DO i = 1, tgt_grid%ngp
      IF (icos_ocn_mask(i)==1) THEN
        v(1) = src_data(tgt_grid%nn(1, i))
        tgt_data(i) = v(1)
      ELSE
        tgt_data(i) = -1.e33
      ENDIF 
    END DO

END SUBROUTINE nn_interp_ocn

! Generic init subroutine for bilinear interpolation
! To be deprecated
SUBROUTINE bilinear_init_fn(grid_file1, n1, grid_file2, n2)
    IMPLICIT NONE
                                                                                                                                            
    INTEGER, INTENT(IN) :: n1, n2
    CHARACTER *(*), INTENT(IN) :: grid_file1, grid_file2
    
    CALL init_intern_fn(grid_file1, n1, grid_file2, n2, 0)

END SUBROUTINE bilinear_init_fn

! Generic init subroutine for nearest neighbor interpolation
! To be deprecated
SUBROUTINE nn_init_fn(grid_file1, n1, grid_file2, n2)
    IMPLICIT NONE
                                                                                                                                            
    INTEGER :: n1, n2
    CHARACTER *(*) :: grid_file1, grid_file2
    
    CALL init_intern_fn(grid_file1, n1, grid_file2, n2, 1)

END SUBROUTINE nn_init_fn

! Special init subroutine for bilinear interpolation
SUBROUTINE bilinear_init(grid_file1, n1, unit2, n2)
    IMPLICIT NONE
    INTEGER :: n1, n2
    CHARACTER *(*),intent(in) :: grid_file1
    integer       ,intent(in) :: unit2
    
    CALL init_intern_fn_unit(grid_file1, n1, unit2, n2, 0)

END SUBROUTINE bilinear_init

! Special init subroutine for nearest neighbor interpolation
SUBROUTINE nn_init(grid_file1, n1, unit2, n2)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n1, n2
    CHARACTER *(*),intent(in) :: grid_file1
    integer       ,intent(in) :: unit2
    
    CALL init_intern_fn_unit(grid_file1, n1, unit2, n2, 1)

END SUBROUTINE nn_init

! Internal subroutines called by those within the module
SUBROUTINE init_intern_fn(grid_file1, n1, grid_file2, n2, nn)
    USE kd, ONLY: init_kd_tree, close_kd_tree
    IMPLICIT NONE
   
    INTEGER :: n1, n2, nn 
    CHARACTER *(*) :: grid_file1, grid_file2

    INTEGER i, j, g_idx, seq
    REAL, ALLOCATABLE :: llpoints(:,:) 
    integer :: ret

    ALLOCATE(llpoints(n1, 2))
!    ret = gptlstart ('read_gridfile')
    OPEN (10,file=grid_file1,status='old',form='unformatted')
!    READ (10)     ! comment out the two read statements for
!    READ (10)     ! a non-icos grid file
    READ (10) llpoints(:, 1), llpoints(:, 2)
    CLOSE(10)
!    ret = gptlstop ('read_gridfile')

!    ret = gptlstart ('init_kd_tree')
    CALL init_kd_tree(llpoints, n1, 1)
!    ret = gptlstop ('init_kd_tree')

    src_grid%type = 1
    src_grid%ngp = n1
    IF (ALLOCATED(src_grid%latlon)) THEN
      DEALLOCATE(src_grid%latlon)
    ENDIF
    IF (ALLOCATED(src_grid%data)) THEN
      DEALLOCATE(src_grid%data)
    ENDIF
    ALLOCATE(src_grid%latlon(2, n1))
    ALLOCATE(src_grid%data(n1))

    DO i = 1, n1
      src_grid%latlon(1,i) = llpoints(i, 1) 
      src_grid%latlon(2,i) = llpoints(i, 2) 
    END DO

    DEALLOCATE(llpoints)

    tgt_grid%type = 1
    tgt_grid%ngp = n2

    IF (ALLOCATED(tgt_grid%latlon)) THEN
      DEALLOCATE(tgt_grid%latlon)
    ENDIF
    IF (ALLOCATED(tgt_grid%data)) THEN
      DEALLOCATE(tgt_grid%data)
    ENDIF
    IF (ALLOCATED(tgt_grid%coeffs)) THEN 
      DEALLOCATE(tgt_grid%coeffs)
    ENDIF
    IF (ALLOCATED(tgt_grid%nn)) THEN
      DEALLOCATE(tgt_grid%nn)
    ENDIF

    ALLOCATE(tgt_grid%latlon(2, n2))
    ALLOCATE(tgt_grid%nn(3, n2))
    ALLOCATE(tgt_grid%coeffs(3, n2))
    ALLOCATE(tgt_grid%data(n2))
    ALLOCATE(llpoints(n2, 2))
      
    OPEN(10,file=grid_file2,status='old',form='unformatted')
    READ (10)    ! comment out the two READ statements 
    READ (10)    ! for non-icos grid fiiles
    READ(10) llpoints(:, 1), llpoints(:, 2)
    CLOSE(10)

    DO i = 1, n2
      tgt_grid%latlon(1,i) = llpoints(i, 1) 
      tgt_grid%latlon(2,i) = llpoints(i, 2) 
    END DO

    CALL coeff_comp(nn)
    CALL close_kd_tree()

END SUBROUTINE init_intern_fn

SUBROUTINE init_intern_fn_unit(grid_file1, n1, unit2, n2, nn)
    USE kd, ONLY: init_kd_tree, close_kd_tree
    IMPLICIT NONE
   
    INTEGER, INTENT(IN) :: n1, n2, nn
    CHARACTER *(*),INTENT(IN) :: grid_file1
    INTEGER       ,INTENT(IN) :: unit2

    REAL*8 lat, lon, d2r, r2d
    INTEGER i, j, g_idx, seq
    REAL, ALLOCATABLE :: llpoints(:,:) 
    INTEGER :: ierr

    d2r = 4.0*ATAN(1.0)/180.0
    r2d = 1 / d2r

    ALLOCATE(llpoints(n1, 2))
    OPEN (30,file=grid_file1,status='old',form='unformatted')
    READ (30,iostat=ierr) llpoints(:, 1), llpoints(:, 2)
    IF (ierr /= 0) then
      WRITE(6,*)'slint.F90:init_intern_fn_unit: Error reading from unit 30'
      WRITE(6,*)'init_interp_fn_unit: read(30) returns ierr=',ierr
    END IF
    CLOSE(30)
    CALL init_kd_tree(llpoints, n1, 1)

    src_grid%type = 1
    src_grid%ngp = n1
    IF (ALLOCATED(src_grid%latlon)) THEN
      DEALLOCATE(src_grid%latlon)
    END IF
    IF (ALLOCATED(src_grid%data)) THEN
      DEALLOCATE(src_grid%data)
    END IF
    ALLOCATE(src_grid%latlon(2, n1))
    ALLOCATE(src_grid%data(n1))

    DO i = 1, n1
      src_grid%latlon(1,i) = llpoints(i, 1) 
      src_grid%latlon(2,i) = llpoints(i, 2) 
    END DO

    DEALLOCATE(llpoints)

    tgt_grid%type = 1
    tgt_grid%ngp = n2

    IF (ALLOCATED(tgt_grid%latlon)) THEN
      DEALLOCATE(tgt_grid%latlon)
    ENDIF
    IF (ALLOCATED(tgt_grid%data)) THEN
      DEALLOCATE(tgt_grid%data)
    ENDIF
    IF (ALLOCATED(tgt_grid%coeffs)) THEN 
      DEALLOCATE(tgt_grid%coeffs)
    ENDIF
    IF (ALLOCATED(tgt_grid%nn)) THEN
      DEALLOCATE(tgt_grid%nn)
    ENDIF

    ALLOCATE(tgt_grid%latlon(2, n2))
    ALLOCATE(tgt_grid%nn(3, n2))
    ALLOCATE(tgt_grid%coeffs(3, n2))
    ALLOCATE(tgt_grid%data(n2))
    ALLOCATE(llpoints(n2, 2))
      
    READ (unit2) llpoints(:, 1)
    READ (unit2) llpoints(:, 2)
    DO i = 1, n2
      tgt_grid%latlon(1,i) = llpoints(i, 1) 
      tgt_grid%latlon(2,i) = llpoints(i, 2) 
    END DO

    DEALLOCATE(llpoints)

    CALL coeff_comp(nn)
    CALL close_kd_tree()

END SUBROUTINE init_intern_fn_unit

! internal routine that does the interpolation
SUBROUTINE interp_intern() 
    IMPLICIT NONE

    INTEGER i, j, mx, my, g_idx
    REAL*4 v(3), c(3)
   
    IF (src_grid%type == 1 .AND. tgt_grid%type == 0) THEN
      mx = tgt_grid%mx
      my = tgt_grid%my
        DO j = 1, my
      DO i = 1, mx
          g_idx = (i + (j - 1) * mx)
          c = tgt_grid%coeffs(1:3,g_idx)
          v(1) = src_grid%data(tgt_grid%nn(1, g_idx))
          v(2) = src_grid%data(tgt_grid%nn(2, g_idx))
          v(3) = src_grid%data(tgt_grid%nn(3, g_idx))
          tgt_grid%data(g_idx) = c(1) * v(1) + c(2) * v(2) + c(3) * v(3) 
        END DO
      END DO
    ELSE IF (src_grid%type == 1 .AND. tgt_grid%type == 1) THEN
      DO i = 1, tgt_grid%ngp
        c = tgt_grid%coeffs(1:3,i)
        v(1) = src_grid%data(tgt_grid%nn(1, i))
        v(2) = src_grid%data(tgt_grid%nn(2, i))
        v(3) = src_grid%data(tgt_grid%nn(3, i))
        tgt_grid%data(i) = c(1) * v(1) + c(2) * v(2) + c(3) * v(3) 
      END DO
    END IF

END SUBROUTINE interp_intern

! Subroutine to compute interpolation coefficients
SUBROUTINE coeff_comp(nn_w)
    IMPLICIT NONE

    INTEGER nn_w

    REAL latlon(2, 3), intsec(2), gcd1, gcd2
    REAL hp1(3), hp2(3), min_dist
    INTEGER i, j,mx, my, g_idx, nn(3), num
    REAL epsilon, r2d , t1, t2
    integer :: ret

    epsilon = 0.00000000001
    r2d = 180.0/ACOS(-1.0)

    IF (src_grid%type == 1 .AND. tgt_grid%type == 0) THEN
      mx = tgt_grid%mx
      my = tgt_grid%my
      DO i = 1, mx
        DO j = 1, my
          g_idx = (i + (j - 1) * mx)
          CALL nsn(tgt_grid%latlon(1,g_idx), nn, num, min_dist)
          tgt_grid%nn(1:3, g_idx) = nn(1:3)
          IF (num == 1 .OR. nn_w == 1) THEN
            tgt_grid%coeffs(1, g_idx) = 1.0  
            tgt_grid%coeffs(2, g_idx) = 0.0  
            tgt_grid%coeffs(3, g_idx) = 0.0  
          ELSE IF (num == 0) THEN
            tgt_grid%coeffs(1, i) = -1.0
            tgt_grid%coeffs(2, i) = 0.0
            tgt_grid%coeffs(3, i) = 0.0
          ELSE IF (num == 2) THEN
            latlon(1:2,1) = src_grid%latlon(1:2,nn(1))
            latlon(1:2,2) = src_grid%latlon(1:2,nn(2))
            CALL gcd_ratio(latlon(1,1), latlon(1,2), tgt_grid%latlon(1,g_idx), gcd1, t1)
            tgt_grid%coeffs(1, g_idx) = (1.0 - t1)   
            tgt_grid%coeffs(2, g_idx) = t1 
            tgt_grid%coeffs(3, g_idx) = 0.0  
          ELSE
            latlon(1:2,1) = src_grid%latlon(1:2,nn(1))
            latlon(1:2,2) = src_grid%latlon(1:2,nn(2))
            latlon(1:2,3) = src_grid%latlon(1:2,nn(3))
            CALL intersection (latlon(1,1), tgt_grid%latlon(1,g_idx), latlon(1,2), latlon(1,3), intsec) 
            CALL gcd_ratio(latlon(1,2), latlon(1,3), intsec, gcd1, t1)
            CALL gcd_ratio(latlon(1,1), intsec, tgt_grid%latlon(1,g_idx), gcd2, t2)

            IF (t1 /= t1 .OR. t2 /= t2) THEN
              PRINT*, 't1 or t2 NaN:', t1, t2
            ENDIF
            
            tgt_grid%coeffs(1, g_idx) = (1.0 - t2)   
            tgt_grid%coeffs(2, g_idx) = t2 * (1.0 - t1)    
            tgt_grid%coeffs(3, g_idx) = t2 * t1 
          END IF
        END DO
      END DO
    ELSE IF (src_grid%type == 1 .AND. tgt_grid%type == 1) THEN
!      ret = gptlstart ('compute_one_intp_coeff')
!$OMP PARALLEL DO PRIVATE (nn, num, min_dist, latlon, gcd1, t1, intsec, gcd2, t2, ret)
      DO i = 1, tgt_grid%ngp
        IF (use_mask) THEN
          IF (tgt_grid%mask(i) /= 1) CYCLE 
        ENDIF
        CALL nsn(tgt_grid%latlon(1,i), nn, num, min_dist)
        tgt_grid%nn(1:3, i) = nn(1:3)
        IF (num == 1 .OR. nn_w == 1) THEN 
          tgt_grid%coeffs(1, i) = 1.0  
          tgt_grid%coeffs(2, i) = 0.0  
          tgt_grid%coeffs(3, i) = 0.0  
        ELSE IF (num == 0) THEN
          tgt_grid%coeffs(1, i) = -1.0
          tgt_grid%coeffs(2, i) = 0.0
          tgt_grid%coeffs(3, i) = 0.0
        ELSE IF (num == 2) THEN
          latlon(1:2,1) = src_grid%latlon(1:2,nn(1))
          latlon(1:2,2) = src_grid%latlon(1:2,nn(2))
          CALL gcd_ratio(latlon(1,1), latlon(1,2), tgt_grid%latlon(1,i), gcd1, t1)
          tgt_grid%coeffs(1, i) = (1.0 - t1)   
          tgt_grid%coeffs(2, i) = t1 
          tgt_grid%coeffs(3, i) = 0.0  
        ELSE
          latlon(1:2,1) = src_grid%latlon(1:2,nn(1))
          latlon(1:2,2) = src_grid%latlon(1:2,nn(2))
          latlon(1:2,3) = src_grid%latlon(1:2,nn(3))
          CALL intersection (latlon(1,1), tgt_grid%latlon(1,i), latlon(1,2), latlon(1,3), intsec) 
          CALL gcd_ratio(latlon(1,2), latlon(1,3), intsec, gcd1, t1)
          CALL gcd_ratio(latlon(1,1), intsec, tgt_grid%latlon(1,i), gcd2, t2)

          IF (t1 /= t1 .OR. t2 /= t2) THEN
            PRINT*, 't1 or t2 NaN:', t1, t2
          ENDIF

          tgt_grid%coeffs(1, i) = (1.0 - t2)   
          tgt_grid%coeffs(2, i) = t2 * (1.0 - t1)    
          tgt_grid%coeffs(3, i) = t2 * t1 
        END IF
      END DO
!      ret = gptlstop ('compute_one_intp_coeff')
    END IF

END SUBROUTINE coeff_comp

! Subroutine to compute interpolation coefficients, distance weight
SUBROUTINE coeff_comp1(nn_w)
    USE kd, ONLY:knn_search
    IMPLICIT NONE

    INTEGER nn_w

    REAL latlon(2, 3), intsec(2), gcd1, gcd2, part_gcd1, part_gcd2
    REAL hp1(3), hp2(3), min_dist
    INTEGER i, j, mx, my, g_idx, nn(3), num
    REAL epsilon, r2d , d1, d2, d3, rd1, rd2, rd3, srds 

    epsilon = 0.00000000001
    r2d = 180.0/(ATAN(1.0) * 4.0)

    IF (src_grid%type == 1 .AND. tgt_grid%type == 0) THEN
      mx = tgt_grid%mx
      my = tgt_grid%my
      DO i = 1, mx
        DO j = 1, my
          g_idx = (i + (j - 1) * mx)
          CALL nsn(tgt_grid%latlon(1,g_idx), nn, num, min_dist)
          tgt_grid%nn(1:3, g_idx) = nn(1:3)
          IF (min_dist < epsilon .OR. nn_w == 1) THEN 
            tgt_grid%coeffs(1, g_idx) = 1.0  
            tgt_grid%coeffs(2, g_idx) = 0.0  
            tgt_grid%coeffs(3, g_idx) = 0.0  
          ELSE 
            latlon(1:2,1) = src_grid%latlon(1:2,nn(1))
            latlon(1:2,2) = src_grid%latlon(1:2,nn(2))
            latlon(1:2,3) = src_grid%latlon(1:2,nn(3))
            d1 =  gc_dist(latlon(1,1), tgt_grid%latlon(1,g_idx))
            d2 =  gc_dist(latlon(1,2), tgt_grid%latlon(1,g_idx))
            d3 =  gc_dist(latlon(1,3), tgt_grid%latlon(1,g_idx))

            rd1 = 1.0/d1
            rd2 = 1.0/d2
            rd3 = 1.0/d3
            srds = rd1 + rd2 + rd3

            tgt_grid%coeffs(1, g_idx) = rd1 / srds    
            tgt_grid%coeffs(2, g_idx) = rd2 / srds    
            tgt_grid%coeffs(3, g_idx) = rd3 / srds
          END IF
        END DO
      END DO
    ELSE IF (src_grid%type == 1 .AND. tgt_grid%type == 1) THEN
      DO i = 1, tgt_grid%ngp
        CALL nsn(tgt_grid%latlon(1,i), nn, num, min_dist)
        tgt_grid%nn(1:3, i) = nn(1:3)
        IF (min_dist < epsilon .OR. nn_w == 1) THEN 
          tgt_grid%coeffs(1, i) = 1.0  
          tgt_grid%coeffs(2, i) = 0.0  
          tgt_grid%coeffs(3, i) = 0.0  
        ELSE
          latlon(1:2,1) = src_grid%latlon(1:2,nn(1))
          latlon(1:2,2) = src_grid%latlon(1:2,nn(2))
          latlon(1:2,3) = src_grid%latlon(1:2,nn(3))
          d1 =  gc_dist(latlon(1,1), tgt_grid%latlon(1,i))
          d2 =  gc_dist(latlon(1,2), tgt_grid%latlon(1,i))
          d3 =  gc_dist(latlon(1,3), tgt_grid%latlon(1,i))

          rd1 = 1.0/d1
          rd2 = 1.0/d2
          rd3 = 1.0/d3
          srds = rd1 + rd2 + rd3

          tgt_grid%coeffs(1, i) = rd1 / srds    
          tgt_grid%coeffs(2, i) = rd2 / srds    
          tgt_grid%coeffs(3, i) = rd3 / srds
        END IF
      END DO
    END IF

END SUBROUTINE coeff_comp1 

SUBROUTINE interp (src_data, tgt_data)
    IMPLICIT NONE

    REAL src_data(*)
    REAL tgt_data(*)

    INTEGER i, n

    n = src_grid%ngp
    DO i = 1, n
      src_grid%data(i) = src_data(i)
    END DO

    CALL interp_intern()

    n = tgt_grid%ngp
    DO i = 1, n
      tgt_data(i) = tgt_grid%data(i)
    END DO
END SUBROUTINE interp

! Two legacy subroutines, keep them for backward compatiablity.
SUBROUTINE nn_int (src_data, tgt_data)
    IMPLICIT NONE

    REAL src_data(*)
    REAL tgt_data(*)
    CALL interp(src_data, tgt_data)

END SUBROUTINE nn_int

SUBROUTINE bl_int (src_data, tgt_data)
    IMPLICIT NONE

    REAL src_data(*)
    REAL tgt_data(*)
    CALL interp(src_data, tgt_data)

END SUBROUTINE bl_int

SUBROUTINE bilinear_interp_uv(src_u, tgt_u, src_v, tgt_v)
    IMPLICIT NONE
! Arguments    
    REAL, intent(in) :: src_u(*), src_v(*)
    REAL, intent(out) :: tgt_u(*), tgt_v(*)

! Local workspace
    REAL u, v, c(3)
    REAL lat, lon, uxyz_src(3, 3), uxyz_tgt(3) 
    REAL mag_src(3), mag_tgt, mag_tmp
    INTEGER i, ipn, nip
    REAL, PARAMETER :: eps = 0.00001

    nip = tgt_grid%ngp

    DO ipn = 1, nip
      DO i = 1, 3
        lat = src_grid%latlon(1, tgt_grid%nn(i, ipn))
        lon = src_grid%latlon(2, tgt_grid%nn(i, ipn))
        u = src_u(tgt_grid%nn(i, ipn))
        v = src_v(tgt_grid%nn(i, ipn))
        CALL uv2xyz(u, v, lat, lon, uxyz_src(1, i), uxyz_src(2, i), uxyz_src(3, i))
        mag_src(i) = sqrt(u*u+v*v)
      ENDDO
      c = tgt_grid%coeffs(1:3,ipn)
      uxyz_tgt(1:3) = c(1) * uxyz_src(1:3, 1) + c(2) * uxyz_src(1:3, 2) &
                    + c(3) * uxyz_src(1:3, 3)
      mag_tgt = c(1) * mag_src(1) + c(2) * mag_src(2) + c(3) * mag_src(3)
      lat = tgt_grid%latlon(1, ipn)
      lon = tgt_grid%latlon(2, ipn)
      CALL xyz2uv(u, v, lat, lon, uxyz_tgt(1), uxyz_tgt(2), uxyz_tgt(3))
      mag_tmp = sqrt(u*u+v*v)
      IF (mag_tmp < eps) THEN
        mag_tgt = 1.0
      ELSE
        mag_tgt = mag_tgt / mag_tmp 
      ENDIF
      tgt_u(ipn) = u * mag_tgt
      tgt_v(ipn) = v * mag_tgt
    END DO
END SUBROUTINE bilinear_interp_uv

SUBROUTINE bilinear_interp_uv2(src_u, tgt_u, src_v, tgt_v)
    IMPLICIT NONE
! Arguments
    REAL, intent(in) :: src_u(*), src_v(*)
    REAL, intent(out) :: tgt_u(*), tgt_v(*)

! Local workspace
    REAL uv_s(2), uv_t(2,3), c(3)
    REAL ll_s(2), ll_t(2)
    INTEGER i, j
    REAL, PARAMETER :: eps = 0.00001
    REAL :: u, v, mag_src(3), mag_tgt, mag_tmp

    DO i = 1, tgt_grid%ngp
      ll_t(1) = tgt_grid%latlon(1, i)
      ll_t(2) = tgt_grid%latlon(2, i)
      DO j = 1, 3
        ll_s(1) = src_grid%latlon(1, tgt_grid%nn(j, i))
        ll_s(2) = src_grid%latlon(2, tgt_grid%nn(j, i))
        uv_s(1) = src_u(tgt_grid%nn(j, i))
        uv_s(2) = src_v(tgt_grid%nn(j, i))
        CALL uvs2uvt(uv_s, ll_s, uv_t(1:2,j), ll_t)
        mag_src(j) = sqrt(uv_s(1)*uv_s(1)+uv_s(2)*uv_s(2))
      ENDDO
      c = tgt_grid%coeffs(1:3,i)
      u = c(1)*uv_t(1,1)+c(2)*uv_t(1,2)+c(3)*uv_t(1,3)
      v = c(1)*uv_t(2,1)+c(2)*uv_t(2,2)+c(3)*uv_t(2,3)
      mag_tmp = sqrt(u*u+v*v)
      mag_tgt = c(1)*mag_src(1)+c(2)*mag_src(2)+c(3)*mag_src(3)
      IF (mag_tmp < eps) THEN
        mag_tgt = 1.0
      ELSE
        mag_tgt = mag_tgt / mag_tmp
      ENDIF
      tgt_u(i) = u * mag_tgt
      tgt_v(i) = v * mag_tgt
    END DO
END SUBROUTINE bilinear_interp_uv2

SUBROUTINE bilinear_init_i2r(mx, my, llpoints, nip)
    USE kd, ONLY:init_kd_tree, close_kd_tree
    IMPLICIT NONE
   
    INTEGER, intent(in) :: mx, my, nip
    REAL   , intent(in) :: llpoints(nip,2)

    INTEGER i, j, g_idx, seq
    REAL pi

    CALL init_kd_tree(llpoints, nip, 1)

    src_grid%type = 1
    src_grid%ngp = nip
    IF (ALLOCATED(src_grid%latlon)) THEN
      DEALLOCATE(src_grid%latlon)
    END IF
    IF (ALLOCATED(src_grid%data)) THEN
      DEALLOCATE(src_grid%data)
    END IF
    ALLOCATE(src_grid%latlon(2, nip))
    ALLOCATE(src_grid%data(nip))

    DO i = 1, nip
      src_grid%latlon(1,i) = llpoints(i, 1) 
      src_grid%latlon(2,i) = llpoints(i, 2) 
    END DO

    tgt_grid%type = 0
    tgt_grid%mx = mx
    tgt_grid%my = my

    IF (ALLOCATED(tgt_grid%latlon)) THEN
      DEALLOCATE(tgt_grid%latlon)
    ENDIF
    IF (ALLOCATED(tgt_grid%data)) THEN
      DEALLOCATE(tgt_grid%data)
    ENDIF
    IF (ALLOCATED(tgt_grid%coeffs)) THEN 
      DEALLOCATE(tgt_grid%coeffs)
    ENDIF
    IF (ALLOCATED(tgt_grid%nn)) THEN
      DEALLOCATE(tgt_grid%nn)
    ENDIF

    ALLOCATE(tgt_grid%latlon(2, mx * my))
    ALLOCATE(tgt_grid%nn(3, mx * my))
    ALLOCATE(tgt_grid%coeffs(3, mx * my))
    ALLOCATE(tgt_grid%data(mx * my))
    pi = 4.0*ATAN(1.0)
    DO i = 1, mx 
      DO j = 1, my
        g_idx = (i + (j - 1) * mx)
        tgt_grid%latlon(1, g_idx) = (REAL(j - 1) - REAL(my - 1) * 0.5)  * pi / REAL(my - 1) 
        tgt_grid%latlon(1, g_idx) = -tgt_grid%latlon(1, g_idx) 
        tgt_grid%latlon(2, g_idx) = REAL(i - 1) * 2.0 * pi / REAL(mx) 
      END DO
    END DO
      
    CALL coeff_comp(0)
    CALL close_kd_tree()

END SUBROUTINE bilinear_init_i2r

SUBROUTINE bilinear_interp_i2r(k, nlevels, vardata, data_xyz) 
    IMPLICIT NONE

    INTEGER k, nlevels
    REAL vardata(*)
    REAL data_xyz(*)
  
    INTEGER i, j, n, mx, my

    n = src_grid%ngp
    DO i = 1, n
      src_grid%data(i) = vardata(k + (i - 1) * nlevels) 
    END DO

    CALL interp_intern()

    mx = tgt_grid%mx
    my = tgt_grid%my
    DO i = 1, mx
      DO j = 1, my
        data_xyz((k - 1) * mx * my + (j - 1) * mx + i) =  &
        tgt_grid%data((j - 1) * mx + i) 
      END DO
    END DO

END SUBROUTINE bilinear_interp_i2r

! Convenient subroutine to init the interpolation from given arrays
SUBROUTINE init_intern_array(grid1, n1, grid2, n2)
    USE kd, ONLY: init_kd_tree, close_kd_tree
    IMPLICIT NONE
   
    INTEGER :: n1, n2
    REAL grid1(n1, 2), grid2(n2, 2)

    INTEGER i, j, g_idx, seq

    CALL init_kd_tree(grid1, n1, 1)

    src_grid%type = 1
    src_grid%ngp = n1
    IF (ALLOCATED(src_grid%latlon)) THEN
      DEALLOCATE(src_grid%latlon)
    ENDIF
    IF (ALLOCATED(src_grid%data)) THEN
      DEALLOCATE(src_grid%data)
    ENDIF
    ALLOCATE(src_grid%latlon(2, n1))
    ALLOCATE(src_grid%data(n1))

    DO i = 1, n1
      src_grid%latlon(1,i) = grid1(i, 1) 
      src_grid%latlon(2,i) = grid1(i, 2) 
    END DO

    tgt_grid%type = 1
    tgt_grid%ngp = n2

    IF (ALLOCATED(tgt_grid%latlon)) THEN
      DEALLOCATE(tgt_grid%latlon)
    ENDIF
    IF (ALLOCATED(tgt_grid%data)) THEN
      DEALLOCATE(tgt_grid%data)
    ENDIF
    IF (ALLOCATED(tgt_grid%coeffs)) THEN
      DEALLOCATE(tgt_grid%coeffs)
    ENDIF
    IF (ALLOCATED(tgt_grid%nn)) THEN
      DEALLOCATE(tgt_grid%nn)
    ENDIF
    IF (ALLOCATED(tgt_grid%mask)) THEN
      DEALLOCATE(tgt_grid%mask)
    ENDIF
    ALLOCATE(tgt_grid%latlon(2, n2))
    ALLOCATE(tgt_grid%data(n2))
    ALLOCATE(tgt_grid%coeffs(3, n2))
    ALLOCATE(tgt_grid%nn(3, n2))
      
    DO i = 1, n2
      tgt_grid%latlon(1,i) = grid2(i, 1) 
      tgt_grid%latlon(2,i) = grid2(i, 2) 
    END DO

    use_mask = .false.
    CALL coeff_comp(0)
    CALL close_kd_tree()

END SUBROUTINE init_intern_array

SUBROUTINE init_intern_array_w_mask(grid1, n1, grid2, n2, mask)
    USE kd, ONLY: init_kd_tree, close_kd_tree
    IMPLICIT NONE
   
    INTEGER, INTENT(IN) :: n1, n2
    REAL, INTENT(IN) :: grid1(n1, 2), grid2(n2, 2)
    INTEGER, INTENT(IN) :: mask(n2)

    INTEGER i, j, g_idx, seq

    CALL init_kd_tree(grid1, n1, 1)

    src_grid%type = 1
    src_grid%ngp = n1
    IF (ALLOCATED(src_grid%latlon)) THEN
      DEALLOCATE(src_grid%latlon)
    ENDIF
    IF (ALLOCATED(src_grid%data)) THEN
      DEALLOCATE(src_grid%data)
    ENDIF
    ALLOCATE(src_grid%latlon(2, n1))
    ALLOCATE(src_grid%data(n1))

    DO i = 1, n1
      src_grid%latlon(1,i) = grid1(i, 1) 
      src_grid%latlon(2,i) = grid1(i, 2) 
    END DO

    tgt_grid%type = 1
    tgt_grid%ngp = n2

    IF (ALLOCATED(tgt_grid%latlon)) THEN
      DEALLOCATE(tgt_grid%latlon)
    ENDIF
    IF (ALLOCATED(tgt_grid%data)) THEN
      DEALLOCATE(tgt_grid%data)
    ENDIF
    IF (ALLOCATED(tgt_grid%coeffs)) THEN
      DEALLOCATE(tgt_grid%coeffs)
    ENDIF
    IF (ALLOCATED(tgt_grid%nn)) THEN
      DEALLOCATE(tgt_grid%nn)
    ENDIF
    IF (ALLOCATED(tgt_grid%mask)) THEN
      DEALLOCATE(tgt_grid%mask)
    ENDIF
    ALLOCATE(tgt_grid%latlon(2, n2))
    ALLOCATE(tgt_grid%data(n2))
    ALLOCATE(tgt_grid%coeffs(3, n2))
    ALLOCATE(tgt_grid%nn(3, n2))
    ALLOCATE(tgt_grid%mask(n2))
      
    DO i = 1, n2
      tgt_grid%latlon(1,i) = grid2(i, 1) 
      tgt_grid%latlon(2,i) = grid2(i, 2) 
      tgt_grid%mask(i) = mask(i)
    END DO

    use_mask = .true.
    CALL coeff_comp(0)
    CALL close_kd_tree()

END SUBROUTINE init_intern_array_w_mask

! subrountines for spherical curve interpolation 
SUBROUTINE gcd_ratio (p1, p2, p, gcd, p_gcd)
    IMPLICIT NONE

    REAL p1(2), p2(2), p(2), gcd, p_gcd 
    REAL gcdp1p2, gcdp1p, gcdp2p
    REAL r2d, eps

    r2d = 180.0 / (atan(1.0)*4.0)
    eps = 1.0E-6

    gcdp1p2 = gc_dist(p1, p2)
    gcdp1p = gc_dist(p1, p) 
    gcdp2p = gc_dist(p2, p) 

    IF (gcdp1p2 <= eps) THEN
      PRINT*, 'nearest neighbor almost overlap!!'
      PRINT*, p1*r2d, p2*r2d, p*r2d
      p_gcd = 0.0
    ELSE IF (gcdp1p <= gcdp1p2 .AND. gcdp2p <= gcdp1p2) THEN ! p inside the p1p2 segment
      p_gcd = gcdp1p / gcdp1p2
    ELSE IF (gcdp1p > gcdp1p2) THEN ! outside of end point p2
      p_gcd = 1.0 ! don't allow ! extrapolation
!      IF (gcdp2p > eps) THEN 
!        PRINT*, 'extrapolation! outside of p2' 
!        PRINT*, p1*r2d, p2*r2d, p*r2d
!        PRINT*, gcdp1p, gcdp1p2, gc_dist(p,p2)
!      ENDIF
    ELSE IF (gcdp2p > gcdp1p2) THEN ! outside of end point p1
      p_gcd = 0.0 ! don't allow ! extrapolation
!      IF (gcdp1p > eps) THEN
!        PRINT*, 'extrapolation! outside of p1'
!        PRINT*, p1*r2d, p2*r2d, p*r2d
!        PRINT*, gcdp2p, gcdp1p2,  gc_dist(p,p1)
!      ENDIF
    ENDIF 
      gcd = gcdp1p2 

END SUBROUTINE gcd_ratio

! Great circle distance calculation, law of cosine formula 
!FUNCTION gc_dist(p1, p2)
!    IMPLICIT NONE
!
!    REAL gc_dist
!    REAL, INTENT(IN) :: p1(2), p2(2)
!
!    gc_dist = ACOS(COS(p1(1)) * COS(p2(1)) * COS(p1(2) - p2(2)) + SIN(p1(1)) * SIN(p2(1))) 
!
!END FUNCTION gc_dist


! Great circle distance calculation, using Haversine formula. 
! It is more accurate to compute small angular distances.
FUNCTION gc_dist(p1, p2)

    IMPLICIT NONE

    REAL gc_dist
    REAL, INTENT(IN) :: p1(2), p2(2)

    REAL dlatov2, dlonov2, a 

    dlatov2 = (p2(1)-p1(1))/2.0
    dlonov2 = (p2(2)-p1(2))/2.0
    a = sin(dlatov2) * sin(dlatov2) + cos(p1(1))*cos(p2(1))*sin(dlonov2)*sin(dlonov2)
    gc_dist = 2.0 * atan2(sqrt(a), sqrt(1.0-a))

END FUNCTION gc_dist

LOGICAL FUNCTION enclosure(p1, p2, p3, p, co_gc)
    IMPLICIT NONE
    REAL, INTENT(IN) :: p1(2), p2(2), p3(2), p(2)
    INTEGER, INTENT(OUT) :: co_gc

    REAL*8 p1_xy(2), p2_xy(2), p3_xy(2), p_xy(2)
    REAL*8 cp1_z, cp2_z, cp3_z, cos_d2c, eps, eps2

    eps  = 0.00000001
    eps2 = 0.00000000001
    eps2 = 0.0000001
    co_gc = 0
    cos_d2c = sin(p(1))*sin(p1(1)) + cos(p(1))*cos(p1(1))*cos(p1(2)-p(2)) 
    p1_xy(1) = (cos(p1(1))*sin(p1(2) - p(2))) / cos_d2c
    p1_xy(2) = (cos(p(1))*sin(p1(1)) - sin(p(1))*cos(p1(1))*cos(p1(2) - p(2))) / cos_d2c

    cos_d2c = sin(p(1))*sin(p2(1)) + cos(p(1))*cos(p2(1))*cos(p2(2)-p(2)) 
    p2_xy(1) = (cos(p2(1))*sin(p2(2) - p(2))) / cos_d2c
    p2_xy(2) = (cos(p(1))*sin(p2(1)) - sin(p(1))*cos(p2(1))*cos(p2(2) - p(2))) / cos_d2c

    cos_d2c = sin(p(1))*sin(p3(1)) + cos(p(1))*cos(p3(1))*cos(p3(2)-p(2)) 
    p3_xy(1) = (cos(p3(1))*sin(p3(2) - p(2))) / cos_d2c
    p3_xy(2) = (cos(p(1))*sin(p3(1)) - sin(p(1))*cos(p3(1))*cos(p3(2) - p(2))) / cos_d2c

    cp1_z = p1_xy(1)*p2_xy(2) - p1_xy(2)*p2_xy(1)
    cp2_z = p2_xy(1)*p3_xy(2) - p2_xy(2)*p3_xy(1)
    cp3_z = p3_xy(1)*p1_xy(2) - p3_xy(2)*p1_xy(1)
    
    IF (abs(cp1_z) < eps2) co_gc = 1
    IF (abs(cp2_z) < eps2) co_gc = 2
    IF (abs(cp3_z) < eps2) co_gc = 3

    IF (cp1_z*cp2_z .LT. -eps2 .OR. cp1_z*cp3_z .LT. -eps2) THEN
      enclosure = .false.
      RETURN
    ENDIF

    enclosure = .true.
    RETURN

END FUNCTION enclosure

LOGICAL FUNCTION co_gc(p1, p2, p3)
    IMPLICIT NONE
    REAL, INTENT(IN) :: p1(2), p2(2), p3(2)

    co_gc = .true.
END FUNCTION co_gc


SUBROUTINE nsn(q_ll, nn, num, min_dist)
    USE kd, ONLY: num_k, knn_search
    USE kd, ONLY: knn_search_ts
    IMPLICIT NONE 
    REAL, INTENT(IN) :: q_ll(2) 
    INTEGER, INTENT(out) :: nn(3), num
    REAL,INTENT(OUT) :: min_dist
    
    REAL nn_ll(2, 3), hp1(3), hp2(3)
    REAL qxyz(3), nnxyz(3), min_d(num_k)
    REAL eps

    INTEGER nni(3), co_gc, nn_swp, num_nn

    eps = 0.00000000001

    hp1 = 0.0
    hp2 = 0.0
!    CALL knn_search(q_ll, nni, min_d, hp1, hp2, 1.0, num_nn)
    CALL knn_search_ts(q_ll, nni, min_d, hp1, hp2, 1.0, 1, num_nn)
!   print*, 'q_ll, nni', q_ll*180/acos(-1.0), nni(1)
    nn(1) = nni(1)
    nn_ll(1:2, 1) = src_grid%latlon(1:2,nn(1))   ! the first vertex  
    min_dist = min_d(1)

    IF (min_dist < eps) THEN   ! if the nearest neighbor is too close
      nn(2) = nn(1)
      nn(3) = nn(1)
      num  = 1
      RETURN
    ENDIF
    
    CALL ll2xyz(q_ll, qxyz)
    CALL ll2xyz(nn_ll(1:2,1), nnxyz)
    hp1 = qxyz - nnxyz / inner_product(qxyz, nnxyz) 
    CALL cross_product2(qxyz, nnxyz, hp2)
    hp1 = hp1 / sqrt(inner_product(hp1, hp1))
    hp2 = hp2 / sqrt(inner_product(hp2, hp2))

!    CALL knn_search(q_ll, nni, min_d, hp1, hp2, 1.0, num_nn)
    CALL knn_search_ts(q_ll, nni, min_d, hp1, hp2, 1.0, 1, num_nn)
    nn(2) = nni(1)
    nn_ll(1:2, 2) = src_grid%latlon(1:2,nn(2))   ! the second vertex  

    hp2 = -hp2
    CALL ll2xyz(nn_ll(1:2,2), nnxyz)
    CALL cross_product2(qxyz, nnxyz, hp1)
!    CALL knn_search(q_ll, nni, min_d, hp1, hp2, 1.0, num_nn)
    CALL knn_search_ts(q_ll, nni, min_d, hp1, hp2, 1.0, 1, num_nn)
    nn(3) = nni(1)
    nn_ll(1:2, 3) = src_grid%latlon(1:2,nn(3))   ! the third vertex  

    IF (enclosure(nn_ll(1:2, 1), nn_ll(1:2, 2), nn_ll(1:2, 3), q_ll, co_gc)) THEN
      num = 3
    ELSE 
      PRINT*, 'inside test fails'
    END IF

    IF (co_gc /= 0) THEN
      num = 2
      IF (nn(1) == nn(3)) THEN
        RETURN
      ELSE IF (co_gc == 1) THEN
        RETURN
      ELSE IF (co_gc == 2) THEN
        nn_swp = nn(1)
        nn(1) = nn(2)
        nn(2) = nn(3)
        nn(3) = nn_swp
        RETURN
      ELSE IF (co_gc == 3) THEN
        nn_swp = nn(2)
        nn(2) = nn(3)
        nn(3) = nn_swp
      RETURN
      ENDIF
    ENDIF
END SUBROUTINE nsn

SUBROUTINE intersection (p1, p2, p3, p4, p)
    IMPLICIT NONE
    REAL p1(2), p2(2), p3(2), p4(2), p(2)
    REAL gc1(3), gc2(3), e(3)
    REAL pi

    pi = ATAN(1.0) * 4.0

    CALL cross_product1(p1, p2, gc1)
    CALL cross_product1(p4, p3, gc2)
    CALL cross_product2(gc1, gc2, e)

    CALL xyz2ll(e, p)

    IF (gc_dist(p2, p) > pi / 4.0) THEN
      p(2) = p(2) + pi
      p(1) = -p(1)
    END IF
 
    IF (p(2) < 0.0) THEN
      p(2) = p(2) + 2.0 * pi
    END IF

END SUBROUTINE intersection

FUNCTION inner_product(x1, x2)
    IMPLICIT NONE

    REAL x1(3), x2(3), inner_product

    inner_product = x1(1)*x2(1) + x1(2)*x2(2) + x1(3)*x2(3)

END FUNCTION inner_product

SUBROUTINE cross_product1(p1, p2, gc)
    IMPLICIT NONE

    REAL p1(2), p2(2), gc(3)
    REAL a, b, c, d, e, f, g

    a = SIN(p1(1) + p2(1))
    b = SIN(p1(1) - p2(1))
    c = SIN((p1(2) + p2(2))/ 2.0)
    d = SIN((p1(2) - p2(2))/ 2.0)
    e = COS((p1(2) + p2(2))/ 2.0)
    f = COS((p1(2) - p2(2))/ 2.0)
    g = COS(p1(1)) * COS(p2(1)) 

    gc(1) = b * c * f  - a * e * d
    gc(2) = b * e * f  + a * c * d
    gc(3) = 2.0 * g * d * f

END SUBROUTINE cross_product1

SUBROUTINE cross_product2(e1, e2, e)
    IMPLICIT NONE
    REAL e1(3), e2(3), e(3)
 
    e(1) = e1(2) * e2(3) - e2(2) * e1(3)
    e(2) = e1(3) * e2(1) - e2(3) * e1(1)
    e(3) = e1(1) * e2(2) - e2(1) * e1(2)

END SUBROUTINE cross_product2

SUBROUTINE xyz2ll(e, p)
    IMPLICIT NONE
    REAL e(3), p(2)
    REAL, PARAMETER :: eps = 1.E-8

    p(1) = atan2(e(3), SQRT(e(1) * e(1) + e(2) * e(2)))
    IF (abs(e(1)) < eps .AND. abs(e(2)) < eps) THEN
      p(2) = 0.0
    ELSE
      p(2) = atan2(-e(2), e(1))
    ENDIF

END SUBROUTINE xyz2ll

SUBROUTINE ll2xyz(p, e)
    IMPLICIT NONE
    REAL p(2)
    REAL e(3)

    e(1) = cos(p(1)) * cos(p(2))
    e(2) = cos(p(1)) * sin(p(2))
    e(3) = sin(p(1))

END SUBROUTINE ll2xyz

END MODULE slint
