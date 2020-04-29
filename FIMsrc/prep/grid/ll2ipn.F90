!=============================================================
!
! Ning Wang, March 2007 
!
!=============================================================
   PROGRAM ll2ipn
#ifdef NAG
   USE F90_UNIX_ENV, ONLY : getarg, iargc
#endif
   USE kd, ONLY:init_kd_tree,knn_search 
   IMPLICIT NONE
 
! command line argument count
#ifndef NAG
   INTEGER :: iargc
#endif
   CHARACTER(len=128) :: lat_s, lon_s
   CHARACTER(len=128) :: grid_file
   CHARACTER(len=16)  :: header
   INTEGER :: glvl, nip, curve, num_nn
   REAL lat, lon, d2r

   REAL ll(2), min_dist(3), hp1(3), hp2(3)
   REAL, ALLOCATABLE :: llpoints(:,:)
   INTEGER nn(3)

   IF (iargc() .NE. 2 ) THEN
     WRITE(0,*) 'Usage: ll2ipn [lat] [lon]'
     STOP
   END IF
   
   CALL getarg(1,lat_s)
   CALL getarg(2,lon_s)
   READ(lat_s, *) lat
   READ(lon_s, *) lon

   d2r = acos(-1.0) / 180.0

   grid_file = "icos_grid_info_level.dat" 
   OPEN (10,file=grid_file,form='unformatted')
   READ (10) header
   READ (header,"('glvl =',I2)") glvl
   READ (10) header
   nip = 10 * (2 ** (2 * glvl)) + 2 
   ALLOCATE(llpoints(nip,2))
   READ(10) llpoints(1:nip, 1), llpoints(1:nip,2)
   CALL init_kd_tree(llpoints, nip, 3)

   ll(1) = lat * d2r
   ll(2) = lon * d2r
   min_dist = 1.0
   hp1 = 0.0
   hp2 = 0.0
   
   CALL knn_search(ll, nn, min_dist, hp1, hp2, 1.0, num_nn)

   PRINT*, 'The three nearest neighbors :', nn(1), nn(2), nn(3)

   END PROGRAM ll2ipn
       

