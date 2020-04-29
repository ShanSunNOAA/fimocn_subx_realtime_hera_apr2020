!=============================================================
! GetTgtGrid subroutine 
!
! This subroutine computes the target grids of various 
! resolutions and projections. These grids data are used 
! in slint remapping subroutines.
!   
!
! Ning Wang, Aug. 2012,    original version. 
!
!
!=============================================================
     SUBROUTINE GetTgtGrid(grid_id, ll_tgt, csrot, mx, my)
       IMPLICIT NONE
        
       INTEGER grid_id, mx, my
       REAL ll_tgt(mx*my,2)
       REAL csrot(mx*my,2)
        
       INTEGER i, j, g_idx, ierr
       REAL d2r, pi
       REAL bl_lat, bl_lon, ur_lat, ur_lon, dx, elov, alatan, alat, alon, elon, orient
       REAL xi, xj, xi_lolf, xj_lolf, xi_org, xj_org, lat_range, lon_range
       REAL*8, ALLOCATABLE :: gl(:)
       REAL*4, ALLOCATABLE :: gl4(:)

       INTEGER :: gds(200), npts, nret
       REAL, ALLOCATABLE :: xpts(:),ypts(:)
       REAL, ALLOCATABLE :: dxdln(:),dxdlt(:),dydln(:),dydlt(:),area(:)

       PRINT*, "Target grid : grid#", grid_id
       pi = acos(-1.0)
       d2r = pi / 180.0
       csrot(:,1) = 1.0
       csrot(:,2) = 0.0

       SELECT CASE(grid_id)
         CASE(228, 45, 4, 3, 173, 174)   ! Global latitude/longitude grids
           IF (grid_id /= 173 .AND. grid_id /= 174) THEN
             xi_org = 0.0
             xj_org = 0.0
           ELSE IF (grid_id == 173) THEN 
             xi_org = 0.042 * d2r
             xj_org = 0.042 * d2r 
           ELSE IF (grid_id == 174) THEN 
             xi_org = 0.062 * d2r
             xj_org = 0.062 * d2r 
           ENDIF
           lat_range = pi - 2.0 * xj_org
           lon_range = 2.0 * pi
           DO i = 1, mx
             DO j = 1, my
               g_idx = (i + (j - 1) * mx)
               ll_tgt(g_idx, 1) = (REAL(j - 1) - REAL(my - 1) * 0.5)  * lat_range / REAL(my - 1)
               IF (grid_id /= 228) THEN
                 ll_tgt(g_idx, 1) = -ll_tgt(g_idx, 1)
               ENDIF
               ll_tgt(g_idx, 2) = xi_org + REAL(i - 1) * lon_range / REAL(mx)
             END DO
           END DO

         CASE(83)             ! Rotated lat/lon, 13 km Rapid Refresh grid
           npts = mx*my
           ALLOCATE(xpts(npts),ypts(npts))
           ALLOCATE(dxdln(1),dxdlt(1),dydln(1),dydlt(1),area(1))
!           gds(1) = 202
!           gds(2) = mx; gds(3) = my
!           gds(4) = 2228; gds(5) = -140481  
!           gds(6) = 136
!           gds(7) = mx; gds(8) = my 
!           gds(9) = 121; gds(10) = 121 
!           gds(11) = 64

           gds(1) = 203
           gds(2) = mx; gds(3) = my
           gds(4) = 2228; gds(5) = -140481  
           gds(6) = 136
           gds(7) = 47500; gds(8) = -104000 
           gds(9) = 121; gds(10) = 121 
           gds(11) = 64
           
           CALL GDSWZD(gds,0,npts,-9999.0,xpts,ypts,ll_tgt(:,2),ll_tgt(:,1), &
                       nret,1,csrot(:,1),csrot(:,2),0,dxdln,dxdlt,dydln,dydlt,area)
           ll_tgt = ll_tgt * d2r
           DEALLOCATE(xpts,ypts,dxdln,dxdlt,dydln,dydlt,area)

         CASE(244)             ! North Atlantic hurricane lat/lon grid
           bl_lat = -0.205 * d2r
           bl_lon = 261.75 * d2r
           ur_lat = 50.75 * d2r
           ur_lon = 330.25 * d2r
           DO i = 1, mx
             DO j = 1, my
               g_idx = (i + (j - 1) * mx)
               ll_tgt(g_idx, 1) =  ur_lat + REAL(j - 1)  * (bl_lat - ur_lat ) / REAL(my - 1)
               ll_tgt(g_idx, 2) = bl_lon + REAL(i - 1) * (ur_lon - bl_lon) / REAL(mx - 1)
             END DO
           END DO

         CASE(236)             ! Regional - CONUS (Lambert conformal)
           bl_lat = 16.281 
           bl_lon = 233.862  ! 233.862 = 360.0 - 126.138
           dx = 40633.25
           elov = 265.0 
           alatan = 25.0

           DO i = 1, mx
             DO j = 1, my
               CALL W3FB12(REAL(i),REAL(j),bl_lat,bl_lon,dx,elov,alatan,alat,elon,ierr)
               g_idx = (i + (j - 1) * mx)
               ll_tgt(g_idx, 1) = alat * d2r
               ll_tgt(g_idx, 2) = elon * d2r 
             END DO
           END DO

         CASE(130)             ! Regional - CONUS (Lambert conformal)
           bl_lat = 16.281 
           bl_lon = 233.862  ! 233.862 = 360.0 - 126.138
           dx = 13545.087
           elov = 265.0 
           alatan = 25.0

           DO i = 1, mx
             DO j = 1, my
               CALL W3FB12(REAL(i),REAL(j),bl_lat,bl_lon,dx,elov,alatan,alat,elon,ierr)
               g_idx = (i + (j - 1) * mx)
               ll_tgt(g_idx, 1) = alat * d2r
               ll_tgt(g_idx, 2) = elon * d2r 
             END DO
           END DO

         CASE(129)             ! Global - (Gaussian grid for T574)

           ALLOCATE(gl(my),stat=ierr)
           if (ierr.ne.0) then
             write (*,'(a)') 'GetTgtGrid: allocation of gl failed'
             call flush(6)
             stop
           endif

           ALLOCATE(gl4(my),stat=ierr)
           if (ierr.ne.0) then
             write (*,'(a)') 'GetTgtGrid: allocation of gl4 failed'
             call flush(6)
             stop
           endif

           CALL gaussian_lats(gl, my)
           gl4=gl
           xi_org = 0.0
           lon_range = 2.0 * pi 
           DO i = 1, mx
             DO j = 1, my
               g_idx = (i + (j - 1) * mx)
               ll_tgt(g_idx, 1) = gl(my - j + 1) * d2r
               ll_tgt(g_idx, 2) = xi_org + REAL(i - 1) * lon_range / REAL(mx)
             END DO
           END DO
           DEALLOCATE(gl,gl4)

         CASE(201)              ! Northern Hemisphere (polar stereographics)
           dx = 94.512 ! km, true at 60N
           orient = 105.0 

           alat = -20.826
           alon = 150.0
           CALL W3FB04(alat,alon,dx,orient,xi_lolf,xj_lolf)
           xi_org = (REAL(mx) + 1.0) / 2.0 
           xj_org = (REAL(my) + 1.0) / 2.0 
           DO i = 1, mx 
             DO j = 1, my
               xi = xi_lolf - (REAL(i) - 1.0)/(xi_org - 1.0) * xi_lolf
               xj = xj_lolf - (REAL(j) - 1.0)/(xj_org - 1.0) * xj_lolf
               CALL W3FB05 (xi, xj, dx, orient, alat, alon)
               g_idx = (i + (j - 1) * mx)
               ll_tgt(g_idx, 1) = alat * d2r
               ll_tgt(g_idx, 2) = -alon * d2r 
             END DO
           END DO

         CASE(219)              ! Northern Hemisphere (polar stereographics)
           dx = 25.4 ! km, true at 60N
           orient = 80.0 !(in west long) 

           alat = 25.032
           alon = 119.56 !(in west long)
           CALL W3FB04(alat,alon,dx,orient,xi_lolf,xj_lolf)
           xi_org = 191 
           xj_org = 231 
           DO i = 1, mx 
             DO j = 1, my
               xi = xi_lolf - (REAL(i) - 1.0)/(xi_org - 1.0) * xi_lolf
               xj = xj_lolf - (REAL(j) - 1.0)/(xj_org - 1.0) * xj_lolf
               CALL W3FB05 (xi, xj, dx, orient, alat, alon)
               g_idx = (i + (j - 1) * mx)
               ll_tgt(g_idx, 1) = alat * d2r
               ll_tgt(g_idx, 2) = -alon * d2r 
             END DO
           END DO
print*, "lltgt"
print*, '1st corner', ll_tgt(1,:)/d2r
print*, '2nd corner', ll_tgt(mx,:)/d2r
print*, '3rd corner', ll_tgt((my-1)*mx+1,:)/d2r
print*, '4th corner', ll_tgt(my*mx,:)/d2r
print*, 'xi_lolf, xj_lolf', xi_lolf, xj_lolf

         CASE(224)              ! Southern Hemisphere (polar stereographics)
           dx = 95.25 ! km, true at 60S
!           dx = 381.0 ! km, true at 60S
           orient = 285.0 

           alat = 20.826
           alon = 240.0
           CALL W3FB04(alat,alon,-dx,orient,xi_lolf,xj_lolf)
           xi_org = 129 
           xj_org = 129 
           DO i = 1, mx 
             DO j = 1, my
               xi = xi_lolf - (REAL(i) - 1.0)/(xi_org - 1.0) * xi_lolf
               xj = xj_lolf - (REAL(j) - 1.0)/(xj_org - 1.0) * xj_lolf
               CALL W3FB05 (xi, xj, -dx, orient, alat, alon)
               g_idx = (i + (j - 1) * mx)
               ll_tgt(g_idx, 1) = alat * d2r
               ll_tgt(g_idx, 2) = -alon * d2r 
             END DO
           END DO
print*, "lltgt"
print*, '1st corner', ll_tgt(1,:)/d2r
print*, '2nd corner', ll_tgt(mx,:)/d2r
print*, '3rd corner', ll_tgt((my-1)*mx+1,:)/d2r
print*, '4th corner', ll_tgt(my*mx,:)/d2r

         CASE(242)              ! Alaska domain (polar stereographics)
           dx = 11.25 ! km, true at 60S
           orient = 135.0  !(in west long)

           alat = 30.0
           alon = 173.0    !(in west long)
           CALL W3FB04(alat,alon,dx,orient,xi_lolf,xj_lolf)
           xi_org = 376.636 
           xj_org = 481.792  
           DO i = 1, mx 
             DO j = 1, my
               xi = xi_lolf - (REAL(i) - 1.0)/(xi_org - 1.0) * xi_lolf
               xj = xj_lolf - (REAL(j) - 1.0)/(xj_org - 1.0) * xj_lolf
               CALL W3FB05 (xi, xj, dx, orient, alat, alon)
               g_idx = (i + (j - 1) * mx)
               ll_tgt(g_idx, 1) = alat * d2r
               ll_tgt(g_idx, 2) = -alon * d2r 
             END DO
           END DO
print*, "lltgt"
print*, '1st corner', ll_tgt(1,:)/d2r
print*, '2nd corner', ll_tgt(mx,:)/d2r
print*, '3rd corner', ll_tgt((my-1)*mx+1,:)/d2r
print*, '4th corner', ll_tgt(my*mx,:)/d2r

          CASE DEFAULT
            WRITE(6,*) 'GetTgtGrid: invalid grid id'
            STOP
       END SELECT 
        
     END SUBROUTINE GetTgtGrid
        
SUBROUTINE gridid2mxmy(gridid, mx, my)
      IMPLICIT NONE
      INTEGER gridid, mx, my

      IF (gridid == 228) THEN
        mx = 144
        my = 73
      ELSEIF (gridid == 45) THEN
        mx = 288
        my = 145
      ELSEIF (gridid == 3) THEN
        mx = 360
        my = 181
      ELSEIF (gridid == 4) THEN
        mx = 720
        my = 361
      ELSEIF (gridid == 244) THEN
        mx = 275
        my = 203
      ELSEIF (gridid == 236) THEN
        mx = 151
        my = 113
      ELSEIF (gridid == 130) THEN
        mx = 451
        my = 337
      ELSEIF (gridid == 219) THEN
        mx = 385
        my = 465
      ELSEIF (gridid == 201) THEN
        mx = 259
        my = 259
      ELSEIF (gridid == 224) THEN
        mx = 257
        my = 257
      ELSEIF (gridid == 242) THEN
        mx = 553
        my = 425
      ELSEIF (gridid == 83) THEN
        mx = 758
        my = 567
      ELSEIF (gridid == 173) THEN
        mx = 4320
        my = 2160
      ELSEIF (gridid == 174) THEN
        mx = 2880
        my = 1440
      ELSEIF (gridid == 129) THEN
        mx = 1760
        my = 880
      ELSE
        WRITE(6,*) 'gridid2mxmy: invalid grid id'
      ENDIF

END SUBROUTINE gridid2mxmy
