!=================================================================
!  This utility program reads in the icosahedren grid information 
!  and a specified raw icosahedral model data file, and plots the
!  data, voronoi sell, grid point (optional), and grid mesh 
!  (optional) in orthographic or cylindrical equidistant 
!  projection.  
!
!  Ning Wang, April. 2008, original version
!  
!  N. Wang, July 2012, Added cylindrical equidistant projection.            
!==================================================================
#define FIMDATA
PROGRAM plot_icos
!     USE fimnamelist, ONLY: readnl,ginfofile,datafile,grid_level,&
!                            var_name,nvlvls,level,proj,latc,lonc,&
!                            extent, map_vis,cell_vis,ll_vis,ipn_label, &
!                            print_version

     IMPLICIT NONE

     ! variables for all icos grid information arrays
     REAL, ALLOCATABLE :: grid(:,:) ! grid location (ll)
     REAL, ALLOCATABLE :: grid_ll(:,:) ! grid location (ll)
     REAL, ALLOCATABLE :: edge(:,:,:,:) ! end points of edges (ll)
     REAL, ALLOCATABLE :: edge_tmp(:,:,:,:) ! end points of edges (ll)
     REAL, ALLOCATABLE :: real_temp(:)  
     INTEGER, ALLOCATABLE :: prox(:,:) ! neighbors' index
     INTEGER, ALLOCATABLE :: prox_tmp(:,:) ! neighbors' index
     INTEGER, ALLOCATABLE ::  nprox(:)
     INTEGER, ALLOCATABLE ::  nprox_tmp(:)
     INTEGER, ALLOCATABLE ::  perm(:), inv_perm(:)
     INTEGER, ALLOCATABLE ::  int_temp(:)
     INTEGER LIND(14)
     CHARACTER*10 llbs(14) 
     CHARACTER*7 :: label
     CHARACTER(len=6) stfmt
                                                                                
     ! variables to hold data
     CHARACTER(80) :: header(10)
     CHARACTER(16) :: head
     REAL, ALLOCATABLE :: vardata(:), vardata_level(:)

     ! variables for namelist
     CHARACTER(len=256) ginfofile,datafile
     CHARACTER(len=2) proj
     CHARACTER(len=8) var_name
     REAL latc, lonc, extent(2)
     INTEGER grid_level,nvlvls,level,map_vis,cell_vis,ll_vis,ipn_label,ierr
     LOGICAL print_version

     ! temp variables
     REAL*8 r2d, edge_i(6,2,2)
     REAL u, v, w, z, u_edge, v_edge, latlon(2), latlon1(2), latlon2(2), zero, label_sz
     REAL x(6), y(6), x1(8), y1(8), x2(8), y2(8), minv, maxv, value
     REAL left_bound, right_bound, upper_bound, lower_bound, bounds(4)
     INTEGER i, j, k, seqnum, lunout, color, draw,l, p1_egct, p2_egct, p_egct
     INTEGER inside, uv_inside, wz_inside, vtx_inside
     INTEGER digits, deciml 
     INTEGER nip, retv;
     CHARACTER ginfoheader(16)
     INTEGER isn, i2, i3
                                                                                
     DATA LIND / 2,3,4,5,6,7,8,9,10,11,12,13,14,15 /

     NAMELIST  /PlotIcos/ ginfofile, datafile, grid_level, var_name, &
       nvlvls, level, proj, latc, lonc, extent, map_vis, cell_vis, ll_vis, &
       ipn_label, print_version

      
     OPEN (10,file="plticos.nl",status='old',action='read',iostat=ierr)
     IF (ierr.ne.0) THEN
       WRITE (0,'(a)') 'ERROR: Could not open PlotIcos namelsit file for read.'
       STOP
     ENDIF
     READ (10,NML=PlotIcos,iostat=ierr)
     IF (ierr /= 0) THEN
       WRITE (0,'(a)') 'ERROR: Could not read PlotIcos namelist.'
       STOP
     ENDIF
     CLOSE(10)

     nip = 10*(2**grid_level)**2+2
     r2d = 180.0 / (4.0*ATAN(1.0))
     zero = 0.0
     IF (extent(1) == 180.0) THEN
       left_bound = -180. 
       right_bound = 180. 
       upper_bound = 90.  
       lower_bound = -90. 
     ELSE
       left_bound = max(-180., lonc-extent(2))
       right_bound = min(180., lonc+extent(2))
       upper_bound = min(90.,  latc+extent(1))
       lower_bound = max(-90., latc-extent(1))
     ENDIF
     IF (proj == 'CE') THEN
       PRINT*, 'The domain of plot is:'
       PRINT*, '(',lower_bound,left_bound,') --> (',upper_bound,right_bound ,')'
     ENDIF

     label_sz = .005 / REAL(grid_level)

! allocate memory for various arrays
     ALLOCATE(grid(2,nip))
     ALLOCATE(edge(6,2,2,nip))
     ALLOCATE(prox(6,nip))
     ALLOCATE(nprox(nip)) 
     ALLOCATE(perm(nip))
     ALLOCATE(inv_perm(nip))
     ALLOCATE(real_temp(nip * 2))
     ALLOCATE(int_temp(nip))

     ALLOCATE(vardata(max(4, nvlvls)*nip))
     ALLOCATE(vardata_level(nip))
     ALLOCATE(grid_ll(nip,2))	

! init and set up NCAR Graphics
     CALL opngks
     CALL dfclrs
     CALL mapstc('OU','PS')
     CALL MAPPOS (.05,.84,.01,.99)
     IF (proj == 'CE') THEN
       latc = 0.0
       lonc = 0.0
       CALL MAPPOS (.03,.97,.10,.95)
     ENDIF
     CALL maproj(proj,latc,lonc, zero)
     IF (proj == 'CE') THEN
       CALL mapset('CO', lower_bound,left_bound,upper_bound,right_bound)  
     ENDIF

     CALL mapint
     !CALL mapdrw

! reads in info file
     OPEN(10,file=ginfofile,form='unformatted')
#ifdef NIMDATA
     READ(10) grid_ll(1:nip,1),grid_ll(1:nip,2),prox(1:6,1:nip), &
              nprox(1:nip),edge(1:6,1:2,1:2,1:nip)
     CLOSE(10)
#endif
#ifdef FIMDATA
     READ(10) ginfoheader
     print*, ginfoheader
     READ(10) ginfoheader
     print*, ginfoheader
     READ(10) grid_ll(1:nip,1), grid_ll(1:nip,2)
     DO isn = 1,size(prox,1)
         READ(10) int_temp(1:nip)
         prox(isn,1:nip) = int_temp(1:nip)
     ENDDO
     READ(10) nprox(1:nip)
     DO i3 = 1,size(edge,3)
       DO i2 = 1,size(edge,2)
         DO isn = 1,size(edge,1)
           READ(10) real_temp(1:nip)
           edge(isn,i2,i3,1:nip) = real_temp(1:nip)
         ENDDO
       ENDDO
     ENDDO
#endif

! convert 'info' file's special curve order back to standard ij order 
     DO i = 1, nip
       inv_perm(i) = i 
     ENDDO

     DO i = 1, nip
       perm(inv_perm(i)) = i 
     ENDDO

     DO i = 1, nip
       grid(1:2,i) = grid_ll(inv_perm(i),1:2)
     END DO
     DEALLOCATE(grid_ll)

     ALLOCATE(prox_tmp(6,nip))
     ALLOCATE(nprox_tmp(nip))
     DO i = 1, nip
       prox_tmp(1:6, i) = prox(1:6,inv_perm(i))
       nprox_tmp(i) = nprox(inv_perm(i))
       DO j = 1, nprox_tmp(i)
         prox_tmp(j, i) = perm(prox_tmp(j, i))
       END DO
     END DO
     prox = prox_tmp
     nprox = nprox_tmp
     DEALLOCATE(prox_tmp)
     DEALLOCATE(nprox_tmp)

     ALLOCATE(edge_tmp(6,2,2,nip))	
     DO i = 1, nip
       edge_tmp(1:6, 1:2, 1:2, i) = edge(1:6,1:2, 1:2, inv_perm(i))
     END DO
     edge = edge_tmp
     DEALLOCATE(edge_tmp)

! reads in data file
     CALL  read_1var(datafile, vardata, var_name, nip, nvlvls)

! starts plotting .... first fill in colors in the plot with a data set
     DO i = 1, nip 
       vardata_level(i) = vardata((i-1) * nvlvls + level)
     END DO

     DEALLOCATE(vardata)
     minv = minval(vardata_level)
     maxv = maxval(vardata_level)

     PRINT *, 'Minmum and maximum values of the data set =',minv, maxv
     digits = INT(ALOG10(ABS(maxv)) + 1.0)
     IF (digits > 4) THEN
       deciml = 0
     ELSE
       deciml = 2
     ENDIF
     WRITE(stfmt,'("(F", i0, ".", i0, ")")') digits+2+deciml, deciml 
     IF (ABS(minv - INT(minv)) < 0.00001 .AND. ABS(maxv - INT(maxv)) < 0.00001) THEN
       DO i = 1, 14
         WRITE(UNIT=llbs(i), FMT=stfmt) (minv + REAL(i-1) * (maxv - minv) / (14.0-1.0))
       ENDDO
     ELSE 
       DO i = 1, 14
         WRITE(UNIT=llbs(i), FMT=stfmt) (minv + REAL(i) * (maxv - minv) / 14.0)
       ENDDO
     ENDIF

!     CALL GSFAIS(1)

     IF (proj == 'OR') THEN
     DO i = 1, nip 
       value = vardata_level(i)
       color = (value - minv) / (maxv - minv) * 14.0 + 2.0 
       CALL GSFACI(color) 
       draw = 1
       edge_i = edge(1:6,1:2,1:2,i) * r2d
       DO j = 1, nprox(i)
         latlon1(1) = edge_i(j, 1, 1)
         latlon1(2) = edge_i(j, 1, 2)
         CALL maptrn(latlon1(1), latlon1(2), u,v)
         latlon2(1) = edge_i(j, 2, 1)
         latlon2(2) = edge_i(j, 2, 2)
         CALL maptrn(latlon2(1), latlon2(2), w, z)
         IF (u < 10000 .AND. v < 10000 .AND. w < 10000 .AND. z < 10000) THEN
           x(j) = u
           y(j) = v
         ELSE
           draw = 0
         ENDIF
       ENDDO
       IF (draw == 1) THEN
         CALL GFA(nprox(i),x,y)
       END IF
       IF (cell_vis == 1) THEN
         CALL GSPLCI(1) ! set voronoi cell boundary color
         IF (draw == 1) THEN
           CALL polygon(nprox(i),x,y)
         END IF
       ENDIF
       IF (ipn_label==1) THEN
         WRITE(label,'(i7)') i
         j=log10(float(i))+1
         CALL pcloqu(REAL(grid(2,i)*r2d), REAL(grid(1,i)*r2d), label(8-j:7),label_sz,0.,0.)
       ENDIF
     END DO

     ELSE IF (proj == 'CE' .AND. left_bound == -180.0 .AND. right_bound == 180.0) THEN

!    CE projection, for whole Earth.
!    rendering the specified field
!    when a cell is splitted into two parts at 180E (-180W)
!              -180.0         0.0          180.0 
!                |             |             |
!                | __2__       |             |
!                |/3    \ 1    |             |
!                |       \     |           3/|
!                |       /     |           4\|
!                |\4____/6     |             |
!                |   5         |             |
!
     DO i = 1, nip 
       value = vardata_level(i)
       color = (value - minv) / (maxv - minv) * 14.0 + 2.0 
       CALL GSFACI(color) 
       p1_egct = 0
       p2_egct = 0
       draw = 1
       edge_i = edge(1:6,1:2,1:2,i) * r2d
       DO j = 1, nprox(i)
         latlon1(1) = edge_i(j, 1, 1)
         latlon1(2) = edge_i(j, 1, 2)
         CALL maptrn(latlon1(1), latlon1(2), u,v)
         latlon2(1) = edge_i(j, 2, 1)
         latlon2(2) = edge_i(j, 2, 2)
         CALL maptrn(latlon2(1), latlon2(2), w, z)
         IF (u < 10000 .AND. v < 10000 .AND. w < 10000 .AND. z < 10000) THEN
           IF (abs(u - w) < 90.0) THEN
             x(j) = u
             y(j) = v
             IF (p1_egct /= 0 .AND. j /= nprox(i) .AND.  u < 0.0) THEN
               x1(p1_egct+1) = w
               y1(p1_egct+1) = z
               p1_egct = p1_egct + 1
             ENDIF
             IF (p2_egct /= 0 .AND. j /= nprox(i) .AND.  u >= 0.0) THEN
               x2(p2_egct+1) = w
               y2(p2_egct+1) = z
               p2_egct = p2_egct + 1
             ENDIF
           ELSE    ! the edge corss the 180E (or 180W) longitude line
             x(j) = u
             y(j) = v
             IF (u < 0.0) THEN   ! from > -180 to < 180  (figure above)
               u_edge = -180.0
               v_edge = v + (z - v) / ((w - 360.0) - u) * (u_edge - u)  
               IF (p1_egct == 0) THEN ! gather all the vertices P1,.. before edge
                 DO k = 1, j 
                   x1(k) = x(k)
                   y1(k) = y(k)
                 ENDDO
                 p1_egct = j
               ENDIF
               x1(p1_egct+1) = u_edge
               y1(p1_egct+1) = v_edge
               p1_egct = p1_egct + 1
               x2(p2_egct+1) = -u_edge
               y2(p2_egct+1) = v_edge
               x2(p2_egct+2) = w
               y2(p2_egct+2) = z
               p2_egct = p2_egct + 2
             ELSE        ! from < 180 to > -180
               u_edge = 180.0
               v_edge = v + (z - v) / ((w + 360.0) - u) * (u_edge - u)  
               IF (p2_egct == 0) THEN ! gather all the vertices P1,.. before edge 
                 DO k = 1, j 
                   x2(k) = x(k)
                   y2(k) = y(k)
                 ENDDO
                 p2_egct = j
               ENDIF
               x2(p2_egct+1) = u_edge
               y2(p2_egct+1) = v_edge
               p2_egct = p2_egct + 1
               x1(p1_egct+1) = -u_edge
               y1(p1_egct+1) = v_edge
               x1(p1_egct+2) = w
               y1(p1_egct+2) = z
               p1_egct = p1_egct + 2
             ENDIF
             draw = 2
           ENDIF
         ELSE
!           PRINT*, latlon2(1), latlon2(2), u,v
!           PRINT*, latlon1(1), latlon1(2), w,z
            draw = 0
         ENDIF
       END DO

       IF (draw == 1) THEN
         CALL GFA(nprox(i),x,y)
       END IF
       IF (draw == 2) THEN
         IF (i == 1 .OR. i == nip) THEN
           CALL polar_cell_fill(y1(1), i)
         ELSE
           CALL GFA(p1_egct,x1,y1)
           CALL GFA(p2_egct,x2,y2)
         ENDIF
       END IF
       
       IF (cell_vis == 1) THEN
         CALL GSPLCI(1) ! set voronoi cell boundary color
         IF (draw == 1) THEN
           CALL polygon(nprox(i),x,y)
         ELSE IF (draw == 2) THEN
           IF (i == 1 .OR. i == nip) THEN
             CALL polar_cell(y1(1), i)
           ELSE
             CALL polygon(p1_egct,x1,y1)
             CALL polygon(p2_egct,x2,y2)
           ENDIF
         ENDIF
       ENDIF 

       IF (ipn_label==1) THEN
         WRITE(label,'(i7)') i
         j=log10(float(i))+1
         latlon(1) = REAL(grid(1,i)*r2d)
         latlon(2) = REAL(grid(2,i)*r2d)
         IF (latlon(2) > 180.0) latlon(2) = latlon(2) - 360.0
         CALL pcloqu(latlon(2), latlon(1), label(8-j:7),label_sz,0.,0.)
       ENDIF
     END DO

     ELSE IF (proj == 'CE' .AND. (left_bound /= -180.0 .OR. right_bound /= 180.0)) THEN

! Note: the current clipping algorithm is intended for processing multiple domains within 
! a defined map. It contains a konwn bug at the present implementation, that the four corner 
! cells could be clipped incorrectly, which needs to be fixed in a later version.  
! If only a single domain is required, one can also use mapset() subroutine to clip 
! the desired domain. 

     bounds(1) = left_bound
     bounds(2) = right_bound
     bounds(3) = lower_bound
     bounds(4) = upper_bound
     DO i = 1, nip 
       value = vardata_level(i)
       color = (value - minv) / (maxv - minv) * 14.0 + 2.0 
       CALL GSFACI(color) 
       p_egct = 0
       draw = 1
       vtx_inside = 0
       edge_i = edge(1:6,1:2,1:2,i) * r2d
       DO j = 1, nprox(i)
         latlon1(1) = edge_i(j, 1, 1)
         latlon1(2) = edge_i(j, 1, 2)
!         CALL maptrn(latlon1(1), latlon1(2), u,v)
         CALL maptrn(latlon1(1), latlon1(2), w,z)
         latlon2(1) = edge_i(j, 2, 1)
         latlon2(2) = edge_i(j, 2, 2)
!         CALL maptrn(latlon2(1), latlon2(2), w, z)
         CALL maptrn(latlon2(1), latlon2(2), u, v)
         IF (u < 10000 .AND. v < 10000 .AND. w < 10000 .AND. z < 10000) THEN
           uv_inside = inside(u,v,bounds)
           wz_inside = inside(w,z,bounds)
           IF (uv_inside == 1 .AND. wz_inside == 1) THEN !  both points in side  
             x(j) = u
             y(j) = v
             vtx_inside = 1
             IF (p_egct /= 0 ) THEN
               x1(p_egct+1) = u
               y1(p_egct+1) = v
               p_egct = p_egct + 1
             ENDIF
           ELSE IF (uv_inside < 0 .AND. wz_inside < 0) THEN ! both points outside 
             IF (draw == 1) draw = 0
             CYCLE   ! discard the line segment 
           ELSE IF (uv_inside < 0) THEN  ! (u,v) exits, (w, z) in domain
             draw = 3
             x(j) = w
             y(j) = z
             IF (uv_inside == -1 .OR. uv_inside == -2) THEN
               IF (uv_inside == -1) u_edge = left_bound
               IF (uv_inside == -2) u_edge = right_bound
               v_edge = v + (z-v)/(w-u)*(u_edge-u)
             ELSE IF (uv_inside == -3 .OR. uv_inside == -4) THEN
               IF (uv_inside == -3) v_edge = lower_bound 
               IF (uv_inside == -4) v_edge = upper_bound 
               u_edge = u + (w-u)/(z-v)*(v_edge-v)
             ENDIF
             IF (p_egct == 0 .AND. vtx_inside == 1) THEN  
               DO k = 1, j-1 
                 x1(k) = x(k)
                 y1(k) = y(k)
               ENDDO
               p_egct = j-1
             ENDIF
             x1(p_egct+1) = u_edge
             y1(p_egct+1) = v_edge
             p_egct = p_egct + 1
           ELSE IF (wz_inside < 0) THEN  ! (w,z) outside, (u,v) enter domain
             draw = 3
             IF (wz_inside == -1 .OR. wz_inside == -2) THEN
               IF (wz_inside == -1) u_edge = left_bound
               IF (wz_inside == -2) u_edge = right_bound
               v_edge = v + (z-v)/(w-u)*(u_edge-u)
             ELSE IF (wz_inside == -3 .OR. wz_inside == -4) THEN
               IF (wz_inside == -3) v_edge = lower_bound 
               IF (wz_inside == -4) v_edge = upper_bound 
               u_edge = u + (w-u)/(z-v)*(v_edge-v)
             ENDIF
             x1(p_egct+1) = u_edge
             y1(p_egct+1) = v_edge
             p_egct = p_egct + 1
             x1(p_egct+1) = u
             y1(p_egct+1) = v
             p_egct = p_egct + 1
           ENDIF
         ELSE
!           PRINT*, latlon2(1), latlon2(2), u,v
!           PRINT*, latlon1(1), latlon1(2), w,z
            draw = 0
         ENDIF
       END DO

       IF (draw == 1) THEN
         CALL GFA(nprox(i),x,y)
       ELSE IF (draw == 2) THEN
         IF (i == 1 .OR. i == nip) THEN
           CALL polar_cell_fill(y1(1), i)
         ELSE
           CALL GFA(p1_egct,x1,y1)
           CALL GFA(p2_egct,x2,y2)
         ENDIF
       ELSE IF (draw == 3) THEN
         IF (p_egct > 2) CALL GFA(p_egct,x1,y1)
       END IF
       
       IF (cell_vis == 1) THEN
         CALL GSPLCI(1) ! set voronoi cell boundary color
         IF (draw == 1) THEN
           CALL polygon(nprox(i),x,y)
         ELSE IF (draw == 2) THEN
           IF (i == 1 .OR. i == nip) THEN
             CALL polar_cell(y1(1), i)
           ELSE
             CALL polygon(p1_egct,x1,y1)
             CALL polygon(p2_egct,x2,y2)
           ENDIF
         ELSE IF (draw == 3) THEN
           CALL polygon(p_egct,x1,y1)
         ENDIF
       ENDIF 

       IF (ipn_label==1) THEN
         IF (draw /= 0 .AND. draw /= 3) THEN
           WRITE(label,'(i7)') i
           j=log10(float(i))+1
           latlon(1) = REAL(grid(1,i)*r2d)
           latlon(2) = REAL(grid(2,i)*r2d)
           IF (latlon(2) > 180.0) latlon(2) = latlon(2) - 360.0
           CALL pcloqu(latlon(2), latlon(1), label(8-j:7),label_sz,0.,0.)
         ENDIF
       ENDIF
     END DO
     ENDIF
     
! draw map and lat/lon lines as requested
     CALL GSPLCI(1) ! set color for map
     IF (map_vis == 1 .AND. ll_vis == 1) THEN
       CALL mapdrw
     ELSE IF (map_vis == 1 .AND. ll_vis == 0) THEN
       CALL maplot
     ELSE IF (map_vis == 0 .AND. ll_vis == 1) THEN
       CALL mapgrd
       CALL maplbl
     ENDIF

! draw color bars (with the appropriate colors for print version)
     IF (print_version) THEN
       CALL LBSETI ('CBL - COLOR OF BOX LINES',0) ! black color boxes
       CALL LBSETI ('CLB - COLOR OF LABELS',0) ! black color labels
     ELSE
       CALL LBSETI ('CBL - COLOR OF BOX LINES',1) ! white color boxes
       CALL LBSETI ('CLB - COLOR OF LABELS',1) ! white color labels
     ENDIF
     IF (proj /= 'CE') THEN
       CALL LBLBAR (1,.90,.99,.13,.87,14,.5,1.,LIND,0,LLBS,14,1)
     ELSE
       CALL LBLBAR (0,.13,.87,.01,.06,14,1.,.5,LIND,0,LLBS,14,1)
     ENDIF
     CALL clsgks

     DEALLOCATE(vardata_level)
     DEALLOCATE(nprox)
     DEALLOCATE(grid)
     DEALLOCATE(prox)
     DEALLOCATE(edge)

END PROGRAM plot_icos
 
SUBROUTINE polar_cell_fill(last_lat, i)
      IMPLICIT NONE
      
      REAL last_lat
      INTEGER i
      REAL x(4), y(4)
      
      x(1) = -180.0
      x(2) = 180.0
      x(3) = 180.0
      x(4) = -180.0

      IF (i == 1) THEN
        y(1) = 90.0 
        y(2) = 90.0 
      ELSE
        y(1) = -90.0 
        y(2) = -90.0 
      ENDIF
      y(3) = last_lat 
      y(4) = last_lat
      CALL GFA(4,x,y)

END SUBROUTINE polar_cell_fill

SUBROUTINE polar_cell(last_lat, i)
      IMPLICIT NONE
      
      REAL last_lat
      INTEGER i
      REAL x(4), y(4)
      
      x(1) = -180.0
      x(2) = 180.0
      x(3) = 180.0
      x(4) = -180.0

      IF (i == 1) THEN
        y(1) = 90.0 
        y(2) = 90.0 
      ELSE
        y(1) = -90.0 
        y(2) = -90.0 
      ENDIF
      y(3) = last_lat 
      y(4) = last_lat
      CALL polygon(4,x,y)

END SUBROUTINE polar_cell

SUBROUTINE polygon(n,x,y)
      IMPLICIT NONE
      
      INTEGER n
      REAL x(n), y(n)
      INTEGER i

      DO i = 1, n-1
        CALL line (x(i), y(i), x(i+1), y(i+1))
      ENDDO
      CALL line(x(n),y(n),x(1),y(1))

END SUBROUTINE polygon 

INTEGER FUNCTION inside (x,y,bounds)
      IMPLICIT NONE
    
      REAL x, y
      REAL bounds(4)

      IF (x < bounds(1)) THEN
        inside = -1
        RETURN
      ELSE IF ( x > bounds(2)) THEN
        inside = -2
        RETURN
      ELSE IF (y < bounds(3)) THEN
        inside = -3
        RETURN
      ELSE IF (y > bounds(4)) THEN 
        inside = -4
        RETURN
      ELSE
        inside = 1
      ENDIF

END FUNCTION inside

SUBROUTINE read_1var(filename, vardata, var_name, nip, nlevels)
      IMPLICIT NONE
      
      INTEGER ::  nip, nlevels
      CHARACTER(256) :: filename
      CHARACTER(8) :: var_name
      REAL vardata(nip * nlevels)

      CHARACTER(80) :: header(10)

      INTEGER i, j, its, time, lunout, is2Dvar

! when all 2d vars, sm3d, and st3d are in one file 
      IF (is2Dvar(var_name) == 1) THEN
        CALL read_2Dvar(filename, var_name, nip, nlevels, vardata)
        RETURN
      ENDIF

      OPEN (29,file=filename,form="unformatted")
      READ(29) header
      READ(29) vardata(1:nip*nlevels) 
      CLOSE(29)

END SUBROUTINE read_1var

SUBROUTINE read_2Dvar(filename, var_name, nip, nlevels, vardata)
      IMPLICIT NONE

      CHARACTER(len=256) :: filename
      CHARACTER(len= 8) :: var_name
      INTEGER  nip, nlevels
      REAL vardata(nip * 4)
      CHARACTER(len=80) :: header(10)
      CHARACTER(len=80) :: varname

      OPEN(29,file=filename,form="unformatted")
      READ(29) header
      READ(header,FMT="(4X A4)") varname
      CALL tolowercase(varname)

      DO WHILE (var_name(1:4) /= varname(1:4))
          IF (varname(1:4) == "sm3d" .OR. varname(1:4) == "st3d") THEN
            READ(29) vardata(1:nip*4)
          ELSE
            READ(29) vardata(1:nip)
          ENDIF
          READ(29) header
          READ(header,FMT="(4X A4)") varname
          CALL tolowercase(varname)
      END DO

      IF (var_name /= varname) THEN
        PRINT*, var_name, ' does not exist in the file ', filename
        STOP
      ENDIF

      READ(29) vardata(1:nip*nlevels)
      CLOSE(29)

END SUBROUTINE read_2Dvar


INTEGER FUNCTION is2Dvar(var_name)
      IMPLICIT NONE

      CHARACTER(len=8) :: var_name

      IF(var_name(3:4) == "2D" .OR. var_name(3:4) == "2d" .OR. & 
         var_name(1:4) == "vtyp" .OR. &
         var_name(1:4) == "sm3d" .OR. var_name(1:4) == "st3d" .OR. &
         var_name(1:4) == "shfx" .OR. var_name(1:4) == "lhfx" ) THEN
        is2Dvar = 1
      ELSE
        is2Dvar = 0
      ENDIF
END FUNCTION is2Dvar


SUBROUTINE tolowercase(str)
      IMPLICIT NONE

      CHARACTER*80 str
      INTEGER len, i

      len = LEN_TRIM(str)
      DO i = 1, len
        IF (ichar(str(i:i)) >= ichar('A') .AND. ichar(str(i:i)) <=  ichar('Z')) THEN
          str(i:i) = char(ichar(str(i:i)) + 32)
        END IF
      END DO
END SUBROUTINE tolowercase
