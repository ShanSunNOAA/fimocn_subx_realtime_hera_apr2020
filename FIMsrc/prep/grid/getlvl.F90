!******************************************************************************
   program getlvl
!  Loads the initial variables and constants to start sgm
!  Alexander E. MacDonald  11/27/05
!  J. Lee                  September, 2005
!  N.W.     Changed to dynamic memory allocation for big arrays.
!           Generalize the computation of nip to work with the icos grid that 
!           has mixed bi-section and tri-section subdivisions. 
!           Changed to a more accurate computation of the middle points of the 
!           Voronoi-cell edges.
!  N.W.     Added additional pre-processing for higher order flux computation. 
!            
!******************************************************************************
use headers,only: testcurveheader,testglvlheader,writecurveheader,writeglvlheader
use fimnamelist   ,only: readnl_serial, glvl,curve, SubdivNum
use LSfitting, only: initLSmod, closeLSmod, pseudoInverse
implicit none

integer, parameter :: npp=6
integer, parameter ::  nd=2
real, parameter :: pi = 3.1415926535897
real, parameter :: ae = 6371220.   !earth radius

!
real*4 map
real*4 conr_xy(npp,2,2)
real*4 prox_xy(npp,2)    ! holds x and y locs for prox pts (m)

real*4 loc_xy(2,6), pinv(5,6)

!
real*4 eltp(4),elnp(4) ! 4 lat/lon surrounding a particular edge
real*4 conr_tmp(1:6,1:2) 

!
real*4, allocatable :: conr_ll(:,:,:,:), lle(:,:,:)
real*4, allocatable :: corner_xy(:,:,:), hfcoef(:,:,:)
integer, allocatable :: nprox(:), prox(:, :), proxs(:, :), inv_perm(:)
real*4, allocatable :: lat(:),lon(:),area(:)
real*4, allocatable :: sideln(:, :),rprox_ln(:, :)
real*4, allocatable :: sidevec_c(:, :, :),sidevec_e(:, :, :)
real*4, allocatable :: pxy(:, :, :)
real*4, allocatable :: cs(:,:,:),sn(:,:,:)

real*4 ap(5), tri_area_r, xim1, xi, yim1, yi

integer :: ipn, isn, ism, ixy, ipt, iprox, ip1, im1,  j, idx, i2, i3
integer :: ibf, ico
integer :: i, sl, nip 
real :: xlat, xlon, xxp, yyp, xxm, yym, xx, yy, xltc, xlnc, rf, sm
real :: p1(2), p2(2), pm(2)

real*4 dist

     CALL readnl_serial()
     sl = 1
     DO i = 1, glvl 
       IF (SubdivNum(i) == 2 .OR. SubdivNUm(i) == 3) THEN
         sl = sl * SubdivNum(i)
       ELSE
         PRINT*, "Namelist variable 'SubdivNum' specification is incomplete."
         PRINT*, "SubdivNum(i) must be either 2 or 3."
         STOP
       ENDIF
     ENDDO 

! num of icos grid points
     nip = 10 * sl * sl + 2

! allocate memory for arrays
ALLOCATE(conr_ll(npp,2,2,nip))
ALLOCATE(corner_xy(npp,2,nip))
ALLOCATE(hfcoef(6,npp,nip))
ALLOCATE(lle(npp,2,nip))
ALLOCATE(nprox(nip))
ALLOCATE(prox(npp, nip))
ALLOCATE(proxs(npp, nip))
ALLOCATE(lat(nip))
ALLOCATE(lon(nip))
ALLOCATE(area(nip))
ALLOCATE(sideln(npp, nip))
ALLOCATE(rprox_ln(npp, nip))
ALLOCATE(sidevec_c(nd, npp, nip))
ALLOCATE(sidevec_e(nd, npp, nip))
ALLOCATE(pxy(5, npp, nip))
ALLOCATE(cs(4,npp,nip))
ALLOCATE(sn(4,npp,nip))
ALLOCATE(inv_perm(nip))

print*, 'start getlvl ... '
!...................................................................

cs       (:,6,:) = 0.0
sn       (:,6,:) = 0.0
sidevec_c(:,6,:) = 0.0
sidevec_e(:,6,:) = 0.0
sideln   (  6,:) = 0.0
rprox_ln (  6,:) = 0.0
corner_xy(6,:,:) = 0.0
pxy      (:,6,:) = 0.0
hfcoef   (:,6,:) = 0.0

OPEN                (10,file='icos_grid_info_level.dat',form='unformatted')
call TestGlvlHeader (10,     'icos_grid_info_level.dat','getlvl',glvl)
call TestCurveHeader(10,     'icos_grid_info_level.dat','getlvl',curve)
READ(10) lat,lon
do isn = 1,size(prox,1)
  READ(10) prox(isn,:)
enddo
READ(10) nprox
do i3 = 1,size(conr_ll,3)
  do i2 = 1,size(conr_ll,2)
    do isn = 1,size(conr_ll,1)
      READ(10) conr_ll(isn,i2,i3,:)
    enddo
  enddo
enddo
READ(10) inv_perm
CLOSE(10)
!
do ipn=1,nip
if(lon(ipn).lt.0.) lon(ipn)=lon(ipn)+2.*pi
end do
!
do ipn=1,nip
!
! shift the corner index back by 1, conr(i) = conr(i-1)
do isn=1,nprox(ipn)
do ixy=1,2
conr_tmp(isn,ixy)=conr_ll(isn,1,ixy,ipn)
end do
end do
!
do isn=1,nprox(ipn)
ism=isn-1
if(isn.eq.1) ism=nprox(ipn)
do ixy=1,2
conr_ll(isn,1,ixy,ipn)=conr_tmp(ism,ixy)
conr_ll(isn,2,ixy,ipn)=conr_tmp(isn,ixy)
end do
end do
!
end do
!
! calculate the edge indexes of the neighboring cells
proxs=-99
do ipn=1,nip
do 20 isn=1,nprox(ipn)
iprox=prox(isn,ipn)
do j=1,nprox(iprox)
if(prox(j,iprox).eq.ipn) then
proxs(isn,ipn)=j
goto 20
end if
end do
20 continue
end do
!
!...................................................................
!
do ipn=1,nip
do isn=1,nprox(ipn)
p1(1:2) = conr_ll(isn,1,1:2,ipn)
p2(1:2) = conr_ll(isn,2,1:2,ipn)
call middle_r(p1, p2, pm)
lle(isn, 1:2, ipn) = pm(1:2)
end do
end do
!
!........................................................................
!Caculate sidevec (in mid-edge and cell-center projections)
!
do ipn=1,nip
do isn=1,nprox(ipn)
xlon=lon(prox(isn,ipn))
xlat=lat(prox(isn,ipn))
call ll2xy(lon(ipn),lat(ipn),xlon,xlat,prox_xy(isn,1),prox_xy(isn,2))
rprox_ln(isn,ipn)=1./(ae*sqrt(prox_xy(isn,1)**2+prox_xy(isn,2)**2))
do ipt=1,2
xlon=conr_ll(isn,ipt,2,ipn)
xlat=conr_ll(isn,ipt,1,ipn)
call ll2xy(lon(ipn),lat(ipn),xlon,xlat,conr_xy(isn,ipt,1),conr_xy(isn,ipt,2))
end do
map=2./(1.+sin(lle(isn,1,ipn))*sin(lat(ipn)) &
+cos(lle(isn,1,ipn))*cos(lat(ipn)) &
       *cos(lle(isn,2,ipn)-lon(ipn)))
do ixy=1,nd
sidevec_c(ixy,isn,ipn)=ae*( conr_xy(isn,2,ixy) &
                         -conr_xy(isn,1,ixy)) *map
end do

call ll2xy(lle(isn,2,ipn),lle(isn,1,ipn),conr_ll(isn,2,2,ipn) &
,conr_ll(isn,2,1,ipn),xxp,yyp)
call ll2xy(lle(isn,2,ipn),lle(isn,1,ipn),conr_ll(isn,1,2,ipn) &
,conr_ll(isn,1,1,ipn),xxm,yym)
sidevec_e(1,isn,ipn)= ae*(xxp-xxm)
sidevec_e(2,isn,ipn)= ae*(yyp-yym)
sideln(isn,ipn)=sqrt(sidevec_e(1,isn,ipn)**2+sidevec_e(2,isn,ipn)**2)
end do ! isn loop

area(ipn)=0.   
do isn=1,nprox(ipn)
xx=ae*.5*(conr_xy(isn,2,1)+conr_xy(isn,1,1))
yy=ae*.5*(conr_xy(isn,2,2)+conr_xy(isn,1,2))
area(ipn)=area(ipn)+.5*(xx*sidevec_c(2,isn,ipn)-yy*sidevec_c(1,isn,ipn))
end do

do ixy = 1, 2
  do isn=1,nprox(ipn)
    corner_xy(isn,ixy,ipn) = ae*conr_xy(isn,1,ixy)
  enddo
enddo

end do ! ipn loop
!
!...................................................................
! Pre-processing for higher order flux computation 
!
pxy(:,:,:) = 0.0
do ipn = 1, nip
  do isn = 1, nprox(ipn)
    do ipt=1,2
      xlon=conr_ll(isn,ipt,2,ipn)
      xlat=conr_ll(isn,ipt,1,ipn)
      call ll2xy(lon(ipn),lat(ipn),xlon,xlat,conr_xy(isn,ipt,1),conr_xy(isn,ipt,2))
     end do
  enddo
  ap(:) = 0.0
  do isn = 1, nprox(ipn)
    xim1 = conr_xy(isn,1,1); xi = conr_xy(isn,2,1)
    yim1 = conr_xy(isn,1,2); yi = conr_xy(isn,2,2)
    xx=ae*.5*(xim1+xi)
    yy=ae*.5*(yim1+yi)
    tri_area_r=.5*(xx*sidevec_c(2,isn,ipn)-yy*sidevec_c(1,isn,ipn))/area(ipn)
    ap(1) =  ap(1) + tri_area_r/3.*(xim1 + xi) 
    ap(2) =  ap(2) + tri_area_r/3.*(yim1 + yi)
    ap(3) =  ap(3) + tri_area_r/6.*(xim1*xim1+xim1*xi+xi*xi)
    ap(4) =  ap(4) + tri_area_r/12.*(2.*xim1*yim1+2.*xi*yi+xim1*yi+xi*yim1)
    ap(5) =  ap(5) + tri_area_r/6.*(yim1*yim1+yim1*yi+yi*yi)
  enddo
  do isn = 1, nprox(ipn)
    xim1 = conr_xy(isn,1,1); xi = conr_xy(isn,2,1)
    yim1 = conr_xy(isn,1,2); yi = conr_xy(isn,2,2)
    pxy(1,isn,ipn) = .5*(xim1+xi) - ap(1) 
    pxy(2,isn,ipn) = .5*(yim1+yi) - ap(2) 
    pxy(3,isn,ipn) = 1./3.*(xim1*xim1+xi*xi+xim1*xi) - ap(3) 
    pxy(4,isn,ipn) = 1./6.*(2.*xim1*yim1+2.*xi*yi+xim1*yi+xi*yim1) - ap(4) 
    pxy(5,isn,ipn) = 1./3.*(yim1*yim1+yi*yi+yim1*yi) - ap(5) 
  enddo
enddo

call initLSmod(6,5,nip)
do ipn = 1, nip
!  loc_xy(:,nprox(ipn)+1) = 0.0   ! origin
  do isn=1,nprox(ipn)
    xlon=lon(prox(isn,ipn))
    xlat=lat(prox(isn,ipn))
    call ll2xy(lon(ipn),lat(ipn),xlon,xlat,loc_xy(1, isn),loc_xy(2, isn))
  enddo
  CALL pseudoInverse(loc_xy, nprox(ipn), 5, pinv)
  do isn=1,nprox(ipn)
    hfcoef(:,isn,ipn) = 0.0
    do ico = 1,nprox(ipn)
      do ibf = 1, 5
        hfcoef(ico,isn,ipn) = hfcoef(ico,isn,ipn) + pxy(ibf,isn,ipn) * pinv(ibf,ico) 
      enddo
    enddo 
  enddo
enddo 

!...................................................................
!
do ipn=1,nip
  do isn=1,nprox(ipn) 
    xltc=lle(isn,1,ipn)
    xlnc=lle(isn,2,ipn)
    ip1=mod(isn,nprox(ipn))+1
    im1=isn-1
    if(im1.eq.0) im1=nprox(ipn)
    eltp(1)=lat(ipn)
    elnp(1)=lon(ipn)
    eltp(2)=lat(prox(isn,ipn))
    elnp(2)=lon(prox(isn,ipn))
    eltp(3)=lat(prox(im1,ipn))
    elnp(3)=lon(prox(im1,ipn))
    eltp(4)=lat(prox(ip1,ipn))
    elnp(4)=lon(prox(ip1,ipn))
    do ipt=1,4
      rf=1.0/(1.0+sin(xltc)*sin(eltp(ipt))+cos(xltc)*cos(eltp(ipt))*cos(elnp(ipt)-xlnc))
      cs(ipt,isn,ipn)=rf*( cos(xltc)*cos(eltp(ipt))+(1.0+sin(xltc)*sin(eltp(ipt)))*cos(elnp(ipt)-xlnc))
      sn(ipt,isn,ipn)=-rf*sin(elnp(ipt)-xlnc)*(sin(xltc)+sin(eltp(ipt)))
    end do
  enddo
enddo

open(unit=28,file="glvl.dat", form="unformatted")
call WriteGlvlHeader (28,glvl )
call WriteCurveHeader(28,curve)
write(28) lat
write(28) lon
write(28) nprox
do isn=1,size(proxs,1)
  write(28) proxs(isn,:)
enddo
do isn=1,size(prox,1)
  write(28) prox(isn,:)
enddo
write(28) area
do isn=1,size(cs,2)
  do idx=1,size(cs,1)
    write(28) cs(idx,isn,:)
  enddo
enddo
do isn=1,size(sn,2)
  do idx=1,size(sn,1)
    write(28) sn(idx,isn,:)
  enddo
enddo
do isn=1,size(sidevec_c,2)
  do idx=1,size(sidevec_c,1)
    write(28) sidevec_c(idx,isn,:)
  enddo
enddo
do isn=1,size(sidevec_e,2)
  do idx=1,size(sidevec_e,1)
    write(28) sidevec_e(idx,isn,:)
  enddo
enddo
do isn=1,size(sideln,1)
  write(28) sideln(isn,:)
enddo
do isn=1,size(rprox_ln,1)
  write(28) rprox_ln(isn,:)
enddo
write(28) inv_perm
do ixy=1,size(corner_xy,2)
  do isn=1,size(corner_xy,1)
    write(28) corner_xy(isn,ixy,:)
  enddo
enddo
do isn=1,size(hfcoef,2)
  do idx=1,size(hfcoef,1)
    write(28) hfcoef(idx,isn,:)
  enddo
enddo
close(28)

print*, 'done saving glvl.dat'
!...................................................................
!
call closeLSmod()
DEALLOCATE(pxy)
DEALLOCATE(conr_ll)
DEALLOCATE(corner_xy)
DEALLOCATE(hfcoef)
DEALLOCATE(lle)
DEALLOCATE(nprox)
DEALLOCATE(prox)
DEALLOCATE(proxs)
DEALLOCATE(lat)
DEALLOCATE(lon)
DEALLOCATE(area)
DEALLOCATE(sideln)
DEALLOCATE(rprox_ln)
DEALLOCATE(sidevec_c)
DEALLOCATE(sidevec_e)
DEALLOCATE(cs)
DEALLOCATE(sn)
DEALLOCATE(inv_perm)
stop  
end program getlvl
!
!!!
!
!#############################################################
!    ll2xy.f
!    Convert lat/lon to (x,y) on General Stereographic Coordinate (GSTC).
!    Original program:  J.Lee - 2004
!    Program testing:   J.Lee - 2004
!    Modified for Non-Structure Grid:  J.Lee - 2004
!############################################################

!    Purpose:  Given latitude and longitude on Spherical coordinate,
!              this subroutine computes X and Y coordinates on GSTC.
!    Reference: J.Lee, G. Browning, and Y. Xie: 
!               TELLUS (1995), p.892-910.
!               
! Input Variables : Angles are assumed in unit of "radian" 
! 
!     (latc,lonc) : the GSTC projected point.
!     ( lat, lon) : Input lat/lon in radians.
!
! OUTPUT Variables: 
!
!             xm  : X-Coordinate values on GSTC.  
!                   positive to East of central longitude
!             ym  : Y-Coordinate values on GSTC.  
!                   positive to North of central latitude.
! Note:  Output variables of xm and ym are 
!               nondimensionalized with "ae", the radius of earth.
!
subroutine ll2xy(lonc,latc,lon,lat,xm,ym)
!
implicit none

integer i
real*4 lonc,latc,lon,lat,mf
real*4 xm,ym
!
mf=2.0/(1.0+sin(lat)*sin(latc)+cos(lat)*cos(latc) &
       *cos(lon-lonc))
xm=mf*(cos(lat)*sin(lon-lonc))
ym=mf*((sin(lat)*cos(latc)-cos(lat) &
        *sin(latc)*cos(lon-lonc)) )
!
return
end

!===============================================================================
!  This subroutine computes the latitude and longitude 
!  of the middle point between two given ponits. The 
!  subroutine is similar to what is in the middle.F90,
!  except that its input and output variables have 
!  single precision and use radians for lat-lon values.
!
!   Ning Wang,   March, 2006
!===============================================================================

SUBROUTINE middle_r(p1,p2,p)

     IMPLICIT NONE

     REAL :: pi, d2r

     ! Two given points in lat/lon:
     REAL :: p1(2),p2(2),p(2)

     REAL :: xyz1(3),xyz2(3),xyz(3)

     pi = atan(1.0) * 4.0
     ! Convert them into Cardesian coor:
     xyz1(1) = cos(p1(1)) * cos(p1(2))
     xyz1(2) = cos(p1(1)) * sin(p1(2))
     xyz1(3) = sin(p1(1))

     xyz2(1) = cos(p2(1)) * cos(p2(2))
     xyz2(2) = cos(p2(1)) * sin(p2(2))
     xyz2(3) = sin(p2(1))

     ! middle point:

     xyz = 0.5 * (xyz1 + xyz2)
     xyz = xyz / sqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2) + xyz(3)*xyz(3))

     ! Convert the middle point to lat/lon coor:
     p(1) = atan2(xyz(3), sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2)))
     p(2) = atan2(xyz(2), xyz(1))

     IF (p(2) < 0.0) THEN
       p(2) = p(2) + 2 * pi
     END IF

END SUBROUTINE middle_r


REAL*4 FUNCTION dist(lat1, lon1, lat2, lon2)

IMPLICIT NONE

REAL*4, INTENT(IN) :: lat1, lon1, lat2, lon2
REAL*4             :: x

x = COS(lat1) * COS(lat2) * COS(lon1 - lon2) + SIN(lat1) * SIN(lat2)
if(abs(x) >= 1.0) then
  dist = 0.0
else
  dist = 6371.220 * ACOS(x)
endif

END FUNCTION dist
