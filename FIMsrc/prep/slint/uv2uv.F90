! =================================================================================!
! This file contains subroutines that parallelly translate a vector (u, v) from
! one specified lat/lon location (source) to another specified lat/lon location 
! (target).  
!
! SUBROUTINE uvs2uvt(uv_src, src_loc, uv_tgt, tgt_loc)
!   Input:
!     uv_src(2) - (u,v) components at source location
!     src_loc(2), tgt_loc(2) - source and target locations in lat/lon 
!   Output:
!     uv_tgt(2) - (u,v) components at source location
!   
! SUBROUTINE rot_coeff(src, tgt, cs, sn)
!   Input:
!     src(2), tgt(2) - source and target locations in lat/lon 
!   Output:
!     cs, sn - components of the transform matrix
!
! Ning Wang,  Mar. 2015, Original version
! ================================================================================!
SUBROUTINE uvs2uvt(uv_src, src_loc, uv_tgt, tgt_loc)
     IMPLICIT NONE

     REAL, INTENT(IN) :: uv_src(2), src_loc(2), tgt_loc(2)
     REAL, INTENT(OUT) :: uv_tgt(2)

     REAL :: cs, sn

     CALL rot_coeff(src_loc, tgt_loc, cs, sn)
  
     uv_tgt(1) = cs * uv_src(1) + sn * uv_src(2)
     uv_tgt(2) = -sn * uv_src(1) + cs * uv_src(2) 

END SUBROUTINE uvs2uvt
 
SUBROUTINE rot_coeff(src, tgt, cs, sn)
     IMPLICIT NONE

     REAL, INTENT(IN) :: src(2), tgt(2)
     REAL, INTENT(OUT) :: cs, sn
    
     REAL :: mf
  
     mf = 1.0+sin(tgt(1))*sin(src(1))+cos(tgt(1))*cos(src(1))*cos(src(2)-tgt(2))
     cs = (cos(tgt(1))*cos(src(1))+(1.0+sin(tgt(1))*sin(src(1)))*cos(src(2)-tgt(2)))/mf
     sn = -(sin(tgt(1))+sin(src(1)))*sin(src(2)-tgt(2))/mf 

END SUBROUTINE rot_coeff     
     
