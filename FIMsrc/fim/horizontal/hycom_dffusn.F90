module hycom_dffusn
use findmaxmin2
use stencilprint
contains
!*********************************************************************
!     dffusn_lyr
!	Diffuse layer variable (thickness-weighted for conservation)
!	S. Sun                           September 2009
!*********************************************************************

  subroutine dffusn_lyr (fld, delp, kdm, kfrst, klast, dffveldt, vrbos)

  use module_control  ,only: npp,nip
  use module_constants,only: prox,rarea,sideln,nedge,perm
  use hycom_constants ,only: wet
  use hycom_control   ,only: kgap
  use fimnamelist     ,only: itest

  implicit none
!  External type and dimension:
  integer,intent (IN)    :: kdm			! vert.dimension of fld,delp
  integer,intent (IN)    :: kfrst,klast		! range of layers to operate on
  real   ,intent (IN)    :: dffveldt		! diffn.velocity x time step
  logical,intent (IN),optional :: vrbos
!SMS$DISTRIBUTE (dh,2) BEGIN
  real   ,intent (INOUT) :: fld (kdm,nip)	! field to be diffused
  real   ,intent (IN)    :: delp(kdm,nip)	! layer thickness
! Local variables:
  real    :: flxdv(kdm,nip)	! line integral of flux along the edge
!SMS$DISTRIBUTE END

! Local variables
  integer   :: k		! layer index
  integer   :: ipn		! icos point index
  integer   :: ipnx             ! neighbor across joint edge
  integer   :: edg              ! icos edge index
  real      :: factor
  character :: string*8
  real,parameter :: thshld = 1.e-11
  real      :: hfharm,a,b
  hfharm(a,b) = a*b/(a+b)	! harmonic average x 0.5

  if (present(vrbos)) then
   if (vrbos) then
    print *,'entering ocn dffusn_lyr, dffveldt =',dffveldt
    do k=kfrst,klast
     if (kfrst.eq.klast .or. mod(k,kgap).eq.1) then
      write (string,'(a,i2)') ' k=',k
      call findmxmn2(fld,kdm,nip,k,'ocn dff_lyr  IN'//trim(string),wet)
     end if
    end do
    call stencl(fld,kdm,10.,'(ocn dffusn_lyr) field x 10 before smoothing')
   end if
  end if

!SMS$PARALLEL (dh,ipn) BEGIN
!SMS$EXCHANGE(delp,fld)
!$OMP PARALLEL DO PRIVATE (ipnx)
  do ipn = 1,nip			! horizontal loop
   flxdv(:,ipn) = 0.
   if (wet(ipn) > 0 ) then
    do edg = 1,nedge(ipn)
     ipnx = prox(edg,ipn)		! neighbor across shared edge
     if (wet(ipnx) > 0 ) then
      do k = kfrst,klast
       flxdv(k,ipn) = flxdv(k,ipn)+(fld(k,ipn)-fld(k,ipnx))		&
         *sideln(edg,ipn)*2.*hfharm(max(delp(k,ipn ),thshld)		&
                                   ,max(delp(k,ipnx),thshld))
       if (present(vrbos)) then
        if (vrbos .and. ipn.eq.itest .and. mod(k,kgap).eq.1) then
         write (*,'(a,2i8,i3,a,2es12.4)') 'ocn dffusn_lay ipn,ipnx,k=',	&
          perm(ipn),perm(ipnx),k,' fld=',fld(k,ipn),fld(k,ipnx)
        end if
       end if
      end do				! loop through layers
     end if				! open ocean (non-coastal) edge
    end do				! loop through edges
   end if				! ocean point
  end do				! horizontal loop
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE (factor)
  do ipn = 1,nip			! horizontal loop
   if (wet(ipn) > 0 ) then
    do k = kfrst,klast			! loop through layers
     factor = -dffveldt*rarea(ipn)/max(delp(k,ipn),thshld)
     fld(k,ipn) = fld(k,ipn) + flxdv(k,ipn)*factor

     if (present(vrbos)) then
      if (vrbos .and. ipn.eq.itest .and. mod(k,kgap).eq.1) then
       write (*,'(i8,i3,a,3es12.4)') perm(ipn),k,			&
         ' ocn_dffusn_lay flxdv,tdcy,field=',flxdv(k,ipn),		&
         flxdv(k,ipn)*factor,fld(k,ipn)
      end if
     end if

    end do				! loop through layers
   end if				! ocean point
  end do				! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  if (present(vrbos)) then
   if (vrbos) then
    do k=kfrst,klast
     if (kfrst.eq.klast .or. mod(k,kgap).eq.1) then
      write (string,'(a,i2)') ' k=',k
      call findmxmn2(fld,kdm,nip,k,'ocn dff_lyr OUT'//trim(string),wet)
     end if
    end do
    call stencl(fld,kdm,10.,'(ocn dffusn_lyr) field x 10 after smoothing')
   end if
  end if

  return
  end subroutine dffusn_lyr


!*********************************************************************
!     dffusn_lev
!	Diffuse level variable (no thickness-weighting)
!	S. Sun                           September 2009
!*********************************************************************

  subroutine dffusn_lev (fld, kdm, kfrst, klast, dffveldt, vrbos,	&
                         special_mask)

  use module_control  ,only: npp,nip
  use module_constants,only: prox,rarea,sideln,nedge,perm
  use hycom_constants ,only: wet
  use hycom_control   ,only: kgap
  use fimnamelist     ,only: itest

  implicit none
!  External type and dimension:
  integer,intent (IN)    :: kdm			! vert.dimension of fld
  integer,intent (IN)    :: kfrst,klast		! range of layers to operate on
  real   ,intent (IN)    :: dffveldt		! diffn.velocity x time step
  logical,intent (IN)    :: vrbos
!SMS$DISTRIBUTE (dh,1) BEGIN
  integer,intent (IN),optional :: special_mask(nip)	! mask replacing 'wet'
  integer                :: mask(nip)		! mask used in diffusing 'fld'
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
  real   ,intent (INOUT) :: fld(kdm,nip)	! field to be diffused
! Local variables:
  real    :: flxdv(kdm,nip)	! line integral of flux along the edge
!SMS$DISTRIBUTE END

! Local variables
  integer   :: k		! layer index
  integer   :: ipn		! icos point index
  integer   :: ipnx             ! neighbor across joint edge
  integer   :: edg              ! icos edge index
  real      :: factor
  character :: string*8

  if (present(special_mask)) then	! special mask different from 'wet'
    mask(:)=special_mask(:)
  else
    mask(:)=wet(:)		! default mask is ocean points 
  end if

  if (vrbos) then
    print *,'entering ocn dffusn_lev, dffveldt =',dffveldt
    do k=kfrst,klast
    if (kfrst.eq.klast .or. mod(k,kgap).eq.1) then
      write (string,'(a,i2)') ' k=',k
      call findmxmn2(fld,kdm,nip,k,'ocn dff_lev  IN'//trim(string),mask)
    end if
    end do
    call stencl(fld,kdm,100.,'(ocn dffusn_lev) field x 100 before smoothing')
  end if

!SMS$PARALLEL (dh,ipn) BEGIN
!SMS$EXCHANGE(fld,mask)
!$OMP PARALLEL DO PRIVATE (ipnx)
  do ipn = 1,nip			! horizontal loop
   flxdv(:,ipn) = 0.
   if (mask(ipn) > 0 ) then
    do edg = 1,nedge(ipn)
     ipnx = prox(edg,ipn)		! neighbor across shared edge
     if (mask(ipnx) > 0 ) then
      do k = kfrst,klast		! loop through levels
       flxdv(k,ipn) = flxdv(k,ipn)+(fld(k,ipn)-fld(k,ipnx))		&
        *sideln(edg,ipn)

	if (vrbos .and. ipn.eq.itest .and. mod(k,kgap).eq.1) then
         write (*,'(a,2i8,i3,a,2es12.4)') 'ocn dffusn_lev ipn,ipnx,k=',	&
          perm(ipn),perm(ipnx),k,' fld=',fld(k,ipn),fld(k,ipnx)
        end if

      end do				! loop through levels
     end if				! open ocean (non-coastal) edge
    end do				! loop through edges
   end if				! ocean point
  end do				! horizontal loop
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE (factor)
  do ipn = 1,nip			! horizontal loop
   if (mask(ipn) > 0 ) then
    do k = kfrst,klast			! loop through levels
     factor = -dffveldt*rarea(ipn)
     fld(k,ipn) = fld(k,ipn) + flxdv(k,ipn)*factor

      if (vrbos .and. ipn.eq.itest .and. mod(k,kgap).eq.1) then
       write (*,'(i8,i3,a,3es12.4)') perm(ipn),k,' flxdv,tdcy,field=',	&
        flxdv(k,ipn),flxdv(k,ipn)*factor,fld(k,ipn)
      end if

    end do				! loop through levels
   end if				! ocean point
  end do				! horizontal loop
!$OMP END PARALLEL DO
!SMS$PARALLEL END

  if (vrbos) then
   print *,'entering ocn dffusn_lev, dffveldt =',dffveldt
   do k=kfrst,klast
    if (kfrst.eq.klast .or. mod(k,kgap).eq.1) then
     write (string,'(a,i2)') ' k=',k
     call findmxmn2(fld,kdm,nip,k,'ocn dff_lev OUT'//trim(string),mask)
    end if
   end do
   call stencl(fld,kdm,100.,'(ocn dffusn_lev) field x 100 after smoothing')
  end if

  return
  end subroutine dffusn_lev
end module hycom_dffusn
