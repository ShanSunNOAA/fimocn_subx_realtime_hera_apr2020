module mdul_smoofld
use module_constants,only: prox,nprox
use module_control  ,only: nip
use findmaxmin2
public :: smolight,smoheavy
contains

!*********************************************************************
!     smolight
!       L i g h t  smoother for data on icos grid
!	R.Bleck				July 2014
!*********************************************************************

! --- Method: Each icos point forms the center of a cluster of 6 (occasionally
! --- 5) triangles formed by triplets of icos points.  Averaging the 3 vertex
! --- values of each triangle yields an approximate areal average of 'field'
! --- in that triangle. This routine replaces the 'field' value at each grid
! --- point by the average of the 6 (or 5) surrounding triangular averages.

  subroutine smolight(field,nvl,npass,what)

  integer, intent(IN) :: nvl		! number of vertical levels
  integer, intent(IN) :: npass		! number of smoothing passes
  character,intent(IN),optional :: what*(*)
!SMS$DISTRIBUTE (dh,2) BEGIN
  real,intent(INOUT) :: field(nvl,nip)	! field to be smoothed
  real    :: work(nvl,nip)
!SMS$DISTRIBUTE END
  real    :: sum(nvl)
  integer :: ipn,edg,iter
  character :: string*13

  do iter=1,npass
   if (present(what)) print *,'smoothing ',what,' field, pass',iter
   work(:,:)=field(:,:)

#ifdef DEBUGPRINT
   do k=1,nvl
    write (string,'(a,i2)') 'smolight k=',k
    call findmxmn2(work,nvl,nip,k,string)
   end do
#endif

!SMS$EXCHANGE(work)
!SMS$PARALLEL (dh,ipn) BEGIN
   do ipn=1,nip
    do k=1,nvl
     sum(k)=work(k,ipn)*.5*nprox(ipn)
     do edg=1,nprox(ipn)
      sum(k)=sum(k)+work(k,prox(edg,ipn))
     end do
     field(k,ipn)=sum(k)/(1.5*nprox(ipn))
    end do
   end do
!SMS$PARALLEL END

  end do		! iter
  return
  end subroutine smolight

!*********************************************************************
!     smoheavy
!       H e a v y   smoother for data on icos grid
!	R.Bleck				June 2016
!*********************************************************************

! --- Method: Replace each icos value by the average of its neighbors.
! --- Repeat once, resulting in a 21-pt stencil.

  subroutine smoheavy(field,nvl,npass,what)
  integer, intent(IN) :: nvl		! number of vertical levels
  integer, intent(IN) :: npass		! number of smoothing passes
  character,intent(IN),optional :: what*(*)
!SMS$DISTRIBUTE (dh,2) BEGIN
  real,intent(INOUT) :: field(nvl,nip)	! field to be smoothed
  real    :: work(nvl,nip)
!SMS$DISTRIBUTE END
  real    :: sum(nvl)
  integer :: ipn,edg,iter
  character :: string*13

  do iter=1,2*npass
   work(:,:)=field(:,:)
   if (present(what)) print *,'smooothing ',what,' field'

#ifdef DEBUGPRINT
   do k=1,nvl
    write (string,'(a,i2)') 'smoheavy k=',k
    call findmxmn2(work,nvl,nip,k,string)
   end do
#endif

!SMS$EXCHANGE(work)
!SMS$PARALLEL (dh,ipn) BEGIN
   do ipn=1,nip
    do k=1,nvl
     sum(k)=0.
     do edg=1,nprox(ipn)
      sum(k)=sum(k)+work(k,prox(edg,ipn))
     end do
     field(k,ipn)=sum(k)/nprox(ipn)
    end do
   end do
!SMS$PARALLEL END

  end do		! iter
  return
  end subroutine smoheavy
end module mdul_smoofld
