! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! this package contains routines for finding maxima/minima in
! distributed multidimensional arrays. in all cases the search extends
! over the  s e c o n d  array index associated with the distributed
! dimension. (obvious exception: 1-D arrays).
! additionally, there are routines for computing the mean value in the
! distributed dimension.
! if 'mask' is supplied, only points where mask > 0 will be searched.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


module findmaxmin1
use module_constants,only: perm
contains
   subroutine findmxmn1(array,dim,what,mask)

! --- find location of maximum and minimum in 1-D distributed array

   implicit none
   integer,intent(IN)   :: dim			! array dimension
   character,intent(IN) :: what*(*)
!sms$distribute (dh,1) begin
   real,intent(IN)      :: array(dim)
   integer,intent(IN),optional :: mask(dim)
!sms$distribute end
   integer n,mini,maxi,i,scrutn
   real locmin,locmax,globmn,globmx
   logical foundmax,foundmin
 

!SMS$PARALLEL (dh,n) BEGIN
   mini=-999
   maxi=-999
   locmin= 1.e33
   locmax=-1.e33
   do n=1,dim
    scrutn=1
    if (present(mask)) scrutn=mask(n)
    if (scrutn > 0) then
     if (array(n).lt.locmin) then
      locmin=array(n)
      mini=n
     end if
     if (array(n).gt.locmax) then
      locmax=array(n)
      maxi=n
     end if
    end if
   end do

   globmn=locmin
!SMS$REDUCE(globmn,MIN)
   if (locmin.eq.globmn) then
    i=mini
    call GetIpnGlobal(i,mini,foundmin)
    mini=perm(mini)
   else
    mini=-1
   end if

   globmx=locmax
!SMS$REDUCE(globmx,MAX)
   if (locmax.eq.globmx) then
    i=maxi
    call GetIpnGlobal(i,maxi,foundmax)
    maxi=perm(maxi)
   else
    maxi=-1
   end if
!SMS$REDUCE(mini,maxi,MAX)
!SMS$PARALLEL END

   write (*,'(a,2(a,es13.6,a,i6))') trim(what),				&
       ':  min=',globmn,' at ',mini,					&
        '  max=',globmx,' at ',maxi

   return
   end subroutine findmxmn1
end module findmaxmin1


module findmaxmin2
use module_constants,only: perm
contains
   subroutine findmxmn2(array,dim1,dim2,n1,what,mask)

! --- find location of max and min in 2nd dimension of 2-D distributed array

   implicit none
   integer,intent(IN)   :: dim1,dim2		! array dimensions
   integer,intent(IN)   :: n1			! array index assoc. with dim1
   character,intent(IN) :: what*(*)
!sms$distribute (dh,1) begin
   integer,intent(IN),optional :: mask(dim2)
!sms$distribute end
!sms$distribute (dh,2) begin
   real,intent(IN)      :: array(dim1,dim2)
!sms$distribute end
   integer n2,mini,maxi,i,scrutn
   real locmin,locmax,globmn,globmx
   logical foundmax,foundmin
 
!SMS$PARALLEL (dh,n2) BEGIN
   mini=-999
   maxi=-999
   locmin= 1.e33
   locmax=-1.e33
   do n2=1,dim2
    scrutn=1
    if (present(mask)) scrutn=mask(n2)
    if (scrutn > 0 ) then
     if (array(n1,n2).lt.locmin) then
      locmin=array(n1,n2)
      mini=n2
     end if
     if (array(n1,n2).gt.locmax) then
      locmax=array(n1,n2)
      maxi=n2
     end if
    end if
   end do

   globmn=locmin
!SMS$REDUCE(globmn,MIN)
   if (locmin.eq.globmn) then
    i=mini
    call GetIpnGlobal(i,mini,foundmin)
    mini=perm(mini)
   else
    mini=-1
   end if

   globmx=locmax
!SMS$REDUCE(globmx,MAX)
   if (locmax.eq.globmx) then
    i=maxi
    call GetIpnGlobal(i,maxi,foundmax)
    maxi=perm(maxi)
   else
    maxi=-1
   end if
!SMS$REDUCE(mini,maxi,MAX)
!SMS$PARALLEL END

   write (*,'(a,2(a,es13.6,a,i6))') trim(what),				&
       ':  min=',globmn,' at ',mini,					&
        '  max=',globmx,' at ',maxi

   return
   end subroutine findmxmn2
end module findmaxmin2


module findmaxmin3
use module_constants,only: perm
contains
   subroutine findmxmn3(array,dim1,dim2,dim3,n1,n3,what,mask)

! --- find location of max and min in 2nd dimension of 3-D distributed array

   implicit none
   integer,intent(IN)   :: dim1,dim2,dim3	! array dimensions
   integer,intent(IN)   :: n1,n3		! indices assoc. with dim1/3
   character,intent(IN) :: what*(*)
!sms$distribute (dh,1) begin
   integer,intent(IN),optional :: mask(dim2)
!sms$distribute end
!sms$distribute (dh,2) begin
   real,intent(IN)      :: array(dim1,dim2,dim3)
!sms$distribute end
   integer n2,mini,maxi,scrutn
   real locmin,locmax,globmn,globmx
   logical foundmax,foundmin
 
!SMS$PARALLEL (dh,n2) BEGIN
   mini=-999
   maxi=-999
   locmin= 1.e33
   locmax=-1.e33
   do n2=1,dim2
    scrutn=1
    if (present(mask)) scrutn=mask(n2)
    if (scrutn > 0 ) then
     if (array(n1,n2,n3).lt.locmin) then
       locmin=array(n1,n2,n3)
       mini=n2
     end if
     if (array(n1,n2,n3).gt.locmax) then
      locmax=array(n1,n2,n3)
      maxi=n2
     end if
    end if
   end do

   globmn=locmin
!SMS$REDUCE(globmn,MIN)
   if (locmin.eq.globmn) then
    call GetIpnGlobal(mini,mini,foundmin)
    mini=perm(mini)
   else
    mini=-1
   end if

   globmx=locmax
!SMS$REDUCE(globmx,MAX)
   if (locmax.eq.globmx) then
    call GetIpnGlobal(maxi,maxi,foundmax)
    maxi=perm(maxi)
   else
    maxi=-1
   end if
!SMS$REDUCE(mini,maxi,MAX)
!SMS$PARALLEL END

   write (*,'(a,2(a,es13.6,a,i6))') trim(what),				&
       ':  min=',globmn,' at ',mini,					&
        '  max=',globmx,' at ',maxi

   return
   end subroutine findmxmn3
end module findmaxmin3


module findmaxmin4
use module_constants,only: perm
contains
   subroutine findmxmn4(array,dim1,dim2,dim3,dim4,n1,n3,n4,what,mask)

! --- find location of max and min in 2nd dimension of 4-D distributed array

   implicit none
   integer,intent(IN)   :: dim1,dim2,dim3,dim4	! array dimensions
   integer,intent(IN)   :: n1,n3,n4		! indices assoc.with dim1/3/4
   character,intent(IN) :: what*(*)
!sms$distribute (dh,1) begin
   integer,intent(IN),optional :: mask(dim2)
!sms$distribute end
!sms$distribute (dh,2) begin
   real,intent(IN)      :: array(dim1,dim2,dim3,dim4)
!sms$distribute end
   integer n2,mini,maxi,scrutn
   real locmin,locmax,globmn,globmx
   logical foundmax,foundmin
 
!SMS$PARALLEL (dh,n2) BEGIN
   mini=-999
   maxi=-999
   locmin= 1.e33
   locmax=-1.e33
   do n2=1,dim2
    scrutn=1
    if (present(mask)) scrutn=mask(n2)
    if (scrutn > 0 ) then
    if (array(n1,n2,n3,n4).lt.locmin) then
      locmin=array(n1,n2,n3,n4)
      mini=n2
     end if
     if (array(n1,n2,n3,n4).gt.locmax) then
      locmax=array(n1,n2,n3,n4)
      maxi=n2
     end if
    end if
   end do

   globmn=locmin
!SMS$REDUCE(globmn,MIN)
   if (locmin.eq.globmn) then
    call GetIpnGlobal(mini,mini,foundmin)
    mini=perm(mini)
   else
    mini=-1
   end if

   globmx=locmax
!SMS$REDUCE(globmx,MAX)
   if (locmax.eq.globmx) then
    call GetIpnGlobal(maxi,maxi,foundmax)
    maxi=perm(maxi)
   else
    maxi=-1
   end if
!SMS$REDUCE(mini,maxi,MAX)
!SMS$PARALLEL END

   write (*,'(a,2(a,es13.6,a,i6))') trim(what),				&
       ':  min=',globmn,' at ',mini,					&
        '  max=',globmx,' at ',maxi

   return
   end subroutine findmxmn4
end module findmaxmin4


   subroutine findmean1(array,dim,what,mask)

! --- find mean value in 1-D distributed array

   implicit none
   integer,intent(IN)   :: dim			! array dimensions
   character,intent(IN) :: what*(*)
!sms$distribute (dh,1) begin
   real,intent(IN)      :: array(dim)
   logical,intent(IN),optional :: mask(dim)
!sms$distribute end
   integer n
   logical scrutn
   real globsum,count
   real*8 locsum
 
!SMS$PARALLEL (dh,n) BEGIN
   locsum=0.
   count=0.
   do n=1,dim
    scrutn=.false.
    if (present(mask)) scrutn=mask(n)
    if (scrutn) then
     locsum=locsum+array(n)
     count=count+1.
    end if
   end do
   globsum=locsum/count
!SMS$REDUCE(globsum,SUM)
!SMS$PARALLEL END

   write (*,'(2a,es15.7)') trim(what),'  mean value:',globsum
   return
   end subroutine findmean1


   subroutine findmean2(array,dim1,dim2,n1,what,mask)

! --- find mean value in 2nd (distributed) dimension of 2-D array

   implicit none
   integer,intent(IN)   :: dim1,dim2		! array dimensions
   integer,intent(IN)   :: n1			! array index assoc. with dim1
   character,intent(IN) :: what*(*)
!sms$distribute (dh,1) begin
   integer,intent(IN),optional :: mask(dim2)
!sms$distribute end
!sms$distribute (dh,2) begin
   real,intent(IN)      :: array(dim1,dim2)
!sms$distribute end
   integer n2, scrutn
   real globsum,count
   real*8 locsum
 
!SMS$PARALLEL (dh,n2) BEGIN
   locsum=0.
   count=0.
   do n2=1,dim2
    scrutn=1
    if (present(mask)) scrutn=mask(n2)
    if (scrutn > 0 ) then
     locsum=locsum+array(n1,n2)
     count=count+1
    end if
   end do
   globsum=locsum/count
!SMS$REDUCE(globsum,SUM)
!SMS$PARALLEL END

   write (*,'(2a,es15.7)') trim(what),'  mean value:',globsum
   return
   end subroutine findmean2


   subroutine findmean3(array,dim1,dim2,dim3,n1,n3,what,mask)

! --- find mean value in 2nd (distributed) dimension of 3-D array

   implicit none
   integer,intent(IN)   :: dim1,dim2,dim3	! array dimensions
   integer,intent(IN)   :: n1,n3		! indices assoc. with dim1/3
   character,intent(IN) :: what*(*)
!sms$distribute (dh,1) begin
   integer,intent(IN),optional :: mask(dim2)
!sms$distribute end
!sms$distribute (dh,2) begin
   real,intent(IN)      :: array(dim1,dim2,dim3)
!sms$distribute end
   integer n2, scrutn
   real globsum,count
   real*8 locsum
 
!SMS$PARALLEL (dh,n2) BEGIN
   locsum=0.
   count=0.
   do n2=1,dim2
    scrutn=1
    if (present(mask)) scrutn=mask(n2)
    if (scrutn > 0 ) then
     locsum=locsum+array(n1,n2,n3)
     count=count+1.
    end if
   end do
   globsum=locsum/count
!SMS$REDUCE(globsum,SUM)
!SMS$PARALLEL END

   write (*,'(2a,es15.7)') trim(what),'  mean value:',globsum
   return
   end subroutine findmean3


   subroutine findmean4(array,dim1,dim2,dim3,dim4,n1,n3,n4,what,mask)

! --- find mean value in 2nd (distributed) dimension of 4-D array

   implicit none
   integer,intent(IN)   :: dim1,dim2,dim3,dim4	! array dimensions
   integer,intent(IN)   :: n1,n3,n4		! indices assoc.with dim1/3/4
   character,intent(IN) :: what*(*)
!sms$distribute (dh,1) begin
   integer,intent(IN),optional :: mask(dim2)
!sms$distribute end
!sms$distribute (dh,2) begin
   real,intent(IN)      :: array(dim1,dim2,dim3,dim4)
!sms$distribute end
   integer n2, scrutn
   real globsum,count
   real*8 locsum
 
!SMS$PARALLEL (dh,n2) BEGIN
   locsum=0.
   count=0.
   do n2=1,dim2
    scrutn=1
    if (present(mask)) scrutn=mask(n2)
    if (scrutn > 0 ) then
     locsum=locsum+array(n1,n2,n3,n4)
     count=count+1.
    end if
   end do
   globsum=locsum/count
!SMS$REDUCE(globsum,SUM)
!SMS$PARALLEL END

   write (*,'(2a,es15.7)') trim(what),'  mean value:',globsum
   return
   end subroutine findmean4


module edgmaxmin
use module_constants,only: nedge,permedge,perm
contains
   subroutine edgmxmn(array,kdm,what,mask)

! --- find min/max in an edge array, i.e. an array dimensioned (kdm,npp,nip)

   use module_control  ,only: nip,npp
   implicit none
   integer  ,intent(IN) :: kdm
!sms$distribute (dh,1) begin
   integer,intent(IN),optional :: mask(nip)
!sms$distribute end
!sms$distribute (dh,3) begin
   real     ,intent(IN) :: array(kdm,npp,nip)
!sms$distribute end
   character,intent(IN) :: what*(*)
   integer ipn,edgcount,edg,k,minedg,maxedg,mink,maxk,minipn,maxipn,	&
           edgmin,edgmax,kmin,kmax,ipnmin,ipnmax,scrutn
   real locmin,locmax,globmin,globmax,valmin,valmax
   character*20 :: empt = '                    '

! !SMS$PARALLEL (dh,ipn) BEGIN
!    valmin=minval(array)
!    valmax=maxval(array)
! !SMS$REDUCE(valmin,MIN)
! !SMS$REDUCE(valmax,MAX)
! !SMS$PARALLEL END
!    print '(2a,2es15.7)',what,' min,max:',valmin,valmax

!SMS$PARALLEL (dh,ipn) BEGIN
! !SMS$HALO_COMP(<1,1>) BEGIN
   locmin= 1.e33
   locmax=-1.e33
   do ipn = 1,nip			! horizontal loop
    scrutn=1
    if (present(mask)) scrutn=mask(ipn)
    if (scrutn > 0 ) then
     do edgcount = 1,nedge(ipn)		! loop through edges
      edg = permedge(edgcount,ipn)
      do k=1,kdm
       if (array(k,edg,ipn).gt.locmax) then
        locmax=array(k,edg,ipn)
        maxedg=edg
        maxk=k
        maxipn=ipn
       end if
       if (array(k,edg,ipn).lt.locmin) then
        locmin=array(k,edg,ipn)
        minedg=edg
        mink=k
        minipn=ipn
       end if
      end do
     end do
    end if
   end do
   globmin=locmin
   globmax=locmax
!SMS$REDUCE(globmin,MIN)
!SMS$REDUCE(globmax,MAX)
   if (globmax.eq.locmax) then
    ipnmax=perm(maxipn)
    kmax=maxk
    edgmax=maxedg
   else
    ipnmax=-1
    kmax=-1
    edgmax=-1
   end if
   if (globmin.eq.locmin) then
    ipnmin=perm(minipn)
    kmin=mink
    edgmin=minedg
   else
    ipnmin=-1
    kmin=-1
    edgmin=-1
   end if
!SMS$REDUCE(ipnmin,ipnmax,kmin,kmax,edgmin,edgmax,MAX)
! !SMS$HALO_COMP END
!SMS$PARALLEL END

   k=len_trim(what)
   print 100,what(1:k),' min',globmin,' is at ipn,edg,k=',ipnmin,edgmin,kmin
   print 100,empt(1:k),' max',globmax,' is at ipn,edg,k=',ipnmax,edgmax,kmax
100 format (2a,es15.7,a,i7,2i3)
   return
   end subroutine edgmxmn
end module edgmaxmin

module findmaxmin8
use module_constants,only: perm
contains
   subroutine findmxmn8(array,dim1,dim2,n1,what,mask)

! --- find location of max and min in 2nd dimension of 2-D distributed array

   implicit none
   integer,intent(IN)   :: dim1,dim2            ! array dimensions
   integer,intent(IN)   :: n1                   ! array index assoc. with dim1
   character,intent(IN) :: what*(*)
!sms$distribute (dh,1) begin
   logical,intent(IN),optional :: mask(dim2)
!sms$distribute end
!sms$distribute (dh,2) begin
   real*8,intent(IN)      :: array(dim1,dim2)
!sms$distribute end
   integer n2,mini,maxi,i
   real locmin,locmax,globmn,globmx
   logical scrutn,foundmax,foundmin

!SMS$PARALLEL (dh,n2) BEGIN
   mini=-999
   maxi=-999
   locmin= 1.e33
   locmax=-1.e33
   do n2=1,dim2
    scrutn=.true.
    if (present(mask)) scrutn=mask(n2)
    if (scrutn) then
     if (array(n1,n2).lt.locmin) then
      locmin=array(n1,n2)
      mini=n2
     end if
     if (array(n1,n2).gt.locmax) then
      locmax=array(n1,n2)
      maxi=n2
     end if
    end if
   end do

   globmn=locmin
!SMS$REDUCE(globmn,MIN)
   if (locmin.eq.globmn) then
    i=mini
    call GetIpnGlobal(i,mini,foundmin)
    mini=perm(mini)
   else
    mini=-1
   end if

   globmx=locmax
!SMS$REDUCE(globmx,MAX)
   if (locmax.eq.globmx) then
    i=maxi
    call GetIpnGlobal(i,maxi,foundmax)
    maxi=perm(maxi)
   else
    maxi=-1
   end if
!SMS$REDUCE(mini,maxi,MAX)
!SMS$PARALLEL END

   write (*,'(a,2(a,es13.6,a,i6))') trim(what),                         &
       ':  min=',globmn,' at ',mini,                                    &
        '  max=',globmx,' at ',maxi

   return
   end subroutine findmxmn8
end module findmaxmin8

