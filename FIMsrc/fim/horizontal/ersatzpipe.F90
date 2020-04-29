module ersatzpipe

!***********************************************************************
! --- serialize distributed array, map to global index space, and
! --- output as binary file for diagnosing distributed memory problems.

! --- approach:
! --- insert calls to 'pipe' in strategic places in the source code, compile, 
! --- and run with different processor counts and/or 'curve' parameters.
! --- 'diff' the output files in the two runtime directories. iterate.
!***********************************************************************

use module_constants,only: inv_perm	! inv_perm(global_index)=local_index
use module_control  ,only: nip, numphr
use hycom_constants ,only: wet		! needed for zeroing out land points
contains
   subroutine pipe2(its,field,kdm,what)

! --- use this entry to output arrays whose  S E C O N D  index is distributed
! --- (set kdm=1 when calling with a 1-D distributed array)

   use module_control  ,only: nip
   use hycom_constants,only: wet	! needed for zeroing out land points
   implicit none
   integer  ,intent(IN) :: its,kdm
   character,intent(IN) :: what*(*)
!SMS$DISTRIBUTE (dh,2) BEGIN
   real     ,intent(IN) :: field(kdm,nip)
   real             masked_field(kdm,nip)
!SMS$DISTRIBUTE END
   real serial_field(kdm,nip)
   integer i
   real valmin,valmax
   character flnm*40

   if (its.gt.numphr*240+1) return	! reached 10 day limit in pipe
   masked_field(:,:)=field(:,:)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- optional: zero out land points

!!SMS$PARALLEL (dh,i) BEGIN
!   do i=1,nip
!    if (wet(i) == 0 ) masked_field(:,i)=0.
!   end do
!!SMS$PARALLEL END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
   write (flnm,'(a,"_",i4.4,".bin")') what,its
!SMS$SERIAL (<masked_field,inv_perm,IN> : default=ignore)  BEGIN
   do i=1,nip
    serial_field(:,i)=masked_field(:,inv_perm(i))
   end do
   open (unit=11,file=flnm,form='unformatted')
   write (11) serial_field
   close (11)
!SMS$SERIAL END

!SMS$PARALLEL (dh,i) BEGIN
   valmin=minval(masked_field(1:kdm,1:nip))
   valmax=maxval(masked_field(1:kdm,1:nip))
!SMS$REDUCE(valmin,MIN)
!SMS$REDUCE(valmax,MAX)
!SMS$PARALLEL END
   print '(3a,2es15.7)','created file ',flnm,'  min,max =',valmin,valmax

   return
   end subroutine pipe2


   subroutine pipe3(its,field,kdm,what)

! --- use this entry to output arrays whose  T H I R D  index is distributed

   implicit none
   integer  ,intent(IN) :: its,kdm
   character,intent(IN) :: what*(*)
!SMS$DISTRIBUTE (dh,3) BEGIN
   real     ,intent(IN) :: field(kdm,6,nip)
   real             masked_field(kdm,6,nip)
!SMS$DISTRIBUTE END
   real serial_field(kdm,6,nip)
   integer i
   real valmin,valmax
   character flnm*40

   masked_field(:,:,:)=field(:,:,:)
   
   write (flnm,'(a,"_",i4.4,".bin")') what,its
!SMS$SERIAL (<masked_field,inv_perm,IN> : default=ignore)  BEGIN
   do i=1,nip
    serial_field(:,:,i)=masked_field(:,:,inv_perm(i))
   end do
   open (unit=11,file=flnm,form='unformatted')
   write (11) serial_field
   close (11)
!SMS$SERIAL END

!SMS$PARALLEL (dh,i) BEGIN
   valmin=minval(masked_field(1:kdm,:,1:nip))
   valmax=maxval(masked_field(1:kdm,:,1:nip))
!SMS$REDUCE(valmin,MIN)
!SMS$REDUCE(valmax,MAX)
!SMS$PARALLEL END
   print '(3a,2es15.7)','created file ',flnm,'  min,max =',valmin,valmax

   return
   end subroutine pipe3


   subroutine maskpipe(its,mask,field,kdm,what)

! --- use this entry to mask out portions of array before outputting

   use module_control,only: nip
   implicit none
   integer  ,intent(IN) :: its,kdm
   character,intent(IN) :: what*(*)
!SMS$DISTRIBUTE (dh,1) BEGIN
   logical  ,intent(IN) :: mask(nip)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
   real     ,intent(IN) :: field(kdm,nip)
   real             masked_field(kdm,nip)
!SMS$DISTRIBUTE END
   real serial_field(kdm,nip)
   integer i,k
   real valmin,valmax
   character flnm*20

   print *,'entering maskpipe'
   call flush(6)

!SMS$PARALLEL (dh,i) BEGIN
   do i=1,nip
    if (mask(i)) then
     masked_field(:,i)=field(:,i)
    else
     masked_field(:,i)=0.
    end if
   end do
!SMS$PARALLEL END
   
   print *,'what,its =',what,its
   call flush(6)
   write (flnm,'(a,"_",i4.4,".bin")') what,its
   print *,'about to write file ',flnm
   call flush(6)
!SMS$SERIAL (<masked_field,inv_perm,IN> : default=ignore)  BEGIN
   do i=1,nip
    serial_field(:,i)=masked_field(:,inv_perm(i))
   end do
   open (unit=11,file=flnm,form='unformatted')
   write (11) serial_field
   close (11)
!SMS$SERIAL END

!SMS$PARALLEL (dh,i) BEGIN
   valmin=minval(masked_field(1:kdm,1:nip))
   valmax=maxval(masked_field(1:kdm,1:nip))
!SMS$REDUCE(valmin,MIN)
!SMS$REDUCE(valmax,MAX)
!SMS$PARALLEL END
   print '(3a,2es15.7)','created file ',flnm,'  min,max =',valmin,valmax

   return
   end subroutine maskpipe
end module ersatzpipe
