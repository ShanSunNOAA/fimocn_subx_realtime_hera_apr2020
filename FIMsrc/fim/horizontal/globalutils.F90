!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! globalutils: contains public routines to compute and print max and min values of distributed
! arrays. These distributed arrays can be dimensioned (k,ipn), or (ipn)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module globalutils
  use global_bounds, only: myrank, npes, ims, ime, ips, ipe
  implicit none

  private
  public :: maxinfo, mininfo

  interface maxinfo
     module procedure maxinfo_2d
     module procedure maxinfo_1d
  end interface

  interface mininfo
     module procedure mininfo_2d
     module procedure mininfo_1d
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! maxinfo_2d: Compute and print global max for distributed array dimensioned (k,nip)
!             NOTE: k-dim can be anything, but dope vector is 1-based even if the input
!             array is something else (e.g. 0:nz).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine maxinfo_2d (arr, descr)
!SMS$DISTRIBUTE (dh,2) BEGIN
    real, intent(in) :: arr(:,:)              ! Model array for which to find max info
!SMS$DISTRIBUTE END
    character(len=*), intent(in) :: descr     ! Description of arr

    real :: arrmaxk(ims:ime)   ! on-rank max value over k
    real :: arrmax             ! on-rank max value
    integer :: ksavek(ims:ime) ! saved value of k for max
    integer :: k, ipn          ! loop indices
    integer :: numk            ! size of k dimension
    integer :: ipnsave         ! horizontal location of max
    integer :: ksave           ! vertical location of max
    integer :: rank = 0        ! MPI rank responsible for max (0 in serial case)

! NOTE: Use of Fortran intrinsic maxloc() here could cut down the number of lines quite a bit.
! Good luck figuring out how to code around the idiocy of maxloc returning an *array*!

    if (size(arr,2) /= ime-ims+1) then
      write(6,*)'maxinfo_2d:', 'Problem when called with descr=', descr
      write(6,*)'Expected size(arr,2)=',ime-ims+1,' got ',size(arr,2)
      return
    end if

    numk = size(arr,1)
    ksavek(:) = -1
    ksave = -1
    ipnsave = -1

!$OMP PARALLEL DO PRIVATE(k) schedule (static)
    do ipn=ips,ipe
      arrmaxk(ipn) = -huge(arrmax)
      do k=1,numk
        if (arr(k,ipn) > arrmaxk(ipn)) then
          arrmaxk(ipn) = arr(k,ipn)
          ksavek(ipn) = k
        end if
      end do
    end do
!$OMP END PARALLEL DO

    arrmax = -huge(arrmax)
    do ipn=ips,ipe
      if (arrmaxk(ipn) > arrmax) then
        arrmax = arrmaxk(ipn)
        ksave = ksavek(ipn)
        ipnsave = ipn
      end if
    end do

! Now have per-rank max and ipn,k responsible.
! Next: Use MPI to expand this info across all compute tasks.
! Algorithm comes from pr_summary.c in GPTL library
    if (npes > 1) then
      call reduce ('max', ipnsave, ksave, arrmax, rank)
    end if
!
! Now rank 0 contains the global information: print it
!
    if (myrank == 0) then
      write(6,100) descr, arrmax, ksave, ipnsave, rank
100   format ('max ',a,'=',1p,e14.7,' at ONE-BASED k=', i3,' ipn=',i8,' MPI rank=', i6)
    end if

    return
  end subroutine maxinfo_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! maxinfo_2d: Compute and print global max for distributed array dimensioned (nip)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine maxinfo_1d (arr, descr)
!SMS$DISTRIBUTE (dh,1) BEGIN
    real, intent(in) :: arr(:)              ! Model array for which to find max info
!SMS$DISTRIBUTE END
    character(len=*), intent(in) :: descr   ! Description of arr

    real :: arrmax                  ! on-rank max value
    integer :: ipn                  ! loop index
    integer :: ipnsave              ! horizontal location of max
    integer :: ksave = 1            ! vertical location of max (needed by reduce)
    integer :: rank = 0             ! MPI rank responsible for max (0 in serial case)

    if (size(arr,1) /= ime-ims+1) then
      write(6,*)'maxinfo_1d:', 'Problem when called with descr=', descr
      write(6,*)'Expected size(arr,1)=',ime-ims+1,' got ',size(arr,1)
      return
    end if

! NOTE: Use of Fortran intrinsic maxloc() here could cut down the number of lines quite a bit.
! Good luck figuring out how to code around the idiocy of maxloc returning an *array*!

    ipnsave = -1

    arrmax = -huge(arrmax)
    do ipn=ips,ipe
      if (arr(ipn) > arrmax) then
        arrmax = arr(ipn)
        ipnsave = ipn
      end if
    end do

! Now have per-rank max and ipn responsible.
! Next: Use MPI to expand this info across all compute tasks.
! Algorithm comes from pr_summary.c in GPTL library
    if (npes > 1) then
      call reduce ('max', ipnsave, ksave, arrmax, rank)
    end if
!
! Now rank 0 contains the global information: print it
!
    if (myrank == 0) then
      write(6,100) descr, arrmax, ipnsave, rank
100   format ('max ',a,'=',1p,e14.7,' ipn=',i8,' MPI rank=', i6)
    end if

    return
  end subroutine maxinfo_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! mininfo_2d: Compute and print global min for distributed array dimensioned (k,nip)
!             NOTE: k-dim can be anything, but dope vector is 1-based even if the input
!             array is something else (e.g. 0:nz).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mininfo_2d (arr, descr)
!SMS$DISTRIBUTE (dh,2) BEGIN
    real, intent(in) :: arr(:,:)              ! Model array for which to find min info
!SMS$DISTRIBUTE END
    character(len=*), intent(in) :: descr     ! Description of arr

    real :: arrmink(ims:ime)   ! on-rank min value over k
    real :: arrmin             ! on-rank min value
    integer :: ksavek(ims:ime) ! saved value of k for min
    integer :: k, ipn          ! loop indices
    integer :: numk            ! size of k dimension
    integer :: ipnsave         ! horizontal location of min
    integer :: ksave           ! vertical location of min
    integer :: rank = 0        ! MPI rank responsible for min (0 in serial case)

    if (size(arr,2) /= ime-ims+1) then
      write(6,*)'maxinfo_2d:', 'Problem when called with descr=', descr
      write(6,*)'Expected size(arr,2)=',ime-ims+1,' got ',size(arr,2)
      return
    end if

! NOTE: Use of Fortran intrinsic minloc() here could cut down the number of lines quite a bit.
! Good luck figuring out how to code around the idiocy of minloc returning an *array*!

    numk = size(arr,1)
    ksavek(:) = -1
    ksave = -1
    ipnsave = -1

!$OMP PARALLEL DO PRIVATE(k) schedule (static)
    do ipn=ips,ipe
      arrmink(ipn) = +huge(arrmin)
      do k=1,numk
        if (arr(k,ipn) < arrmink(ipn)) then
          arrmink(ipn) = arr(k,ipn)
          ksavek(ipn) = k
        end if
      end do
    end do
!$OMP END PARALLEL DO

    arrmin = +huge(arrmin)
    do ipn=ips,ipe
      if (arrmink(ipn) < arrmin) then
        arrmin = arrmink(ipn)
        ksave = ksavek(ipn)
        ipnsave = ipn
      end if
    end do

! Now have per-rank min and ipn,k responsible.
! Next: Use MPI to expand this info across all compute tasks.
! Algorithm comes from pr_summary.c in GPTL library
    if (npes > 1) then
      call reduce ('min', ipnsave, ksave, arrmin, rank)
    end if
!
! Now rank 0 contains the global information: print it
!
    if (myrank == 0) then
      write(6,100) descr, arrmin, ksave, ipnsave, rank
100   format ('min ',a,'=',1p,e14.7,' at ONE-BASED k=', i3,' ipn=',i8,' MPI rank=', i6)
    end if

    return
  end subroutine mininfo_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! mininfo_1d: Compute and print global min for distributed array dimensioned (nip)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mininfo_1d (arr, descr)
!SMS$DISTRIBUTE (dh,1) BEGIN
    real, intent(in) :: arr(:)              ! Model array for which to find min info
!SMS$DISTRIBUTE END
    character(len=*), intent(in) :: descr   ! Description of arr

    real :: arrmin                  ! on-rank min value
    integer :: ipn                  ! loop index
    integer :: ipnsave              ! horizontal location of min
    integer :: ksave = 1            ! vertical location of min (needed by reduce)
    integer :: rank = 0             ! MPI rank responsible for min (0 in serial case)

    if (size(arr,1) /= ime-ims+1) then
      write(6,*)'mininfo_1d:', 'Problem when called with descr=', descr
      write(6,*)'Expected size(arr,1)=',ime-ims+1,' got ',size(arr,1)
      return
    end if

! NOTE: Use of Fortran intrinsic minloc() here could cut down the number of lines quite a bit.
! Good luck figuring out how to code around the idiocy of minloc returning an *array*!

    ipnsave = -1

    arrmin = +huge(arrmin)
    do ipn=ips,ipe
      if (arr(ipn) < arrmin) then
        arrmin = arr(ipn)
        ipnsave = ipn
      end if
    end do

! Now have per-rank min and ipn responsible.
! Next: Use MPI to expand this info across all compute tasks.
! Algorithm comes from pr_summary.c in GPTL library
    if (npes > 1) then
      call reduce ('min', ipnsave, ksave, arrmin, rank)
    end if
!
! Now rank 0 contains the global information: print it
!
    if (myrank == 0) then
      write(6,100) descr, arrmin, ipnsave, rank
100   format ('min ',a,'=',1p,e14.7,' ipn=',i8,' MPI rank=', i6)
    end if

    return
  end subroutine mininfo_1d

#ifndef SERIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reduce: Given task-local values for a distributed array and the ipn and k values associated,
!         reduce to global value, and retain associated ipn and k on rank 0.
!         NOTE: Eventually this belongs in SMS and can be removed from model code.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine reduce (flag, ipnsave, ksave, arrval, rank)

#ifndef SERIAL
#ifdef USEMPIMOD
    use mpi
#endif
#endif

#ifndef SERIAL
#ifndef USEMPIMOD
    include 'mpif.h'
#endif
#endif

    character(len=3),intent(in) :: flag      ! whether we are doing max or min
    integer                     :: comm      ! communicator
    integer                     :: iarr(2)   ! put 2 integers together to minimize number of MPI msgs
    integer                     :: ierr      ! error return from MPI
    integer                     :: incr      ! rank increment in MPI reduction algorithm
    integer                     :: ipnsave_p ! ipnsave received from neighbor
    integer                     :: ksave_p   ! ksave received from neighbor
    integer                     :: p         ! rank to receive from
    integer                     :: sendto    ! rank to send to
    integer                     :: twoincr   ! 2*incr
    integer,intent(inout)       :: ipnsave   ! on input: on-rank ipn. on output: global ipn
    integer,intent(inout)       :: ksave     ! on input: on-rank k. on output: k for global ipn
    integer,intent(out)         :: rank      ! output rank associate with arrval
    integer,parameter           :: tag = 477 ! tag for MPI communication
    logical                     :: dorecv    ! whether to recv data this step
    logical                     :: dosend    ! whether to send data this step
    real                        :: arrval_p  ! val value received from neighbor
    real,intent(inout)          :: arrval    ! on input: on-rank max or min value, on ouput: domain max or min value on rank 0

#ifndef SERIAL
    integer :: status(MPI_STATUS_SIZE) ! status returned from MPI
#endif

!sms$get_communicator(comm)
    incr = 1
    rank = 0
    do while (incr < npes)
      twoincr = 2*incr
      sendto = myrank - incr
      p = myrank + incr
!
! The .and. part of the next 2 stmts prevents sending to or receiving from
! outside communicator bounds when npes is not a power of 2
!
      dorecv = mod((myrank + twoincr), twoincr) == 0 .and. p < npes
      dosend = mod((myrank + incr), twoincr) == 0 .and. sendto > -1
      if (dosend) then
        if (dorecv) then
!sms$ignore begin
          write(6,*) 'WARNING: myrank=',myrank,': dosend and dorecv both true: possible hang?'
          call flush(6)
!sms$ignore end
        end if
! Put the 2 integers together to minimize number of MPI msgs
        iarr(1) = ipnsave
        iarr(2) = ksave
#ifndef SERIAL
        call mpi_send (iarr, 2, MPI_INTEGER, sendto, tag, comm, ierr)
        call mpi_send (arrval, 1, MPI_REAL, sendto, tag, comm, ierr)
#endif
      end if

      if (dorecv) then
        if (dosend) then
!sms$ignore begin
          write(6,*) 'WARNING: myrank=',myrank,': dosend and dorecv both true: possible hang?'
          call flush(6)
!sms$ignore end
        end if

#ifndef SERIAL
        call mpi_recv (iarr, 2, MPI_INTEGER, p, tag, comm, status, ierr)
        call mpi_recv (arrval_p, 1, MPI_REAL, p, tag, comm, status, ierr)
#endif
        ipnsave_p = iarr(1)
        ksave_p = iarr(2)
!
! Merge stats just received with those currently on this MPI task.
! NOTE: ">" test in one case vs. "<" in the other (as opposed to >= and <=) guarantees the same
! result regardless of number of ranks: When there are ties, always favor the lower rank.
! The above ASSUMES that Fortran intrinsics minloc and maxloc behave similarly, which has been
! verified with ifort, pgf90 and gfortran.
!
        if (flag == 'max') then
          if (arrval_p > arrval) then
            arrval = arrval_p
            ipnsave = ipnsave_p
            ksave = ksave_p
            rank = p
          end if
        else if (flag == 'min') then
          if (arrval_p < arrval) then
            arrval = arrval_p
            ipnsave = ipnsave_p
            ksave = ksave_p
            rank = p
          end if
        end if
      end if
      incr = twoincr
    end do
    return
  end subroutine reduce
#endif /* parallel */
end module globalutils
