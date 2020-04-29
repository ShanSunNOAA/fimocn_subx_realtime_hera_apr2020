program spe2nc
  use netcdf
  use sigio_module
  use sptsnv, only: init_sp, end_sp, sptez_mod, sptezmv_mod
  use nc, only: initnc, handle_err, ncid, latid, lonid, vertid, &
                hsid, psid, tvid, qid, uid, vid, ozid, clwid

  implicit none

  integer, parameter :: lusig = 82   ! Big endian!!!
  integer, parameter :: idrt = 4     ! This value flags the NCEP grid as gaussian
  integer :: imax, jmax, nvp         ! gaussian grid dimensions
  integer :: maxwv                   ! max number of waves
  integer :: i, j, k                 ! gaussian grid indices
  integer :: iret                    ! Return code from sigio routines
  real :: r2d                        ! Convert radians to degrees
  real, allocatable :: arr(:,:)      ! 2D slab of information on gaussian grid
  real, allocatable :: lat(:,:)      ! latitude information
  real, allocatable :: lon(:,:)      ! longitude information
  real, allocatable :: lev(:)        ! vertical level information for netcdf file
  type(sigio_head) :: head           ! NCEP spectral data header
  type(sigio_data) :: data           ! NCEP spectral data
  character(len=120) :: sanlfile    = 'sanlfile not set'     ! spectral coefficient file
  character(len=80) :: gfsltln_file = 'gfsltln_file not set' ! lats/lons
  character(len=80) :: outfile      = 'outfile.nc'           ! output netcdf file
  namelist /spe2nc_nl/ sanlfile, gfsltln_file, outfile

  read (5, nml=spe2nc_nl, iostat=iret)
  if (iret /= 0) then
    write(6,*)'Error reading spec2nc_nl namelist input'
    stop 999
  end if

  r2d = 180./acos(-1.)
  write(6,*)'Reading spectal data...'
  call sigio_srohdc (lusig, sanlfile, head, data, iret)
  if (iret /= 0) then
    write(6,*)'Error from sigio_srohdc'
    stop 999
  end if

  nvp  = head%levs
  imax = head%lonb
  jmax = head%latb
  maxwv = head%jcap
  write(6,*)'nvp,imax,jmax=',nvp,imax,jmax

  call init_sp (0, maxwv, 4, imax, jmax) ! 4 means gaussian grid

  allocate (arr(imax,jmax))
  allocate (lat(imax,jmax))
  allocate (lon(imax,jmax))
  allocate (lev(nvp))

  open (file=gfsltln_file, unit=1, status='old', form='unformatted')
  read (1) lat(:,:), lon(:,:)
  close (1)

! Verify lon and lat can be 1-d
  lat(:,:) = lat(:,:)*r2d
  lon(:,:) = lon(:,:)*r2d

  write(6,*)'Checking lat/lon values...'
  do j=1,jmax
    do i=2,imax
      if (lat(i,j) /= lat(i-1,j)) then
        write(6,*)'bad lat=',lat(i,j),' at i,j=',i,j,'. lat(i-1,j)=',lat(i-1,j)
      end if
    end do
  end do

  do i=1,imax
    do j=2,jmax
      if (lon(i,j) /= lon(i,j-1)) then
        write(6,*)'bad lon=',lon(i,j),' at i,j=',i,j,'. lon(i,j-1)=',lon(i,j-1)
      end if
    end do
  end do

  write(6,*)'Initializing the netcdf file...'
! Initialize netcdf file, write metadata
  call initnc (sanlfile, outfile, imax, jmax, nvp)

! Write dimension variables
  call handle_err (nf90_put_var (ncid, lonid, lon(:,jmax/2)))
  call handle_err (nf90_put_var (ncid, latid, lat(imax/2,:)))
  do k=1,nvp
    lev(k) = k
  end do
  call handle_err (nf90_put_var (ncid, vertid, lev))

! xt1
  write(6,*)'Transforming hs...'
  call sptez_mod (0, maxwv, idrt, imax, jmax, data%hs, arr, 1)
  call handle_err (nf90_put_var (ncid, hsid, arr))  ! sfc hgt

  write(6,*)'Transforming ps...'
  call sptez_mod (0, maxwv, idrt, imax, jmax, data%ps, arr, 1)
  arr(:,:) = exp(arr(:,:))*1.e3
  call handle_err (nf90_put_var (ncid, psid, arr))  ! ps
  deallocate (arr)

!$OMP PARALLEL DO
  do k=1,nvp
    call xt2 (imax,jmax,k)
  end do
!$OMP END PARALLEL DO


!$OMP PARALLEL DO
  do k=1,nvp
    call xt4 (imax,jmax,k)
  end do
!$OMP END PARALLEL DO


!$OMP PARALLEL DO
  do k=1,nvp
    call xt5 (imax,jmax,k)
  end do
!$OMP END PARALLEL DO


!$OMP PARALLEL DO
  do k=1,nvp
    call xt6 (imax,jmax,k)
  end do
!$OMP END PARALLEL DO

  write(6,*)'Done'
  call handle_err (nf90_close (ncid))

  stop

  contains

    subroutine xt2(imax,jmax,k)
      integer, intent(in) :: imax, jmax, k

      integer :: start(3)                ! netcdf slice starting locations
      integer :: kount(3)                ! netcdf slice count information
      real :: arr(imax,jmax)

      start(:) = (/1,1,k/)
      kount(:) = (/imax,jmax,1/)
      write(6,*)'Transforming Q k=',k
      call sptez_mod (0, maxwv, idrt, imax, jmax, data%q(:,k,1), arr, 1)
!$OMP CRITICAL
      write(6,*)'Writing Q k=',k
      call handle_err (nf90_put_var (ncid, qid, arr, start=start, count=kount))
!$OMP END CRITICAL
    end subroutine xt2

    subroutine xt4(imax,jmax,k)
      integer, intent(in) :: imax, jmax, k

      integer :: start(3)                ! netcdf slice starting locations
      integer :: kount(3)                ! netcdf slice count information
      real :: arr(imax,jmax)

      start(:) = (/1,1,k/)
      kount(:) = (/imax,jmax,1/)
      write(6,*)'Transforming T k=',k
      call sptez_mod (0, maxwv, idrt, imax, jmax, data%t(:,k), arr, 1)
!$OMP CRITICAL
      write(6,*)'Writing T k=',k
      call handle_err (nf90_put_var (ncid, tvid, arr, start=start, count=kount))
!$OMP END CRITICAL
    end subroutine xt4

    subroutine xt5(imax,jmax,k)
      integer, intent(in) :: imax, jmax, k

      integer :: start(3)                ! netcdf slice starting locations
      integer :: kount(3)                ! netcdf slice count information
      real :: arr1(imax,jmax)
      real :: arr2(imax,jmax)

      start(:) = (/1,1,k/)
      kount(:) = (/imax,jmax,1/)

      write(6,*)'Transforming D,Z -> U,V k=',k
      call sptezmv_mod (0, maxwv, idrt, imax, jmax, 1, data%d(:,k), data%z(:,k), arr1, arr2, 1)
!$OMP CRITICAL
      write(6,*)'Writing U k=',k
      call handle_err (nf90_put_var (ncid, uid, arr1, start=start, count=kount))
      write(6,*)'Writing V k=',k
      call handle_err (nf90_put_var (ncid, vid, arr2, start=start, count=kount))
!$OMP END CRITICAL
    end subroutine xt5


    subroutine xt6(imax,jmax,k)
      integer, intent(in) :: imax, jmax, k

      integer :: start(3)                ! netcdf slice starting locations
      integer :: kount(3)                ! netcdf slice count information
      real :: arr1(imax,jmax)
      real :: arr2(imax,jmax)

      start(:) = (/1,1,k/)
      kount(:) = (/imax,jmax,1/)

      write(6,*)'Transforming O3 k=',k
      call sptez_mod (0, maxwv, idrt, imax, jmax, data%q(:,k,2), arr1, 1)
      write(6,*)'Transforming CLW k=',k
      call sptez_mod (0, maxwv, idrt, imax, jmax, data%q(:,k,3), arr2, 1)
!$OMP CRITICAL
      write(6,*)'Writing O3 k=',k
      call handle_err (nf90_put_var (ncid, ozid, arr1, start=start, count=kount))
      write(6,*)'Writing CLW k=',k
      call handle_err (nf90_put_var (ncid, clwid, arr2, start=start, count=kount))
!$OMP END CRITICAL
    end subroutine xt6
end program spe2nc
