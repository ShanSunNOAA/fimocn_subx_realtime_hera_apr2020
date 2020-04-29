!*********************************************************************
! cnvMtnvar.F90
! This program converts T574 mtnvar file T1534 mtnvar file.
!*********************************************************************
!!#define F2C_INTERP

program cnvPattern
  use slint, only:slint_init, slint_init_f2c, bilinear_interp
  implicit none

  integer, parameter :: ntime = 241 
  real, allocatable:: gg_ll(:,:), icos_ll(:,:)
  character(len=128) :: gfsltln_file, icosltln_file
  integer :: imax, jmax, nip
  integer :: l, ierr

  imax=1760; jmax=880
  nip = 655362

  allocate(gg_ll(imax*jmax,2), icos_ll(nip,2))

  gfsltln_file = "gfsltln_t574.dat"
  OPEN (10,file=gfsltln_file, status='old',form='unformatted')
  READ (10,iostat=ierr) gg_ll(:, 1), gg_ll(:, 2)
  IF (ierr /= 0) then
    WRITE(6,*)'cnvMtnvar: Error reading from gfsltln_t574.dat, ierr=',ierr
  END IF
  CLOSE(10)
  
  icosltln_file = "icos_grid_info_level.dat"
  OPEN (10,file=icosltln_file, status='old',form='unformatted')
  READ (10,iostat=ierr) icos_ll(:, 1), icos_ll(:, 2)
  IF (ierr /= 0) then
    WRITE(6,*)'cnvPattern: Error reading from', icosltln_file, 'err=',ierr
  END IF
  CLOSE(10)

  PRINT*, 'Start slint_init() ...'  
#ifdef F2C_INTERP
  CALL slint_init_f2c(gg_ll, imax1*jmax1, icos_ll, nip)
  DEALLOCATE(gg_ll, icos_ll)
#else
  CALL slint_init(gg_ll, imax1*jmax1, icos_ll, nip)
  DEALLOCATE(gg_ll, icos_ll)
#endif
  PRINT*, 'Done slint_init().'  

  allocate(pattern_gg(imax,jmax,ntime), pattern_icos(nip,ntime))

  pat_dat_file = "patterns.dat" 
  call read_pattern(pattern_gg,imax,jmax,ntime,pat_dat_file)

  do l=1,ntime
    CALL bilinear_interp(pattern_gg(:,:,l), pattern_icos(:,l))
  end do

  pat_dat_file = "patterns_icos.dat" 
  call write_mtnvar(pattern_icos,nip,1,ntime,pat_dat_file)
  deallocate(pattern_gg, pattern_icos)

end program cnvPattern

