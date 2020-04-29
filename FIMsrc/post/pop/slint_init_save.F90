PROGRAM slint_init_save_prog
#ifdef NAG
    USE F90_UNIX_ENV, ONLY : getarg, iargc
#endif
    USE slint, ONLY: slint_init_save
    IMPLICIT NONE

#ifndef NAG
    INTEGER :: iargc
#endif
    INTEGER :: grid_id
    CHARACTER(len=3) :: grid_id_asc
    CHARACTER(len=128) :: src_grid_file, tgt_grid_file, coeffs_file
    REAL, ALLOCATABLE :: ll_src(:,:), ll_tgt(:,:)
    INTEGER :: ngp_src, ngp_tgt
    INTEGER :: cmd_arg_ct, num, ioerr
    integer :: ierr

    cmd_arg_ct = iargc()

    IF (cmd_arg_ct == 0) THEN
      CALL slint_init_save_default()
      STOP
    ELSE IF (cmd_arg_ct == 3) THEN
      CALL getarg(1, src_grid_file)
      CALL getarg(2, tgt_grid_file)
      CALL getarg(3, coeffs_file)
    ELSE
      WRITE(6,*)'slint_init_save_prog: Incorrect argument list!'
      STOP
    ENDIF

    OPEN (10, file=src_grid_file, action='read', form='unformatted', iostat=ioerr)
    IF (ioerr == 0) THEN
      WRITE(6,*)'slint init program: successfully opened init file=', trim(src_grid_file)
    ELSE
      WRITE(6,*)'slint init program : failed to open init file=', trim(src_grid_file)
      STOP
    ENDIF 

    READ(10) ngp_src  ! the number of grid points for source grid

    ALLOCATE(ll_src(ngp_src, 2),stat=ierr)
    IF (ierr.ne.0) THEN
      WRITE (*,'(a)') 'slint_init_save_prog: allocation of ll_src failed'
      CALL flush(6)
      STOP
    ENDIF

    READ (10, iostat=ioerr) ll_src(:, 1), ll_src(:, 2)
    IF (ioerr == 0) THEN
      WRITE(6,*)'slint init program: successfully read ll_src #elem =',ubound(ll_src,1)
    ELSE
      WRITE(6,*)'slint init program: bad attempt to read ', trim(src_grid_file), '#elem=', &
                ubound(ll_src,1), ' iostat=', ioerr
      STOP
    ENDIF 
    CLOSE(10)

    OPEN (10, file=tgt_grid_file, action='read', form='unformatted', iostat=ioerr)
    IF (ioerr == 0) THEN
      WRITE(6,*)'slint init program: successfully opened init file=', trim(tgt_grid_file)
    ELSE
      WRITE(6,*)'slint init program : failed to open init file=', trim(tgt_grid_file)
      STOP
    ENDIF 

    READ(10) ngp_tgt  ! the number of grid points for target grid

    ALLOCATE(ll_tgt(ngp_tgt,2), stat=ierr)
    IF (ierr.ne.0) THEN
      WRITE (*,'(a)') 'slint_init_save_prog: allocation of ll_tgt and failed'
      CALL flush(6)
      STOP
    ENDIF

    READ (10, iostat=ioerr) ll_tgt(:, 1), ll_tgt(:, 2)
    IF (ioerr == 0) THEN
      WRITE(6,*)'slint init program: successfully read ll_tgt #elem =',ubound(ll_tgt,1)
    ELSE
      WRITE(6,*)'slint init program: bad attempt to read ', trim(tgt_grid_file), '#elem=', &
                ubound(ll_tgt,1), ' iostat=', ioerr
      STOP
    ENDIF 
    CLOSE(10)

    CALL slint_init_save(coeffs_file, ll_src, ngp_src, ll_tgt, ngp_tgt) 
    DEALLOCATE(ll_src, ll_tgt)

END PROGRAM slint_init_save_prog 

SUBROUTINE slint_init_save_default()
    use headers,only: testcurveheader,testglvlheader
    USE fimnamelist, ONLY: readnl_serial, glvl, SubdivNum, curve, &
        NumPostGrids, PostGridIds
    USE slint, ONLY: slint_init_save
    IMPLICIT NONE

    INTEGER :: i, sl, nip
    INTEGER :: mx, my, ioerr, ret
    CHARACTER(len=256) :: init_file, coeff_file
    REAL, ALLOCATABLE :: ll_src(:,:), ll_tgt(:,:)
    REAL, ALLOCATABLE :: cs_rot(:,:)
    INTEGER grid_ids(10), n_tgts
    integer::ierr

! Initialize SubdivNum, grid_ids to 0...
    SubdivNum = 0
    grid_ids  = 0

    call readnl_serial
    n_tgts = NumPostGrids

! ...then fill in with values from namelist 
    do i = 1, n_tgts
      grid_ids(i)  = PostGridIds(i)
    enddo

    sl = 1
    i = 1
    DO WHILE (SubdivNum(i) /= 0 .AND. i <= glvl) 
      sl = sl * SubdivNum(i)
      i = i + 1
    ENDDO 

    IF (i /= glvl + 1) THEN
      PRINT*, "Namelist variable 'SubdivNum' specification is incomplete."
      STOP
    ENDIF

    nip = 10 * sl * sl + 2

    init_file = './icos_grid_info_level.dat'

    ALLOCATE(ll_src(nip, 2),stat=ierr)
    if (ierr.ne.0) then
      write (*,'(a)') 'slint_init_save_default: allocation of ll_src failed'
      call flush(6)
      stop
    endif

    OPEN (10, file=init_file, action='read', form='unformatted', iostat=ioerr)
    IF (ioerr == 0) THEN
      WRITE(6,*)'slint init program : successfully opened init_file=', trim(init_file)
    ELSE
      WRITE(6,*)'slint init program : failed to open init_file=', trim(init_file)
      STOP
    ENDIF 

    CALL TestGlvlHeader (10,init_file,'slint init',glvl )
    CALL TestCurveHeader(10,init_file,'slint init',curve)
    READ (10, iostat=ioerr) ll_src(:, 1), ll_src(:, 2)
    IF (ioerr == 0) THEN
      WRITE(6,*)'slint init program: successfully read ll_src #elem =',ubound(ll_src,1)
    ELSE
      WRITE(6,*)'slint init program: bad attempt to read ', trim (init_file), '#elem=', &
                ubound(ll_src,1), ' iostat=', ioerr
      STOP
    ENDIF 
    CLOSE(10)

    DO i = 1, n_tgts
      CALL gridid2mxmy(grid_ids(i), mx, my)

      ALLOCATE(ll_tgt(mx*my, 2),cs_rot(mx*my, 2), stat=ierr)
      IF (ierr.ne.0) THEN
        WRITE (*,'(a)') 'slint_init_save_default: allocation of ll_tgt cnd cs_rot  failed'
        CALL flush(6)
        STOP
      ENDIF

      CALL getTgtGrid(grid_ids(i), ll_tgt, cs_rot, mx,my)
      WRITE(coeff_file, '("grid_",I3.3,"_coeffs")') grid_ids(i)
      CALL slint_init_save(coeff_file, ll_src, nip, ll_tgt, mx * my) 
      PRINT*, 'Coefficient file ', trim(coeff_file), ' is saved.' 
      DEALLOCATE(ll_tgt, cs_rot)
    ENDDO

    DEALLOCATE(ll_src)

END SUBROUTINE slint_init_save_default 

