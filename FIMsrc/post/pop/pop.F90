!=============================================================
! Post processing utility
!
!
! Main packages: slint (spherical linear interpolation)
!                vlint (vertical linear interpolation)
!                gribio (FIM grib utility routines)
!                varinfo (variable information)
!
! Ning Wang, March 2007 
!
!
!=============================================================

      program pop
#ifdef NAG
      USE F90_UNIX_ENV, ONLY : iargc
#endif

      use headers,only: testcurveheader,testglvlheader
      USE varinfo
      use fimnamelist,only: TotalTime, ArchvIntvl, curve,           &
          NumCacheBlocksPerPE, PrintMAXMINtimes, FixedGridOrder,    &
          TimingBarriers, glvl, nvl, yyyymmddhhmm, ArchvTimeUnit,   &
          post_datadir, outputdir, input, output, output_fmt,       &
          max_vars, var_list, multiple_output_files, gribtable,     &
          grid_id,  mx, my, latlonfld, is, vres, mode, nsmooth_var, &
          t1, t2, delta_t,multiHrTu
      use module_control, only : control, nvlp, pres_hpa, nip, numphr
      USE slint, ONLY: slint_init, bilinear_interp,  bilinear_init_i2r, &
                       bilinear_interp_i2r, bilinear_interp_uv

      IMPLICIT NONE

! define a derived type for variable data
      TYPE array_pointer
         REAL, DIMENSION(:,:), POINTER :: p
      END TYPE array_pointer

      INTEGER                 :: nverlvs  !The number of vertical levels defined in the Makefile
#ifndef NAG
      INTEGER                 :: iargc    ! command line argument count
#endif
      INTEGER                 :: nvars, nfct
      CHARACTER(len=char_len) :: init_file
      CHARACTER(len=6)        :: ahr 
      CHARACTER(len=8)        :: FMT='(I3.3)'
      CHARACTER(len=80)       :: FMT2
      INTEGER                 :: file_handle 
      INTEGER                 :: var_num, nlevels
      INTEGER                 :: year, month, day, hour, minute, jday, IW3JDN, is2Dvar
      INTEGER                 :: nt, i, j, k , idx, ierr 
      REAL, ALLOCATABLE       :: data_xyz(:,:,:), vardata(:), vardata_n(:), src_data(:), tgt_data(:)
      REAL, ALLOCATABLE       :: data_xyz2(:,:,:), vardata2(:), src_data2(:), tgt_data2(:)
      REAL, ALLOCATABLE       :: data_xyz_var(:,:,:)
      REAL, ALLOCATABLE       :: data_xyz_pr(:,:,:)
      REAL, ALLOCATABLE       :: ll_src(:,:)
      REAL, ALLOCATABLE       :: ll_tgt(:,:)
      REAL, ALLOCATABLE       :: cs_rot(:,:)
      CHARACTER(len=char_len) :: var_name, var_description, units, var_name2
      CHARACTER(len=19      ) :: date_str, date_str2
      CHARACTER(len=10      ) :: jdate 
      CHARACTER(len=char_len) :: gribfile
      REAL                    :: missing_value, r2d, pi
      integer :: ioerr, iret_rd1var
      integer :: ret      ! returned from subroutine call
      integer :: maxlevs  ! max number of levels for array allocation

! check and get the commandline arguments
      IF (iargc() .NE. 0) THEN
        WRITE(0,*) 'Usage: pop'
        WRITE(0,*) 'Note: Make sure that namelist file pop.nl is in the'// &
                    ' current directory'
        STOP
      END IF

! set up control variables

      call control ()
      nverlvs = nvl

      PRINT*, t1, t2, delta_t

      nvars = 0
      do while (var_list(nvars+1) /= ' ' .and. nvars < max_vars)
        nvars = nvars + 1
      end do

!JR is=1 in default FIM => horiz. interp and no vert. interp
      IF (is == 3 .AND. nvars > 3) THEN
        WRITE(0,*)  &
        'Only allow maximum 3 variables for vertical cross  sections'
        STOP
      ENDIF

! get date info from the date string
      READ(UNIT=yyyymmddhhmm(1:4), FMT='(I4)') year
      READ(UNIT=yyyymmddhhmm(5:6), FMT='(I2)') month
      READ(UNIT=yyyymmddhhmm(7:8), FMT='(I2)') day
      READ(UNIT=yyyymmddhhmm(9:10), FMT='(I2)') hour
      READ(UNIT=yyyymmddhhmm(11:12), FMT='(I2)') minute

! create a year 'month-date-hour-minute' date string
      date_str = yyyymmddhhmm(1:4) // "-" // yyyymmddhhmm(5:6) //  &
                 "-" // yyyymmddhhmm(7:8) // "-" // yyyymmddhhmm(9:10) & 
                 // ":" // yyyymmddhhmm(11:12) // ":00"
      date_str2 = date_str
     
! create the jdate string
      jday = IW3JDN(year,month,day) - IW3JDN(year,1, 1) + 1
      WRITE(UNIT=jdate(1:2), FMT='(I2.2)') MOD (year, 100) 
      WRITE(UNIT=jdate(3:5), FMT='(I3.3)') jday 
      WRITE(UNIT=jdate(6:7), FMT='(I2.2)') hour 
      jdate = jdate(1:7) // '000'

! compute the number of forecast time and number of icosahedral grid
      nfct = (t2 - t1) / delta_t 

      IF (output_fmt == "grib") THEN
        CALL gridid2mxmy(grid_id, mx, my)
      ENDIF

      IF (is == 0) THEN
        mx = 1024
        my = nip / 1024 + 1
      ENDIF

      IF (is == 0) THEN
        ALLOCATE(vardata(nip * (nverlvs+1)),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of vardata 0 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(vardata_n(mx * my * (nverlvs+1)),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of vardata_n 0 failed'
          call flush(6)
          stop
        endif

      ELSE IF (is == 1) THEN
!JR The following is a HACK to get pop to work when number of pressure levels exceeds number of model levels
        maxlevs = max (nverlvs+1,nvlp)

        ALLOCATE(vardata(nip * maxlevs),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of vardata 1 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(vardata2(nip * maxlevs),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of vardata2 1 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(data_xyz(mx, my, maxlevs),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of data_xyz 1 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(data_xyz2(mx, my, maxlevs),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of data_xyz2 1 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(src_data(nip), tgt_data(mx*my),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of src_data 1 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(src_data2(nip), tgt_data2(mx*my),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of src_data2 1 failed'
          call flush(6)
          stop
        endif

      ELSE IF (is == 2) THEN
        ALLOCATE(vardata(nip * (nverlvs+1)),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of vardata 2 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(data_xyz_var(mx,my,nverlvs+1),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of data_xyz_var 2 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(data_xyz_pr(mx,my,nverlvs+1),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of data_xyz_pr 2 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(data_xyz(mx, my, nvlp),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of data_xyz 2 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(src_data(nip), tgt_data(mx*my),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of src_data 2 failed'
          call flush(6)
          stop
        endif

      ELSE IF (is == 3) THEN

        ALLOCATE(vardata(nip * (nverlvs+1)),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of vardata 3 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(data_xyz_var(mx,my,nverlvs+1),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of data_xyz_var 3 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(data_xyz_pr(mx,my,nverlvs+1),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of data_xyz_pr 3 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(data_xyz(mx, my, vres),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of data_xyz 3 failed'
          call flush(6)
          stop
        endif

        ALLOCATE(src_data(nip), tgt_data(mx*my),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of src_data 3 failed'
          call flush(6)
          stop
        endif

      ENDIF

      missing_value = -99.0

! open and init a GRIB file
      IF(is == 0) THEN
        CALL set_model_nlevels(nverlvs)
        CALL initgrib(gribtable)
      ELSE IF(is == 1) THEN
        IF (output_fmt == "grib") THEN
          CALL set_model_nlevels(nverlvs)
          CALL initgrib(gribtable)
          IF (.not. multiple_output_files) CALL opengrib(output)
        END IF
      ELSE IF (is == 2) THEN
        IF (output_fmt == "grib") THEN
          CALL set_model_nlevels(nverlvs)
          CALL initgrib(gribtable)
          IF (.not. multiple_output_files) CALL opengrib(output)
        END IF
      ELSE IF (is == 3) THEN
        IF (output_fmt == "grib") THEN
          CALL set_model_nlevels(nverlvs)
          CALL initgrib(gribtable)
          IF (.not. multiple_output_files) CALL opengrib(output)
        END IF
      ENDIF
    
      post_datadir = post_datadir(1:LEN_TRIM(post_datadir)) // '/'

! if interpolation scheme is not 0, init the horizontal interpolation.
      IF (is /= 0) THEN
        IF (FixedGridOrder) THEN
          FMT2 = '(a,"latlonIJ.dat")'
          write(init_file,FMT2) post_datadir(1:LEN_TRIM(post_datadir))

          ALLOCATE(ll_src(nip, 2),stat=ierr)
          if (ierr.ne.0) then
            write (*,'(a)') 'pop: allocation of ll_src 0 failed'
            call flush(6)
            stop
          endif

          OPEN (10, file=init_file, action='read', form='unformatted', iostat=ioerr)
          if (ioerr == 0) then
            write(6,*)'pop: successfully opened init_file=', trim(init_file)
          else
            write(6,*)'pop: failed to open init_file=', trim(init_file)
          end if
          call TestGlvlHeader(10,init_file,'pop',glvl)
          READ (10, iostat=ioerr) ll_src(:, 1), ll_src(:, 2)
          if (ioerr /= 0) then
            write(6,*)'pop: bad attempt to read ', trim (init_file), &
                      'nelem=', ubound(ll_src,1), ' iostat=', ioerr
          end if
          CLOSE(10)
        ELSE
          init_file = './icos_grid_info_level.dat'

          ALLOCATE(ll_src(nip, 2),stat=ierr)
          if (ierr.ne.0) then
            write (*,'(a)') 'pop: allocation of ll_src 1 failed'
            call flush(6)
            stop
          endif

          OPEN (10, file=init_file, action='read', form='unformatted', iostat=ioerr)
          if (ioerr == 0) then
            write(6,*)'pop: successfully opened init_file=', trim(init_file)
          else
            write(6,*)'pop: failed to open init_file=', trim(init_file)
          end if
          call TestGlvlHeader (10,init_file,'pop',glvl )
          call TestCurveHeader(10,init_file,'pop',curve)
          READ (10, iostat=ioerr) ll_src(:, 1), ll_src(:, 2)
          if (ioerr == 0) then
!            write(6,*)'pop: successfully read ll_src nelem=',ubound(ll_src,1)
          else
            write(6,*)'pop: bad attempt to read ', trim (init_file), & 
                      'nelem=', ubound(ll_src,1), ' iostat=', ioerr
          end if
          CLOSE(10)
        ENDIF
!        write(6,*)'JR pop: calling bilinear_init_i2r ll_src=', ll_src(:,1)

        ALLOCATE(ll_tgt(mx*my,2),cs_rot(mx*my,2), stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocations of ll_tgt and cs_rot failed'
          call flush(6)
          stop
        endif

        CALL getTgtGrid(grid_id, ll_tgt, cs_rot, mx,my)
        CALL slint_init(ll_src, nip, ll_tgt, mx * my) 

        DEALLOCATE(ll_src)
        DEALLOCATE(ll_tgt)

      ENDIF
    
      IF (latlonfld) THEN
        pi = 4.0*ATAN(1.0)
        r2d = 180.0 / pi 

        ALLOCATE(ll_src(nip, 2),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'pop: allocation of ll_src 2 failed'
          call flush(6)
          stop
        endif

        FMT2 = '(a,"latlonIJ.dat")'
        write(init_file,FMT2) post_datadir(1:LEN_TRIM(post_datadir))
        OPEN(10,file=init_file,status='old',form='unformatted')
        call TestGlvlHeader(10,'latlonIJ.dat','pop',glvl)
        READ(10) ll_src(:, 1), ll_src(:, 2)
        CLOSE(10)
        ll_src(1:nip, 1) = ll_src(1:nip, 1) * r2d
        ll_src(1:nip, 2) = (ll_src(1:nip, 2) - pi) * r2d
      END IF

      IF (is == 0) THEN  ! no interpolation, native grid can only be saved in GRIB file format at this point.
        DO nt = t1, t2, delta_t
          IF (output_fmt.eq."grib" .and. multiple_output_files) THEN
            WRITE(ahr, FMT) nt
            gribfile = outputdir(1:LEN_TRIM(outputdir)) // '/' // jdate  // ahr
            CALL opengrib(gribfile)
          ELSE
             PRINT*,'Native grid can only be saved in GRIB file format.'
             STOP
          ENDIF
          DO var_num = 1, nvars ! for each variable
            var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
            CALL read_1var(nt, delta_t, post_datadir, vardata, var_name, nip, nverlvs + 1,ArchvTimeUnit, nvlp, iret_rd1var)
            IF (iret_rd1var == -1) THEN
              CYCLE
            ENDIF
            CALL var_info(var_list(var_num), var_description, units, nlevels, nvlp)
            DO i = 1, nlevels  
              DO j = 1, nip
                vardata_n(mx*my*(i-1)+j) = vardata(i + (j - 1) * (nverlvs + 1))
              END DO
              Do j = nip + 1, mx * my
                vardata_n(mx*my*(i-1)+j) = vardata(i + (nip - 1) * (nverlvs + 1))
              END DO
            ENDDO
            IF (nlevels > 2) THEN
              var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num))) // '_B'
            ELSE
              var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
            ENDIF
            CALL writegrib(var_name, grid_id, nlevels, glvl, curve, vardata_n, jdate, nt, &
                           tba(nt,t1,delta_t,var_name), delta_t, nvlp, pres_hpa,  &
                           ArchvTimeUnit, multiHrTu)
          ENDDO
          IF (latlonfld) THEN
            vardata_n(1:nip) = ll_src(1:nip,1) 
            vardata_n(nip+1:mx * my) = 0.0 
            var_name = "LAT"
            CALL writegrib(var_name, grid_id, 1, glvl, curve, vardata_n, jdate, t1, 0, &
                           delta_t, nvlp, pres_hpa, ArchvTimeUnit, multiHrTu)
            vardata_n(1:nip) = ll_src(1:nip,2) 
            vardata_n(nip+1:mx * my) = 0.0 
            var_name = "LON"
            CALL writegrib(var_name, grid_id, 1, glvl, curve, vardata_n, jdate, t1, 0, &
                           delta_t, nvlp, pres_hpa, ArchvTimeUnit, multiHrTu)
            DEALLOCATE(ll_src)
          ENDIF
          IF (output_fmt.eq."grib" .and. multiple_output_files) CALL closegrib()
#ifdef NAG
          write (*, '(A1)', ADVANCE = "NO" ) '*'
#else
          PRINT "('*'$)"  ! progress '*'
#endif
        ENDDO
      ELSE IF (is == 1) THEN ! Horizontal interpolation
        DO nt = t1, t2, delta_t
          IF (output_fmt.eq."grib" .and. multiple_output_files) THEN
            WRITE(ahr, FMT) nt
            gribfile = outputdir(1:LEN_TRIM(outputdir)) // '/' // jdate // ahr
            CALL opengrib(gribfile)
          ENDIF
          DO var_num = 1, nvars ! for each variable
            IF (input /= "") THEN
              CALL read_1var_direct(input, nip, nverlvs, vardata)
              IF (iret_rd1var == -1) THEN
                CYCLE
              ENDIF
              nlevels = nverlvs
            ELSE
              var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
              IF ((var_name == 'r42D' .OR. var_name == 'r52D' .OR.  &
                   var_name == 'r62D' .OR. var_name == 's62D')   &
                   .AND. mod(nt, 6) /= 0) THEN
                CYCLE
              ENDIF 
              CALL read_1var(nt, delta_t, post_datadir, vardata, var_name, nip, nverlvs + 1,ArchvTimeUnit, nvlp, iret_rd1var)
              IF (iret_rd1var == -1) THEN
                CYCLE
              ENDIF
              CALL var_info(var_list(var_num), var_description, units, nlevels, nvlp)
              IF (var_name == 'us3D' .OR. var_name == 'up3P' .OR. var_name == 'u12D' ) THEN
                var_name2 = 'v' // var_name(2:LEN_TRIM(var_name))
                CALL read_1var(nt, delta_t, post_datadir, vardata2, var_name2, nip, nverlvs + 1,ArchvTimeUnit, nvlp, iret_rd1var)
                IF (iret_rd1var == -1) THEN
                  CYCLE
                ENDIF
              ENDIF
            ENDIF

            IF (var_name == 'us3D' .OR. var_name == 'up3P' .OR. var_name == 'u12D' ) THEN
              DO k=1,nlevels 
                DO i = 1, nip
                  src_data(i) = vardata(k + (i - 1) * nlevels) 
                  src_data2(i) = vardata2(k + (i - 1) * nlevels) 
                END DO
                CALL bilinear_interp_uv(src_data, tgt_data, src_data2, tgt_data2) 
                DO i = 1, mx
                  DO j = 1, my
                    idx = (j - 1) * mx + i
                    data_xyz(i,j,k) = cs_rot(idx,1)*tgt_data(idx)-cs_rot(idx,2)*tgt_data2(idx)
                    data_xyz2(i,j,k) = cs_rot(idx,2)*tgt_data(idx)+cs_rot(idx,1)*tgt_data2(idx) 
                  END DO
                END DO
              ENDDO
            ELSE IF ( var_name == 'vs3D' .OR. var_name == 'vp3P' .OR. var_name == 'v12D') THEN
              CYCLE
            ELSE
              DO k=1,nlevels 
                DO i = 1, nip
                  src_data(i) = vardata(k + (i - 1) * nlevels) 
                END DO
                CALL bilinear_interp(src_data, tgt_data) 
                DO i = 1, mx
                  DO j = 1, my
                    data_xyz(i, j, k) = tgt_data((j - 1) * mx + i) 
                  END DO
                END DO
              END DO 
            ENDIF

            DO i = 1, nsmooth_var(var_num)
              CALL smooth(data_xyz,mx,my,nlevels,0.2)
            END DO
            
            IF (var_name == 'us3D' .OR. var_name == 'up3P' .OR. var_name == 'u12D' ) THEN
              DO i = 1, nsmooth_var(var_num)
                CALL smooth(data_xyz2,mx,my,nlevels,0.2)
              END DO
            ENDIF

            IF (nlevels < nverlvs + 1) THEN
              data_xyz(:,:,nlevels + 1:) = 0.0
            END IF
            IF (output_fmt == "grib") THEN
              IF (nlevels > 2) THEN
                var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
!JR The "_B" means the data are on native levels rather than pressure levels
                IF(var_name /= "hgtP" .AND. var_name /= "tmpP" .AND. &
                  var_name /= "up3P" .AND. var_name /= "vp3P" .AND. &
                  var_name /= "oc1P" .AND. var_name /= "oc2P" .AND. &
                  var_name /= "bc1P" .AND. var_name /= "bc2P" .AND. &
                  var_name /= "so2P" .AND. var_name /= "slfP" .AND. &
                  var_name /= "d1sP" .AND. var_name /= "d2sP" .AND. &
                  var_name /= "d3sP" .AND. var_name /= "d4sP" .AND. &
                  var_name /= "d5sP" .AND. var_name /= "s1sP" .AND. &
                  var_name /= "s2sP" .AND. var_name /= "s3sP" .AND. &
                  var_name /= "s4sP" .AND. var_name /= "dmsP" .AND. &
                  var_name /= "msaP" .AND. var_name /= "p25P" .AND. &
                  var_name /= "rh3P" .AND. var_name /= "p10P" .AND. &
                  var_name /= "qc3P" .AND.                          &
                  var_name /= "oz3P" .AND. var_name /= "vv3P") THEN
                  var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num))) // '_B'
                ENDIF
              ELSE
                var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
              ENDIF
              CALL writegrib(var_name, grid_id, nlevels, 0, 0, data_xyz, jdate, nt, &
                             tba(nt,t1,delta_t,var_name), delta_t, nvlp, pres_hpa, &
                             ArchvTimeUnit, multiHrTu)
              IF (var_name == 'us3D_B' .OR. var_name == 'up3P' .OR. var_name == 'u12D' ) THEN
                var_name2 = 'v' // var_name(2:LEN_TRIM(var_name))
                CALL writegrib(var_name2, grid_id, nlevels, 0, 0, data_xyz2, jdate, nt, &
                               tba(nt,t1,delta_t,var_name), delta_t, nvlp, pres_hpa, &
                               ArchvTimeUnit, multiHrTu)
              ENDIF
            ENDIF
          END DO 
          CALL mmdddateadj(date_str, delta_t)  
          IF (output_fmt.eq."grib" .and. multiple_output_files) CALL closegrib()
#ifdef NAG
          write (*, '(A1)', ADVANCE = "NO" ) '*'
#else
          PRINT "('*'$)"  ! progress '*'
#endif
        END DO 
      ELSE ! vertical interpolation
        var_list(nvars+1) = "pr3D"
        DO nt = t1, t2, delta_t
          IF (output_fmt.eq."grib" .and. multiple_output_files) THEN
            WRITE(ahr, FMT) nt
            gribfile = outputdir(1:LEN_TRIM(outputdir)) // '/' // jdate // ahr
            CALL opengrib(gribfile)
          ENDIF

          ! interpolate pressure to latlon grid 
          var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
          CALL read_1var(nt, delta_t, post_datadir, vardata, var_name, nip, nverlvs + 1,ArchvTimeUnit,nvlp, iret_rd1var)
          CALL var_info(var_list(nvars+1), var_description, units, nlevels,nvlp)
          DO k = 1, nlevels  
            CALL bilinear_interp_i2r (k, nlevels, vardata, data_xyz_pr) 
          END DO

          ! interpolate the specified variables to latlon grid 
          DO var_num = 1, nvars ! for each variable
            var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
            CALL read_1var(nt, delta_t, post_datadir, vardata, var_name, nip, nverlvs + 1,ArchvTimeUnit, nvlp, iret_rd1var)
            CALL var_info(var_list(var_num), var_description, units, nlevels, nvlp)
            DO k = 1, nlevels 
              CALL bilinear_interp_i2r (k, nlevels, vardata, data_xyz_var) 
            END DO 
            ! interpolate to vertical plane
            data_xyz = 0.0
            IF (is == 2) THEN
              IF (is2Dvar(var_list(var_num)) == 1) THEN
                var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
                DO i = 1, nsmooth_var(var_num)
                  CALL smooth(data_xyz_var,mx,my,1,0.2)
                END DO
                CALL writegrib(var_name, grid_id, 1, 0, 0, data_xyz_var, jdate, nt, &
                               tba(nt,t1,delta_t,var_name), delta_t,nvlp,pres_hpa, &
                               ArchvTimeUnit, multiHrTu) 
                CYCLE
              ENDIF
              CALL vlint2coor(mx, my, nlevels, nverlvs + 1, data_xyz_pr, data_xyz_var, data_xyz, pres_hpa, nvlp)  
            ELSE IF (nlevels == nverlvs) THEN
              IF (mode == "step") THEN
                CALL v_interp(mx, my, nlevels, data_xyz_pr, data_xyz_var, vres, data_xyz, missing_value, 1)  
              ELSE IF (mode == "linear") THEN
                CALL v_interp(mx, my, nlevels, data_xyz_pr, data_xyz_var, vres, data_xyz, missing_value, 0)  
              ENDIF
            ELSEIF (nlevels == nverlvs + 1) THEN ! variables defined on interface
              CALL v_interp_lvlvar(mx, my, nlevels, data_xyz_pr, data_xyz_var, vres, data_xyz, missing_value, 0)
            ENDIF
            DO i = 1, nsmooth_var(var_num)
              IF (is == 2) THEN
                CALL smooth(data_xyz,mx,my,nvlp,0.2)
              ELSE if (is == 3) THEN
                CALL smooth(data_xyz,mx,my,vres,0.2)
              END IF
            END DO

! write data to the GRIB file, the variable interpolated
            IF (output_fmt == "grib") THEN
              var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num))) // '_P'
              IF (is == 2) THEN
                CALL writegrib(var_name, grid_id, nvlp, 0, 0, data_xyz, jdate, nt, &
                               tba(nt,t1,delta_t,var_name), delta_t,nvlp,pres_hpa, &
                               ArchvTimeUnit, multiHrTu)
              ELSE
                CALL writegrib(var_name, grid_id, vres, 0, 0, data_xyz, jdate, nt, &
                               tba(nt,t1,delta_t,var_name), delta_t,nvlp,pres_hpa, &
                               ArchvTimeUnit, multiHrTu)
              ENDIF
            ENDIF
          END DO
          IF (is == 3) THEN
            DO k=1,vres 
              IF (k <= nverlvs + 1) THEN
                data_xyz(1:mx, 1:my, k) = data_xyz_pr(1:mx, 1:my, k)
              ELSE
                data_xyz(1:mx, 1:my, k) = 0
              ENDIF
            END DO !level loop
            CALL var_info(var_list(nvars+1), var_description, units, nlevels, nvlp)
            IF (output_fmt == "grib") THEN
              var_name = var_list(nvars+1)(1:LEN_TRIM(var_list(nvars+1))) // '_B'
              CALL writegrib(var_name, grid_id, nlevels, 0, 0, data_xyz, jdate, nt, &
                             tba(nt,t1,delta_t,var_name), delta_t, nvlp, pres_hpa, &
                             ArchvTimeUnit, multiHrTu)
            ENDIF
          ENDIF
          CALL mmdddateadj(date_str, delta_t)  
          IF (output_fmt.eq."grib" .and. multiple_output_files) CALL closegrib()
#ifdef NAG
          write (*, '(A1)', ADVANCE = "NO" ) '*'
#else
          PRINT "('*'$)"  ! progress '*'
#endif
! time iteration
        END DO 
      ENDIF

      IF (output_fmt == "grib") THEN
        CALL endgrib()
      ENDIF

      PRINT*, '  '

      IF (is == 0) THEN
        DEALLOCATE(vardata, vardata_n)
      ELSE IF (is == 1) THEN
        DEALLOCATE(vardata, data_xyz)
      ELSE
        DEALLOCATE(data_xyz_var, data_xyz_pr, vardata, data_xyz)
      ENDIF

      write (*,'(a)') 'pop completed successfully'

      STOP

    contains

      integer function tba (nt, t1, delta_t, varname)
        implicit none

        integer, intent(in) :: nt, t1, delta_t
        character(len=*), intent(in) :: varname
        if (nt >= delta_t .and. (varname == 'r12D' .or.  &
            varname == 'r22D' .or. varname == 'r32D')) then
          tba = nt - delta_t
        else if (nt >= 6 .and. (varname == 'r42D' .or.  &
            varname == 'r52D' .or. varname == 'r62D')) then
          tba = nt - 6
        else if (nt >= delta_t .and. (varname == 's12D' .or.  &
            varname == 's22D' .or. varname == 's32D')) then
          tba = nt - delta_t
        else if (nt >= 6 .and. (varname == 's42D' .or.  &
            varname == 's52D' .or. varname == 's62D')) then
          tba = nt - 6
        else
          tba = 0
        end if
        return
      end function tba
    end program pop

SUBROUTINE read_1var(nt, delta_t, post_datadir, vardata, var_name, nip, nverlvs,ArchvTimeUnit, nvlp, iret)
      USE varinfo, ONLY: var_info,char_len
      IMPLICIT NONE
      
      INTEGER nt, delta_t, nip, nverlvs, nvlp
      CHARACTER *(*) :: post_datadir
      CHARACTER *(*) :: var_name
      CHARACTER(80) :: header(10)
!JR Don't know the actual dimension of vardata
      REAL vardata(*)
      character(2),intent(in)::ArchvTimeUnit

      CHARACTER (len=char_len)  :: var_description, units
      INTEGER :: nlevels

      INTEGER time, lunout, is2Dvar, ioerr, iret, iret_2d
      CHARACTER (len=char_len) :: filename

      iret = 0

! when all 2d vars, sm3d, and st3d are in one file, comment the following statements
! if these variables are in seperate files.
      IF (is2Dvar(var_name) == 1) THEN
        CALL read_2Dvar(nt, post_datadir, var_name, nip, nverlvs, vardata,ArchvTimeUnit, nvlp, iret_2d)
        iret = iret_2d
        RETURN
      ENDIF

      lunout = 29
      time = nt
      WRITE(filename,"('fim_out_',a4,i6.6,a2)") var_name,time,ArchvTimeUnit
      filename = post_datadir(1:LEN_TRIM(post_datadir)) // filename
      CALL var_info(var_name, var_description, units, nlevels, nvlp)
      OPEN (lunout, file=filename, form="unformatted",status='old', &
            action='read', iostat=ioerr)

      if (ioerr /= 0) then
        write(6,*)'read_1var: Failed to open file ', & 
                   trim(filename), '.'
        iret = -1
        return
      end if

      READ (lunout, iostat=ioerr) header
      if (ioerr /= 0) then
        write(6,*)'read_1var: Failed to read header from file ', &
                   trim(filename), '. Stopping'
        iret = -1
        return
      end if

      READ (lunout, iostat=ioerr) vardata(1:nip*nlevels)
      if (ioerr /= 0) then
        write(6,*)'read_1var: Failed to read vardata from file ', &
                   trim(filename), '. ioerr=', ioerr
        write(6,*)'var_name=', var_name, 'nlevels=', nlevels
        write(6,*)'Stopping.'
        iret = -1
        return
      end if

      CLOSE(lunout)
      ! special treament for qv3d and qw3d
      IF (var_name == "qv3D" .OR. var_name == "qw3D") THEN
        vardata(1:nip*nlevels) = vardata(1:nip*nlevels) * 1000.0
      ENDIF
      ! special treament for ph3d 
      IF (var_name == "ph3D") THEN
!JR Changed to multiply by reciprocal to match scalefactor passed in when gribout=true
        vardata(1:nip*nlevels)  = vardata(1:nip*nlevels) * (1./9.8)
      ENDIF
      ! special treament for oz3d
      IF (var_name == "oz3D") THEN
          vardata(1:nip*nlevels)  = vardata(1:nip*nlevels) * 1000.0
      ENDIF
!print*, 'var ', var_name, minval(vardata(1:nip*nlevels)), maxval(vardata(1:nip*nlevels))
END SUBROUTINE read_1var

SUBROUTINE read_2Dvar(nt, post_datadir, var_name, nip, nverlvs, vardata,ArchvTimeUnit, nvlp, iret)
      USE varinfo, ONLY: var_info, char_len
      IMPLICIT NONE
      
      CHARACTER *(*) :: post_datadir, var_name
      INTEGER  nt, nip, nverlvs, ioerr, iret

!JR Don't know the actual dimension of vardata
      REAL vardata(*)

      INTEGER i, j, lunout, nlevels, nvlp
      CHARACTER(len=char_len) :: filename

      CHARACTER(len=80) :: header(10)
      CHARACTER(len=char_len) :: varname, varname_uc
      CHARACTER(len=char_len)  :: var_description, units
      character(2),intent(in)::ArchvTimeUnit
      
      iret = 0

      lunout = 29
      varname_uc = var_name(1:LEN_TRIM(var_name))
      CALL tolowercase(var_name)
      WRITE(filename,"('fim_out_2D__',i6.6,a2)") nt,ArchvTimeUnit
      filename = post_datadir(1:LEN_TRIM(post_datadir)) // filename
      OPEN(lunout,file=filename,status='old',form="unformatted",  &
            action='read', iostat=ioerr)

      IF (ioerr /= 0) THEN
        WRITE(6,*)'read_2Dvar: Failed to open file ', & 
                   trim(filename), '.'
        iret = -1
        RETURN
      ENDIF

      READ(lunout, iostat=ioerr) header
      READ(header,FMT="(4X,A4)") varname  ! for now 10:00 am
      CALL tolowercase(varname)

      DO WHILE (var_name /= varname) 
          IF (varname == "sm3d" .OR. varname == "st3d") THEN
            READ(lunout, iostat=ioerr) vardata(1:4*nip)
          ELSE
            READ(lunout, iostat=ioerr) vardata(1:nip)
          ENDIF
          READ(lunout, iostat=ioerr) header
          READ(header,FMT="(4X,A4)") varname 
          CALL tolowercase(varname)
      END DO

      IF (ioerr /= 0) THEN
        iret = -1
        RETURN
      END IF 

      var_name = varname_uc(1:LEN_TRIM(varname_uc))
      CALL var_info(var_name, var_description, units, nlevels, nvlp)

      READ (lunout, iostat=ioerr) vardata(1:nip*nlevels)
      CLOSE(lunout)

      IF (ioerr /= 0) THEN
        iret = -1
        RETURN
      END IF 

END SUBROUTINE read_2Dvar

SUBROUTINE read_1var_direct(filename, nip, nverlvs, vardata)
      USE varinfo
      IMPLICIT NONE
      
      CHARACTER *(*) :: filename
      INTEGER  nip, nverlvs
      REAL vardata(nip * nverlvs)

      INTEGER :: lunout=30

      OPEN (lunout,file=filename,form="unformatted")
      READ(lunout) vardata(1:nip*nverlvs)
      CLOSE(lunout)
END SUBROUTINE read_1var_direct

INTEGER FUNCTION is2Dvar(var_name)
      IMPLICIT NONE

      CHARACTER *(*) :: var_name

      IF(var_name(3:4) == "2D"   .or. var_name(3:4) == "2d"   .or. var_name(1:4) == "iash" .or.	&
         var_name(1:4) == "hfss" .or. var_name(1:4) == "hfls" .or. var_name(1:4) == "rsds" .or.	&
         var_name(1:4) == "rlds" .or. var_name(1:4) == "rlut" .or. var_name(1:4) == "prpv" .or. &
         var_name(1:4) == "rsus" .or. var_name(1:4) == "rlus" .or. var_name(1:4) == "w080" .or. &
         var_name(1:4) == "ustr" .or. var_name(1:4) == "vstr" .or. var_name(1:4) == "tmax" .or. &
         var_name(1:4) == "ssta" .or. var_name(1:4) == "cice" .or. var_name(1:4) == "tmin" .or. &
         var_name(1:4) == "mslp" .or. var_name(1:4) == "cltt" .or. var_name(1:4) == "runo" .or. &
         var_name(1:4) == "thpv" .or. var_name(1:4) == "uspv" .or. var_name(1:4) == "vspv" ) THEN
        is2Dvar = 1
      ELSE
        is2Dvar = 0
      ENDIF

      RETURN
END FUNCTION is2Dvar


SUBROUTINE mmdddateadj(mmdddate, delta_t)
      IMPLICIT NONE
      CHARACTER*19 mmdddate
      INTEGER delta_t

      INTEGER year, month, day, hour, minute,second,ly 
      INTEGER month_ny(12), month_ly(12)
      CHARACTER dash, col
      DATA month_ny/31,28,31,30,31,30,31,31,30,31,30,31/
      DATA month_ly/31,29,31,30,31,30,31,31,30,31,30,31/
      dash = '-'
      col = ':'

      READ(mmdddate,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') year,month,day,hour,minute,second
!print*, year,month,day,hour,minute,second
      IF (mod(year, 4) == 0 .AND. mod(year, 100) /= 0) THEN
        ly = 1
      ELSE
        ly = 0 
      ENDIF
      IF (hour + delta_t .GE. 24) THEN
        day = day + 1
        IF ((ly == 0 .AND. day > month_ny(month)) .OR. (ly == 1 .AND. day > month_ly(month))) THEN
          day = 1
          month = month + 1
        ENDIF 
        IF (month > 12) THEN
          month = 1
          year = year + 1
        ENDIF
      ENDIF
      hour = MOD(hour + delta_t, 24)
!print*, year,month,day,hour,minute,second
      WRITE(mmdddate,'(i4.4,A1,i2.2,A1,i2.2,A1,i2.2,A1,i2.2,A1,i2.2)') year,dash,month,dash,day,dash,hour,col,minute,col,second 
     
END SUBROUTINE mmdddateadj

SUBROUTINE jdateadj(jdate, delta_t)
      IMPLICIT NONE
      CHARACTER*7 jdate
      INTEGER delta_t

      INTEGER yr, jday, hr, ly 

      READ(jdate,'(i2,i3,i2)') yr,jday,hr 
      IF (hr + delta_t .GE. 24) THEN
        jday = jday + 1
        IF (mod(yr, 4) == 0 .AND. mod(yr, 100) /= 0) THEN
          ly = 1
        ELSE
          ly = 0 
        ENDIF
        IF (jday > 365 + ly) THEN
          jday = 1
          yr = yr + 1
        ENDIF
      ENDIF
      hr = MOD(hr + delta_t, 24)
      WRITE(jdate,'(i2.2,i3.3,i2.2)') yr,jday,hr
     
END SUBROUTINE jdateadj

