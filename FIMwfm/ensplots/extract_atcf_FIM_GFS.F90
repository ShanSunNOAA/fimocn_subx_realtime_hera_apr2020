PROGRAM extract_atcf

! called by plot_ellipses.pro, this should extract the desired GFS ensemble forecast 
! information and write it to idl_fcst.dat, which idl will read back in.

! NOTE: will need to change directory name of output file, possibly.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HARDWIRED CHANGES TO MAKE !!!!! - need to change cmem in 2 places!!!
! PARAMETER (nmembers = 4)
! DATA cmem /'01','02','03','04'/
! CHARACTER*2, DIMENSION(4) :: cmem
! build_filename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PARAMETER (nmembers = 20)
PARAMETER (nleads   = 29)  ! 28

TYPE ens_fcst
  CHARACTER(LEN=2)  :: basin
  INTEGER           :: id
  INTEGER           :: ibasetime
  INTEGER           :: iflead
  REAL              :: ctr_lat(nmembers)
  REAL              :: ctr_lon(nmembers)
  REAL              :: centralpressure(nmembers)
  REAL              :: windspeed(nmembers)
END TYPE

TYPE (ens_fcst), DIMENSION(nleads) :: hurricane_ens_fcst
TYPE (ens_fcst), DIMENSION(nmembers, nleads) :: GEFS

CHARACTER*150 :: infile, outfile, outfile_ens
CHARACTER*2, DIMENSION(nmembers)  :: cmem
CHARACTER*2   :: stormnumber

CHARACTER*10  :: cyyyymmddhh_in
CHARACTER*2   :: cbasin_in
CHARACTER*2   :: cstormno_in
CHARACTER*150  :: rundir_in
CHARACTER*150  :: gfs_rundir_in
CHARACTER*150  :: plotPath
CHARACTER*2   :: glvl_in
CHARACTER*3   :: fcst_len_in
CHARACTER*3   :: fcst_output_int_in

INTEGER       :: fcst_len
INTEGER       :: fcst_output_int
INTEGER       :: ios

LOGICAL iex

DATA cmem /'01','02','03','04','05','06','07','08','09','10','11','12', '13',&
'14', '15', '16', '17', '18', '19', '20'/


! ----------------------------------------------------------------
! input beginning and end date to process, and then find indices of
! dates to process in cmmdd array
! ----------------------------------------------------------------

CALL getarg (1, cyyyymmddhh_in)
CALL getarg (2, cstormno_in)
CALL getarg (3, cbasin_in)   
CALL getarg (4, rundir_in)   
CALL getarg (5, gfs_rundir_in)   
CALL getarg (6, plotPath)   
CALL getarg (7, glvl_in)   
CALL getarg (8, fcst_len_in)   
CALL getarg (9, fcst_output_int_in)   
read (fcst_len_in,*)fcst_len  
read (fcst_output_int_in,*)fcst_output_int

! ----------------------------------------------------------------
! set up forecast structure, and write a warning to the output file
! that will be overwritten if the program later runs successfully
! to conclusion
! ----------------------------------------------------------------
fcst_output_int=6

 print *, 'in extract: fcst_len: ',fcst_len,' fcst_output_int: ',fcst_output_int
 print *, 'in extract: rundir_in: ',rundir_in,' plotPath: ',plotPath
 print *, 'in extract: gfs_rundir_in: ',gfs_rundir_in
DO i = 1, nleads
  CALL init_ens_fcst_structure (hurricane_ens_fcst(i))
END DO

! print *, 'in extract: before outfile '
! outfile = trim(rundir_in) // '/fimens/idl_fcst.dat'
outfile = trim(plotPath) // '/idl_fcst.dat'
print *, 'in extract: outfile: ',outfile
OPEN (UNIT=1, FILE=outfile, STATUS='replace',FORM='formatted')
WRITE (1, 460) 666   ! this will indicate an error unless replaced later
CLOSE (1)

! ----------------------------------------------------------------
! read in the ensemble forecast data for this lead
! ----------------------------------------------------------------
DO ilead = 0,fcst_len,fcst_output_int   
  PRINT *,'calling get_ens_fcstinfo, lead = ',ilead,cyyyymmddhh_in
  CALL get_ens_fcstinfo(ilead, cyyyymmddhh_in, cbasin_in, &
     cstormno_in,rundir_in, glvl_in, fcst_len_in ,hurricane_ens_fcst,iex)
  PRINT *,'returning from get_ens_fcstinfo'
END DO    

! --------------------------------------------------------------
! write to file
! --------------------------------------------------------------

! outfile = '/lfs1/projects/rtfim/FIMYENS/FIMwfm/ensplots/idl_fcst.dat'
PRINT *,'Program extract_atcf_FIM_GFS writing to ',TRIM(outfile)
OPEN (UNIT=1, FILE=outfile, STATUS='replace',FORM='formatted', IOSTAT=ios)
if (ios .ne. 0) then
   print *, 'error: ',ios,' opening ',outfile
endif

DO ilead = 0, fcst_len, fcst_output_int
  ileadidx = 1 + ilead/fcst_output_int
  WRITE (1, 460, iostat=ios) ileadidx
  if (ios .ne. 0) then
   print *, 'error: ',ios,' writing ',ileadidx
  endif
  460 FORMAT(i3)
  DO imem = 1, nmembers
    WRITE (1,461,iostat=ios) imem, hurricane_ens_fcst(ileadidx)%iflead, &
      hurricane_ens_fcst(ileadidx)%ctr_lat(imem),&
      hurricane_ens_fcst(ileadidx)%ctr_lon(imem),&
      hurricane_ens_fcst(ileadidx)%centralpressure(imem),&
      hurricane_ens_fcst(ileadidx)%windspeed(imem)
    461 FORMAT (i2,1x,i3,1x,4(f8.2,1x))
    if (ios .ne. 0) then
     print *, 'error: ',ios,' writing values'
    endif
!   print *, imem, hurricane_ens_fcst(ileadidx)%iflead, &
!     hurricane_ens_fcst(ileadidx)%ctr_lat(imem),&
!     hurricane_ens_fcst(ileadidx)%ctr_lon(imem),&
!     hurricane_ens_fcst(ileadidx)%centralpressure(imem),&
!     hurricane_ens_fcst(ileadidx)%windspeed(imem)
  END DO
END DO

! PRINT *,'end program extract_atcf.x'

CONTAINS

! ===================================================================

SUBROUTINE get_ens_fcstinfo(ilead_in, cyyyymmddhh_in, cbasin_in, &
  cstormno_in, rundir_in, glvl_in, fcst_len_in, hurricane_ens_fcst, iex)

! --------------------------------------------------------------------
! read ATCF hurricane ens forecast files.
! the chosen lead in hours (ilead).  Return the forecast information 
! for all in the storms that have been tracked for this lead time.
! --------------------------------------------------------------------

INTEGER, INTENT(IN)      :: ilead_in ! forecast lead, h
CHARACTER*10, INTENT(IN) :: cyyyymmddhh_in
CHARACTER*2, INTENT(IN)  :: cbasin_in
CHARACTER*2, INTENT(IN)  :: cstormno_in
CHARACTER*150, INTENT(IN) :: rundir_in
CHARACTER*2, INTENT(IN)  :: glvl_in
CHARACTER*3, INTENT(IN)  :: fcst_len_in

TYPE (ens_fcst), INTENT(OUT), DIMENSION(29)  :: hurricane_ens_fcst
LOGICAL, INTENT(OUT)     :: iex ! did the forecast files exist?

CHARACTER*150 :: infile
CHARACTER*112 cline
CHARACTER*2 cbasin
CHARACTER*2 cstormno
CHARACTER*2, DIMENSION(nmembers) :: cmem

INTEGER imemktr, ileady

DATA cmem /'01','02','03','04','05','06','07','08','09','10','11','12', '13',&
'14', '15', '16', '17', '18', '19', '20'/
           

! ---------------------------------------------------------------------
! check to make sure that all ensemble forecast members were computed;
! if not all 10, don't make plot...
! ---------------------------------------------------------------------

! PRINT *,' get_ens_fcstinfo: ------ ilead = ',ilead_in  
iex = .TRUE.
imemktr = 0
DO imem = 1,10
  CALL build_filename_fim(cyyyymmddhh_in, cmem(imem), rundir_in, glvl_in, fcst_len_in, infile)
  INQUIRE (file=infile,exist=iex)
  ! PRINT *,'infile: ',TRIM(infile)
  IF (iex) THEN
    imemktr = imemktr + 1
  !  PRINT *,TRIM(infile), " FOUND"
  ENDIF
END DO
DO imem = 11,20
    CALL build_filename_gfs(cyyyymmddhh_in, cmem(imem), gfs_rundir_in, infile)
    ! PRINT *,'infile: ',TRIM(infile)
    INQUIRE (file=infile,exist=iex)
    IF (iex) THEN
      imemktr = imemktr + 1
    !  PRINT *,TRIM(infile), " FOUND"
    ENDIF
END DO
  PRINT *,imemktr, ' member forecast files found'

ileadidx = 1 + ilead/fcst_output_int
!print *,ileadidx,ilead,trim(infile)
IF (imemktr .eq. nmembers) THEN   ! only process if all members available
  print *,'ALL MEMBERS AVAIL'

  ifound = 0
  DO imem = 1, nmembers
   ! CALL build_filename(cyyyymmddhh_in, cmem(imem), infile)
   ! CALL build_filename(cyyyymmddhh_in, cmem(imem), rundir_in, glvl_in, fcst_len_in, infile)
    IF (imem .LE. 10) THEN
       CALL build_filename_fim(cyyyymmddhh_in, cmem(imem), rundir_in, glvl_in, fcst_len_in, infile)
    ENDIF
    IF (imem .GT. 10) THEN
       CALL build_filename_gfs(cyyyymmddhh_in, cmem(imem), gfs_rundir_in, infile)
    ENDIF
    !PRINT *, TRIM(infile)
    OPEN (UNIT=1, FILE=infile, STATUS='old', FORM='formatted')
    DO
      READ (1,'(a112)', END=2000) cline
      ! PRINT *,'cline: ',cline
      READ (cline(31:33), '(i3)', ERR=234) ileady
      cstormno = cline(5:6)
      READ (cline(5:6), '(i2)') id
      cbasin = cline(1:2)

      ifoundtc = 0
      IF (cbasin .eq. cbasin_in .and. cstormno .eq. cstormno_in .and. &
      ileady .eq. ilead_in) THEN

        ! --------------------------------------------------------
        ! This member is tracking the storm.  Get the central pressure, 
        ! max wind speed, center's lat/lon for this member
        ! --------------------------------------------------------

        ifoundtc = 1
        READ (cline(36:38), '(i3)') ilat
          print *, " *******ifoundtc: ",ifoundtc, 'ilat: ',ilat
        IF (ilat .ne. 0) THEN
          IF (cline(39:39) .eq. 'N' .or. cline(39:39) .eq. 'n') THEN
            hurricane_ens_fcst(ileadidx)%ctr_lat(imem) = REAL(ilat)/10.
          ELSE
           hurricane_ens_fcst(ileadidx)%ctr_lat(imem) = - REAL(ilat)/10.
          END IF

          READ (cline(42:45), '(i4)') ilon
          IF (cline(46:46) .eq. 'E' .or. cline(46:46) .eq. 'e') THEN
            hurricane_ens_fcst(ileadidx)%ctr_lon(imem) = REAL(ilon)/10.
          ELSE
            hurricane_ens_fcst(ileadidx)%ctr_lon(imem) = 360. - REAL(ilon)/10.
          END IF

          IF (cline(54:57) .NE. ' -99') THEN
            READ (cline(54:57), '(i4)') imslp
            hurricane_ens_fcst(ileadidx)%centralpressure(imem) = REAL(imslp)
          ENDIF

          IF (cline(49:51) .NE. '***') THEN
            READ (cline(49:51), '(i3)') iwindkt  ! in knots
            hurricane_ens_fcst(ileadidx)%windspeed(imem) = REAL(iwindkt)*.514444
          ENDIF
          GOTO 2000
        END IF
      END IF
      234 CONTINUE
    END DO
    2000 CLOSE (1)
! print *, " mem: ",imem," *******windspeed: ",hurricane_ens_fcst(ileadidx)%windspeed(imem)
  END DO   ! imem = 1, 10

  hurricane_ens_fcst (ileadidx)%basin      = cline(1:2)
  hurricane_ens_fcst (ileadidx)%id         = idsv
  READ (cline(9:18), '(i10)') ibasetime
  hurricane_ens_fcst (ileadidx)%ibasetime  = ibasetime
  hurricane_ens_fcst (ileadidx)%iflead     = ilead

ENDIF ! imemktr = 10
print *,  imemktr

RETURN
END SUBROUTINE get_ens_fcstinfo
  
! ======================================================================

! SUBROUTINE build_filename(cyyyymmddhh, cmem, infile)
SUBROUTINE build_filename_fim(cyyyymmddhh, cmem, rundir_in, glvl_in, fcst_len_in, infile)

CHARACTER*10, INTENT(IN) :: cyyyymmddhh
CHARACTER*2, INTENT(IN) :: cmem     ! member number
CHARACTER*150, INTENT(IN) :: rundir_in
CHARACTER*2, INTENT(IN)  :: glvl_in
CHARACTER*3, INTENT(IN)  :: fcst_len_in
CHARACTER*150, INTENT(OUT) :: infile

CHARACTER*150 :: cdir
LOGICAL iex

cdir = TRIM(rundir_in) // '/tracker_0' // TRIM(cmem) // '/' // TRIM(fcst_len_in)
infile = TRIM(cdir)//'/track.'//TRIM(cyyyymmddhh)//'00.FE'// TRIM(cmem)

print*,'in build_filename_fim: reading ',infile
INQUIRE (file=infile,exist=iex)
! print*,'reading iex: ',iex

RETURN
END SUBROUTINE build_filename_fim


! ======================================================================
SUBROUTINE build_filename_gfs(cyyyymmddhh, cmem, rundir_in, infile)

CHARACTER*10, INTENT(IN) :: cyyyymmddhh
CHARACTER*2, INTENT(IN) :: cmem     ! member number
CHARACTER*150, INTENT(IN) :: rundir_in
CHARACTER*150, INTENT(OUT) :: infile

CHARACTER*150 :: cdir
LOGICAL iex

cdir = TRIM(rundir_in) 
infile = TRIM(cdir)//'/track.'//cyyyymmddhh//'.GE'//cmem
print*,'in build_filename_gfs: reading ',infile

INQUIRE (file=infile,exist=iex)

RETURN
END SUBROUTINE build_filename_gfs

! =============================================================

SUBROUTINE init_ens_fcst_structure(eforecast)
TYPE (ens_fcst), INTENT(OUT) :: eforecast

eforecast%basin = 'ZZ'
eforecast%id = -999

eforecast%ibasetime            = -999
eforecast%iflead               = -999

eforecast%ctr_lat(:)           = -999.99
eforecast%ctr_lon(:)           = -999.99
eforecast%centralpressure(:)   = -999.99
eforecast%windspeed(:)         = -999.99

RETURN
END SUBROUTINE init_ens_fcst_structure

! ===============================================================

END PROGRAM extract_atcf

