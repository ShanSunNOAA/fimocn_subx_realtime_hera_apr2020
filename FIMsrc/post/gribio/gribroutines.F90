!********************************************************************
!   SUBROUTINE initgrib
!
!   This subroutine initializes the grib table from the text file
!   and opens a grib output file for writing.
!
!   Syntax:
!     CALL initgrib(gribtable)
!     
! 
!    N. Wang, May 2007, initial verision.
!********************************************************************
      subroutine initgrib (gribtable)
        use grib_datastru
        implicit none

        character(len=*), intent(in) :: gribtable  ! rosetta stone of grib data
        character(len=128) dummy_line
        integer i
        integer :: ioerr

        open (unit=10, file=gribtable, form='formatted', action='read', iostat=ioerr)
        if (ioerr /= 0) then
          write(6,*)'initgrib: bad attempt to open ', trim(gribtable)
!TODO: pass back a bad return code
        end if

! skip header lines 

        do i=1,5
          read (unit=10, fmt='(a128)', iostat=ioerr) dummy_line
          if (ioerr /= 0) then
            write(6,*)'initgrib: bad attempt to skip line ', i, & 
                      ' of file ', trim(gribtable)
!TODO: pass back a bad return code
          end if
        end do    
 
! Read in the gribtable 
        do i=1,tbl_sz 
          read (unit=10, fmt=50, end=100, iostat=ioerr) &
               varnames(i),parm(i),ztype(i), iz1(i),iz2(i),itime_range(i),dscal(i),varabv(i)
 50       format(a60,      i3,3x,  i3,2x,    i4,2x, i4,2x, i3,5x,         i2,3x,   a7)

          if (ioerr /= 0) then
            write(6,*)'initgrib: bad attempt to read line ', i, & 
                      ' of file ', trim(gribtable)
!TODO: pass back a bad return code
          end if
!          print *, 'initgrib: varname: ',trim(varnames(i)),' parm: ',parm(i),' ztype: ',ztype(i), ' iz1: ',iz1(i),' iz2: ',iz2(i),' itime: ', itime_range(i),' dscal: ',dscal(i),'var: ',trim(varabv(i))
        end do 
 100    nvars_in_tbl = i - 1
        close (10)
          
        grib_lun = 0
!        CALL opengrib(gribfile)
           
      END SUBROUTINE initgrib



!************************************************************************************
!   SUBROUTINE writegrib
!
!   This subroutine writes out (a variable of) grid data in 3d volume.
!   It is called after the call to subroutine initgrib()
!
!   Syntax:
!     CALL writegrib(varname,igrid,nz,glvl,ct,data,date)
!       varname - variable name;
!       igrid -- grid id number  
!       nz -- sizes for each dimension when it is rectangular grid;
!       glvl -- grid refinement levels, only for icosahedral grid;
!       ct -- curve type, only for icosahedral grid;
!       data -- data array to be coded;
!       date -- year, julian day and hours for the model. 
!     
!   Note:
!     nx -- the total number of grid points, when ny = 0, indicating icosahedral grid;
! 
!    N. Wang, May 2007, initial verision.
!    N. Wang, Aug. 2007, added icosahedral model native grid.
!    N. Wang, Apr. 2014, changed argument list to pass-in grid number - igrid
!************************************************************************************
      SUBROUTINE writegrib(var_name,igrid,nz,glvl,ct,data,date, nt, tba, delta_t, nvlp, pres_hpa, atu, mh_tu)
        USE grib_datastru
        IMPLICIT NONE

        CHARACTER*(*) var_name, date
        INTEGER, INTENT(IN):: igrid, nz, nt, delta_t, glvl, ct, nvlp
        INTEGER, INTENT(IN)::tba
        REAL data(*)
        INTEGER, INTENT(IN):: pres_hpa(*)
        CHARACTER(2), INTENT(IN):: atu
        LOGICAL, INTENT(IN):: mh_tu

        INTEGER pds_sz, gds_sz, bms_sz, max_vlevels
        INTEGER :: nx, ny, dt_sc
        PARAMETER(pds_sz=28) 
        PARAMETER(gds_sz=50) 
        PARAMETER(max_vlevels=150) 

        CHARACTER*80 varname

! declare variables for grib table entries
        CHARACTER da_char(80)
        INTEGER id(pds_sz)
        INTEGER idx
        INTEGER yr,jday,hr, mo, day
        INTEGER type, bitl, pflag, gflag, comp, bflag, blen, da_int(1)
        INTEGER tot_len, err, npts, i, c_write,j
        INTEGER levels(max_vlevels)

        INTEGER bdsfl(9)
        DATA bdsfl/9*0/

        character, ALLOCATABLE :: grib_buf(:)
        REAL, ALLOCATABLE :: bms(:)
        INTEGER, ALLOCATABLE :: gds(:)
        INTEGER :: ierr, tu
  
        CALL gridid2mxmy(igrid, nx, ny)
!JR If the size of grib_buf becomes too big, Ning says the 10 could probably
!JR be made as small as 4. The whole key is compression, specified in the
!JR gribtable.
!TODO: Reduce the 10 to something smaller for high-resolution runs.
        ALLOCATE(grib_buf(nx * ny * nz * 10),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'writegrib: allocation of grib_buf failed'
          call flush(6)
          stop
        endif

        ALLOCATE(gds(gds_sz),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'writegrib: allocation of gds failed'
          call flush(6)
          stop
        endif

        ALLOCATE(bms(nx*ny),stat=ierr)
        if (ierr.ne.0) then
          write (*,'(a)') 'writegrib: allocation of bms failed'
          call flush(6)
          stop
        endif

! First search the grib table for variable entry
        varname = var_name
!        PRINT *,'writegrib:  var_name: ',trim(var_name),' varname: ',trim(varname)
        CALL touppercase(varname)
        idx = 0
        DO i = 1,nvars_in_tbl
          IF (varname .EQ. varabv(i)) THEN
            idx  = i
            EXIT 
          ENDIF
        ENDDO
        IF (idx .EQ. 0) THEN 
          PRINT *,'ERROR! Variable ',varname,' does not match any variables in GRIB table'
          PRINT *,'writegrib terminates!!'
          DEALLOCATE(grib_buf, gds, bms)
          RETURN
        ENDIF 

        READ(date,'(i2,i3,i2)') yr,jday,hr !for later
!        PRINT *,'Variable found ',varnames(idx),parm(idx),ztype(idx), iz1(idx),iz2(idx),itime_range(idx),dscal(idx)

! Second, assign values to id array, which is used to create PDS by w3fi68, 
! which in turn called bt w3fi72.
        id(1) =  28        ! number of bytes in PDS
        id(2) =  2         ! parm_table version
        id(3) =  59        ! id_center:  GSD = 59, NCEP = 7
        id(4) =  106       ! id_model: 105 = RUC, 106 = FIM
!        igrid = 255        ! 255 -- unknown grid, be defined in GDS.
!        IF (nx == 144) igrid = 228 !(144, 73)
!        IF (nx == 275) igrid = 244 !(275, 203) North Atlantic Hurricane scale
!        IF (nx == 151) igrid = 236 !(151, 113) Regional - CONUS scale 
!        IF (nx == 451) igrid = 130 !(451, 337) Regional - CONUS scale Hi-Res
!        IF (nx == 288) igrid = 45  !(288, 145)
!        IF (nx == 360) igrid = 3 !(360, 181)
!        IF (nx == 720) igrid = 4 !(720, 361)
!        IF (nx == 4320) igrid = 173 !(4320, 2610)
!        IF (nx == 2880) igrid = 174 !(2880, 1440)
!        IF (nx == 385) igrid = 219 !(385, 465)
!        IF (nx == 259) igrid = 201 !(259, 259)
!        IF (nx == 1760) igrid = 129 !(1760, 880)
!        IF (nx == 758) igrid = 83 !(758, 576)
        id(5) =  igrid     ! predefined grid and proj.
        id(6) = 1          ! gds_flag
        id(7) = 0          ! bms_flag 
      
        id(8) = parm(idx) ! indicator of param. and units. (Table2)
        id(9) = ztype(idx) ! indicator of level type.
        id(10) = iz1(idx)  ! value 1 of level
        id(11) = iz2(idx)  ! value 2 of level
  
        id(12) = yr        ! year of century
        CALL jday2moday(yr, jday, mo, day)
        id(13) = mo        ! month of year
        id(14) = day        ! day of the month
        id(15) = hr        ! hour of the day
        id(16) = 0         ! minute of hour 
        
        dt_sc = 1
        tu = 1
        IF (atu == 'mi') THEN
          tu = 0
        ELSE IF (atu == 'hr') THEN
          tu = 1
          IF (mh_tu) THEN
            IF (mod(delta_t,12) == 0) THEN
              tu = 12
              dt_sc = 12
            ELSE IF (mod(delta_t, 6) == 0) THEN 
              tu = 11
              dt_sc = 6
            ELSE IF (mod(delta_t, 3) == 0) THEN  
              tu =  10
              dt_sc = 3
            ENDIF
          ENDIF
        ELSE IF (atu == 'dy') THEN
          tu = 2
        ELSE IF (atu == 'mo') THEN
          tu = 3
        ENDIF

        id(17) = tu         ! forecast time unit: 0 minute, 1 hour, 2 day, 3 month
        IF (itime_range(idx) == 4) THEN
          id(18) = tba/dt_sc! time for the beginning of accumulation
          id(19) = nt/dt_sc ! time for the end of the accumulation
        ELSE
          id(18) = nt/dt_sc ! p1 period of time, 0 for initial analysis.
          id(19) = 0        ! p2 period of time, time interval between successive
        ENDIF               ! analyses, or forecasts undergoing averaging.
        id(20) = itime_range(idx)
        IF (id(20) /= 4) THEN
          id(20) = 10
        ENDIF

        id(21) = 0         ! number included in average
        id(22) = 0         ! number missing from average
        IF(yr.gt.95.or.yr.eq.0) then
          id(23) = 20
        ELSE
          id(23) = 21
        ENDIF
        id(24) = 0         ! sub_center
        id(25) = dscal(idx)! decimal scale factor

        type = 0           ! 0 floating point number, 1 integer number
        bitl = 0           ! computer determines the length for packing data
        pflag = 0          ! pds flag, 0:make pds from caller supplied array (id)
        gflag = 0          ! gds flag, 0: make gds based on igrid value
!        comp = 1           ! 0 earth oriented wind, 1 grid oriented wind 
        comp = 0           ! 0 earth oriented wind, 1 grid oriented wind 
        bflag = 1          ! bitmap flag, 0: make bitmap from caller supplied data
        blen = 0           ! length of bit map, 

! if data is icosahedral model grid      
        IF (ny == 1) THEN
          gds(1) = 0
          gds(2) = 255
          gds(3) = 12       ! grid type we give to our icosahedral hexagnal grid
          gds(4) = 1024     ! the "two dimension sizes for the grid"
          gds(5) = nx / 1024 + 1
          gds(6) = 26565    ! millidegrees for the lat. of the first anchor point 
          gds(7) = 10000    ! millidegrees for the lon. of the first anchor point 
          gds(8) = 0        ! N/A for our grid
          gds(9) = -26565    ! millidegrees for the lat. of the last anchor point
          gds(10) = 334000  ! millidegrees for the lon. of the last anchor point 
          gds(11) = 0       ! N/A for our grid
          gds(12) = 0       ! N/A for our grid
          gds(13) = ct      ! curve type
          gds(14) = glvl    ! grid refinement levels
          gflag = 1
        ENDIF
!       print *,'in writegrib before create_levels ztype: ',ztype(idx),' id(11): ', id(11) 
        CALL create_levels(ztype(idx),nz, levels, nvlp, pres_hpa, id(11)) 
        DO i = 1, nz
          id(11) = levels(i)
!         print *,'in writegrib after create_levels: id(11): ', id(11) 
          if(ny==1) then
            j = (i-1)*gds(4)*gds(5) + 1
          else
            j = (i-1)*nx*ny + 1
          endif
          CALL w3fi72(type,data(j),da_int,bitl,pflag,id,da_char,gflag,igrid,gds,comp,bflag, &
                      da_int,blen,bdsfl,npts,grib_buf,tot_len,err)
          IF (err /= 0) THEN
             CALL errormsg(err)
             EXIT
          ELSE
!JR Need an error check here?
            err = c_write(0, tot_len, grib_buf, grib_lun)
          END IF
        END DO

1111    continue
        DEALLOCATE(grib_buf, gds, bms)
        RETURN
      END SUBROUTINE writegrib

!********************************************************************
!   SUBROUTINE endgrib
!
!   This subroutine ends the writing of grib file.
!
!   Syntax:
!     CALL endgrib()
!     
! 
!    N. Wang, May 2007, initial verision.
!********************************************************************
      SUBROUTINE endgrib()
        USE grib_datastru
        IMPLICIT NONE

!        CALL closegrib()
 
      END SUBROUTINE endgrib


! opens an grib file using 'C' style open statement
      SUBROUTINE opengrib(gribfile)
        USE grib_datastru, ONLY: grib_lun 
        IMPLICIT NONE
        CHARACTER*80 gribfile
        INTEGER c_open

!JR Should there be an error return code check here?
        grib_lun = c_open(gribfile(1:LEN_TRIM(gribfile))//CHAR(0), 'w'//CHAR(0))
!        write(6,*) 'opengrib: opened gribfile ', trim(gribfile), ' to grib_lun=', grib_lun
         
      END SUBROUTINE opengrib
          

! closes an grib file using 'C' style close statement
      SUBROUTINE closegrib()
     
        USE grib_datastru, ONLY: grib_lun 
        IMPLICIT NONE
        INTEGER ret, c_close

        ret = c_close(grib_lun)

       END SUBROUTINE closegrib

      SUBROUTINE create_levels(ztype, nz, levels, nvlp, pres_hpa, iz2) 
        IMPLICIT NONE
        INTEGER ztype, nz, levels(*), nvlp, iz2
        INTEGER pres_hpa(*)
        
        INTEGER n_std_pres_levels, n_reg_pres_levels, n_reg_pres_levels2
        PARAMETER(n_std_pres_levels=17) 
        PARAMETER(n_reg_pres_levels2=111) 
        INTEGER i
        INTEGER std_pres_levels(n_std_pres_levels)
        DATA std_pres_levels/1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10/ 
        n_reg_pres_levels = nvlp

!       print *, 'in create_levels: ztype: ',ztype,' nz: ',nz,' iz2: ',iz2
        
        ! pressure levels
        IF (ztype == 100) THEN
          IF (nz == n_std_pres_levels) THEN
            DO i = 1, n_std_pres_levels
              levels(i) = std_pres_levels(i)
            END DO
          ELSE IF(nz == n_reg_pres_levels) THEN
            DO i = 1, n_reg_pres_levels
              levels(i) = pres_hpa(i)
            END DO
          ELSE
            DO i = 1, n_reg_pres_levels2
              levels(i) = (n_reg_pres_levels2 - i) * 10 
            END DO
          ENDIF
        ! m above ground - nz should be 1
        ELSE IF (ztype == 105) THEN
          DO i = 1, nz
            levels(i) = iz2
          END DO
        ! hybrid levels
        ELSE
          DO i = 1, nz
            levels(i) = i
          END DO
        END IF
            
      END SUBROUTINE create_levels

! Misc subroutines
       SUBROUTINE jday2moday(yr, jday, mo, day)
      
         IMPLICIT NONE
         INTEGER yr, jday, mo, day
         INTEGER jdate(13), jdate_ly(13), i
         LOGICAL leap_year

         DATA jdate/0, 31,59,90,120,151,181,212,243,273,304,334,365/
         DATA jdate_ly/0, 31,60,91,121,152,182,213,244,274,305,335,366/

         IF (mod(yr, 4) == 0 .AND. mod(yr, 100) /= 0) THEN
           leap_year = .TRUE.
         ELSE
           leap_year = .FALSE.
         ENDIF

         DO i = 2, 13
           IF (.NOT. leap_year .AND. jday <= jdate(i))  THEN
             mo = i - 1
             day = jday - jdate(mo) 
             EXIT
           END IF
           IF (leap_year .AND. jday <= jdate_ly(i))  THEN
	     mo = i - 1
             day = jday - jdate_ly(mo)
             EXIT
           END IF
         END DO
       END SUBROUTINE jday2moday


       SUBROUTINE touppercase(str)
         IMPLICIT NONE

         CHARACTER*80 str
         INTEGER len, i

         len = LEN_TRIM(str)
         DO i = 1, len
           IF (ichar(str(i:i)) >= ichar('a') .AND. ichar(str(i:i)) <=  ichar('z')) THEN
             str(i:i) = char(ichar(str(i:i)) - 32)
           END IF
         END DO

       END SUBROUTINE touppercase 

       SUBROUTINE tolowercase(str)
         IMPLICIT NONE

         CHARACTER*80 str
         INTEGER len, i

         len = LEN_TRIM(str)
         DO i = 1, len
           IF (ichar(str(i:i)) >= ichar('A') .AND. ichar(str(i:i)) <=  ichar('Z')) THEN
             str(i:i) = char(ichar(str(i:i)) + 32)
           END IF
         END DO
       END SUBROUTINE tolowercase 

  
       SUBROUTINE errormsg(err)
         IMPLICIT NONE

         INTEGER err
         IF (err == 1) THEN
           PRINT*, 'PDS flag error (It should be 1 or 0)'
         ELSEIF (err == 2) THEN
           PRINT*, 'GDS flag error (It should be 1 or 0)'
         ELSEIF (err == 3) THEN
           PRINT*, 'Error converting IEEE floating point number to IBM 370 floating point number'
         ELSEIF (err == 4) THEN
           PRINT*, 'Grid id not defined'
         ELSEIF (err == 5) THEN
           PRINT*, 'W3fi74 error: grid representation type not valid'
         ELSEIF (err == 6) THEN
           PRINT*, 'Grid too large for packer dimension arrays'
         ELSEIF (err == 7) THEN
           PRINT*, 'Length of bitmap not equal to size of the filed'
         ELSEIF (err == 8) THEN
           PRINT*, 'W3fi73 error:  bitmap values all zero'
         ELSEIF (err == 9) THEN
           PRINT*, 'W3fi75(58) error:  pack routine dynamic range overflow'
         ELSE
           PRINT*, 'Error code is ', err
         END IF
           PRINT*, 'Gribbing failed'
       
       END SUBROUTINE errormsg
