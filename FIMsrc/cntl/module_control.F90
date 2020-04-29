module module_control

  !***********************************************************
  ! This module specifies control variables for the global fim
  !   A. E. MacDonald  October 11, 2004
  !   J. LEE           September,  2005
  !***********************************************************

  use fimnamelist, only : readnl, nvl, glvl, SubdivNum,           &
                          max_vars, var_list, t1, t2, delta_t,    &
                          TotalTime, ArchvIntvl, ArchvTimeUnit,   &
                          isobaric_levels_file, curve,            &
                          nts, dophysics, enkfio_in, digifilt,    &
                          chem_opt, cu_physics, mp_physics,       &
                          ra_sw_physics, readrestart,yyyymmddhhmm,&
                          restart_freq, dt_reducer_numerator,     &
                          dt_reducer_denominator, FixedGridOrder, &
                          PrintDiagProgVars, PrintDiagNoise
  
  implicit none

  save

  !  MODEL GRID LEVEL: glvl
  !
  !  Grid levels for the globe are:
  !
  ! glvl Number of grid points Linear Scale (km)
  ! ----- --------------------- -----------------
  ! 0    12                    7,071
  ! 1    42                    3,779
  ! 2    162                   1,174
  ! 3    642                   891
  ! 4    2562                  446
  ! 5    10,242                223
  ! 6    40,962                111
  ! 7    163,842               56
  ! 8    655,362               28
  ! 9    2,621,442             14
  ! 10   10,485,762            7
  ! 11   41,943,042            3.5
  ! 12   167,772,162           1.75

  ! Define variables & set initial values. (Paramters set in following section.)

  character(len =80)  :: sst_dat              = 'NO_SUCH_FILE'               ! Name of SST file

  integer :: ArchvStep           =  1      ! archive interval in time steps
  integer :: ArchvStep0          =  0      ! steps prior to official start due to time lag
  integer :: ArchvStep6h         =  1      ! 6h archive interval in time steps
  integer :: CallPhysics                   ! Timestep interval to call physics
  integer :: CallRadiation                 ! Timestep interval to call radiation
  integer :: dospptdt            =  4      ! Timestep interval to apply SPPT 
  integer :: CallSST                       ! Timestep interval to call radiation
  integer :: filename_len                  ! max length of output filenames
  integer :: i                   =  0      !
  integer :: ipnDiagLocal        =  0      ! Local ipn value at which to print diagnostics
  integer :: ipnDiagPE           = -1      ! Processor on which ipnDiag resides.
  integer :: itsDFI              =  1      ! index of starting time step for DFI
  integer :: kbl                           ! # of thin layers in surface boundary lyr
  integer :: nabl                          ! # of adams bashforth levels
  integer :: nd                            ! # of directions (x,y)
  integer :: nip                           ! # of icosahedral points
  integer :: npp                           ! # of proximity points(max)
  integer :: nr                            ! # of rhombi
  integer :: ntra                =  4      ! # of tracers advected on small dt: 1=theta 2=qv 3=qw 4=O3
  integer :: ntrb                =  0      ! # of tracers advected on large dt: will include chemistry
  integer :: numphr                        ! # of time steps/hr
  integer :: nvar2d                        ! # of extra 2d diagnostic variables for output
  integer :: nvarp               =  7      ! # of isobaric variables: 1=height 2=temp 3=RH 4=U 5=V 6=vert.velocity,
                                           !                          7=cloud/hydrometeor condensate (g/g)
  integer :: nvarsig                       ! # of sigma variables: 1=pres 2=dpres 3=virtual temp 4=U 5=V 6=specific humidity
                                           !                        7=condensate (g/g) 8=ozone mixing ratio
  integer :: nvartracersig                 ! # of tracers for sigma (6,7,8)
  integer :: nvlp                =  0      ! # of isobaric vertical levels - ex. 1000-25 hPa
  integer :: nvlp1                         ! # of vertical levels ( = layers+1)
  integer :: nvlsig                        ! # of GFS sigma levels for direct output of FIM gridded data on GFS native levels
  integer :: nx                            ! rhombus x dimension
  integer :: ny                            ! rhombus y dimension
  integer :: pres_hpa(:)                   ! Holds isobaric vertical levels (hPa)
  integer :: pres_pa(:)                    ! Holds isobaric vertical levels (Pa)
  integer :: prev_date(4)                  ! time information for in-memory SST & ICE FRACTION slices
  integer :: next_date(4)                  ! time information for in-memory SST & ICE FRACTION slices
  integer :: step_per_unit                 ! time steps per ArchvTimeUnit
  integer :: TestDiagNoise                 ! Calculated test value for printing diagnostic gravity wave noise
  integer :: TestDiagProgVars              ! Calculated test value for printing diagnostic prognosis variables
  logical :: have_next_sst       = .false. ! Day=mid month & SSTs have been updated to next month?
  real :: dt                               ! model time step (seconds)
  real :: hrs_in_month           = 730.    ! length of month in hrs (= 24*365/12)
  integer, parameter :: numsecs=10         ! max.# of massflx diagno.transects
  character*40 :: thruflsecs(numsecs)=' '  ! transects for thruflow diagnostics

  ! Set parameter values

  parameter (            &
    filename_len  =  80, &
    kbl           =  2,  &
    nabl          =  3,  &
    nd            =  2,  &
    npp           =  6,  &
    nr            =  10, &
    nvar2d        =  7,  &
    nvarsig       =  8,  &
    nvartracersig =  3   &
    )

  ! Declare allocatables

  allocatable :: pres_hpa,pres_pa

contains

  subroutine control(quiet_arg)

    use module_wrf_control ,only : num_chem,num_moist,num_emis_ant,num_emis_vol,wrf_control ! to re-compute ntra if needed
    use units              ,only : getunit, returnunit

    ! arguments

    logical,optional,intent(in)   :: quiet_arg

    ! local variables

    logical                       :: quiet
    logical                       :: wrf_flag
    integer                       :: j, unitno, rsl, hr


    quiet = .false.
    if (present(quiet_arg)) then
      quiet = quiet_arg
    end if

    unitno = getunit ()
    if (unitno < 0) then
      print*,'module_control: getunit failed for input files. Stopping'
      stop
    end if

    call readnl

    ! Ensure against old-style input fields in POSTnamelist by checking
    ! for embedded spaces in the field name
      do j=1,max_vars
        if (index (trim (var_list(j)), ' ') /= 0) then
          write(*,'(a,i0,a,a)') 'module_control: var_list element ',j, &
                                ' contains an embedded blank character ', &
                                ' in POSTnamelist which is not allowed.'
          write(*,'(a,a)') 'The first bad input field is:', trim (var_list(j))
          write(*,'(a,a)') 'Perhaps you are using the old-style syntax of "var1 var2 var3 ..."', &
                           ' instead of the new syntax of "var1", "var2", "var3", ...?'
          STOP
        end if
      end do

    ! Set t1, t2, and delta_t that control the pop post-processing range
      if ( ( t1 < 0 ).or.( t2 < 0 ).or.( delta_t < 0 ) ) then
        write (*,'(a)') 'module_control: t1, t2 or delta_t not set in POSTnamelist'
        write (*,'(a)') 'module_control: t1, t2, and delta_t are being set to: '
        t1      = 0
        t2      = TotalTime
        delta_t = ArchvIntvl
      else
        write (*,'(a)') 'module_control: using t1, t2, delta_t set in POSTnamelist'
      end if
      write (6,*) 'module_control: t1, t2, delta_t = ', t1, t2, delta_t

    ! need to read file and do initialization here because this is needed by both fim and pop
    open (unitno, file="./"//isobaric_levels_file, form="formatted", status="old", err=80)
    print *,'module_control: opened: ', isobaric_levels_file
    read (unitno, *, err=81) nvlp

    if (.not.dophysics) then
      print *,'control: dophysics=.F. which means physics will NOT be called ever'
    end if

    IF (readrestart) THEN
       IF (enkfio_in) THEN
          readrestart=.FALSE.
          PRINT *,'For enkfio_in=.true. readrestart need to be .false.'
          PRINT *,'stopping'
          STOP
       END IF
    ENDIF

    if (readrestart) then
      if (digifilt) then
        print *,'control: restart does not yet work with digifilt=.true.'
        call flush(6)
        stop
      end if
    end if


    allocate (pres_hpa(nvlp))
    allocate (pres_pa(nvlp))
    read (unitno, *, err=82) pres_hpa
    close (unitno)
    call returnunit (unitno)

    do i=1,nvlp
      pres_pa(i) = pres_hpa(i) * 100
    end do

    print *,'*** nvlp: ', nvlp
    print '(a/(10i4))','pres (hPa):',pres_hpa
    print '(a/(10i7))','pres (pa):',pres_pa

    nip = 1
    do i = 1, glvl
      nip = nip * SubdivNum(i)
    enddo
    nip = 10 * nip * nip + 2

    rsl = 1
    do i = 1, glvl
      rsl = rsl * SubdivNum(i)
    enddo
    rsl =  INT(REAL(rsl) / 2.0)
    dt = 1.2*4800./REAL(rsl)

! returns .true. iff any WRF physics or chemistry is turned on
    wrf_flag = ( chem_opt /= 0 ).or.( cu_physics /= 0 ).or. &
               ( mp_physics /= 0 ).or.( ra_sw_physics /= 0 )

    if ((glvl>7 .or. nip==163842) .and.ArchvTimeUnit/='ts'.and.ArchvTimeUnit/='mi'.and.ArchvTimeUnit/='mo') then
	dt =  dt * dt_reducer_numerator / dt_reducer_denominator
    elseif ((nip==368642) .and.ArchvTimeUnit/='ts'.and.ArchvTimeUnit/='mi'.and.ArchvTimeUnit/='mo') then
	dt =  50.
    end if
    nvlp1             = nvl+1                 ! # of vertical levels (= layers+1)
    !dt               = 5760./2**(glvl-1)     ! model time step (seconds)
    numphr            = nint(3600./dt)        ! # of time steps/hr
    !dt               = 3600./float(numphr)
    nx                = 2**glvl               ! rhombus x dimension
    ny                = 2**glvl               ! rhombus y dimension
    num_chem = 0
    num_moist =3
    num_emis_ant =0
    num_emis_vol =0
    if (wrf_flag) then
      ! add WRF variables to ntr-dimensioned arrays, if needed
!      call wrf_control(nvarp,ntrb,num_chem,num_moist,num_emis_ant,num_emis_vol)
      call wrf_control(nvarp,ntrb)
    end if
    !ntra = ntr

    if (curve == 0) then ! The grid order is already IJ for curve=0.
      FixedGridOrder=.false.
    end if

    if (glvl>=5.and.abs(numphr-3600./dt).gt.0.02) then
      print *,' number of time step in one hour is ',3600./dt, ' needs to be an integer'
      stop 'wrong: number of time step in one hour needs to be an integer'
    end if

    read(yyyymmddhhmm(9:10),'(i2)') hr

    if (ArchvTimeUnit == 'ts') then
      step_per_unit = 1
    else if (ArchvTimeUnit == 'mi') then
      step_per_unit = nint(60./dt+0.5)
      dt = 60./step_per_unit
    else if (ArchvTimeUnit == 'hr') then
      step_per_unit = nint(3600./dt)
      ArchvStep6h  = 6 * step_per_unit
      if (ArchvIntvl==24) then	! 24hr mean goes from 0z to 0z
        if (hr .ge. 12) then
          ArchvStep0 = nint((24-hr)*3600./dt)
        else
          ArchvStep0 = nint(-hr*3600./dt)
        end if
      end if
      print *,'hr =',hr,' ArchvStep0=',ArchvStep0
    else if (ArchvTimeUnit == 'dy') then
      step_per_unit = nint(86400./dt)
      ArchvStep6h  = step_per_unit/4.
    else if (ArchvTimeUnit == 'mo') then
      step_per_unit = nint(hrs_in_month*3600./dt)
      ArchvStep6h  = step_per_unit/hrs_in_month
    else
      write (*,'(a,a)') 'ERROR in module_control unrecognized output time unit: ',ArchvTimeUnit
      stop
    end if
!   nts          = TotalTime  * step_per_unit              ! convert from time to time step
    nts          = TotalTime  * step_per_unit + ArchvStep0 ! integrate TotalTime after ArchvStep0
    ArchvStep    = ArchvIntvl * step_per_unit
    restart_freq = restart_freq*step_per_unit

    if      (PrintDiagProgVars == 0) then
      TestDiagProgVars = 1
    else if (PrintDiagProgVars <  0) then
      TestDiagProgVars = 2*nts
    else
      TestDiagProgVars = PrintDiagProgVars*numphr
    end if

    if      (PrintDiagNoise == 0) then
      TestDiagNoise = 1
    else if (PrintDiagNoise <  0) then
      TestDiagNoise = 2*nts
    else
      TestDiagNoise = PrintDiagNoise*numphr
    end if

    return

80  print *,'module_control: error reading isobaric_levels_file: COULD NOT OPEN: ',isobaric_levels_file,' program aborted'
    call flush(6)
    stop

81  print *,'module_control: error reading isobaric_levels_file - nvlp -  program aborted'
    call flush(6)
    stop

82  print *,'module_control: error reading isobaric_levels_file - pres_ha - program aborted'
    call flush(6)
    stop

  end subroutine control

end module module_control
