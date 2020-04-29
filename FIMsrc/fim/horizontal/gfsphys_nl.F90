module gfsphys_nl
! Define a namelist for all things related to the GFS physics. All variables
! come from namelist_def.f in the GFS physics directory.
! Contained routines initialize reals and integers to a bad value that will
! subsequently be tested, read the namelist, and print its contents.
! Logicals from namelist_def are not part of the namelist, but are set to
! their appropriate value.
  use machine, only: kind_phys
  use namelist_def
  use fimnamelist, only: ictm, isol, ico2, iaer, iaer_mdl, ialb, iems,  &
                          iovr_sw, iovr_lw, isubc_sw, isubc_lw,          &
                          n2dzhaocld, n3dzhaocld, n2dcldpdf, n3dcldpdf,  &
                          n3dflxtvd, fhswr, fhlwr

  implicit none

  private
  save
  public :: set_gfsphys_vars

contains

  subroutine set_gfsphys_vars ()
! Initialize namelist variables to "huge" so that we can detect those
! that have not been set in the namelist, but which must be.
! Then read the namelist.
! Also: Set logical variables to what FIM requires, since one can not
! test whether logical variables have been set by the user.

    use units, only: getunit, returnunit

    integer :: unitno
    integer :: ierr

! Set logicals to what they need to be. These are defined in namelist_def.f
! and NOT put in the namelist. But they could easily be added to the namelist
! if needed.
    crick_proof  = .false.   
    ccnorm       = .false. 
    norad_precip = .false. 
    ras          = .false.  ! False => call SAS
!   sashal       = .true.   ! jbao from gfs namelist
    lsswr        = .true.   ! jbao
    lslwr        = .true.   ! jbao
    ldiag3d      = .false.
    lssav        = .true.   !JFM false means FLUXR is not calculated in grrad,

    call check_gfsphys_vars ()

  end subroutine set_gfsphys_vars

  subroutine check_gfsphys_vars ()
! Check that user set all variables in namelist /gfsphys/
! Integers
    call checkvari (ictm, 'ictm')
    call checkvari (isol, 'isol')
    call checkvari (ico2, 'ico2')
    call checkvari (iaer, 'iaer')
    call checkvari (iaer_mdl, 'iaer_mdl')
    call checkvari (ialb, 'ialb')
    call checkvari (iems, 'iems')
    call checkvari (iovr_sw, 'iovr_sw')
    call checkvari (iovr_lw, 'iovr_lw')
    call checkvari (isubc_sw, 'isubc_sw')
    call checkvari (isubc_lw, 'isubc_lw')
    call checkvari (n2dzhaocld, 'n2dzhaocld')
    call checkvari (n3dzhaocld, 'n3dzhaocld')
    call checkvari (n2dcldpdf, 'n2dcldpdf')
    call checkvari (n3dcldpdf, 'n3dcldpdf')
    if (n3dcldpdf /= 0) then
! jbao, jr: In MPI mode, n3cldpdf /= 0 will not work properly since ipn=1
! jbao, jr: is used to initialize si_loc
      print*,'n3dcldpdf/=0 requires code mods to properly set si_loc in MPI mode'
      stop
    end if
      
    call checkvari (n3dflxtvd, 'n3dflxtvd')

! Reals
    call checkvarr (fhswr, 'fhswr')
    call checkvarr (fhlwr, 'fhlwr')
  end subroutine check_gfsphys_vars

  subroutine checkvari (var, name)
    integer, intent(in) :: var
    character(len=*), intent(in) :: name

    if (var == huge(var)) then
      write(6,*) 'check_gfsphys_vars:', name,' must be set in FIMnamelist'
      call flush(6)
      stop
    end if
  end subroutine checkvari

  subroutine checkvarr (var, name)
    real(kind_phys), intent(in) :: var
    character(len=*), intent(in) :: name

    if (var == huge(var)) then
      write(6,*) 'check_gfsphys_vars:', name,' must be set in FIMnamelist'
      call flush(6)
      stop
    end if
  end subroutine checkvarr
end module gfsphys_nl
