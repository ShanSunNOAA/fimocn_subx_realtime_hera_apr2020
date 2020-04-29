module nc
  use netcdf
  implicit none

! dimensions
  integer :: ncid = -1
  integer :: lon_dimid = -1
  integer :: lat_dimid = -1
  integer :: vert_dimid = -1
! dimension variables
  integer :: lonid = -1
  integer :: latid = -1
  integer :: vertid = -1
! variables
  integer :: hsid = -1
  integer :: psid = -1
  integer :: tvid = -1
  integer :: qid = -1
  integer :: uid = -1
  integer :: vid = -1
  integer :: ozid = -1
  integer :: clwid = -1

contains

  subroutine initnc (sanlfile, outfile, imax, jmax, nvp)
    character(len=*), intent(in) :: sanlfile
    character(len=*), intent(in) :: outfile
    integer, intent(in) :: imax, jmax, nvp

    integer :: dimids(3)
    integer :: lev(nvp)
    integer :: old_mode      ! returned from nf90_set_fill (unused)

    call handle_err (nf90_create (trim(outfile), &
                                  cmode=or(NF90_CLOBBER, NF90_64BIT_OFFSET), &  ! netcdf3
                                  ncid=ncid))
! Set nofill for efficiency
    call handle_err (nf90_set_fill (ncid, NF90_NOFILL, old_mode))

! Define dimensions
    call handle_err (nf90_def_dim (ncid, 'lon', imax, lon_dimid))
    call handle_err (nf90_def_dim (ncid, 'lat', jmax, lat_dimid))
    call handle_err (nf90_def_dim (ncid, 'lev', nvp, vert_dimid))
! Dimension variables
    call handle_err (nf90_def_var (ncid, 'lon',  NF90_DOUBLE, lon_dimid, lonid))
    call handle_err (nf90_def_var (ncid, 'lat',  NF90_DOUBLE, lat_dimid, latid))

    call handle_err (nf90_def_var (ncid, 'lev', NF90_DOUBLE, vert_dimid, vertid))
    call handle_err (nf90_put_att (ncid, vertid, 'long_name', 'layer index'))
    call handle_err (nf90_put_att (ncid, vertid, 'units',     '1'))

! Attributes of dimension variables
    call handle_err (nf90_put_att (ncid, lonid, 'long_name', 'longitude'))
    call handle_err (nf90_put_att (ncid, lonid, 'units', 'degrees_east'))

    call handle_err (nf90_put_att (ncid, latid, 'long_name', 'latitude'))
    call handle_err (nf90_put_att (ncid, latid, 'units', 'degrees_north'))

! Variables from spectral data
    dimids(1) = lon_dimid
    dimids(2) = lat_dimid
    dimids(3) = vert_dimid

    call handle_err (nf90_def_var (ncid, 'hs', NF90_FLOAT, dimids(1:2), hsid))
    call handle_err (nf90_put_att (ncid, hsid, 'long_name', 'sfc height'))
    call handle_err (nf90_put_att (ncid, hsid, 'units',     'UNK'))

    call handle_err (nf90_def_var (ncid, 'ps', NF90_FLOAT, dimids(1:2), psid))
    call handle_err (nf90_put_att (ncid, psid, 'long_name', 'sfc pressure'))
    call handle_err (nf90_put_att (ncid, psid, 'units',     'MB'))

    call handle_err (nf90_def_var (ncid, 'tv', NF90_FLOAT, dimids, tvid))
    call handle_err (nf90_put_att (ncid, tvid, 'long_name', 'virtual temperature'))
    call handle_err (nf90_put_att (ncid, tvid, 'units',     'degrees K'))

    call handle_err (nf90_def_var (ncid, 'q', NF90_FLOAT, dimids, qid))
    call handle_err (nf90_put_att (ncid, qid, 'long_name', 'specific humidity'))
    call handle_err (nf90_put_att (ncid, qid, 'units',     'g/kg?'))

    call handle_err (nf90_def_var (ncid, 'u', NF90_FLOAT, dimids, uid))
    call handle_err (nf90_put_att (ncid, uid, 'long_name', 'eastward wind'))
    call handle_err (nf90_put_att (ncid, uid, 'units',     'm/s'))

    call handle_err (nf90_def_var (ncid, 'v', NF90_FLOAT, dimids, vid))
    call handle_err (nf90_put_att (ncid, vid, 'long_name', 'northward wind'))
    call handle_err (nf90_put_att (ncid, vid, 'units',     'm/s'))

    call handle_err (nf90_def_var (ncid, 'ozone', NF90_FLOAT, dimids, ozid))
    call handle_err (nf90_put_att (ncid, ozid, 'long_name', 'ozone conc'))
    call handle_err (nf90_put_att (ncid, ozid, 'units',     'UNK'))

    call handle_err (nf90_def_var (ncid, 'clw', NF90_FLOAT, dimids, clwid))
    call handle_err (nf90_put_att (ncid, clwid, 'long_name', 'cloud condensate'))
    call handle_err (nf90_put_att (ncid, clwid, 'units',     'UNK'))

! Done defining metada. End define mode so we can now write data
    call handle_err (nf90_enddef (ncid))
  end subroutine initnc

  subroutine handle_err (ret)
    integer, intent(in) :: ret

    if (ret /= NF90_NOERR) then
      write(6,*) nf90_strerror (ret)
      stop 999
    end if
  end subroutine handle_err
end module nc

