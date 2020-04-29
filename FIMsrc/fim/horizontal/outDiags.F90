module module_outDiags
  use module_constants,        only: deg_lat,deg_lon,thetac
  USE module_sfc_variables,    only: slmsk2d
  use module_outqv,            only: outqv
  use module_outqv_wsp,        only: outqv_wsp
  use module_outqv_mn,         only: outqv_mn
  use module_outqv_mn_lat,     only: outqv_mn_lat
  use module_outqv_mn_lat_abs, only: outqv_mn_lat_abs
  use module_outqv_mn_theta_abs, only: outqv_mn_theta_abs
  use module_outqv_mn_lat_land,only: outqv_mn_lat_land
  use module_out4d_mn,         only: out4d_mn
  use module_control,          only: nip,dt
  use fimnamelist,             only: ArchvTimeUnit

  implicit none

contains
  subroutine outDiags (its, nvl, nvlp1, ntra, pr3d, ph3d, tr3d, rn2d, rc2d, &
                       sdot, dp3d, us3d, vs3d, rh3d, rlds, rsds, hf2d, qf2d)
! External variable declarations:
    integer,intent(IN) :: its,nvl,nvlp1,ntra
!sms$distribute (dh,1) begin
    real   ,intent(IN) :: rsds(nip),rlds(nip)
    real   ,intent(IN) :: rn2d(nip),rc2d(nip)
    real   ,intent(IN) :: hf2d(nip),qf2d(nip)
!sms$distribute end
!sms$distribute (dh,2) begin
    real   ,intent(IN) :: us3d(nvl  ,nip),vs3d(nvl,nip),dp3d(nvl,nip)
    real   ,intent(IN) :: pr3d(nvlp1,nip)
    real   ,intent(IN) :: ph3d(nvlp1,nip)
    real   ,intent(IN) :: tr3d(nvl    ,nip,ntra)
    real   ,intent(IN) :: rh3d(nvl  ,nip)
    real   ,intent(IN) :: sdot(nvlp1,nip)	! mass flux across interfaces, sdot*(dp/ds)
!sms$distribute end

!factors for changing units for outqv range, max/min
!water vapor - g/g to g/kg
    real, parameter :: fact_qv = 1.e3
    real, parameter :: fact_qc = 1.e3
!pressure - Pa to hPa
    real, parameter :: fact_pres = 1.e-2
    real, parameter :: fact    = 1.
    real :: time
    integer :: LB
    real, external :: its2time

    time = its2time(its)

    write (6,*) 'thetac - outDiag',thetac(1)
    write (6,*) 'thetac - outDiag',thetac(2)

    write(6,100) its, time, ArchvTimeUnit
100 format ('outDiags: outqv* for its=', i8, ' time=', f9.2, ' ArchvTimeUnit=', a2)

    call outqv_mn     (pr3d, nvlp1,                   fact_pres, 'Pressure')
    call outqv_mn     (ph3d, nvlp1,                   1.,        'Geopotential height')
    call outqv        (ph3d, nvlp1, deg_lat, deg_lon, fact,      'Geopotential height')

    call out4d_mn     (tr3d, nvl,   ntra, 1., 1, 'THETA')
    call out4d_mn     (rn2d, 1,     1,    1., 1, 'Precip')

    call outqv_mn_lat (rn2d, 1,     deg_lat, deg_lon, 30., fact, 'Precip')
    call outqv        (rn2d, 1,     deg_lat, deg_lon,      fact, 'Precip')
    call out4d_mn     (rc2d, 1,     1,    1., 1,                 'Precip-conv')
    call outqv_mn_lat (rc2d, 1,     deg_lat, deg_lon, 30., fact, 'Precip-conv')
    call outqv        (rc2d, 1,     deg_lat, deg_lon,      fact, 'Precip-conv')
    call out4d_mn     (sdot, nvlp1, 1,    1., 1,                 'Vertical velocity')
    call outqv_mn_lat (sdot, nvlp1, deg_lat, deg_lon, 30., fact, 'Vertical velocity')
    call outqv        (sdot, nvlp1, deg_lat, deg_lon,      fact, 'Vertical velocity')
    call outqv_mn_lat_abs (sdot, nvlp1, deg_lat, deg_lon, 30.,fact, 'Abs Vertical velocity')

!sms$ignore begin
    LB = LBOUND(tr3d,2)
!sms$ignore end
    call outqv_mn_theta_abs &
        (sdot,tr3d(1,LB,1),ntra,thetac,nvlp1, nvl, deg_lat, deg_lon, 0.2,fact, 'Abs Vertical velocity')
    write(6,*)'Water vapor - tr(2) fact=', fact_qv
    call outqv        (tr3d(1,LB,2), nvl, deg_lat, deg_lon, fact_qv,   'Water vapor - tr(2)')

    write(6,*)'Cloud water - tr(3) fact=', fact_qc
    call outqv        (tr3d(1,LB,3), nvl, deg_lat, deg_lon, fact_qc,   'Cloud water - tr(3)')
    call outqv        (dp3d,         nvl, deg_lat, deg_lon, fact_pres, 'DP3d')
    call outqv        (tr3d(1,LB,1), nvl, deg_lat, deg_lon, fact,      'pot temp')

    write(6,*)'Water vapor - tr(2) fact=', fact_qv
    call outqv_mn     (tr3d(1,LB,2), nvl, fact_qv, 'Water vapor - tr(2)')

    write(6,*)'Cloud water - tr(3) fact=', fact_qc
    call outqv_mn     (tr3d(1,LB,3), nvl, fact_qc, 'Cloud water - tr(3)')

    write(6,*)'Cloud water - tr(3) fact=', fact_qc
    call outqv_mn_lat (tr3d(1,LB,3), nvl,   deg_lat, deg_lon, 30., fact_qc, 'Cloud water - tr(3)')
    call outqv_mn     (tr3d(1,LB,1), nvl, fact,                             'Pot temp')
    call outqv_mn     (us3d,         nvl, 1.,                               'Zonal wind')
    call outqv        (us3d,         nvl, deg_lat, deg_lon, fact,           'Zonal wind')
    call outqv_mn     (vs3d,         nvl, 1.,                               'Meridional wind')
    call outqv        (vs3d,         nvl, deg_lat, deg_lon, fact,           'Meridional wind')

!  Call for wind speed
    call outqv_wsp    (us3d, vs3d,   nvl, deg_lat, deg_lon, fact,           'Wind speed')
    call outqv_mn     (rh3d,         nvl, 1.,                               'Relative humidity')
    call outqv_mn     (rlds,         1,   1.,                               'Longwave')
    call outqv_mn     (rsds,         1,   1.,                               'Shortwave')
    call outqv_mn     (hf2d,         1,   1.,                               'Sensible heat flux')
    call outqv        (hf2d,         1,   deg_lat, deg_lon, 1.,             'Sensible heat flux')
    call outqv_mn_lat_land (hf2d, 1, deg_lat, deg_lon, 30., slmsk2d, 1, 1., 'Sensible heat flux')
    call outqv_mn     (qf2d,         1,   1.,                               'Latent heat flux')
    call outqv        (qf2d,         1,   deg_lat, deg_lon, 1.,             'Latent heat flux')
    call outqv_mn_lat_land (qf2d, 1, deg_lat, deg_lon, 30., slmsk2d, 1, 1., 'Latent heat flux')
    
    return
  end subroutine outDiags
end module module_outDiags
