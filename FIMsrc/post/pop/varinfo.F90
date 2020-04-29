!-------------------------------------------------------------------------
! This file contsins the subroutines to retrieve variable information
!
! Adapted from module fimnc, N.W. 02-19-2013.
!-------------------------------------------------------------------------
!#define VARIABLE_LEVELS
#define ERROR_CHECK(a) call IO_DEBUG(a, __LINE__, __FILE__)
      MODULE varinfo
        IMPLICIT NONE 

        INTEGER, PARAMETER :: dummy_com = 999 ! Dummy MPI communicator
        INTEGER, PARAMETER :: WRF_REAL = 104
        INTEGER, PARAMETER :: char_len = 512

        INTEGER :: nvls, STATUS
        CHARACTER(len=char_len), dimension(3) :: dim_names
        INTEGER, DIMENSION(4) :: domain_start = (/1, 1, 1, 1/)
        INTEGER, DIMENSION(4) :: domain_end
        CHARACTER(len=char_len)  :: sysdepinfo = " "

        INTEGER :: dim_allo_set = 0
        REAL, ALLOCATABLE :: lons(:)
        REAL, ALLOCATABLE :: lats(:)
        REAL, ALLOCATABLE :: levels(:)
        REAL, ALLOCATABLE :: times(:)

      CONTAINS

! subroutine to obtain the variable info through variable number
        SUBROUTINE var_info(var_name, var_description, units, nlevels, nvlp)
          CHARACTER(len=*)  :: var_name, var_description, units
          INTEGER :: nlevels, nvlp

          SELECT CASE (var_name)
      
! 3D variables for dynamics 
          CASE ("us3D")
          var_description = "U wind"
          units = "meter/second"
          nlevels = nvls 

          CASE ("vs3D")
          var_description = "V wind"
          units = "meter/second"
          nlevels = nvls

          CASE ("dp3D")
          var_description = "3-D Delta Pressure"
          units = "pascals"
          nlevels = nvls

          CASE ("pr3D")
          var_description = "3-D Pressure"
          units = "pascals"
          nlevels = nvls + 1

          CASE ("mp3D")
          var_description = "Montgomery Potential"
          units = "meter2/second2"
          nlevels = nvls

          CASE ("th3D")
          var_description = "Potential Temperature"
          units = "kelvin"
          nlevels = nvls

          CASE ("ph3D")
          var_description = "Geo Potential"
          units = "meter2/second2"
          nlevels = nvls + 1

          CASE ("qv3D")
          var_description = "Specific humidity"
          units = "non-dimensional"
          nlevels = nvls

          CASE ("rh3D")
          var_description = "Relative humidity"
          units = "percentage"
          nlevels = nvls

          CASE ("vr3D")
          var_description = "Vorticity"
          units = "1/second"
          nlevels = nvls

          CASE ("ws3D")
          var_description = "Omega"
          units = "Micro-bar/second"
          nlevels = nvls

          CASE ("tk3D")
          var_description = "Temperature"
          units = "Kelvin"
          nlevels = nvls

          CASE ("td3D")
          var_description = "Dew-point temperature"
          units = "Kelvin"
          nlevels = nvls

! 3D variables for chemistry
          CASE ("d5st")
          var_description = "dust particles bin 3"
          units = "ug/kg"
          nlevels = nvls

          CASE ("d4st")
          var_description = "dust particles bin 3"
          units = "ug/kg"
          nlevels = nvls

          CASE ("d3st")
          var_description = "dust particles bin 3"
          units = "ug/kg"
          nlevels = nvls

          CASE ("d2st")
          var_description = "dust particles bin 2"
          units = "ug/kg"
          nlevels = nvls

          CASE ("d1st")
          var_description = "dust particles bin 1"
          units = "ug/kg"
          nlevels = nvls

          CASE ("s1ea")
          var_description = "seasalt particles bin 1"
          units = "ug/kg"
          nlevels = nvls

          CASE ("s2ea")
          var_description = "seasalt particles bin 2"
          units = "ug/kg"
          nlevels = nvls

          CASE ("s3ea")
          var_description = "seasalt particles bin 3"
          units = "ug/kg"
          nlevels = nvls

          CASE ("s4ea")
          var_description = "seasalt particles bin 4"
          units = "ug/kg"
          nlevels = nvls

          CASE ("pmsa")
          var_description = "msa"
          units = "ppm"
          nlevels = nvls

          CASE ("dms1")
          var_description = "dms"
          units = "ppm"
          nlevels = nvls

          CASE ("pso2")
          var_description = "so2"
          units = "ppm"
          nlevels = nvls

          CASE ("sulf")
          var_description = "sulfate"
          units = "ppm"
          nlevels = nvls

          CASE ("pp25")
          var_description = "other primary pm25 +volcanic ash "
          units = "ug/Kg"
          nlevels = nvls

          CASE ("pp10")
          var_description = "other primary pm10 +volcanic ash"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("obc1")
          var_description = "hydrophobic organic carbon"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("obc2")
          var_description = "hydrophillic organic carbon"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("pbc1")
          var_description = "hydrophobic black carbon"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("pbc2")
          var_description = "hydrophillic black carbon"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("ash1")
          var_description = "volcanic ash size bin 1"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("ash2")
          var_description = "volcanic ash size bin 2"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("ash3")
          var_description = "volcanic ash size bin 3"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("ash4")
          var_description = "volcanic ash size bin 4"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("c13D")
          var_description = "radioactive tracer 1, explosive emissions with height"
          units = "?"
          nlevels = nvls

          CASE ("c23D")
          var_description = "radioactive tracer 2, linear emissions with height"
          units = "?"
          nlevels = nvls

          CASE ("oc1P")
          var_description = "hydrophobic organic carbon"
          units = "ug/Kg"
          nlevels = nvlp

          CASE ("oc2P")
          var_description = "hydrophobic organic carbon"
          units = "ug/Kg"
          nlevels = nvlp

          CASE ("bc1P")
          var_description = "hydrophobic black carbon"
          units = "ug/Kg"
          nlevels = nvlp

          CASE ("bc2P")
          var_description = "hydrophobic black carbon"
          units = "ug/Kg"
          nlevels = nvlp

          CASE ("so2P")
          var_description = "so2"
          units = "ppm"
          nlevels = nvlp

          CASE ("slfP")
          var_description = "sulfate"
          units = "ppm"
          nlevels = nvlp

          CASE ("d1sP")
          var_description = "dust particles"
          units = "ppm"
          nlevels = nvlp

          CASE ("d2sP")
          var_description = "dust particles"
          units = "ppm"
          nlevels = nvlp

          CASE ("d3sP")
          var_description = "dust particles"
          units = "ppm"
          nlevels = nvlp

          CASE ("d4sP")
          var_description = "dust particles"
          units = "ppm"
          nlevels = nvlp

          CASE ("d5sP")
          var_description = "dust particles"
          units = "ppm"
          nlevels = nvlp

          CASE ("s1sP")
          var_description = "seasalt"
          units = "ppm"
          nlevels = nvlp

          CASE ("s2sP")
          var_description = "seasalt"
          units = "ppm"
          nlevels = nvlp

          CASE ("s3sP")
          var_description = "seasalt"
          units = "ppm"
          nlevels = nvlp

          CASE ("s4sP")
          var_description = "seasalt"
          units = "ppm"
          nlevels = nvlp

          CASE ("dmsP")
          var_description = "dms"
          units = "ppm"
          nlevels = nvlp

          CASE ("msaP")
          var_description = "msa"
          units = "ppm"
          nlevels = nvlp

          CASE ("p25P")
          var_description = "other primary pm25"
          units = "ug/Kg"
          nlevels = nvlp

          CASE ("p10P")
          var_description = "other primary pm25"
          units = "ug/Kg"
          nlevels = nvlp

! 3D variables for physics 
          CASE ("qw3D")
          var_description = "liquid cloud mixing ratio"
          units = "Kg/Kg"
          nlevels = nvls

          CASE ("hl3D")
          var_description = "long-wave heating rate"
          units = "Kelvin/second"
          nlevels = nvls

          CASE ("hs3D")
          var_description = "short-wave heating rate"
          units = "Kelvin/second"
          nlevels = nvls

          CASE ("oz3D")
          var_description = "ozone"
          units = "Kg/Kg"
          nlevels = nvls

          CASE ("ar3D")
          var_description = "aerosol"
          units = "Kg/Kg"
          nlevels = nvls

          CASE ("cf3D")
          var_description = "cloud fraction"
          units = "percent(%)"
          nlevels = nvls

          CASE ("st3D")
          var_description = "st"
          units = "kg/meter^2"
          nlevels = 4

          CASE ("sm3D")
          var_description = "sm"
          units = "kg/meter^2"
          nlevels = 4

!2D variables for physics
          CASE ("rn2D")
          var_description = "rainfall(accumulated total)"
          units = "millimeter"
          nlevels = 1

          CASE ("rc2D")
          var_description = "rainfall(accumulated conv.)"
          units = "millimeter"
          nlevels = 1

          CASE ("rg2D")
          var_description = "rainfall(accumulated large-scale)"
          units = "millimeter"
          nlevels = 1

          CASE ("sn2D")
          var_description = "snowfall (run accumulated total)"
          units = "millimeter"
          nlevels = 1

          CASE ("sa2D")
          var_description = "snow water equivalent on ground"
          units = "meter"
          nlevels = 1

          CASE ("r12D")
          var_description = "precipitation (total, since last output)"
          units = "millimeter"
          nlevels = 1

          CASE ("s12D")
          var_description = "snow precip (water eq, total, since last output)"
          units = "millimeter"
          nlevels = 1

          CASE ("r22D")
          var_description = "precipitation (conv, since last output)"
          units = "millimeter"
          nlevels = 1

          CASE ("r32D")
          var_description = "precipitation (large-scale, since last output)"
          units = "millimeter"
          nlevels = 1

          CASE ("r42D")
          var_description = "rainfall(6-hour accumulated total)"
          units = "millimeter"
          nlevels = 1

          CASE ("r52D")
          var_description = "rainfall(6-hour accumulated conv.)"
          units = "millimeter"
          nlevels = 1

          CASE ("r62D")
          var_description = "rainfall(6-hour accumulated large-scale)"
          units = "millimeter"
          nlevels = 1

          CASE ("s62D")
          var_description = "snow precip (water eq, 6-hour accumulated)"
          units = "millimeter"
          nlevels = 1

          CASE ("pw2D")
          var_description = "precipitable water"
          units = "millimeter"
          nlevels = 1

          CASE ("pq2D")
          var_description = "vertically integrated cloud condensate"
          units = "millimeter"
          nlevels = 1

          CASE ("ts2D")
          var_description = "skin temperature"
          units = "deg. Kelvin"
          nlevels = 1

          CASE ("us2D")
          var_description = "friction velocity"
          units = "meter/sec."
          nlevels = 1

          CASE ("u12D")
          var_description = "u-component of 10m wind"
          units = "meter/sec."
          nlevels = 1

          CASE ("v12D")
          var_description = "v-component of 10m wind"
          units = "meter/sec."
          nlevels = 1

          CASE ("t22D")
          var_description = "2m temperature"
          units = "Kelvin"
          nlevels = 1

          CASE ("q22D")
          var_description = "2m specific humidity"
          units = "kg/kg"
          nlevels = 1

          CASE ("hfss")
          var_description = "sensible heat flux"
          units = "W/m**2"
          nlevels = 1

          CASE ("hfls")
          var_description = "latent heat flux"
          units = "W/m**2"
          nlevels = 1

          CASE ("rsds")
          var_description = "radiation shortwave downward at surface"
          units = "W/m**2"
          nlevels = 1

          CASE ("rlds")
          var_description = "radiation longwave downward at surface"
          units = "W/m**2"
          nlevels = 1

          CASE ("rsus")
          var_description = "radiation shortwave upward at surface"
          units = "W/m**2"
          nlevels = 1

          CASE ("rlus")
          var_description = "radiation longwave upward at surface"
          units = "W/m**2"
          nlevels = 1

          CASE ("ms2D")
          var_description = "mean sea level pressure"
          units = "Pa"
          nlevels = 1

          CASE ("ct2D")
          var_description = "cloud top height"
          units = "m"
          nlevels = 1

          CASE ("cb2D")
          var_description = "cloud base height"
          units = "m"
          nlevels = 1

          CASE ("io2D")
          var_description = "integrated organic carbon"
          units = "ug/kg"
          nlevels = 1

          CASE ("ib2D")
          var_description = "integrated black carbon"
          units = "ug/kg"
          nlevels = 1

          CASE ("id2D")
          var_description = "integrated fine dust"
          units = "ug/kg"
          nlevels = 1

          CASE ("is2D")
          var_description = "integrated sulfate"
          units = "ppm"
          nlevels = 1

          CASE ("ia2D")
          var_description = "integrated PM25"
          units = "ug/m3"
          nlevels = 1

          CASE ("ao2D")
          var_description = "Aerosol Optical Depth"
          units = "unitless"
          nlevels = 1

          CASE ("iash")
          var_description = "Vertically integrated volcanic ash"
          units = "ug/kg"
          nlevels = 1

          CASE ("fl2D")
          var_description = "fallout"
          units = "?"
          nlevels = 1

          CASE ("rp2D")
          var_description = "relative humidity with respect to precipitable water"
          units = "%"
          nlevels = 1

          CASE ("rlut")
          var_description = "outgoing LW radiation at top of atmosphere"
          units = "W/m**2"
          nlevels = 1

          CASE ("w080")
          var_description = "wind speed at 80 meters fixed height"
          units = "m/s"
          nlevels = 1

! FIM output variables at standard pressure levels
          CASE ("hgtP")
          var_description = "height at pressure levels"
          units = "meter2/second2"
          nlevels = nvlp

          CASE ("tmpP")
          var_description = "temperature at pressure levels"
          units = "deg"
          nlevels = nvlp

          CASE ("rp3P")
          var_description = "Relative humidity"
          units = "percentage"
          nlevels = nvlp

          CASE ("up3P")
          var_description = "U wind"
          units = "meter/second"
          nlevels = nvlp 

          CASE ("vp3P")
          var_description = "V wind"
          units = "meter/second"
          nlevels = nvlp
      
          CASE ("oz3P")
          var_description = "ozone mixing ratio"
          units = "kg per kg"
          nlevels = nvlp
 
          CASE ("vv3P")
          var_description = "vertical velocity"
          units = "Pa/s"
          nlevels = nvlp

          CASE ("qc3P")
          var_description = "cloud-hydrometeor condensate mixing ratio"
          units = "g/g"
          nlevels = nvlp

!temp variables 
          CASE ("t1xx")
          var_description = "temporary variable 1"
          units = "unit"
          nlevels = nvls

          CASE ("t2xx")
          var_description = "temporary variable 2"
          units = "unit"
          nlevels = nvls

          CASE ("t3xx")
          var_description = "temporary variable 3"
          units = "unit"
          nlevels = nvls + 1

          CASE ("t4xx")
          var_description = "temporary variable 4"
          units = "unit"
          nlevels = 4 

          CASE ("t5xx")
          var_description = "temporary variable 5"
          units = "unit"
          nlevels = 1 

!temp 2D variables
          CASE ("prpv")
          var_description = "pressure at tropopause"
          units = "pascals"
          nlevels = 1 

          CASE ("thpv")
          var_description = "potential temperature at tropopause"
          units = "kelvin"
          nlevels = 1 

          CASE ("uspv")
          var_description = "u-component at tropopause"
          units = "m/s"
          nlevels = 1 

          CASE ("vspv")
          var_description = "v-component at tropopause"
          units = "m/s"
          nlevels = 1 

! subseasonal forecast variables -- begin
          CASE ("td2D")
          var_description = "dewpoint temperature at 2m"
          units = "K"
          nlevels = 1

          CASE ("ssta")
          var_description = "averaged SST"
          units = "K"
          nlevels = 1

          CASE ("cice")
          var_description = "averaged sea ice"
          units = "%"
          nlevels = 1

          CASE ("sm2D")
          var_description = "2m vertically integrated soil moisture"
          units = "%"
          nlevels = 1

          CASE ("tmax")
          var_description = "maximum 2m air temperature"
          units = "K"
          nlevels = 1

          CASE ("tmin")
          var_description = "minimum 2m air temperature"
          units = "K"
          nlevels = 1

          CASE ("ustr")
          var_description = "zonal surface wind stress"
          units = "N/m2"
          nlevels = 1

          CASE ("vstr")
          var_description = "meridional surface wind stress"
          units = "N/m2"
          nlevels = 1

          CASE ("runo")
          var_description = "total runoff"
          units = "kg/m2"
          nlevels = 1

          CASE ("mslp")
          var_description = "sea level pressure"
          units = "Pa"
          nlevels = 1

          CASE ("cltt")
          var_description = "total cloud fraction"
          units = "%"
          nlevels = 1
! subseasonal forecast variables -- end


          CASE default
          write(6,*)'var_info: unknown input variable:', var_name(1:len_trim(var_name))
          var_description = "unknown variable"
          units = "unknown units"
          nlevels = 1 

          END SELECT  

        END SUBROUTINE var_info

! subroutine to set number of levels for the model
        SUBROUTINE set_model_nlevels(nlevels)
          INTEGER nlevels
          nvls = nlevels
        END SUBROUTINE set_model_nlevels
      END MODULE varinfo
