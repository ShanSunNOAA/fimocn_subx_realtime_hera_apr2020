#!/bin/ksh -l
#
############################################################################
# pullVars_FlatFiles.sh
# Created by: Ben Green
# Last updated on: April 17, 2017
# This script takes the DAILY files (one per each time-lagged ensemble member)
#   and then pulls out specific variables into "flat" netcdf files for CPC/IRI
#   Also, lots of variable/attribute renaming and unit conversion is required here
############################################################################

module load intel 
module load nco
##
# Print out value of required environment variables
echo
echo "FIM_RUN = ${FIM_RUN}"
echo "YEAR    = ${YEAR}"
echo "MONTH   = ${MONTH}"
echo "DAY     = ${DAY}"
echo "HOUR    = ${HOUR}"
echo
yyyy=${YEAR}
mm=${MONTH}
dd=${DAY}
hh=${HOUR}

# Set up paths to shell commands
MKDIR=/bin/mkdir
CP=/bin/cp
RM=/bin/rm
UNZIP=/usr/bin/unzip
TAR=/bin/tar

yy2=${yyyy}
mm2=${mm}
dd2=${dd}

if [[ $hh -eq 12 ]] ; then
  em=1
  # Move on to next day
  date1=`date -d "${yyyy}-${mm}-${dd}" +%Y%m%d`
  date2=`date -d "${date1} + 1 day"`
  yy2=`date -d "${date2}" +%Y`
  mm2=`date -d "${date2}" +%m`
  dd2=`date -d "${date2}" +%d`
elif [[ $hh -eq 18 ]] ; then
  em=2
  date1=`date -d "${yyyy}-${mm}-${dd}" +%Y%m%d`
  date2=`date -d "${date1} + 1 day"`
  yy2=`date -d "${date2}" +%Y`
  mm2=`date -d "${date2}" +%m`
  dd2=`date -d "${date2}" +%d`
elif [[ $hh -eq 0 ]] ; then
  em=3
elif [[ $hh -eq 6 ]] ; then
  em=4
fi
# find the dir
dy2=${yy2}${mm2}${dd2}00   # need work to use Wed date
ncfile=${NC3D_DIR}/fim${GLVL}_${yyyy}${mm}${dd}${hh}.nc
tmpdir=${NC3D_DIR}/${yyyy}${mm}${dd}${hh}
outDirStart=${NC_SUBX_DIR}
mkdir $tmpdir
cd $tmpdir

    if   [ $mm2 -eq 1  ] ; then
      mon2="jan"
    elif [ $mm2 -eq 2  ] ; then
      mon2="feb"
    elif [ $mm2 -eq 3  ] ; then
      mon2="mar"
    elif [ $mm2 -eq 4  ] ; then
      mon2="apr"
    elif [ $mm2 -eq 5  ] ; then
      mon2="may"
    elif [ $mm2 -eq 6  ] ; then
      mon2="jun"
    elif [ $mm2 -eq 7  ] ; then
      mon2="jul"
    elif [ $mm2 -eq 8  ] ; then
      mon2="aug"
    elif [ $mm2 -eq 9  ] ; then
      mon2="sep"
    elif [ $mm2 -eq 10 ] ; then
      mon2="oct"
    elif [ $mm2 -eq 11 ] ; then
      mon2="nov"
    elif [ $mm2 -eq 12 ] ; then
      mon2="dec"
    else
      echo "Invalid month!"
      exit
    fi # End of getting 3-character month

outStr=FIM_${dd2}${mon2}${yy2}_00z_d01_d32
echo "do nc_subx for ${yyyy}${mm}${dd}${hh}; string=$outStr"

ncks -h -F -d time,2,33 ${ncfile} tempDaily_${dy2}.nc

# Need to add 4 fields together and multiply by -1 to get positive up net surface radiation
  ncflint -h -O --fix_rec_crd -v rlds -w '-1.0',0.0 tempDaily_${dy2}.nc tempDaily_${dy2}.nc tempRLDS_${dy2}.nc
  ncrename -h -O -v rlds,rad tempRLDS_${dy2}.nc

  ncflint -h -O --fix_rec_crd -v rlus -w '-1.0',0.0 tempDaily_${dy2}.nc tempDaily_${dy2}.nc tempRLUS_${dy2}.nc
  ncrename -h -O -v rlus,rad tempRLUS_${dy2}.nc

  ncflint -h -O --fix_rec_crd -v rsds -w '-1.0',0.0 tempDaily_${dy2}.nc tempDaily_${dy2}.nc tempRSDS_${dy2}.nc
  ncrename -h -O -v rsds,rad tempRSDS_${dy2}.nc

  ncflint -h -O --fix_rec_crd -v rsus -w '-1.0',0.0 tempDaily_${dy2}.nc tempDaily_${dy2}.nc tempRSUS_${dy2}.nc
  ncrename -h -O -v rsus,rad tempRSUS_${dy2}.nc

  outDir=${outDirStart}/rad_sfc/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncbo -h --op_typ=add tempRLDS_${dy2}.nc tempRLUS_${dy2}.nc tempRLTS_${dy2}.nc
  rm tempRLDS_${dy2}.nc tempRLUS_${dy2}.nc
  ncbo -h --op_typ=add tempRSDS_${dy2}.nc tempRSUS_${dy2}.nc tempRSTS_${dy2}.nc
  rm tempRSDS_${dy2}.nc tempRSUS_${dy2}.nc
  ncbo -h --op_typ=add tempRLTS_${dy2}.nc tempRSTS_${dy2}.nc ${outDir}/rad_sfc_${outStr}_m0${em}.nc
  rm tempRLTS_${dy2}.nc tempRSTS_${dy2}.nc

  # Finally rename the field and add attributes
  ncatted -h -O -a units,rad,m,c,"W m^-2" ${outDir}/rad_sfc_${outStr}_m0${em}.nc
  ncatted -h -O -a long_name,rad,m,c,"net_surface_radiation" ${outDir}/rad_sfc_${outStr}_m0${em}.nc


  # Do unit conversion for various fields (extract them here and save to separate flat files)

  # Precipitation extraction, unit conversion, variable renaming, etc.
  outDir=${outDirStart}
  echo "chk outDIR is $outDir $outDirStart "
  ncflint -h -O --fix_rec_crd -v r12D -w '1.157407407e-5',0.0 tempDaily_${dy2}.nc tempDaily_${dy2}.nc ${outDir}/pr_sfc_${outStr}_m0${em}.nc
  ncrename -h -O -v r12D,pr ${outDir}/pr_sfc_${outStr}_m0${em}.nc
  ncatted -h -O -a units,pr,m,c,"kg m^-2 s^-1" ${outDir}/pr_sfc_${outStr}_m0${em}.nc
  ncatted -h -O -a long_name,pr,m,c,"precipitation_flux" ${outDir}/pr_sfc_${outStr}_m0${em}.nc

  # OLR extraction, sign reversal (positive = UP)
  outDir=${outDirStart}/rlut_toa/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncflint -h -O --fix_rec_crd -v rlut -w '-1.0',0.0 tempDaily_${dy2}.nc tempDaily_${dy2}.nc ${outDir}/rlut_toa_${outStr}_m0${em}.nc
  ncatted -h -O -a units,rlut,m,c,"W m^-2" ${outDir}/rlut_toa_${outStr}_m0${em}.nc
  ncatted -h -O -a long_name,rlut,m,c,"toa_outgoing_longwave_flux" ${outDir}/rlut_toa_${outStr}_m0${em}.nc

  # Surface sensible heat flux extraction, sign reversal (positive = UP)
  outDir=${outDirStart}/hfss_sfc/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncflint -h -O --fix_rec_crd -v hfss -w '-1.0',0.0 tempDaily_${dy2}.nc tempDaily_${dy2}.nc ${outDir}/hfss_sfc_${outStr}_m0${em}.nc
  ncatted -h -O -a units,hfss,m,c,"W m^-2" ${outDir}/hfss_sfc_${outStr}_m0${em}.nc
  ncatted -h -O -a long_name,hfss,m,c,"surface_upward_sensible_heat_flux" ${outDir}/hfss_sfc_${outStr}_m0${em}.nc

  # Surface latent heat flux extraction, sign reversal (positive = UP)
  outDir=${outDirStart}/hfls_sfc/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncflint -h -O --fix_rec_crd -v hfls -w '-1.0',0.0 tempDaily_${dy2}.nc tempDaily_${dy2}.nc ${outDir}/hfls_sfc_${outStr}_m0${em}.nc
  ncatted -h -O -a units,hfls,m,c,"W m^-2" ${outDir}/hfls_sfc_${outStr}_m0${em}.nc
  ncatted -h -O -a long_name,hfls,m,c,"surface_upward_latent_heat_flux" ${outDir}/hfls_sfc_${outStr}_m0${em}.nc

  # Surface zonal stress extraction, sign reversal (positive = UP): FIM code has correct units, but fim2nc writes them wrong
  outDir=${outDirStart}/stx_sfc/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncflint -h -O --fix_rec_crd -v ustr -w '-1.0',0.0 tempDaily_${dy2}.nc tempDaily_${dy2}.nc ${outDir}/stx_sfc_${outStr}_m0${em}.nc
  ncrename -h -O -v ustr,stx ${outDir}/stx_sfc_${outStr}_m0${em}.nc
  ncatted -h -O -a units,stx,m,c,"N m^-2" ${outDir}/stx_sfc_${outStr}_m0${em}.nc
  ncatted -h -O -a long_name,stx,m,c,"surface_zonal_stress_positive_to_the_west" ${outDir}/stx_sfc_${outStr}_m0${em}.nc

  # Surface meridional stress extraction, sign reversal (positive = UP): FIM code has correct units, but fim2nc writes them wrong
  outDir=${outDirStart}/sty_sfc/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncflint -h -O --fix_rec_crd -v vstr -w '-1.0',0.0 tempDaily_${dy2}.nc tempDaily_${dy2}.nc ${outDir}/sty_sfc_${outStr}_m0${em}.nc
  ncrename -h -O -v vstr,sty ${outDir}/sty_sfc_${outStr}_m0${em}.nc
  ncatted -h -O -a units,sty,m,c,"N m^-2" ${outDir}/sty_sfc_${outStr}_m0${em}.nc
  ncatted -h -O -a long_name,sty,m,c,"surface_meridional_stress_positive_to_the_south" ${outDir}/sty_sfc_${outStr}_m0${em}.nc
  
  # Vertically integrated soil moisture, initially "volumetric" (unitless!) -- multiply by 1000 kg/m3 then by 2 m (depth) to get kg/m2
  outDir=${outDirStart}/mrso_sfc/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncflint -h -O --fix_rec_crd -v sm2D -w 2000.0,0.0 tempDaily_${dy2}.nc tempDaily_${dy2}.nc ${outDir}/mrso_sfc_${outStr}_m0${em}.nc
  ncrename -h -O -v sm2D,mrso ${outDir}/mrso_sfc_${outStr}_m0${em}.nc
  ncatted -h -O -a units,mrso,m,c,"kg m^-2" ${outDir}/mrso_sfc_${outStr}_m0${em}.nc
  ncatted -h -O -a long_name,mrso,m,c,"soil_moisture_content" ${outDir}/mrso_sfc_${outStr}_m0${em}.nc


  # Now move on to variable/attribute renaming
  ncrename -h -O -v hgtP,zg tempDaily_${dy2}.nc
  ncatted -h -O -a long_name,zg,m,c,"geopotential_height" tempDaily_${dy2}.nc

  ncrename -h -O -v up3P,ua tempDaily_${dy2}.nc
  ncatted -h -O -a units,ua,m,c,"m s^-1" tempDaily_${dy2}.nc
  ncatted -h -O -a long_name,ua,m,c,"eastward_wind" tempDaily_${dy2}.nc

  ncrename -h -O -v vp3P,va tempDaily_${dy2}.nc
  ncatted -h -O -a units,va,m,c,"m s^-1" tempDaily_${dy2}.nc
  ncatted -h -O -a long_name,va,m,c,"northward_wind" tempDaily_${dy2}.nc

  ncrename -h -O -v tmpP,ta tempDaily_${dy2}.nc
  ncatted -h -O -a long_name,ta,m,c,"air_temperature" tempDaily_${dy2}.nc

  ncrename -h -O -v vv3P,wap tempDaily_${dy2}.nc
  ncatted -h -O -a units,wap,m,c,"Pa s^-1" tempDaily_${dy2}.nc
  ncatted -h -O -a long_name,wap,m,c,"lagrangian_tendency_of_air_pressure" tempDaily_${dy2}.nc

  ncrename -h -O -v qv3P,huss tempDaily_${dy2}.nc
  ncatted -h -O -a units,huss,m,c,"1" tempDaily_${dy2}.nc
  ncatted -h -O -a long_name,huss,m,c,"specific_humidity" tempDaily_${dy2}.nc

  ncrename -h -O -v t22D,tas tempDaily_${dy2}.nc # 2-m air temperature
  ncatted -h -O -a units,tas,m,c,"K" tempDaily_${dy2}.nc
  ncatted -h -O -a long_name,tas,m,c,"air_temperature" tempDaily_${dy2}.nc

  ncrename -h -O -v ts2D,ts tempDaily_${dy2}.nc # Skin/surface temperature
  ncatted -h -O -a long_name,ts,m,c,"surface_temperature" tempDaily_${dy2}.nc

  ncrename -h -O -v u12D,uas tempDaily_${dy2}.nc # 10-m zonal wind
  ncatted -h -O -a units,uas,m,c,"m s^-1" tempDaily_${dy2}.nc
  ncatted -h -O -a long_name,uas,m,c,"eastward_wind" tempDaily_${dy2}.nc

  ncrename -h -O -v v12D,vas tempDaily_${dy2}.nc # 10-m meridional wind
  ncatted -h -O -a units,vas,m,c,"m s^-1" tempDaily_${dy2}.nc
  ncatted -h -O -a long_name,vas,m,c,"northward_wind" tempDaily_${dy2}.nc

  ncrename -h -O -v tmax,tasmax tempDaily_${dy2}.nc # 24-hour maximum 2-m air temperature
  ncatted -h -O -a units,tasmax,m,c,"K" tempDaily_${dy2}.nc
  ncatted -h -O -a long_name,tasmax,m,c,"air_temperature" tempDaily_${dy2}.nc

  ncrename -h -O -v tmin,tasmin tempDaily_${dy2}.nc # 24-hour minimum 2-m air temperature
  ncatted -h -O -a units,tasmin,m,c,"K" tempDaily_${dy2}.nc
  ncatted -h -O -a long_name,tasmin,m,c,"air_temperature" tempDaily_${dy2}.nc

  ncrename -h -O -v td2D,tdps tempDaily_${dy2}.nc # 2-m dewpoint temperature
  ncatted -h -O -a units,tdps,m,c,"K" tempDaily_${dy2}.nc
  ncatted -h -O -a long_name,tdps,m,c,"dew_point_temperature" tempDaily_${dy2}.nc

  ncrename -h -O -v ms2D,psl tempDaily_${dy2}.nc # MSLP
  ncatted -h -O -a long_name,psl,m,c,"air_pressure_at_sea_level" tempDaily_${dy2}.nc

  ncrename -h -O -v sa2D,swe tempDaily_${dy2}.nc # Snow-water equivalent
  ncatted -h -O -a units,swe,m,c,"kg m^-2" tempDaily_${dy2}.nc # Multiply by 1000/1000 to convert from mm to kg/m2
  ncatted -h -O -a long_name,swe,m,c,"snow_water_equivalent" tempDaily_${dy2}.nc 

  ncrename -h -O -v fice,sic tempDaily_${dy2}.nc # Sea ice concentration (already in percent)
  ncatted -h -O -a units,sic,m,c,"%" tempDaily_${dy2}.nc
  ncatted -h -O -a long_name,sic,m,c,"sea_ice_area_fraction" tempDaily_${dy2}.nc 


  ##########  
  # Now extract more variables
  # Priority 1 daily averages from snapshots

  #---------
  # 500 hPa geopotential height
  outDir=${outDirStart}/zg_500/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v zg,lat,lon,'time' -d pressure,500. tempDaily_${dy2}.nc zg_500_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure zg_500_${outStr}_m0${em}.nc zg_500_${outStr}_m0${em}.nc
  mv zg_500_${outStr}_m0${em}.nc ${outDir}/zg_500_${outStr}_m0${em}.nc
  #----

  #---------
  # 200 hPa geopotential height
  outDir=${outDirStart}/zg_200/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v zg,lat,lon,'time' -d pressure,200. tempDaily_${dy2}.nc zg_200_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure zg_200_${outStr}_m0${em}.nc zg_200_${outStr}_m0${em}.nc
  mv zg_200_${outStr}_m0${em}.nc ${outDir}/zg_200_${outStr}_m0${em}.nc
  #----

  #---------
  # 850 hPa zonal wind
  outDir=${outDirStart}/ua_850/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v ua,lat,lon,'time' -d pressure,850. tempDaily_${dy2}.nc ua_850_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure ua_850_${outStr}_m0${em}.nc ua_850_${outStr}_m0${em}.nc
  mv ua_850_${outStr}_m0${em}.nc ${outDir}/ua_850_${outStr}_m0${em}.nc
  #----

  #---------
  # 850 hPa meridional wind
  outDir=${outDirStart}/va_850/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v va,lat,lon,'time' -d pressure,850. tempDaily_${dy2}.nc va_850_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure va_850_${outStr}_m0${em}.nc va_850_${outStr}_m0${em}.nc
  mv va_850_${outStr}_m0${em}.nc ${outDir}/va_850_${outStr}_m0${em}.nc
  #----

  #---------
  # 200 hPa zonal wind
  outDir=${outDirStart}/ua_200/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v ua,lat,lon,'time' -d pressure,200. tempDaily_${dy2}.nc ua_200_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure ua_200_${outStr}_m0${em}.nc ua_200_${outStr}_m0${em}.nc
  mv ua_200_${outStr}_m0${em}.nc ${outDir}/ua_200_${outStr}_m0${em}.nc
  #----

  #---------
  # 200 hPa meridional wind
  outDir=${outDirStart}/va_200/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v va,lat,lon,'time' -d pressure,200. tempDaily_${dy2}.nc va_200_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure va_200_${outStr}_m0${em}.nc va_200_${outStr}_m0${em}.nc
  mv va_200_${outStr}_m0${em}.nc ${outDir}/va_200_${outStr}_m0${em}.nc
  #----


  #----------
  # Priority 1 daily averages from accumulated fields
  outDir=${outDirStart}/tas_2m/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -O -v tas tempDaily_${dy2}.nc ${outDir}/tas_2m_${outStr}_m0${em}.nc # 2-m temperature

  outDir=${outDirStart}/ts_sfc/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -O -v ts tempDaily_${dy2}.nc ${outDir}/ts_sfc_${outStr}_m0${em}.nc  # Surface temperature
  #----
  

  ##########
  # Priority 2 daily averages from snapshots

  #----------
  # 850 hPa specific humidity
  outDir=${outDirStart}/huss_850/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v huss,lat,lon,'time' -d pressure,850. tempDaily_${dy2}.nc huss_850_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure huss_850_${outStr}_m0${em}.nc huss_850_${outStr}_m0${em}.nc
  mv huss_850_${outStr}_m0${em}.nc ${outDir}/huss_850_${outStr}_m0${em}.nc
  #----

  #----------
  # 500 hPa vertical velocity
  outDir=${outDirStart}/wap_500/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v wap,lat,lon,'time' -d pressure,500. tempDaily_${dy2}.nc wap_500_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure wap_500_${outStr}_m0${em}.nc wap_500_${outStr}_m0${em}.nc
  mv wap_500_${outStr}_m0${em}.nc ${outDir}/wap_500_${outStr}_m0${em}.nc
  #----

  #----------
  # 100 hPa zonal wind
  outDir=${outDirStart}/ua_100/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v ua,lat,lon,'time' -d pressure,100. tempDaily_${dy2}.nc ua_100_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure ua_100_${outStr}_m0${em}.nc ua_100_${outStr}_m0${em}.nc
  mv ua_100_${outStr}_m0${em}.nc ${outDir}/ua_100_${outStr}_m0${em}.nc
  #----

  #----------
  # 100 hPa meridional wind
  outDir=${outDirStart}/va_100/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v va,lat,lon,'time' -d pressure,100. tempDaily_${dy2}.nc va_100_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure va_100_${outStr}_m0${em}.nc va_100_${outStr}_m0${em}.nc
  mv va_100_${outStr}_m0${em}.nc ${outDir}/va_100_${outStr}_m0${em}.nc
  #----

  #----------
  # MSLP
  outDir=${outDirStart}/psl_sfc/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -O -v psl tempDaily_${dy2}.nc ${outDir}/psl_sfc_${outStr}_m0${em}.nc
  #----


  ##########
  # Priority 2 daily averages from accumulated fields
  outDir=${outDirStart}/uas_10m/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -O -v uas tempDaily_${dy2}.nc ${outDir}/uas_10m_${outStr}_m0${em}.nc   # 10-m U wind

  outDir=${outDirStart}/vas_10m/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -O -v vas tempDaily_${dy2}.nc ${outDir}/vas_10m_${outStr}_m0${em}.nc   # 10-m V wind

  outDir=${outDirStart}/tdps_2m/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -O -v tdps tempDaily_${dy2}.nc ${outDir}/tdps_2m_${outStr}_m0${em}.nc  # 10-m dew point

  outDir=${outDirStart}/swe_sfc/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -O -v swe tempDaily_${dy2}.nc ${outDir}/swe_sfc_${outStr}_m0${em}.nc   # Snow-water equivalent

  outDir=${outDirStart}/sic_sfc/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -O -v sic tempDaily_${dy2}.nc ${outDir}/sic_sfc_${outStr}_m0${em}.nc   # Sea ice fraction

  # Priority 2 max and min
  outDir=${outDirStart}/tasmax_2m/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -O -v tasmax tempDaily_${dy2}.nc ${outDir}/tasmax_2m_${outStr}_m0${em}.nc

  outDir=${outDirStart}/tasmin_2m/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -O -v tasmin tempDaily_${dy2}.nc ${outDir}/tasmin_2m_${outStr}_m0${em}.nc
  ####


  ##########
  # Priority 3 daily averages from snapshots

  #----------
  # 850 hPa geopotential height
  outDir=${outDirStart}/zg_850/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v zg,lat,lon,'time' -d pressure,850. tempDaily_${dy2}.nc zg_850_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure zg_850_${outStr}_m0${em}.nc zg_850_${outStr}_m0${em}.nc
  mv zg_850_${outStr}_m0${em}.nc ${outDir}/zg_850_${outStr}_m0${em}.nc
  #----

  #----------
  # 50 hPa geopotential height
  outDir=${outDirStart}/zg_50/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v zg,lat,lon,'time' -d pressure,50. tempDaily_${dy2}.nc zg_50_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure zg_50_${outStr}_m0${em}.nc zg_50_${outStr}_m0${em}.nc
  mv zg_50_${outStr}_m0${em}.nc ${outDir}/zg_50_${outStr}_m0${em}.nc
  #----

  #----------
  # 30 hPa geopotential height
  outDir=${outDirStart}/zg_30/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v zg,lat,lon,'time' -d pressure,30. tempDaily_${dy2}.nc zg_30_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure zg_30_${outStr}_m0${em}.nc zg_30_${outStr}_m0${em}.nc
  mv zg_30_${outStr}_m0${em}.nc ${outDir}/zg_30_${outStr}_m0${em}.nc
  #----

  #----------
  # 10 hPa geopotential height
  outDir=${outDirStart}/zg_10/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v zg,lat,lon,'time' -d pressure,10. tempDaily_${dy2}.nc zg_10_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure zg_10_${outStr}_m0${em}.nc zg_10_${outStr}_m0${em}.nc
  mv zg_10_${outStr}_m0${em}.nc ${outDir}/zg_10_${outStr}_m0${em}.nc
  #----

  #----------
  # 100 hPa temperature
  outDir=${outDirStart}/ta_100/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v ta,lat,lon,'time' -d pressure,100. tempDaily_${dy2}.nc ta_100_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure ta_100_${outStr}_m0${em}.nc ta_100_${outStr}_m0${em}.nc
  mv ta_100_${outStr}_m0${em}.nc ${outDir}/ta_100_${outStr}_m0${em}.nc
  #----

  #----------
  # 50 hPa temperature
  outDir=${outDirStart}/ta_50/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v ta,lat,lon,'time' -d pressure,50. tempDaily_${dy2}.nc ta_50_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure ta_50_${outStr}_m0${em}.nc ta_50_${outStr}_m0${em}.nc
  mv ta_50_${outStr}_m0${em}.nc ${outDir}/ta_50_${outStr}_m0${em}.nc
  #----

  #----------
  # 30 hPa temperature
  outDir=${outDirStart}/ta_30/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v ta,lat,lon,'time' -d pressure,30. tempDaily_${dy2}.nc ta_30_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure ta_30_${outStr}_m0${em}.nc ta_30_${outStr}_m0${em}.nc
  mv ta_30_${outStr}_m0${em}.nc ${outDir}/ta_30_${outStr}_m0${em}.nc
  #----

  #----------
  # 10 hPa temperature
  outDir=${outDirStart}/ta_10/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v ta,lat,lon,'time' -d pressure,10. tempDaily_${dy2}.nc ta_10_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure ta_10_${outStr}_m0${em}.nc ta_10_${outStr}_m0${em}.nc
  mv ta_10_${outStr}_m0${em}.nc ${outDir}/ta_10_${outStr}_m0${em}.nc
  #----

  #----------
  # 50 hPa zonal wind
  outDir=${outDirStart}/ua_50/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v ua,lat,lon,'time' -d pressure,50. tempDaily_${dy2}.nc ua_50_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure ua_50_${outStr}_m0${em}.nc ua_50_${outStr}_m0${em}.nc
  mv ua_50_${outStr}_m0${em}.nc ${outDir}/ua_50_${outStr}_m0${em}.nc
  #----

  #----------
  # 30 hPa zonal wind
  outDir=${outDirStart}/ua_30/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v ua,lat,lon,'time' -d pressure,30. tempDaily_${dy2}.nc ua_30_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure ua_30_${outStr}_m0${em}.nc ua_30_${outStr}_m0${em}.nc
  mv ua_30_${outStr}_m0${em}.nc ${outDir}/ua_30_${outStr}_m0${em}.nc
  #----

  #----------
  # 10 hPa zonal wind
  outDir=${outDirStart}/ua_10/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v ua,lat,lon,'time' -d pressure,10. tempDaily_${dy2}.nc ua_10_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure ua_10_${outStr}_m0${em}.nc ua_10_${outStr}_m0${em}.nc
  mv ua_10_${outStr}_m0${em}.nc ${outDir}/ua_10_${outStr}_m0${em}.nc
  #----

  #----------
  # 50 hPa meridional wind
  outDir=${outDirStart}/va_50/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v va,lat,lon,'time' -d pressure,50. tempDaily_${dy2}.nc va_50_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure va_50_${outStr}_m0${em}.nc va_50_${outStr}_m0${em}.nc
  mv va_50_${outStr}_m0${em}.nc ${outDir}/va_50_${outStr}_m0${em}.nc
  #----

  #----------
  # 30 hPa meridional wind
  outDir=${outDirStart}/va_30/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v va,lat,lon,'time' -d pressure,30. tempDaily_${dy2}.nc va_30_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure va_30_${outStr}_m0${em}.nc va_30_${outStr}_m0${em}.nc
  mv va_30_${outStr}_m0${em}.nc ${outDir}/va_30_${outStr}_m0${em}.nc
  #----

  #----------
  # 10 hPa meridional wind
  outDir=${outDirStart}/va_10/${yy2}/${mm2}
  outDir=${outDirStart}
  mkdir -p ${outDir}
  ncks -h -C -v va,lat,lon,'time' -d pressure,10. tempDaily_${dy2}.nc va_10_${outStr}_m0${em}.nc
  ncwa -h -O -a pressure va_10_${outStr}_m0${em}.nc va_10_${outStr}_m0${em}.nc
  mv va_10_${outStr}_m0${em}.nc ${outDir}/va_10_${outStr}_m0${em}.nc
  #----
  ##########

  ##########
  # Priority 4
  # Eventually fill in TOA shortwave, and net surface radiation fluxes (longwave and shortwave separately, LW + SW together is in Priority 2)
  ##########  
##  zip all *.nc files
  correctsize=308000  ### for theia
  correctsize=271620  ### for hera
  correctsize=222712  ### for hera
  correctsize=207552  ### for no3
  
  cd ${outDir}
  date2=${dd2}${mon2}${yy2}
  zipfile=fim_${date2}_m0${em}.zip
  echo "zipfile will be in ${zipfile}"
  zip -r $zipfile *FIM_${date2}*_m0${em}.nc
  sleep 60
  
  if [[ -f $zipfile ]]; then
    actualsize=$(du -k $zipfile | cut -f 1)
    if [[ $actualsize -lt $correctsize ]]; then
      echo "chk2 wrong zip size of $actualsize, should be $correctsize for ${zipfile}"
    else
      echo "chk2 $zipfile has size of $actualsize; removing nc files"
      /bin/rm *FIM_${date2}*_m0${em}.nc
    fi
  fi
  
  # Delete temporary file and (local) daily file [local daily file already saved elsewhere]
  if [[ -d $tmpdir ]]; then
    /bin/rm -rf $tmpdir
  fi
exit # Done with shell
