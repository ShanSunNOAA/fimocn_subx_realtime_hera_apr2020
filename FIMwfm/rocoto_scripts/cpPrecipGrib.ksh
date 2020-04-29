#!/bin/ksh --login

## this script copies the GRIB2 file to /rt/fim/precip_gaussian 
##
##   J. Henderson    03/2014

# Print out value of required environment variables
echo entering cpPrecipGrib.ksh....
echo "FIM_RUN    = ${FIM_RUN}"
echo "MEMBER_ID  = ${MEMBER_ID}"
echo "yyyymmddhh = ${yyyymmddhh}"
echo "filename   = ${filename}"
echo

# initialize
CP=/bin/cp
grib2dir=${FIM_RUN}/${yyyymmddhh}/post_${MEMBER_ID}/129/NAT/grib2/
rtDir=/rt/fim/precip_gaussian/

# copy to /rt directory to be transferred to /public
echo "***copying to /rt...."
echo "$CP -p ${grib2dir}/${filename} $rtDir"
$CP -p ${grib2dir}/${filename} $rtDir/${filename}
status=$?
if [ $status != 0 ]; then
  echo "error $status copying ${grib2dir}/${filename} to ${rtDir}"
  return $status
fi
