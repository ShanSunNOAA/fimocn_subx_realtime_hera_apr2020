#!/bin/ksh --login

## driver script to call cpPrecipGrib.ksh, which copies the GRIB2 file to 
##   /rt/fim/precip_gaussian. 
##
##   NOTE:  cpPrecipGrib.ksh needs to run on SERVICE node since compute nodes can not see /rt
##
##   J. Henderson    03/2014

# Print out value of required environment variables
echo entering precipGrib.ksh...
echo "PROJECT          = ${PROJECT}"
echo "SCRIPTS          = ${SCRIPTS}"
echo "FIM_HOME         = ${FIM_HOME}"
echo "FIM_RUN          = ${FIM_RUN}"
echo "MEMBER_ID        = ${MEMBER_ID}"
echo "yyyymmddhh       = ${yyyymmddhh}"
echo "filename         = ${filename}"
echo

# initialize
LOGDIR=${FIM_HOME}/FIMwfm/log/precip/
if [[ ! -d ${LOGDIR} ]]; then
  mkdir $LOGDIR
fi

# call script to copy to /rt directory 
echo "***calling cpPrecipGrib.ksh...."
qsub -V -A $PROJECT -d $LOGDIR -l partition=service -l procs=1 -l walltime=00:05:00 -j oe -N cpPrecip${filename} -e $LOGDIR -e $LOGDIR ${SCRIPTS}/cpPrecipGrib.ksh
