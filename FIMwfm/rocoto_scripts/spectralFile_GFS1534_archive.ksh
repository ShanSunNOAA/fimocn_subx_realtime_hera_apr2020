#!/bin/ksh -l

# Print out value of required environment variables
echo
echo "FIM_RUN      = ${FIM_RUN}"
echo "yyyymmddhh   = ${yyyymmddhh}"
echo "YR           = ${YR}"
echo "HOUR         = ${HOUR}"
echo "DOY          = ${DOY}"
echo

# Set up paths to shell commands
MV=/bin/mv
MKDIR=/bin/mkdir
CP=/bin/cp 
RM=/bin/rm 

# load HPSS module
module load hpss

# initialize
gfs_spectral_dir="/NCEPDEV/hpssuser/g01/wx20rt/WCOSS/prhw14"     #00Z 01 Jan - 00Z 15 Jul 2014
##gfs_spectral_dir="/NCEPDEV/hpssuser/g01/glopara/WCOSS/prhw14"    #18Z 14 Jul - 00Z 09 Oct 2014
##gfs_spectral_dir="/5year/NCEPDEV/emc-global/emc.glopara/WCOSS/prhw14"  #06Z 09 Oct - present
workdir=${FIM_RUN}/${yyyymmddhh}/ensics

# Set up a work directory and cd into it
echo ${workdir}
if [[ -d ${workdir} ]]; then
  ${RM} -rf ${workdir}
fi
${MKDIR} -p ${workdir}
cd ${workdir}

# Pull initial condition files from mass store
htar -xvf ${gfs_spectral_dir}/${yyyymmddhh}gfs.tar sfcanl.gfs.${yyyymmddhh} siganl.gfs.${yyyymmddhh}
status=$?
if [ $status != 0 ]; then
  echo "error $status extracting files from ${yyyymmddhh}gfs.tar in $gfs_spectral_dir"
  return $status
fi

echo $gfs_spectral_dir
echo ${yyyymmddhh}gfs.tar

# Rename files to our naming convention
#    YYDDDHH00.gfs.tHHz.sfcanl
#    YYDDDHH00.gfs.tHHz.sanl
#
mv sfcanl.gfs.${yyyymmddhh} ${YR}${DOY}${HOUR}00.gfs.t${HOUR}z.sfcanl
status=$?
if [ $status != 0 ]; then
  echo "error $status moving sfcanl.gfs.${yyyymmddhh} "
  return $status
fi
mv siganl.gfs.${yyyymmddhh} ${YR}${DOY}${HOUR}00.gfs.t${HOUR}z.sanl
status=$?
if [ $status != 0 ]; then
  echo "error $status moving siganl.gfs.${yyyymmddhh} "
  return $status
fi

