#!/bin/ksh -l

##
## spectralFile_archive_ens.ksh
##    this script extracts GEFS initial condition files from the mass store
##    assumes data is stored in /BMC/fim/5year/GEFS_INIT/YYYY/YYYYMMDD00.tar files
##
##  J.K.Henderson    07-2015
##

# Print out value of required environment variables
echo
echo "FIM_RUN      = ${FIM_RUN}"
echo "yyyymmdd     = ${yyyymmdd}"
echo "YYYY         = ${YYYY}"
echo "MEMBER_ID    = ${MEMBER_ID}"
echo "YR           = ${YR}"
echo "HOUR         = ${HOUR}"
echo "JDY          = ${JDY}"
echo

# Set up paths to shell commands
MV=/bin/mv
MKDIR=/bin/mkdir
RM=/bin/rm 

# load HPSS module
module load hpss

# initialize
gefs_spectral_dir="/BMC/fim/5year/GEFS_INIT/${YYYY}/"     
workdir=${FIM_RUN}/${yyyymmdd}${HOUR}/ensics
echo $gefs_spectral_dir
echo ${yyyymmdd}.tar

# Set up a work directory and cd into it
workdir=${FIM_RUN}/${yyyymmdd}${HOUR}/ensics_${MEMBER_ID}
echo ${workdir}
if [[ -d ${workdir} ]]; then
  ${RM} -rf ${workdir}
fi
${MKDIR} -p ${workdir}
status=$?
if [ $status != 0 ]; then
  echo "error $status creating directory $work_dir"
  return $status
else 
  cd ${workdir}
fi

# Pull initial condition files from mass store
htar -xvf ${gefs_spectral_dir}/${yyyymmdd}.tar gec00.t00z.sfcanl gep${MEMBER_ID}.t00z.sanl
status=$?
if [ $status != 0 ]; then
  echo "error $status extracting files from ${yyyymmdd}.tar in $gefs_spectral_dir"
  return $status
fi

# Rename files to our naming convention
#    YYDDDHH00.gfs.t00z.sfcanl
#    YYDDDHH00.gfs.t00z.sanl
mvcmd="${MV} gec00.t${HOUR}z.sfcanl ${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sfcanl"
echo "cpcmd: ${mvcmd}"
${mvcmd}
mvcmd="${MV} gep${MEMBER_ID}.t${HOUR}z.sanl ${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sanl"
echo "mvcmd: ${mvcmd}"
${mvcmd}
