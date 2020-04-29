#!/bin/ksh -l

##   spectralFile_archive_fdr.ksh
##      extracts GFS T1534 spectral files from mass store
##      /BMC/fdr/Permanent/YYYY/MM/DD/data/grids/gf/spectral_t1534/YYYYMMDD0000.zip
##      renames files, and copies them to $FIM_HOME/FIM_RUN/YYYYMMDDHH/ensics
##
##    J.K.Henderson   10-2015
##

# Print out value of required environment variables
echo
echo "FIM_RUN      = ${FIM_RUN}"
echo "YYYY         = ${YYYY}"
echo "MM           = ${MM}"
echo "DD           = ${DD}"
echo "HOUR         = ${HOUR}"
echo "GLVL         = ${GLVL}"
echo "NVL          = ${NVL}"
echo "PES          = ${PES}"
echo

# Set up paths to shell commands
MKDIR=/bin/mkdir
CP=/bin/cp 
RM=/bin/rm 
UNZIP=/usr/bin/unzip

# load HPSS module
module load hpss

# initialize
gfs_spectral_dir="/BMC/fdr/Permanent/${YYYY}/${MM}/${DD}/data/grids/gfs/spectral_t1534"      ## T1534
#workdir=${FIM_RUN}/${yyyymmdd}${HOUR}/ensics
workdir=${FIM_RUN}/fim_${GLVL}_${NVL}_${PES}_${YYYY}${MM}${DD}${HOUR}00/ensics
echo $gfs_spectral_dir

# Set up a work directory and cd into it
echo ${workdir}
${MKDIR} -p ${workdir}
cd ${workdir}

# Pull zip archive file from mass store
archive=${YYYY}${MM}${DD}0000.zip 
hsi get ${gfs_spectral_dir}/${archive}
status=$?
if [ $status != 0 ]; then
  echo "error $status retrieving zip file $gefs_spectral_dir/${archive}"
  return $status
fi

# extract initial condition files from archive file
${UNZIP} -o ${archive} \*t${HOUR}z\*
status=$?
if [ $status != 0 ]; then
  echo "Can't extract file from ${archive}"
  return $status
fi

# remove archive file
rmcmd="${RM} $archive"
echo "rmcmd:  $rmcmd"
${RM} $archive

exit $?
