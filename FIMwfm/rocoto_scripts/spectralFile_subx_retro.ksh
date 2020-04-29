#!/bin/ksh -l

##
# Print out value of required environment variables
echo
echo "FIM_RUN = ${FIM_RUN}"
echo "YEAR    = ${YEAR}"
echo "MONTH   = ${MONTH}"
echo "DAY     = ${DAY}"
echo "HOUR    = ${HOUR}"
echo "PES     = ${PES}"
echo

# Set up paths to shell commands
MKDIR=/bin/mkdir
CP=/bin/cp
RM=/bin/rm
UNZIP=/usr/bin/unzip
TAR=/bin/tar

# load HPSS module
module load hpss

# initialize
workdir=${FIM_RUN}/fim${GLVL}_${NVL}_${PES}_${YEAR}${MONTH}${DAY}${HOUR}/prep

# Set up a work directory and cd into it
${MKDIR} -p ${workdir}
echo "workdir=$workdir "
cd ${workdir}

# Pull zip archive file from mass store
archive=cfsv2_${YEAR}${MONTH}${DAY}${HOUR}.tar
echo "IC tar is ${MASS_CFS_IC_DIR}/${YEAR}${MONTH}/${archive} "
hsi get ${MASS_CFS_IC_DIR}/${YEAR}${MONTH}/${archive}
status=$?
if [ $status != 0 ]; then
  echo "error $status retrieving tar file $MASS_CFS_IC_DIR/${YEAR}${MONTH}/${archive}"
  return $status
fi

# extract initial condition files from archive file
#${UNZIP} -o ${archive} 
${TAR} -xvf ${archive} 
status=$?
if [ $status != 0 ]; then
  echo "Cannot extract file from ${archive}"
  return $status
fi

# remove archive file
rmcmd="${RM} $archive"
echo "rmcmd:  $rmcmd"
${RM} $archive

exit $?
