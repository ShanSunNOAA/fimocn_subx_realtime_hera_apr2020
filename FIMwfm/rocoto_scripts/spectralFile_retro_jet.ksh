#!/bin/ksh -l

##   spectralFile_retro.ksh
##      reads spectral files from ${DATADR2}/MM/
##      copies to $FIM_RUN/YYYYMMDDHH/ensics
##

# Set the SGE queueing options
#$ -S /bin/ksh
#$ -pe comp 1
#$ -l h_rt=00:30:00
#$ -N SpectralFiles
#$ -j y
#$ -V

# Set up paths to shell commands
MKDIR=/bin/mkdir
CP=/bin/cp 
RM=/bin/rm 
TAR=/bin/tar 
UNZIP=/usr/bin/unzip

# Set up a work directory and cd into it
workdir=${FIM_RUN}/fim_${GLVL}_${NVL}_${PES}_${yyyymmddhhmm}/ensics
echo ${workdir}
${MKDIR} -p ${workdir}
cd ${workdir}

# Print out value of required environment variables
echo
echo "FIM_RUN = ${FIM_RUN}"
echo "YEAR    = ${YEAR}"
echo "MONTH   = ${MONTH}"
echo "DAY     = ${DAY}"
echo "HOUR    = ${HOUR}"
echo "DATADR2 = ${DATADR2}"
echo

DATEDIR=${MONTH}
echo "DATEDIR: ${DATADR2}/${DATEDIR}"

# tar.gz archive file
if [ -e ${DATADR2}/${DATEDIR}/${YEAR}${MONTH}${DAY}${HOUR}00.tar.gz ]
then
  archive=${YEAR}${MONTH}${DAY}${HOUR}00.tar.gz
  ext=.tar.gz
  excmd=$TAR
# zip archive file
elif [ -e ${DATADR2}/${DATEDIR}/${YEAR}${MONTH}${DAY}0000.zip ]
then
  archive=${YEAR}${MONTH}${DAY}0000.zip
  ext=.zip
  excmd=$UNZIP
else
  echo "no archive file found in ${DATADR2/$DATADIR}"
  exit 1
fi

# copy file to local work directory
cpcmd="${CP} ${DATADR2}/${DATEDIR}/${archive} ${workdir}"
echo "cpcmd: ${cpcmd}"
${cpcmd}

# extract files from archive file
cd ${workdir}
if [ -e ${YEAR}${MONTH}${DAY}${HOUR}00.tar.gz ]      ##  tar file
then
  ${excmd} xvfz $archive
elif [ -e ${YEAR}${MONTH}${DAY}0000.zip ]
then
  ${excmd} $archive \*t${HOUR}z\*
else
  echo "Can't extract file from $archive"
fi

# remove archive file
rmcmd="${RM} $archive"
echo "rmcmd:  $rmcmd"
${RM} $archive

exit $?
