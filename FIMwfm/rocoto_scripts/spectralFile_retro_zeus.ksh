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

# XUE comments out 20140408
#DATEDIR=${MONTH}
#echo "DATEDIR: ${DATADR2}/${DATEDIR}"


#cpcmd="${CP} ${DATADR2}/${DATEDIR}/${YEAR}${MONTH}${DAY}${HOUR}00.tar.gz ${workdir}"
#echo "cpcmd: ${cpcmd}"
#${cpcmd}

#cd ${workdir}
#${TAR} xvfz ${YEAR}${MONTH}${DAY}${HOUR}00.tar.gz
#${RM} ${YEAR}${MONTH}${DAY}${HOUR}00.tar.gz

DATEDIR=${YEAR}${MONTH}${DAY}${HOUR}
# for our gfs format
SANLFILE="${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sanl"
SFCANLFILE="${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sfcanl"

# for ncep gfs format
#SANLFILE="gfs.t${HOUR}z.sanl"
#SFCANLFILE="gfs.t${HOUR}z.sfcanl"

#cpcmd="${CP} ${DATADR2}/${DATEDIR}/${SANLFILE} ${workdir}/${SANLFILE}"
#echo "cpcmd: ${cpcmd}"
#${cpcmd}

#cpcmd="${CP} ${DATADR2}/${DATEDIR}/${SFCANLFILE} ${workdir}/${SFCANLFILE}"
#echo "cpcmd: ${cpcmd}"
#${cpcmd}

cpcmd="${CP} ${DATADR2}/${DATEDIR}00/ensics/${SANLFILE} ${workdir}/${SANLFILE}"
echo "cpcmd: ${cpcmd}"
${cpcmd}

cpcmd="${CP} ${DATADR2}/${DATEDIR}00/ensics/${SFCANLFILE} ${workdir}/${SFCANLFILE}"
echo "cpcmd: ${cpcmd}"
${cpcmd}


#cpcmd="${CP} ${DATADR2}/${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sanl ${workdir}/${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sanl"
#echo "cpcmd: ${cpcmd}"
#${cpcmd}
#cpcmd="${CP} ${DATADR2}/${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sfcanl ${workdir}/${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sfcanl"
#echo "cpcmd: ${cpcmd}"
#${cpcmd}

exit $?
