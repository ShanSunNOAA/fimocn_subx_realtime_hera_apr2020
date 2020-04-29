#!/bin/ksh -l

# Set the SGE queueing options
#$ -S /bin/ksh
#$ -pe comp 1
#$ -l h_rt=00:30:00
#$ -N SpectralFiles
#$ -j y
#$ -V

# Set up paths to shell commands
RM=/bin/rm
MKDIR=/bin/mkdir
CP=/bin/cp 

# Set up a work directory and cd into it
workdir=${FIM_RUN}/fim_${GLVL}_${NVL}_${PES}_${yyyymmddhhmm}/ensics
echo ${workdir}
${RM} -rf ${workdir}
${MKDIR} -p ${workdir}
cd ${workdir}

# Print out value of required environment variables
echo
echo "FIM_RUN      = ${FIM_RUN}"
echo "yyyymmddhhmm = ${yyyymmddhhmm}"
echo "YR           = ${YR}"
echo "HOUR         = ${HOUR}"
echo "JDY          = ${JDY}"
echo "DATADR2      = ${DATADR2}"
echo "GLVL         = ${GLVL}"
echo "NVL          = ${NVL}"
echo "PES          = ${PES}"
echo

# /public/data/grids/gfs/spectral 082651200.gfs.t12z.sanl, 082651200.gfs.t12z.sfcanl 
echo "getting from public"
cpcmd="${CP} ${DATADR2}/${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sanl ${workdir}/${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sanl"
echo "cpcmd: ${cpcmd}"
${cpcmd}
cpcmd="${CP} ${DATADR2}/${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sfcanl ${workdir}/${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sfcanl"
echo "cpcmd: ${cpcmd}"
${cpcmd}

exit $?
