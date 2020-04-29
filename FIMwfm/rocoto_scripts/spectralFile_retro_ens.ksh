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

# Print out value of required environment variables
echo
echo "FIM_HOME     = ${FIM_HOME}"
echo "MEMBER_ID    = ${MEMBER_ID}"
echo "yyyymmdd     = ${yyyymmdd}"
echo "YR           = ${YR}"
echo "HOUR         = ${HOUR}"
echo "JDY          = ${JDY}"
echo "DATADR2      = ${DATADR2}"
echo "GLVL         = ${GLVL}"
echo "NVL          = ${NVL}"
echo "PES          = ${PES}"
echo

year=`echo $YR | cut -c3-4`
echo "year           = ${year}"

# Set up a work directory and cd into it
workdir=${FIM_HOME}/FIMrun/fim_${GLVL}_${NVL}_${PES}_${yyyymmdd}${HOUR}00/ensics_${MEMBER_ID}
echo ${workdir}
${RM} -rf ${workdir}
${MKDIR} -p ${workdir}
cd ${workdir}

#
#   /scratch4/BMC/fim/GEFS_INIT/YYYY/YYYYMMDD
#       gec00.tHHz.sfcanl
#       gep01.tHHz.sanl, ..., gep10.tHHz.sanl
echo "getting from retro directory"
cpcmd="${CP} ${DATADR2}/${YR}/${yyyymmdd}/gec00.t${HOUR}z.sfcanl ${workdir}/${year}${JDY}${HOUR}00.gfs.t${HOUR}z.sfcanl"
echo "cpcmd: ${cpcmd}"
${cpcmd}
cpcmd="${CP} ${DATADR2}/${YR}/${yyyymmdd}/gep${MEMBER_ID}.t${HOUR}z.sanl ${workdir}/${year}${JDY}${HOUR}00.gfs.t${HOUR}z.sanl"
echo "cpcmd: ${cpcmd}"
${cpcmd}

exit $?
