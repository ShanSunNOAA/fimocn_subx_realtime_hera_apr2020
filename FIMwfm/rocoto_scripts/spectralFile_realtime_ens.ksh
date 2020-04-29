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
echo

# Set up a work directory and cd into it
workdir=${FIM_HOME}/FIMrun/${yyyymmdd}${HOUR}/ensics_${MEMBER_ID}
echo ${workdir}
${RM} -rf ${workdir}
${MKDIR} -p ${workdir}
cd ${workdir}

#
#   /public/data/grids/gefs/init/YYYYMMDD/HH/
#       gec00.tHHz.sfcanl
#       gep01.tHHz.sanl, ..., gep10.tHHz.sanl
echo "getting from public"
cpcmd="${CP} ${DATADR2}/${yyyymmdd}/${HOUR}/gec00.t${HOUR}z.sfcanl ${workdir}/${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sfcanl"
echo "cpcmd: ${cpcmd}"
${cpcmd}
cpcmd="${CP} ${DATADR2}/${yyyymmdd}/${HOUR}/gep${MEMBER_ID}.t${HOUR}z.sanl ${workdir}/${YR}${JDY}${HOUR}00.gfs.t${HOUR}z.sanl"
echo "cpcmd: ${cpcmd}"
${cpcmd}

exit $?
