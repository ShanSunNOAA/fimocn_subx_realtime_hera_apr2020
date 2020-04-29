#!/bin/ksh -l

# Set the SGE queueing options
#$ -S /bin/ksh
#$ -pe comp 1
#$ -l h_rt=00:30:00
#$ -N GFSfile
#$ -j y
#$ -V

# Set up paths to shell commands
RM=/bin/rm
MKDIR=/bin/mkdir
GUNZIP=/bin/gunzip
GTAR=/bin/gtar
TAR=/bin/tar 
CP=/bin/cp 

# Print out value of required environment variables
echo
echo "FIM_RUN      = ${FIM_RUN}"
echo "CNVGRIB      = ${CNVGRIB}"
echo "yyyymmddhhmm = ${yyyymmddhhmm}"
echo "GLVL         = ${GLVL}"
echo "NVL          = ${NVL}"
echo "PES          = ${PES}"
echo "ANX_MODEL    = ${ANX_MODEL}"
echo

# load modules
module load cnvgrib/1.4.0

# Set up a work directory and cd into it
workDir=${FIM_RUN}/fim_${GLVL}_${NVL}_${PES}_${yyyymmddhhmm}/verif/${ANX_MODEL}
# ${RM} -rf ${workDir}
${MKDIR} -p ${workDir}
cd ${workDir}
echo `pwd`

yyyymmddhhmm=${yyyymmddhhmm}

yr=$(expr substr $yyyymmddhhmm 1 4)
mm=$(expr substr $yyyymmddhhmm 5 2)
dd=$(expr substr $yyyymmddhhmm 7 2)
hh=$(expr substr $yyyymmddhhmm 9 2)

fileYYMMDDHH=`date +%Y%m%d%H -u -d "$mm/$dd/$yr $hh:00" `
fileYYJDYHR=`date +%y%j%H -u -d "$mm/$dd/$yr $hh:00" `
fileDateInSeconds=`date +%s -u -d "$mm/$dd/$yr $hh:00" `

outName="${workDir}/${fileYYJDYHR}000000.grib1"
# call will return 0 or 1
grib1size=`ls -l ${outName} 2>/dev/null | wc -l`
echo outName: ${outName} grib1size: ${grib1size}
if [ ${grib1size} -eq 1 ] ; then
  exit;
fi

gfsDir="${GFS_VERIF_DIR}"
fName=${fileYYJDYHR}000000
file=${fName}.grib1
if [[ ! -s "${workDir}/${file}" ]] ; then
        cmd="${CP} ${gfsDir}/$fName ${workDir}"
        echo "CPCMD: ${cmd}"
        ${cmd}
        `${CP} ${gfsDir}/$fName ${workDir}`
        cmd="${CNVGRIB} -g21 ${workDir}/${fName} ${workDir}/${file} "
        echo "CNVGRIBCMD: ${cmd}"
        ${cmd}
else
  exit 1;
fi
exit $?
