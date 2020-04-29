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
echo "FIM_RUN              = ${FIM_RUN}"
echo "FIM_VERIF_CYCLE      = ${FIM_VERIF_CYCLE}"
echo "CNVGRIB              = ${CNVGRIB}"
echo "yyyymmddhhmm         = ${yyyymmddhhmm}"
echo "GFS_MSS_ROOT_DIR     = ${GFS_MSS_ROOT_DIR}"
echo "GFS_ONE_DEG_DIR      = ${GFS_ONE_DEG_DIR}"
echo "GFS_ONE_HALF_DEG_DIR = ${GFS_ONE_HALF_DEG_DIR}"
echo "GFS_PUBLIC_DIR       = ${GFS_PUBLIC_DIR}"
echo "ANX_MODEL            = ${ANX_MODEL}"
echo "ANX_DIR              = ${ANX_DIR}"
echo "ANX_FILE_NAME        = ${ANX_FILE_NAME}"
echo "T1                   = ${T1}"
echo

# Set up a work directory and cd into it
# workDir=${FIM_RUN}/fim_${FIM_VERIF_CYCLE}_${yyyymmddhhmm}/verif/${ANX_MODEL}
workDir=${ANX_DIR}
# ${RM} -rf ${workDir}
${MKDIR} -p ${workDir}
cd ${workDir}

echo `pwd`

yyyymmddhhmm=${yyyymmddhhmm}

yr=$(expr substr $yyyymmddhhmm 1 4)
mm=$(expr substr $yyyymmddhhmm 5 2)
dd=$(expr substr $yyyymmddhhmm 7 2)
hh=$(expr substr $yyyymmddhhmm 9 2)

fcst=$T1

fileYYYYMMDDHH=`date +%Y%m%d%H -u -d "$mm/$dd/$yr $hh:00 $fcst hours" `
fileYYJDYHR=`date +%y%j%H -u -d "$mm/$dd/$yr $hh:00 $fcst hours" `
fileDateInSeconds=`date +%s -u -d "$mm/$dd/$yr $hh:00 $fcst hours" `

echo fileYYYYMMDDHH: ${fileYYYYMMDDHH}
echo fileYYJDYHR: ${fileYYJDYHR}

outName="${workDir}/${fileYYJDYHR}000000.grib1"
# call will return 0 or 1
grib1size=`ls -l ${outName} 2>/dev/null | wc -l`
echo outName: ${outName} grib1size: ${grib1size}
if [ ${grib1size} -eq 1 ] ; then
  exit;
fi

now=`date +%s`
sevenDays=$((7 * 86400))

# try to get from default directory
verifDir="/lfs1/projects/rtfim/verif/GFS/${fileYYYYMMDDHH}"
fName="${fileYYJDYHR}000000.grib1"
if [ -s "${verifDir}/${fName}" ]; then
   cmd="${CP} ${verifDir}/${fName} ${workDir}/${fName}"
   echo "cmd: $cmd"
   ${cmd}
   error=$?
   if [[ $error -eq 0 ]]; then
     return 0
   fi
fi


diff=$(($now - $fileDateInSeconds))

if [ ${diff} -le ${sevenDays} ]; then
      gfsDir="${GFS_PUBLIC_DIR}"
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
     fi
else 
    yr=$(expr substr $fileYYYYMMDDHH 1 4)
    mm=$(expr substr $fileYYYYMMDDHH 5 2)
    dd=$(expr substr $fileYYYYMMDDHH 7 2)
    hh=$(expr substr $fileYYYYMMDDHH 9 2)
    if [ ${yr} -ge 2007 ]; then
      gfsdir="${GFS_MSS_ROOT_DIR}/${yr}/${mm}/${dd}/${GFS_ONE_HALF_DEG_DIR}"
    else 
      gfsdir="${GFS_MSS_ROOT_DIR}/${yr}/${mm}/${dd}/${GFS_ONE_DEG_DIR}"
    fi
    if [[ ! -s "${workDir}/${fileYYJDYHR}000000.grib1" ]] ; then

       cpcmd="${CP} ${gfsdir}/${fileYYYYMMDDHH}00.tar.gz ${workDir}/${fileYYYYMMDDHH}00.tar.gz"
       echo "cpcmd: ${cpcmd}"
       ${cpcmd}

       unzipcmd="${GUNZIP} ${workDir}/${fileYYYYMMDDHH}00.tar.gz"
       echo "unzipcmd: ${unzipcmd}"
       ${unzipcmd}

       tarcmd="${TAR} -xvf ${workDir}/${fileYYYYMMDDHH}00.tar"
       echo "tarcmd: ${tarcmd}"
       ${tarcmd}

       cnvgribcmd="${CNVGRIB} -g21 ${workDir}/${fileYYJDYHR}000000 ${workDir}/${fileYYJDYHR}000000.grib1"
       echo "cnvgribcmd: ${cnvgribcmd}"
       ${cnvgribcmd}
    fi
fi

if [[ -d $verifDir ]]; then
   echo "$verifDir exists"
else
  cmd="${MKDIR} $verifDir"
  echo "cmd: $cmd"
  $cmd
  error=$?
fi

if [[ -d $verifDir ]]; then
  cpcmd="${CP} ${workDir}/${fileYYJDYHR}000000.grib1 ${verifDir}/${fileYYJDYHR}000000.grib1"
  echo "cpcmd: $cpcmd"
  $cpcmd
  error=$?
  echo $error
else
   echo "cannot $cpcmd $verifDir does not exist"
fi

exit $?
