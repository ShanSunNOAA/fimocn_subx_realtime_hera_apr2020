#!/bin/ksh 

# -Example for FIM global model (with large-scale graphics only). Execute the following in the command line:

# mkdir -p /pan2/projects/hfipprd/${USER}/TEST-HFIP-GRAPHICS-SCRIPTS/
# cd /pan2/projects/hfipprd/${USER}/TEST-HFIP-GRAPHICS-SCRIPTS/
# mkdir -p TEMP FINAL
# cp /pan2/projects/hfipprd/Thiago.Quirino/HFIP-GRAPHICS/SAMPLE/globalfim*.txt ./
# /pan2/projects/hfipprd/Thiago.Quirino/HFIP-GRAPHICS/hfipGlobalModelGraphics.ksh 2012051412 FIM0 "ESRL/AMB" globalfimdata.txt TEMP/ FINAL/
# ls -l FINAL/

# Print run parameters
ECHO=/bin/echo
yyyymmddhh=$(expr substr $yyyymmddhhmm 1 10)
${ECHO}
${ECHO} "hfipPlots.ksh"
${ECHO}
${ECHO} "           T1=${T1}"
${ECHO} "           T2=${T2}"
${ECHO} "           yyyymmddhhmm=${yyyymmddhhmm}"
${ECHO} "           yyyymmddhh=${yyyymmddhh}"
${ECHO} "           FIM_RUN=${FIM_RUN}"
${ECHO} "           FCST_INTERVAL=${FCST_INTERVAL}"
${ECHO} "           HFIP_SCRIPT=${HFIP_SCRIPT}"
${ECHO}

RM=/bin/rm
MKDIR=/bin/mkdir
MV=/bin/mv

# Set up the work directory and cd into it
hfipDir="${FIM_RUN}/${yyyymmddhh}/hfipPlots"
rootDir="${FIM_RUN}/${yyyymmddhh}/hfipPlots_$T1_$T2/"
tempDir="${rootDir}/TEMP/"
echo "hfipDir: $hfipDir"
if [[ ! -d $hfipDir ]]; then
  $MKDIR -p ${hfipDir}
  status=$?
  if [ $status != 0 ]; then
    echo "error $status making $hfipDir"
    return $status
  fi
fi
echo "tempDir: $tempDir"
$RM -rf ${tempDir}
$MKDIR -p ${tempDir}
status=$?
if [ $status != 0 ]; then
  echo "error $status making $tempDir"
  return $status
fi
finalDir="${rootDir}/FINAL/"
echo "finalDir: $finalDir"
$RM -rf ${finalDir}
$MKDIR -p ${finalDir}
status=$?
if [ $status != 0 ]; then
  echo "error $status making $finalDir"
  return $status
fi


# Get yyjjjHHMM
datestr=`echo ${yyyymmddhhmm} | sed 's/^\([0-9]\{4\}\)\([0-9]\{2\}\)\([0-9]\{2\}\)\([0-9]\{2\}\)\([0-9]\{2\}\)/\1\/\2\/\3 \4\:\5/'`
yyjjjhhmm=`date +%y%j%H%M -d "${datestr}"`
yyyymmddhh=$(expr substr $yyyymmddhhmm 1 10)
gribVersion=1

t=$T1


# currently running once since T1=T2
while [ $t -le $T2 ]; do
  txtFile="${rootDir}/globalfimdata_$t.txt"
  echo "txtFile: $txtFile"
  if [[ -e "$txtFile" ]]; then
     $RM -f $txtFile
  fi
  echo "Generating plots for forecast time ${t}"
  FCST_TIME=$t
  typeset -Z3 FCST_TIME
  dir_fcst="$(echo $FCST_TIME | sed 's/^[0]*//')"
  gribName="$FIM_RUN/${yyyymmddhh}/post_${MEMBER_ID}/fim/NAT/grib1/${yyjjjhhmm}0${FCST_TIME}"
  echo "gribName: $gribName"
  trackName="$FIM_RUN/${yyyymmddhh}/tracker_${MEMBER_ID}/$dir_fcst/track.${yyyymmddhhmm}.FIM9"
  echo "trackName: $trackName"
  echo "$t,$gribVersion,$gribName" >  $txtFile
  status=$?
  if [ $status != 0 ]; then
    echo "error $status writing to $txtFile"
    return $status
  fi
  if  [[ -s $trackName ]]; then
     echo "track exists"
     cmd="$HFIP_SCRIPT $yyyymmddhh FIM9 ESRL/GSD $txtFile $tempDir $finalDir $trackName"
  else 
     cmd="$HFIP_SCRIPT $yyyymmddhh FIM9 ESRL/GSD $txtFile $tempDir $finalDir"
  fi
  echo "plot cmd: $cmd"
  $cmd
  status=$?
  if [ $status != 0 ]; then
    echo "error $status running $cmd"
    return $status
  fi
    (( t=t+${FCST_INTERVAL} ))
done


# $RM -rf ${tempDir}
echo "moving $finalDir/*tgz to $hfipDir...."
$MV $finalDir/*.tgz $hfipDir
# $RM -rf ${rootDir}

return 0
