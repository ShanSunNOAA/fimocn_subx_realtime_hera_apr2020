#!/bin/sh 

# check for correct number of parameters
if [ $# -lt 5 ]
then
   script=`basename $0`
   echo '   Usage: $script YYYYMMDDHH model jet|alt "interval" "keep"'
   exit 1
fi

# initialize
yyyymmddhh=$1
model=$2
system=$3
interval=$4
keep=$5
ucmodel=`echo $model | tr "a-z" "A-Z"`
ucsystem=`echo $system | tr "a-z" "A-Z"`
if [[ "$ucmodel" = *ALT* && $ucsystem = "JET" ]]; then
  echo "error - incompatible model and system: ${model} ${system}"
  exit 1
fi
if [[ $ucsystem = "ALT" ]]; then
  if [[ "$ucmodel" = *ALT* ]]; then
    ucmodelsys=$ucmodel
  else
    ucmodelsys=$ucmodel$ucsystem
  fi
else
  ucmodelsys=$ucmodel
fi
echo "ucmodelsys is $ucmodelsys"
FIMHOME=/home/rtfim/$ucmodelsys
script=${FIMHOME}/FIMwfm/rocoto_scripts/fimziprsync.manual.rb
FIMRUN=${FIMHOME}/FIMrun
lcmodel=`echo $model | tr "A-Z" "a-z"`
if [[ "$lcmodel" = *alt* ]]; then
  lcmodel=${lcmodel%alt}
  modelsys=${lcmodel}_${system}
else
  modelsys=${lcmodel}_${system}
fi
server=clank.fsl.noaa.gov
logfile=${FIMHOME}/FIMwfm/log/sync/sync_${model}

echo "$script $FIMRUN $modelsys $server $interval $keep $yyyymmddhh"

$script "$FIMRUN" "$modelsys" "$server" "$interval" "$keep" "$yyyymmddhh" > ${logfile}_manual.${yyyymmddhh}00 2>&1

echo "done!"


