#!/bin/csh

# check for number of arguments
#
if ($#argv != 5) then
  echo "$0 yyyymmddhh fcst PLOTENS PLOTPROB PLOTELL"
  exit
endif

set python=/usr/bin/python

module load intel
module load idl

setenv PYTHONPATH "/pan2/projects/gfsenkf/whitaker/lib64/python"
setenv LD_LIBRARY_PATH "/lfs1/projects/fim/whitaker/lib:${LD_LIBRARY_PATH}"

# set parameters
set yyyymmddhh=$1
set fcst=$2
set PLOTENS=$3
set PLOTPROB=$4
set PLOTELL=$5

setenv FIM_HOME /pan2/projects/fim-njet/FIMENS
setenv PLOT_HOME ${FIM_HOME}/FIMwfm/ensplots/
setenv FIM_RUN ${FIM_HOME}/FIMrun_sjet_p
setenv GLVL 7
setenv NVL 42
setenv PES 160
setenv NUM_ENS_MEMBERS 10
setenv TCVITALS /pan2/projects/fim-njet/tcvitals
setenv FCST_INTERVAL 6
setenv maxfhr 168

set ymd=`echo ${yyyymmddhh} | cut -c-8`
set hh=`echo ${yyyymmddhh} | cut -c9-10`
set yydddhh=`date --date=$ymd +%y%j`$hh
set T2=`printf %03d $fcst`
set maxfhr=$maxfhr
set date=$yyyymmddhh
set gribPath=$FIM_RUN/fim_${GLVL}_${NVL}_${PES}_${yyyymmddhh}00/
set fimensPath=$FIM_RUN/fim_${GLVL}_${NVL}_${PES}_${yyyymmddhh}00/fimens/${yyyymmddhh}/
set nmembers=$NUM_ENS_MEMBERS
set tcvitalsPath=$TCVITALS
set T2=$T2
set fhr=$fcst
set f_int=$FCST_INTERVAL
set nplot=0
set PLOT_HOME=$PLOT_HOME
set glvl_in=$GLVL
set fcst_len_in=$T2
set rundir_in=$fimensPath
set fcst_output_int=$FCST_INTERVAL
set yydddhh=$yydddhh
set FIM_BIN=$FIM_HOME/FIMwfm/ensplots

echo "date: $date"
echo "gribPath: $gribPath"
echo "fimensPath: $fimensPath"
echo "nmembers: $nmembers"
echo "tcvitalsPath: $tcvitalsPath"
echo "T2: $T2"
echo "fhr: $fhr"
echo "maxfhr: $maxfhr"
echo "f_int: $f_int"
echo "glvl_in: $glvl_in"
echo "fcst_len_in: $fcst_len_in"
echo "rundir_in: $rundir_in"
echo "fcst_output_int: $fcst_output_int"
echo "FIM_BIN: $FIM_BIN"

mkdir -p $fimensPath
mkdir -p $fimensPath/Atlantic
mkdir -p $fimensPath/West_Pacific
mkdir -p $fimensPath/East_Pacific

set charfhr=`printf "%03i" $fhr`        ## force to be 3 characters
echo "before call to plotens.py charfhr: $charfhr"

# ** postage stamp plots **
if ($PLOTENS == "T") then
  $python $PLOT_HOME/plotens_FIM.py $date $charfhr $gribPath $fimensPath $tcvitalsPath $nmembers 
endif

# ** probability plots **
if ($PLOTPROB == "T") then
  cd $fimensPath
  $python $PLOT_HOME/plotprob_FIM.py $date $charfhr $gribPath $yydddhh $fimensPath $tcvitalsPath $nmembers $T2 
endif

# ** ellipse Plots **
if ($fhr == $maxfhr && $PLOTELL == "T") then 
  echo "*** csh $PLOT_HOME/plot_ellipse.csh $date $tcvitalsPath $glvl_in $fcst_len_in $nmembers $fcst_output_int $gribPath $fimensPath $PLOT_HOME $FIM_BIN"
  csh $PLOT_HOME/plot_ellipse.csh $date $tcvitalsPath $glvl_in $fcst_len_in $nmembers $fcst_output_int $gribPath $fimensPath $PLOT_HOME $FIM_BIN
endif
