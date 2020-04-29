#!/bin/csh
set python=/usr/bin/python

module load intel
module load idl

# initialize
setenv PYTHONPATH "/pan2/projects/gfsenkf/whitaker/lib64/python"
setenv LD_LIBRARY_PATH "/lfs1/projects/fim/whitaker/lib:${LD_LIBRARY_PATH}"

# set T1=0
# set T2=12
# set fhr=$T1
# set nplot=0
# set date=2011083100
# set gribPath=/lfs1/projects/rtfim/FIMYENS/FIMrun_mvapich_p/fim_5_64_240_201108310000/
# set fimensPath=/lfs1/projects/rtfim/FIMYENS/FIMrun_mvapich_p/fim_5_64_240_201108310000/fimens
# set nmembers=4
# set tcvitalsPath=/lfs1/projects/rtfim/tcvitals

# NOTE:  t1 and t2 are now the same!!!!!

set maxfhr=168
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
set FIM_BIN=$PLOT_HOME                   ## used in ellipses
set glvl_in=$GLVL
set fcst_len_in=$T2
set rundir_in=$fimensPath
set fcst_output_int=$FCST_INTERVAL
set yyjdyhh=$yyjdyhh

echo "date: $date"
echo "gribPath: $gribPath"
echo "fimensPath: $fimensPath"
echo "nmembers: $nmembers"
echo "tcvitalsPath: $tcvitalsPath"
echo "T2: $T2"
echo "fhr: $fcst"
echo "maxfhr: $maxfhr"
echo "f_int: $f_int"
echo "glvl_in: $glvl_in"
echo "fcst_len_in: $fcst_len_in"
echo "rundir_in: $rundir_in"
echo "fcst_output_int: $fcst_output_int"

mkdir -p $fimensPath
mkdir -p $fimensPath/Atlantic
mkdir -p $fimensPath/West_Pacific
mkdir -p $fimensPath/East_Pacific

set charfhr=`printf "%03i" $fhr`        ## force to be 3 characters
echo "before call to plotens.py charfhr: $charfhr"

# ** postage stamp plots **
$python $PLOT_HOME/plotens_FIM.py $date $charfhr $gribPath $fimensPath $tcvitalsPath $nmembers 
if ($status != 0) then
  echo "ERROR in plotens_FIM.py"
  exit $status
endif

# ** probability plots **
echo "******* Probability Plots *********"
cd $fimensPath
$python $PLOT_HOME/plotprob_FIM.py $date $charfhr $gribPath $yyjdyhh $fimensPath $tcvitalsPath $nmembers $T2 
if ($status != 0) then
  echo "ERROR in plotprob_FIM.py"
  exit $status
endif

# ** ellipse Plots **
if ("$fhr" == "$maxfhr") then 
  echo "******* ELLIPSES *********"
  echo "*** csh $PLOT_HOME/plot_ellipse.csh $date $tcvitalsPath $glvl_in $fcst_len_in $nmembers $fcst_output_int $gribPath $fimensPath $PLOT_HOME $FIM_BIN"
  csh $PLOT_HOME/plot_ellipse.csh $date $tcvitalsPath $glvl_in $fcst_len_in $nmembers $fcst_output_int $gribPath $fimensPath $PLOT_HOME $FIM_BIN
else
  echo "fhr:  $fhr,  maxhr:  $maxfhr"
endif
if ($status != 0) then
  echo "ERROR in plot_ellipse.csh"
  exit $status
endif
