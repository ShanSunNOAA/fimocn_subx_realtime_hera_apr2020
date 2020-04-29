#!/bin/sh

## initialize
FIMHOME=/pan2/projects/fim-njet
elevenDays=`date --date='11 days ago' +%s`
fifteenDays=`date --date='15 days ago' +%s`

## loop through FIM models
for m in FIM FIM9
do
  runDir=$FIMHOME/$m/FIMrun_jet_p
  logDir=$FIMHOME/$m/FIMwfm/log
  if [ $m = "FIM" ]; then 
    keep=$fifteenDays 
    ndays=15
  else 
    keep=$elevenDays 
    ndays=11
  fi
  
  echo "processing $runDir...."

  ## remove all directories if greater than $keep days
  echo DIRS over $ndays days old:   
  dirs=`find $runDir/ -maxdepth 1 -mindepth 1 -name 2\*`
  for d in $dirs
  do
    datedir=`basename $d`
    yr=`echo $datedir | cut -c 1-4`
    mm=`echo $datedir | cut -c 5-6`
    dd=`echo $datedir | cut -c 7-8`
    datestr=`echo $yr-$mm-$dd`
    sec=`date --date="$datestr 00" +%s`
    if [ $sec -lt $keep ] ; then
      echo "   removing $d..."
      /bin/rm -rf $d/
    fi
  done
  
  ## remove directories 
  ##    if < 11 days, delete ensics, fim_C, hfipPlots_*, prep_C, soundings, verif
  ##    if > 5 days,  delete all domains under post_C except for fim
  ##    
  echo DIRS over 3 days:
  find $runDir/2*/ -maxdepth 1 -mindepth 1 -type d -mtime +3 ! -name ncl_C ! -name post_C ! -name tracker_C -print
  find $runDir/2*/ -maxdepth 1 -mindepth 1 -type d -mtime +3 ! -name ncl_C ! -name post_C ! -name tracker_C -exec /bin/rm -rf {} \;
  
  echo POST DIRS over 6 days:
  find $runDir/2*/post_C -maxdepth 1 -mindepth 1 -type d -mtime +6 ! -name fim -print 
  find $runDir/2*/post_C -maxdepth 1 -mindepth 1 -type d -mtime +6 ! -name fim -exec /bin/rm -rf {} \;
  
  ## remove log files
  echo "***********************"
  echo "*****  LOG FILES  *****"
  find $logDir -mindepth 2 -maxdepth 2 -type f -mtime +10 ! -name entries -print
  find $logDir -mindepth 2 -maxdepth 2 -type f -mtime +10 ! -name entries -exec /bin/rm -rf {} \;
done
