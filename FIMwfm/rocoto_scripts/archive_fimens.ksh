#!/bin/ksh -l

module load hpss
module list

# initialize
#inpDir=/pan2/projects/fim-njet/FIMENS_2014/FIMrun/
#mssDir=FIMENS_2014
#yyyymmddhh=$1
echo "****************************"
echo RUNDIR:   $runDir
echo mssDIR:   $mssDir
echo DATE:     $yyyymmddhh

# for each directory, archive GEFS IC, GEFS grib, FIM grib, and tracker files to mass store
cd $runDir/$yyyymmddhh
echo "in $runDir/$yyyymmddhh...."
echo "Archiving ${yyyymmddhh} to mss"
cmd="htar -cPvf /BMC/fim/1year/${mssDir}/${yyyymmddhh}.tar ensics* fimens*tgz gefs* post_*/fim/NAT/grib1/1* tracker* "
$cmd
status=$?
if [ $status != 0 ] ; then
  printf "Error : [%d] when executing htar command: '$cmd'" $status
  exit $status
fi
