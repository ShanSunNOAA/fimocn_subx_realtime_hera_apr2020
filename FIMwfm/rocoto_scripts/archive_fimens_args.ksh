#!/bin/ksh -l

module load hpss

# check for correct number of parameters
if [ $# -lt 1 ]; then
  echo "   Usage:  $0  yyyymm[ddhh] "
  exit 1
fi

# initialize
inpDir=/pan2/projects/fim-njet/FIMENS_2014/FIMrun/
mssDir=FIMENS_2014
filter=$1

# for each directory, archive FIM binary files and grib files to mass store
cd $inpDir
for dir in ${filter}*
do
  cd $inpDir/$dir
  echo "Archiving ${dir} to mss"
  htar -cPvf /BMC/fim/1year/${mssDir}/${dir}.tar ensics* fimens*tgz gefs* post_*/fim/NAT/grib1/1* tracker*
  status=$?
  if [ $status != 0 ]; then
    printf "Error : [%d] when executing htar command: '$cmd'" $status
    exit $status
  fi
done
