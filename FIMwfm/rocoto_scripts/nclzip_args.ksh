#!/bin/ksh -l

## this interactive script zips all *png files in ${NCLDIR}_C/domain directories
##
##    ./nclzip_args.ksh /pan2/projects/fim-njet/FIM/FIMrun_mvapich_p ncl|ncldiff  yyyymmddhh
##
##   J. Henderson    03/24/2014

# check for correct number of parameters
if [ $# -lt 3 ]; then
  echo "Usage:  $0 /pan2/projects/fim-njet/FIM/FIMrun_mvapich_p ncl|ncldiff  yyyymmddhh"
  exit 1
fi

# initialize
GRID_NAMES=fimD174D201D130D244D83D129
MEMBER_ID=C
FIM_RUN=$1
NCLDIR=$2
yyyymmddhh=$3
grids=$(echo $GRID_NAMES|sed 's/D/ /g')

# create files.zip file in each domain directory
# domains 130 and 174 have other sub-domains
#     -n  no compression
for GRID_NAME in $grids
do

  echo processing ${GRID_NAME}...

  if [[ "$GRID_NAME" = "130" ]]; then
    for SUB_DIR in 130 t1 t2 t3 t4 t5 t6 t7
    do
      dir=${FIM_RUN}/$yyyymmddhh/${NCLDIR}_${MEMBER_ID}/$SUB_DIR
      echo "dir is $dir" 
      if [[ -d ${dir} ]]; then
        echo "zipping  $GRID_NAME"
        cd ${dir} 
        if [ -f *.png ]; then zip -n .png files.zip * -i \*.png; fi
      else
        echo "$dir not found!"
      fi
    done
  fi

  if [[ "$GRID_NAME" = "174" ]]; then
    for SUB_DIR in africa e_pacific europe floating w_pacific cambodia
    do
      dir=${FIM_RUN}/$yyyymmddhh/${NCLDIR}_${MEMBER_ID}/$SUB_DIR
      echo "dir is $dir" 
      if [[ -d ${dir} ]]; then
        echo "zipping  $GRID_NAME"
        cd ${dir} 
        if [ -f *.png ]; then zip -n .png files.zip * -i \*.png; fi
      else
        echo "$dir not found!"
      fi
    done
  fi

  if [[ "$GRID_NAME" != "130" && "$GRID_NAME" != "174" ]]; then
    dir=${FIM_RUN}/$yyyymmddhh/${NCLDIR}_${MEMBER_ID}/${GRID_NAME}
    echo "dir is $dir" 
    if [[ -d ${dir} ]]; then
      echo "zipping  $GRID_NAME"
      cd ${dir} 
      if [ -f *.png ]; then zip -n .png files.zip * -i \*.png; fi
    else
      echo "$dir not found!"
    fi
  fi

done
