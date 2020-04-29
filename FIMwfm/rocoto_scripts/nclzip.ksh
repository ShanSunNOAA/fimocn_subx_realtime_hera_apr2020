#!/bin/ksh -l

## this script zips all *png files in ${NCLDIR}_C/domain directories
##
##   J. Henderson    03/24/2014

# Print out value of required environment variables
echo entering nclzip.ksh...
echo "FIM_RUN    = ${FIM_RUN}"
echo "NCLDIR     = ${NCLDIR}" 
echo "MEMBER_ID  = ${MEMBER_ID}"
echo "yyyymmddhh = ${yyyymmddhh}"
echo "GRID_NAMES = ${GRID_NAMES}" 
echo

# get grid names  -- replace D with ' '
#     fimD236D201D244  ==>  fim 236 201 244
grids=$(echo $GRID_NAMES|sed 's/D/ /g')

# add kiosk to the list of grids, so those get zipped also
grids="${grids} kiosk"

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
