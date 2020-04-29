#!/bin/ksh --login

## this script submit ncl jobs for each domain listed in GRID_NAMES 
##
##   J. Henderson    07/2014

# Print out value of required environment variables
echo entering ncl_driver.ksh....
date
echo "WFM              = ${WFM}"
echo "PROJECT          = ${PROJECT}"
echo "PARTITION        = ${PARTITION}"
echo "PP_RESERVATION   = ${PP_RESERVATION}"
echo "HOUR             = ${HOUR}"
echo "VMEM             = ${VMEM}"
echo "IS               = ${IS}"
echo "FIM_HOME         = ${FIM_HOME}"
echo "SCRIPTS          = ${SCRIPTS}"
echo "MEMBER_ID        = ${MEMBER_ID}"
echo "ATCFNAME         = ${ATCFNAME}"
echo "yyyymmddhhmm     = ${yyyymmddhhmm}"
echo "PES              = ${PES}"
echo "GLVL             = ${GLVL}"
echo "NVL              = ${NVL}"
echo "T                = ${T}"
echo "T2               = ${T2}"
echo "GRID_NAMES       = ${GRID_NAMES}"

# initialize
LOGDIR=${FIM_HOME}/FIMwfm/log/ncl/

# Export variables
export IS 

#parse out domains
for domain in $(echo $GRID_NAMES | tr "D" " ")
do
  date
  export GRID_NAME=$domain 
  echo "in nclfim_driver, FCST_LENGTH = ${FCST_LENGTH}"
  export FCST_LENGTH=${FCST_LENGTH}
  echo "$jobs: Running ncl: $T2:$GRID_NAME"
  ${FIM_HOME}/FIMrun/batchTemplate-ncl >> ${LOGDIR}/ncl_NAT_${MEMBER_ID}_${GRID_NAME}_${T2}_${yyyymmddhhmm}.log 2>&1 
  date
  status=$?
  if [ ${status} -ne 0 ]; then
    echo "ncl for ${GRID_NAME} failed!  Exit status=${status}"
    echo "See log at  ${LOGDIR}/ncl_NAT_${MEMBER_ID}_${GRID_NAME}_${T2}_${yyyymmddhhmm}.log "
    return ${status}
  fi
done
