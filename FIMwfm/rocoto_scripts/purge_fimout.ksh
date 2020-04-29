#!/bin/ksh --login

## this script deletes fimout* files for each forecast hour
##
##   J. Henderson    06/2014

# Print out value of required environment variables
echo entering purge_fimout.ksh....
echo "FIMRUN         = ${FIMRUN}"
echo "yyyymmddhh     = ${yyyymmddhh}"
echo "T              = ${T}"

# Delete fim_out files for specified cycle
 fimout_files=`printf "fim_out*%04dhr" $T`
 cd ${FIMRUN}/${yyyymmddhh}/fim_C
 cmd=`rm -f $fimout_files`
 echo "Removing $fimout_files for ${FIMRUN}/${yyyymmddhh}/fim_C"
 $cmd
 status=$?
 if [ $status != 0 ]; then
   echo "error $status deleting fimout files for time $T "
   return $status
 fi
