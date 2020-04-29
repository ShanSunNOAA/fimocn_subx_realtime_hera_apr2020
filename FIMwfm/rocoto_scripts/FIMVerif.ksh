#!/bin/ksh -l 

# Set the SGE queueing options
#$ -S /bin/ksh
#$ -pe comp 1
#$ -l h_rt=00:15:00
#$ -N verif
#$ -j y
#$ -V

createFimFileNames()
{
yyyymmddhhmm=${yyyymmddhhmm}

i=0
echo "FCST_LENGTH:  $FCST_LENGTH"

fcstsArr=""
while [[ $i -le ${FCST_LENGTH} ]]; do
   fcst=$(printf "%03d" $i)
   echo "i: $i fcst: $fcst"
   fcstsArr="$fcstsArr $fcst"
   ((i+=12))
done

echo "**** fcstsArr: $fcstsArr"

yr=$(expr substr $yyyymmddhhmm 1 4)
mm=$(expr substr $yyyymmddhhmm 5 2)
dd=$(expr substr $yyyymmddhhmm 7 2)
hh=$(expr substr $yyyymmddhhmm 9 2)
echo in createFimFileNames yyyymmddhhmm: $yyyymmddhhmm

fileNames=""
for fcst in $fcstsArr ; do 
    dirDate=`date +%Y%m%d%H -u -d "$mm/$dd/$yr $hh:00 $fcst hours ago" `
    dir=${FIM_RUN}/fim_${FIM_VERIF_CYCLE}_${dirDate}00/${FIM_POST_DIR_NAME}
    fname=`date +%y%j%H -u -d "$mm/$dd/$yr $hh:00 ${fcst} hours ago" `
    file=${dir}/${fname}000${fcst}
    len=`printf $fileNames | wc -c`
    if [ len -gt 5 ] ; then
      fileNames=$fileNames:$file
    else
      fileNames=$file
    fi
done
}


# Set up paths to shell commands
RM=/bin/rm
MKDIR=/bin/mkdir

# Print out value of required environment variables
echo
echo "FIM_RUN            = ${FIM_RUN}"
echo "FIM_POST_DIR_NAME  = ${FIM_POST_DIR_NAME}"
echo "FIM_VERIF_CYCLE    = ${FIM_VERIF_CYCLE}"
echo "FIMDIR             = ${FIMDIR}"
echo "yyyymmddhhmm       = ${yyyymmddhhmm}"
echo "DIFFGB             = ${DIFFGB}"
echo "MODEL              = ${MODEL}"
echo "ANX_FILE_NAME      = ${ANX_FILE_NAME}"
echo "ANX_DIR            = ${ANX_DIR}"
echo "ANX_MODEL          = ${ANX_MODEL}"
echo "DBI_DSN            = ${DBI_DSN}";
echo "DBI_USER           = ${DBI_USER}";
echo "DBI_PASS           = ${DBI_PASS}";
echo "CLIMATE_FILE       = ${CLIMATE_FILE}";
echo "DIFFGB             = ${DIFFGB}";
echo "SCRIPTS            = ${SCRIPTS}";
echo "FCST_LENGTH        = ${FCST_LENGTH}";
echo

createFimFileNames
echo fileNames: $fileNames

validTime=`date +%Y-%m-%d -u -d "$mm/$dd/$yr $hh:00"`
validHour=$(expr substr $yyyymmddhhmm 9 2)

echo validHour $validHour
echo "after call"
echo IN MAIN fileNames:  ${fileNames}

export MODEL_FILE_NAMES=${fileNames}
export VALID_TIME=${validTime}
export VALID_HOUR=${validHour}

# Run the Perl script that does the verification
${SCRIPTS}/verif.pl
exit $?
