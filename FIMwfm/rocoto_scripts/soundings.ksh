#!/bin/ksh -l

# Set the SGE queueing options
#$ -S /bin/ksh
#$ -pe comp 1
#$ -l h_rt=00:30:00
#$ -N soundings
#$ -j y
#$ -V

# Set up paths to shell commands
RM=/bin/rm
MKDIR=/bin/mkdir



# Print out value of required environment variables
echo
echo "SCRIPTS_DIR   = ${SCRIPTS_DIR}"
echo "yyjjhh        = ${yyjjjhh}"
echo "FILENAME      = ${FILENAME}"
echo "WORK_DIR      = ${WORK_DIR}"
echo "MODEL         = ${MODEL}"
echo "FIM_SOUNDINGS = ${FIM_SOUNDINGS}"
echo "DBI_USER = ${DBI_USER}"
echo "DBI_PASS = ${DBI_PASS}"
echo "DBI_DSN_RUC_UA = ${DBI_DSN_RUC_UA}"
echo "DBI_DSN_SOUNDINGS = ${DBI_DSN_SOUNDINGS}"
echo "FILENAME = ${FILENAME}"
echo

${MKDIR} -p -m 777 ${WORK_DIR}
#cd ${workdir}

# Run the Perl script that does the verification
${SCRIPTS_DIR}/fim_grib2_sdg_1file.pl
exit $?
