#!/bin/ksh -l

# Set the SGE queueing options
#$ -S /bin/ksh
#$ -pe comp 1
#$ -l h_rt=00:30:00
#$ -N surface
#$ -j y
#$ -V

# Set up paths to shell commands
RM=/bin/rm
MKDIR=/bin/mkdir



# Print out value of required environment variables
echo
echo "FIM_RUN         = ${FIM_RUN}"
echo "SCRIPTS_DIR     = ${SCRIPTS_DIR}"
echo "yyjjjhh         = ${yyjjjhh}"
echo "FILENAME        = ${FILENAME}"
echo "WORK_DIR        = ${WORK_DIR}"
echo "MODEL           = ${MODEL}"
echo "FIM_SURFACE     = ${FIM_SURFACE}"
echo "DBI_USER        = ${DBI_USER}"
echo "DBI_PASS        = ${DBI_PASS}"
echo "DBI_DSN_MADIS   = ${DBI_DSN_MADIS}"
echo "DBI_DSN_SURFACE = ${DBI_DSN_SURFACE}"
echo "FILENAME        = ${FILENAME}"
echo

${MKDIR} -p -m 777 ${WORK_DIR}
#cd ${workdir}

# Run the Perl script that does the verification
#./retro.pl FIMRETRO_r4720 1381104000 1381104000
#${SCRIPTS_DIR}/surface.pl $MODEL $yyjjjhh $yyjjjhh

#echo "/home/rtfim/surface_verif/surface.pl $MODEL $yyjjjhh $yyjjjhh"
#/home/rtfim/surface_verif/surface.pl $MODEL $yyjjjhh $yyjjjhh

#${SCRIPTS_DIR}/surface.pl $MODEL $yyjjjhh $yyjjjhh

${SCRIPTS_DIR}/surface.pl $MODEL $yyjjjhh $yyjjjhh $SCRIPTS_DIR $RUN_DIR
exit $?
