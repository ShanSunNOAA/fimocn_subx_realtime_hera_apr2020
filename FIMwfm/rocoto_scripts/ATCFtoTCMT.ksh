#!/bin/ksh --login

# Print out value of required environment variables
echo
echo "yyyymmddhh         = ${yyyymmddhh}"
echo "SCRIPTS            = ${SCRIPTS}";
echo "fcst               = ${fcst}";
echo

# Run the Perl script that does the verification
${SCRIPTS}/cp_ATCFtoTCMT_retro.pl $yyyymmddhh $fcst
exit $?
