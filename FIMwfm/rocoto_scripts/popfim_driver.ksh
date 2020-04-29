#!/bin/ksh --login

CONTEXT="popfim_driver"

## this script submit pop jobs for each domain listed in GRID_NAMES 
##
##   J. Henderson    06/2014

# Print out value of required environment variables
echo entering popfim_driver.ksh....
date
echo "WFM              = ${WFM}"
echo "GRIBOUT          = ${GRIBOUT}"
echo "fimnamelist      = ${fimnamelist}"
echo "FIM_HOME         = ${FIM_HOME}"
echo "FIM_RUN          = ${FIM_RUN}"
echo "SCRIPTS          = ${SCRIPTS}"
echo "MEMBER_ID        = ${MEMBER_ID}"
echo "yyyymmddhhmm     = ${yyyymmddhhmm}"
echo "yyyymmddhh       = ${yyyymmddhh}"
echo "FCST_INTERVAL    = ${FCST_INTERVAL}"
echo "T                = ${T}"
echo "T1               = ${T1}"
echo "T2               = ${T2}"
echo "GRID_NAMES       = ${GRID_NAMES}"
echo "NLFILE           = ${NLFILE}"

# source functions.ksh
cd $FIM_RUN
. $FIM_RUN/functions.ksh
ksh_insist # Ensure that we are running in ksh93 

# Run batchTemplate-setup if it has not already been run.
test -z "$batchTemplate_setup_ran" && xsource ./batchTemplate-setup

# Enter the appropriate run directory (as defined by batchTemplate-setup).

FIMRUN="$PWD"
cd $DIR || fail "Cannot cd to $DIR."

# Make the post directory. For WFM runs, the post directory should already
# exist except for the first incremental batch and mkdir's -p option exits
# with success if the directory already exists.

mkdir -p $POST || fail "Cannot make directory $POST."

# Copy namelist from the appropriate fim directory.
if [[ ! -f $POST/$NLFILE ]]
then
  if [[ -d "$FIMDIR" ]]
  then
    cp $FIMDIR/$NLFILE $POST || fail "Cannot copy $FIMDIR/$NLFILE -> $POST."
  else
    cp $fimnamelist $POST/$NLFILE || \
      fail "Cannot copy $fimnamelist -> $POST/$NLFILE."
  fi
fi

if [[ ! -f $POST/fim_gribtable ]]
then
  cp $fimgribtable $POST/fim_gribtable || \
    fail "Cannot copy $fimgribtable $POST/fim_gribtable."
fi

if [[ ! -f $POST/REDUCEinput ]]
then
  cp $reduceinput $POST/REDUCEinput || \
  fail "Cannot copy $reduceinput $POST/REDUCEinput."
fi

# Enter the post directory.

cd $POST || fail "Cannot cd to $POST."

get_nl_value "$fimnamelist" ISOBARICnamelist isobaric_levels_file ISOBARIC_LEVELS_FILE
test -f $ISOBARIC_LEVELS_FILE || cp $PREP/$ISOBARIC_LEVELS_FILE $POST/$ISOBARIC_LEVELS_FILE 

# Link files.

test -f pop || linksafe $BINDIR/pop
test -f reduce || linksafe $BINDIR/reduce
test -f pop_read_init || linksafe $BINDIR/pop_read_init
test -f "$INFO_FILE" || linksafe $PREP/$INFO_FILE
for f in $PREP/grid_???_coeffs
do
 file=`basename $f`
 test -f $file || linksafe $f
done

# initialize
LOGDIR=${FIM_HOME}/FIMwfm/log/pop/

# Set these variables from namelist file if they are not already in the
# environment.

if [[ -z "$IS" || -z "$SMOOTH_VAR" || -z "$VAR_LIST" ]]
then
  gets="-g postnamelist:is -g postnamelist:nsmooth_var -g postnamelist:var_list"
  vals=($($(findabove nml) -n -i $NLFILE $gets))
  test -z "$IS"         &&         IS="${vals[0]}"
  test -z "$SMOOTH_VAR" && SMOOTH_VAR="${vals[1]}"
  test -z "$VAR_LIST"   &&   VAR_LIST="${vals[2]}"
fi

# Export variables
export IS SMOOTH_VAR VAR_LIST FIM POST PREP BINDIR

#parse out domains
for domain in $(echo $GRID_NAMES | tr "D" " ")
do
  echo "processing domain:  $domain"
  export GRID_NAME=$domain
  export GRID_SPEC=$domain
  if [ "${GRIBOUT}" == "FALSE" ]; then
    echo "$jobs: Running pop: $T:$GRID_NAME:$GRID_SPEC"
    date
    ${FIM_RUN}/batchTemplate-pop >> $LOGDIR/pop_NAT_${MEMBER_ID}_${T}_${yyyymmddhhmm}_${GRID_NAME}.log 2>&1 
    status=$?
    if [ ${status} -ne 0 ]; then
      echo "pop FIM failed!  Exit status=${status}"
      echo "See log at  ${FIM_HOME}/FIMwfm/log/pop/pop_NAT_${MEMBER_ID}_${T}_${yyyymmddhhmm}.log "
      return ${status}
    fi
    date
  fi
done

