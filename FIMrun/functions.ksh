# This code provides shared ksh functions for the run-automation scripts.

function check_nems
{
  # Check that run configuration is NEMS-compatible
  get_build_config
  test -z "$COMPARE_VAR_ON" && fail "$0: COMPARE_VAR_ON undefined."
  if [[ "$BUILD_CONFIG" =~ "nems" ]]
  then
    logically_true "$COMPARE_VAR_ON" && fail "Cannot use NEMS with COMPARE_VAR."
    # NEMS doesnt yet work when restart enabled
    get_nl_value "$fimnamelist" OUTPUTnamelist readrestart READRESTART
    logically_true "$READRESTART" && fail "Cannot use NEMS with READRESTART."
  fi
}

function clean_submit_dir
{
  test -z "$1" && fail "$0: No submit script name specified."
  rm $1 # the batch-submit script
  get_nl_value "$fimnamelist" cntlnamelist nvl NVL
  rm_all_but_in "dpsig.*\.txt" "dpsig${NVL}(_[0-9]+)?.txt" $PWD
  rm_all_but_in "theta_coor.*\.txt" "theta_coor${NVL}.txt" $PWD
  rm -f FIMnamelist.* git.diff
}

function compare_var_setup
{
  # See ../README for for instructions on using SMS COMPARE_VAR.
  typeset ctreq ctsum cvnt1 cvnt2 fimnl fimout gribout nwt ps ron smsnl
  test -n "$1" && fimnl=$1 || fail "$0: No namelist filename supplied."
  test -n "$2" && ctreq=$2 || fail "$0: No PES value supplied."
  test -n "$3" &&    ps=$3 || fail "$0: No parallelism value supplied."
  smsnl="SMSnamelist"
  source get_nl_value "$smsnl" smsnamelist compare_var_on COMPARE_VAR_ON
  if ( logically_true "$COMPARE_VAR_ON" ); then
    test "$ps" = "serial" && fail "Cannot use COMPARE_VAR with serial runs."
    source get_nl_value "$smsnl" smsnamelist compare_var_ntasks_1 cvnt1
    source get_nl_value "$smsnl" smsnamelist compare_var_ntasks_2 cvnt2
    let ctsum=$cvnt1+$cvnt2
    test "$ctsum" -ne "$ctreq" && fail "Cores requested in $fimnl ($ctreq)" \
      "must equal COMPARE_VAR sum ($cvnt1+$cvnt2=$ctsum) in $smsnl."
    source get_nl_value "$fimnl" writetasknamelist root_own_node ron
    logically_true "$ron" && fail "COMPARE_VAR requires root_own_node=.false."
    source get_nl_value "$fimnl" writetasknamelist num_write_tasks nwt
    test "$nwt" -gt 0 && fail "Cannot use COMPARE_VAR with write tasks."
    source get_nl_value "$fimnl" postnamelist fimout fimout
    logically_true "$fimout" && fail "COMPARE_VAR requires fimout=.false."
    source get_nl_value "$fimnl" postnamelist gribout gribout
    logically_true "$gribout" && fail "COMPARE_VAR requires gribout=.false."
    COMPARE_VAR_NTASKS_1=$cvnt1 # set global
    COMPARE_VAR_NTASKS_2=$cvnt2 # set global
  else
    logically_false "$COMPARE_VAR_ON" || fail "COMPARE_VAR_ON = $COMPARE_VAR_ON"
  fi
}

function context_peek
{
  # Print the CONTEXT name most recently pushed onto the stack.
  print $CONSTCK | sed "s/^\([^:][^:]*\).*/\1/"
}

function context_pop
{
  # Pop a value of the context stack and set CONTEXT to its value.
  CONTEXT=$(context_peek)
  CONSTCK=$(print $CONSTCK | sed "s/^$CONTEXT[:]*//")
}

function context_push
{
  # Push a context value onto the stack.
  test -z "$1" && fail "No argument to push supplied."
  test -z "$CONSTCK" && CONSTCK=$1 || CONSTCK="$1:$CONSTCK"
}

function endian_big
{
  # Enable big-endian handling of the space-separated list of logical unit
  # numbers given as the function's argument.
  # sed: replace spaces with commas
  typeset luns1=$(print "$@" | sed 's/ /,/g')
  # sed: insert '-T' at the beginning and ',-T' anywhere a space is found
  typeset luns2=$(print "$@" | sed 's/^/-T/;s/ /,-T/g')
  export F_UFMTENDIAN="big:$luns1" # intel
  export FORT90L="-Wl,$luns2" # lahey
  for i in $@; do
    export FORT_CONVERT${i}=BIG_ENDIAN # nag
  done
  export GFORTRAN_CONVERT_UNIT="big_endian:$luns1" # gfortran
  print "F_UFMTENDIAN=$F_UFMTENDIAN FORT90L=$FORT90L \
GFORTRAN_CONVERT_UNIT=$GFORTRAN_CONVERT_UNIT"
}

function endian_little
{
  # Enable little-endian handling of the space-separated list of logical unit
  # numbers given as the function's argument.
  # sed: replace spaces with commas
  typeset luns1=$(print "$@" | sed 's/ /,/g')
  export XLFRTEOPTS="ufmt_littleendian=$luns1" # ibm
  print "XLFRTEOPTS=$XLFRTEOPTS"
}

function endian_reset
{
  # Disable all endianness control variables.
  unset F_UFMTENDIAN          # intel
  unset FORT90L               # lahey
  unset GFORTRAN_CONVERT_UNIT # gfortran
  unset XLFRTEOPTS            # ibm
  for var in $(env | grep -E "^FORT_CONVERT[0-9]+=.*" | cut -d= -f1); do
    unset $var # nag
  done
}

function errhandler
{
  # Handle trapped errors. The line number where the error occurred is expected
  # as the sole argument. Print an error message including the failed line
  # number, then re-enable xtrace, which presumably was disabled by trap_on().
  # This function is meant to be called by the trap mechanism, not directly by
  # the sourcing script. It is also expected that the sourcing script will exit
  # with an informative error message after this function returns. For example:
  #
  # cd /no/such/directory || fail "Cannot cd to /no/such/directory"
  #
  # will trigger errhandler() if trap_on() has previously been called. The
  # failed line number will be reported and command passed back to the sourcing
  # script, which will then fail with a more-informative message.
  typeset LINE=$1
  print "$CONTEXT[$LINE]: An error occurred, see stdout" >&2
  set -o xtrace
}

function fail
{
  # Print a failure message and terminate.
  test $# -gt 0 && print "ERROR: $@"
  exit 1
}

function fimrun
{
  # The exit status of some mpiruns is unreliable, so disable trapping
  trap_off
  "$@" >> stdout 2>&1
  # Enable trapping
  trap_on
  # Check for completion messages in stdout file.
  fimstatus="fail"
  # normal completion message
  grep 'Program exited normally' stdout && fimstatus="ok"
  # NEMS completion message
  grep 'PROGRAM nems      HAS ENDED' stdout && fimstatus="ok"
  if [[ "$fimstatus" == "fail" ]]
  then
    fail "$FIMEXEBASE failed."
  else
    print "\n$FIMEXEBASE finished\n"
  fi
}

function findabove
{
  # Run in subshell to preserve origin directory.
  (
    until test -e $1 -o $PWD = "/" ; do cd .. ; done
    test -e $1 || fail "$0: Could not find $1"
    print $PWD/$1
  )
}

function get_build_config
{
  # Must be called from FIMrun_*, or a subdirectory thereof.
  BUILD_CONFIG=$(print $PWD | sed -e 's/^.*FIMrun_\([^/]*\).*$/\1/' -e 's/_debug$//')
  test -n "$BUILD_CONFIG" || fail "Could not determine build config."
}

function get_co2file
{
  # If a CO2 data file for the current year exists, use it. Otherwise, use the
  # 2012 data file.

  test -z "$yyyymmddhhmm" && fail "$0: yyyymmddhhmm not set"
  co2file=$DATADIR/co2historicaldata_${yyyymmddhhmm:0:4}.txt
  test -f $co2file || co2file=$DATADIR/co2historicaldata_2012.txt
}

function get_fimnamelist
{
  fimnamelist="$PWD/FIMnamelist"
}

function get_nl_value
{
  # Assigns to the specified variable the value corresponding to the given
  # namelist file, namelist and key.
  #
  # usage: get_nl_value namelist_file namelist key variable

  test -z "$1" && fail "$0: No namelist filename supplied."
  test -r "$1" || fail "$0: Cannot read namelist file '$1'."
  test -z "$2" && fail "$0: No namelist supplied."
  test -z "$3" && fail "$0: No key supplied."
  test -z "$4" && fail "$0: No variable supplied."
  typeset contents="$(cat $1)"
  if [[ "$GET_NL_VALUE_F" != "$1" || "$GET_NL_VALUE_C" != "$contents" ]]
  then
    GET_NL_VALUE_F="$1"
    GET_NL_VALUE_C="$contents"
    typeset nml=$(findabove "nml")
    test -x "$nml" || fail "Did not find 'nml' tool in or above $PWD"
    typeset nmlquery="$($nml -f ksh -i $1)"
    test -n "$nmlquery" || fail "$nml returned no function."
    eval "$nmlquery" || fail "Function returned by $nml did not evaluate."
  fi
  typeset value=$(nmlquery $2 $3 | sed -e s/^[\'\"]// -e s/[\'\"]\$//)
  test -z "$value" && fail "$0: No value found for $2:$3 in $1."
  eval "$4=\"$value\""
}

function get_parallelism
{
  parallelism="parallel"
  get_build_config
  no_build_info=$(print $BUILD_CONFIG | sed -e 's/[^_]*//')
  print $no_build_info | grep -qE "^_s" && parallelism="serial"
}

function get_pes
{
  get_nl_value "$fimnamelist" queuenamelist computetasks COMPUTETASKS
  PES="$COMPUTETASKS"
  print $PES | grep -qE "^[0-9][0-9]*$" || fail "Bad computetasks value: '$PES'"
}

function get_runmaster
{
  # Must be called from FIMrun_*, or a subdirectory thereof.
  RUNMASTER=$(print $PWD | sed -e 's/^\(.*FIMrun_[^/]*\).*$/\1/')
  test -n "$RUNMASTER" || fail "Could not determine run directory."
}

function ksh_check
{
  # Check ksh version. If it is ksh93 (regardless of its name), simply return.
  # Otherwise, see if a 'ksh93' binary is available: If so, we can later call
  # ksh_fix() to use ksh93; if not, fail with an informative message.
  test "$(ksh_version)" == "93" && return 0
  test -z $(whence ksh93) && ksh_insist
}

function ksh_fix
{
  # If 'ksh93' is available on the user's path, modify the copied run-automation
  # scripts to use it.
  ksh_check && return
  typeset KSH93 tmp x
  KSH93=$(whence ksh93)
  tmp=sed.tmp
  for x in $(ls)
  do
    grep -q "^#!.*ksh" $x || continue
    sed "s:^#\!.*ksh\(.*\):#\!$KSH93\1:" $x > $tmp || fail
    mv $tmp $x || fail
    chmod +x $x || fail
  done
}

function ksh_insist
{
  if [[ "$(ksh_version)" == "88" ]]
  then
    print "
FIM run automation requires ksh93. This appears to be ksh88. If you are running
via a queue-submission script, it may be sufficient to make 'ksh93' available on
your path (perhaps via a symbolic link) and re-run this script. If you are
calling run-automation scripts (e.g. batchTemplate-prep, via Rocoto) directly,
you may need to modify the scripts' initial #! lines to specify the path to a
ksh93 binary.
"
    fail # Comment out to permit ksh88 use, at your own risk.
  fi
}

function ksh_version
{
  # SECONDS contains decimal in ksh93, but an integer in ksh88 and pdksh.
  print $SECONDS | grep -q "\." && print 93 || print 88
}

function linksafe
{
  test -z "$1" && fail "usage: $0 target [link]"
  test -z "$2" && typeset link="$PWD" || typeset link=$2
  test -e $1 || fail "$0: $1 not found."
  ln -sf $1 $link || fail "$0: Cannot link ($link -> $1)."
}

function load_modules 
{
  get_build_config
  get_runmaster
  xsource_notrace $RUNMASTER/envset
  module list
}

function logical_false
{
  print "f"
}

function logical_true
{
  print "t"
}

function logically_false
{
  test -z "$1" && fail "$0: No value supplied."
  [[ "$1" == $(logical_false) ]]
}

function logically_true
{
  test -z "$1" && fail "$0: No value supplied."
  [[ "$1" == $(logical_true) ]]
}

function rm_all_but_in
{
  test -z "$1" && fail "$0: No include pattern specified."
  test -z "$2" && fail "$0: No exclude pattern specified."
  test -z "$3" && fail "$0: No directory specified."
  test -d "$3" || fail "$0: Not a directory: $3"
  for x in $(ls $3 | grep -E "$1"); do
    print $x | grep -qE "$2" || rm $x
  done
}

function set_nl_value
{
  # Set the value for the given namelist file, namelist and key.
  #
  # usage: set_nl_value namelist_file namelist key value

  test -f "$1" || fail "$0: Namelist file '$1' not found."
  test -z "$2" && fail "$0: No namelist supplied."
  test -z "$3" && fail "$0: No key supplied."
  test -z "$4" && fail "$0: No value supplied."
  $(findabove nml) -i $1 -o $1 -s $2:$3=$4
}

function syncfiles
{
# Usage:  syncfiles $source_dir $target_dir $file_type
#         $source_dir = source directory with files or links
#         $target_dir = copy files or links here
#         $file_type  = f for files, or l for links
  test -d "$1" || fail "Not a directory: $1"
  test -d "$2" || fail "Not a directory: $2"
  for file in $(find $1 -maxdepth 1 -type $3); do
    rsync -aLu $file $2 || fail "Cannot rsync $file -> $2."
  done
}

function trap_off
{
  # Disable error trapping.
  trap - ERR
}

function trap_on
{
  # Enable error trapping. ksh issues the fake signal ERR when any command
  # returns a non-zero status. When this function is called to enable trapping,
  # subsequent commands in the sourcing script that return a non-zero status
  # will trigger the trap. The trap 1) disables xtrace to prevent verbose
  # tracing inside the error handler, and 2) calls the error handler with the
  # line number of the failed command as its argument.
  trap 'set +o xtrace;errhandler $LINENO' ERR
}

function xsource
{
  # Save the current context, source the argument(s), then restore the previous
  # context. Works recursively due to context stack.
  test -z "$CONTEXT" && CONTEXT="unknown"
  print $CONTEXT | grep -q ':' && fail "Colon not allowed in CONTEXT '$CONTEXT'."
  context_push $CONTEXT
  test -z "$XSOURCE_NOTRACE" && set -o xtrace
  . $* || fail "Problem sourcing $1."
  set +o xtrace
  context_pop
  set -o xtrace
}

function xsource_notrace
{
  # Perform xsource but do not enable tracing of sourced script.
  XSOURCE_NOTRACE=1
  xsource $*
  unset XSOURCE_NOTRACE
}

PS4='$CONTEXT[$LINENO]: '

# Turn on error trapping.

trap_on

# Record that this script has been sourced.

functions_sourced="true"

# Turn on xtrace.

set -o xtrace
