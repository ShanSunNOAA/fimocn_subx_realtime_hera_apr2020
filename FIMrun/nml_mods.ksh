#!/bin/ksh
#---------------------------------------------------------------
# Description: Modify a copy of FIMnamelist in place for specific 
#              machines/architectures.  This script is called 
#              from makefim and the namelists are modified in 
#              directory ${RUNMASTER} after the files are copied
#              to ${RUNMASTER} and before MAKE.
#
# NOTE:  These modifications are made only once for the first build.
#        Once FIMnamelist exists in ${RUNMASTER}, it will not be 
#        modified again by the next build, to preserve namelist
#        changes made by the user.
#
# Usage: . nml_mods.ksh $RUNMASTER $MACHINE $PARTITION $OMP $SERIAL
#
# Input:
#     $RUNMASTER = FIMrun_${BUILD}_${P}${PARTITION}{OMP_SUFFIX}${DEBUG_SUFFIX}
#     $MACHINE   = jet, theia and stampede currently supported
#     $PARTITION = vjet (or null, default) tjet ujet sjet xjet 
#     $OMP       = TRUE if using openmp compiler flags, FALSE otherwise
#     $SERIAL    = s for serial builds, p for parallel builds
#---------------------------------------------------------------
RUNMASTER=$1
MACHINE=$2
PARTITION=$3
OMP=$4
SERIAL=$5

nml_file_name="${RUNMASTER}/FIMnamelist"

if [[ ${MACHINE} == "jet" || ${MACHINE} == "theia" || ${MACHINE} == "stampede" ]] ; then
  echo "Modifying ${nml_file_name} for ${MACHINE}" 
  echo "Partition is ${PARTITION}"
  echo "OMP is ${OMP}"
  echo "SERIAL is ${SERIAL}"

# root_own_node is generally false, except for higher resolution cases
# like t1534, g8 and g9, and for some test suite cases
# For now, using default from namelist

  if [[ ${MACHINE} == "stampede" ]] ; then
    root_own_node=.false.
  else
    root_own_node=.true.
  fi

  if [[ ${MACHINE} == "jet" ]] ; then
# Default from namelist
    compute_tasks=10
    if [[ ${PARTITION} == "vjet" || ${PARTITION} == "null" ]] ; then
      cpn=16
    elif [[ ${PARTITION} == "tjet" ]] ; then
      cpn=12
    elif [[ ${PARTITION} == "ujet" ]] ; then
      cpn=12
    elif [[ ${PARTITION} == "sjet" ]] ; then
      cpn=16
    elif [[ ${PARTITION} == "xjet" ]] ; then
      cpn=24
    else
      echo "ERROR:  Unrecognized partition on jet ${PARTITION}"
    fi
  fi

  if [[ ${MACHINE} == "theia" ]] ; then
# Default from namelist
    compute_tasks=4
    if [[ ${OMP} == "TRUE" ]] ; then
      cpn=48
    else
      cpn=24
    fi
  fi

  if [[ ${MACHINE} == "stampede" ]] ; then
# Default from namelist
    compute_tasks=8
    cpn=16
  fi

# Set mpipn 
  if [[ ${OMP} == "TRUE" ]] ; then
    mpipn=2
  else
    mpipn=${cpn}
  fi

  if [[ ${SERIAL} == "s" ]] ; then
    compute_tasks=1
  fi

  (( nthreads = cpn / mpipn ))

  set -A VarName ComputeTasks cpn mpipn nthreads root_own_node

  set -A VarValue "'${compute_tasks}' ! Number of compute tasks for FIM (1 for serial) " \
                  "${cpn}      ! Number of cores per node " \
                  "${mpipn}    ! Number of MPI tasks per node " \
                  "${nthreads} ! Number of threads per task if OMP threading enabled " \
                  "${root_own_node} ! whether root process has node to itself "

  let i=0
  for var_name in ${VarName[@]} ; do
    echo "In FIMnamelist, setting ${var_name} to ${VarValue[$i]}"
    sed -i "s,^.*${var_name}.*,   ${var_name} \= ${VarValue[$i]}," ${nml_file_name}
    let i+=1
  done
else
  echo "Namelists are not modified for ${MACHINE}"
fi   # jet, theia and stampede
