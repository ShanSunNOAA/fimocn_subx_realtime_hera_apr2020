#!/bin/ksh

# TODO Move some or all of this into Makefile.

HW=$1

fail() { test -n "$1" && print "ERROR: $@"; exit 1; }

link() { ln -sf $1 $2           || fail; }
sync() { rsync -a --delete $@ . || fail; }

ncep_root="../framework/nems"
rm="rm -fv"
topfile="FIM_HORIZONTAL_OBJS_TOP"

# Sync and link files in preparation for build.

sync ../../cntl/*.F90
sync ../../cntl/fimnamelist.exec
sync ../../prep/ss2icos/*.F90
sync ../../utils/module_initial_chem_namelists.F90
sync ../../utils/wtinfo.F90
sync ../column_chem/module_chemvars.F90
sync ../column_chem/module_initial_chem_namelist_defaults.F90
sync ../wrfphys/module_wrfphysvars.F

link module_wrfphysvars.F module_wrfphysvars.F90

# Prep for GPU build.

if [[ "$HW" = "gpu" ]]; then
  sync cuda/*
  for x in cuda/*; do $rm ${$(basename $x)%.*}.F90; done
fi

# Optionally prep for NEMS build, and link FIM_HORIZONTAL_OBJS_TOP.

if [[ -n "$NEMS" ]]; then
  test "$HW" = "GPU" && fail "Cannot build for both GPU and NEMS."
  for x in $(grep OBJS_TOP $topfile.fim | cut -d= -f2); do $rm ${x%.*}.F90; done
  for x in $ncep_root/*; do sync $x; done
  link $topfile.ncep $topfile
else
  link $topfile.fim $topfile
fi
