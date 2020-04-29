#!/bin/ksh -l

#############################################################
## split_grib1_files.ksh
##   thie script splits the full GRIB1 file into smaller
##   GRIB1 files to avoid 2GB limitation:  hybrid levels,
##   isobaric levels, fields used for Mike's post-processing,
##   surface or 2D fields, fields needed for NCL graphics,
##   fields needed for HIWPP.
##   
## JK Henderson                                  April 2014
##    added HIWPP inventory                     Sept 2014
#############################################################


#############################################################
##   create_grib1_inv fcst_hour INV_HYB1 INV_HYB2 INV_NCL INV_PRS INV_SFC INV_W2 INV_HIWPP
##
##     this function creates multiple inventory GRIB1 files
##       cloud.inv.$FCST  -- 2D fields:  CWAT, OLR
##       hyb1.inv.$FCST   -- hybrid levels:  HGT, TMP, RH, UGRD, VGRD  
##       hyb2.inv.$FCST   -- hybrid levels:  PRES, DPT, VVEL, O3MR
##       ncl.inv.$FCST    -- needed for NCL graphics (HIWPP is subset)
##       prs.inv.$FCST    -- pressure levels:
##       sfc.inv.$FCST    -- 
##       w2flds.inv.$FCST -- needed for Mike's post-processing
##
##    inputs
##      fcst      - forecast hour
##      INV_HYB1  - generate inventory for HGT, TMP, RH, UGRD, VGRD if set to T
##                  (hybrid levels)
##      INV_HYB2  - generate inventory for PRES, DPT, VVEL, O3MR  if set to T
##                  (hybrid levels)
##      INV_NCL   - generate inventory for NCL fields if set to T
##      INV_PRS   - generate inventory for isobaric levels if set to T
##      INV_SFC   - generate inventory for everything except hybrid or
##                  isobaric levels if set to T
##      INV_W2    - generate inventory for Mike's post-processing if set to T
##      INV_HIWPP - generate inventory for HIWPP if set to T
#############################################################
create_grib1_inv()
{

  echo "*************************************************************"
  echo "entering create_grib1_inv..."

  # initialize
  fcst=$1
  inv_hyb1=$2
  inv_hyb2=$3
  inv_ncl=$4
  inv_prs=$5
  inv_sfc=$6
  inv_w2=$7
  inv_hiwpp=$8
  RM=/bin/rm
  inventory=full_inventory.$fcst                  ## full inventory of GRIB1 file
  echo "Creating Inventory files for $fcst"
 
  ###  hybrid levels  ###
  ###    split into 2 files to keep under 2GB 
  hybrid1_fields="HGT|TMP|RH|UGRD|VGRD"
  egrep "($hybrid1_fields)"  $inventory | grep "hybrid" > hyb1.inv.$fcst
  hybrid2_fields="PRES|DPT|VVEL|O3MR"
  egrep "($hybrid2_fields)" $inventory | grep "hybrid" > hyb2.inv.$fcst
  
  ###  isobaric levels  ###
  egrep "mb:" $inventory > prs.inv.$fcst
  
  ###  everything except hybrid or isobaric levels  ###
  cat $inventory | grep -v hybrid | grep -v mb > sfc.inv.$fcst
  
  ### cloud fields ###
  egrep "(NLWRT|CWAT)" sfc.inv.$fcst > cloud.inv.$fcst

  ###  w2 fields (Mike's post-processing)  ###
  if [ "$inv_w2" == "T" ]; then
    mandatory_levels="1000 mb|925 mb|850 mb|700 mb|500 mb|400 mb|300 mb|250 mb|200 mb|150 mb|100 mb"
    w2flds="PWAT|MSLMA"
    egrep "($mandatory_levels)" prs.inv.$fcst | grep -v VVEL > w2flds.inv.$fcst
    egrep "TMP|RH" hyb1.inv.$fcst | grep "hybrid lev 1:" >> w2flds.inv.$fcst
    egrep "APCP|ACPCP" sfc.inv.$fcst >> w2flds.inv.$fcst
    egrep "TMP" sfc.inv.$fcst | grep "sfc" >> w2flds.inv.$fcst
    egrep "($w2flds)" sfc.inv.$fcst >> w2flds.inv.$fcst
    egrep "UGRD|VGRD" sfc.inv.$fcst | grep "10 m" >> w2flds.inv.$fcst
  fi
  
  ###  NCL fields  ###
  if [ "$inv_ncl" == "T" ]; then
    other_levels=":25 mb|:20 mb|:10 mb|:5 mb"
    nclflds="SHTFL|LHTFL|DSWRF|DLWRF|WEASD|CNWAT|PRES|POT|UGRD|VGRD|USWRF|ULWRF|TMP|NLWRT|HGT|CWAT|WIND"
    fimx_prsflds="GFLUX|CIN"
    fimx_sfcflds="MSLET|4LFTX|LFTX|KX|SX"
    egrep "($other_levels)" prs.inv.$fcst | grep -v VVEL > ncl.inv.$fcst
    egrep "($nclflds)" sfc.inv.$fcst >> ncl.inv.$fcst
    egrep "HGT" hyb1.inv.$fcst  | egrep "(hybrid lev 50|hybrid lev 1:)" >> ncl.inv.$fcst
    egrep "DPT" hyb2.inv.$fcst  | grep "hybrid lev 1:" >> ncl.inv.$fcst
    egrep "PRES" hyb2.inv.$fcst  | grep "hybrid lev 1:" >> ncl.inv.$fcst
    egrep "VVEL" prs.inv.$fcst | grep "700 mb" >> ncl.inv.$fcst
    egrep "O3MR" hyb2.inv.$fcst | grep "hybrid lev 50" >> ncl.inv.$fcst
    egrep "($fimx_prsflds)" prs.inv.$fcst >> ncl.inv.$fcst
    egrep "($fimx_sfcflds)" sfc.inv.$fcst >> ncl.inv.$fcst
  fi
  
  ###  HIWPP fields  ###
  if [ "$inv_hiwpp" == "T" ]; then
    #  currently need to remove USWRF and ULWRF since having these variables in the GRIB2 files
    #     causes Java2NetCDF to crash in NEIS  (grib table center 59 issue)
    sfcflds="PCP|SHTFL|LHTFL|DSWRF|DLWRF|MSLMA|WEASD|CNWAT|SPFH|USWRF|ULWRF|TMP|NLWRT|HGT|CWAT|PWAT|WIND"
    cp prs.inv.$fcst hiwpp.inv.$fcst
    egrep "($sfcflds)" sfc.inv.$fcst >> hiwpp.inv.$fcst
    egrep "UGRD|VGRD" sfc.inv.$fcst | grep "10 m" >> hiwpp.inv.$fcst
    egrep "HGT" hyb1.inv.$fcst  | egrep "(hybrid lev 1:)" >> hiwpp.inv.$fcst
    egrep "DPT" hyb2.inv.$fcst  | grep "hybrid lev 1:" >> hiwpp.inv.$fcst
    egrep "PRES" hyb2.inv.$fcst  | grep "hybrid lev 1:" >> hiwpp.inv.$fcst
 fi

  ## delete inventories if no longer needed
  if [ "$inv_hyb1" == "F" ]; then $RM hyb1.inv.$fcst; fi
  if [ "$inv_hyb2" == "F" ]; then $RM hyb2.inv.$fcst; fi
  if [ "$inv_prs" == "F" ]; then $RM prs.inv.$fcst; fi
  if [ "$inv_sfc" == "F" ]; then $RM sfc.inv.$fcst; fi
  
}

#############################################################
##   gen_grib1_files model rev forecast_hour gribFile GRB_NCL
##
##     this function generates multiple GRIB1 files based on
##       the different inventory files
##
##      naming convention:
##        MODEL.REV.TYPE.YYYYDDDHH0000hh.grb1
##         
##          where MODEL = FIM9, FIM8
##                REV   = revision number in repository
##                TYPE  = hyb1, hyb2, ncl, prs, sfc, w2flds
##                YYYY  = year
##                DDD   = day of year
##                HH    = runtime
##                hh    = forecast hour
##
##      wgrib file > inventory
##      wgrib -i -grib -o file.hyb file < hyb.inv
##
##    inputs
##      model 
##      rev 
##      fcst 
##      gribFile 
##      GRB_NCL
#############################################################
gen_grib1_files()
{

  echo "*************************************************************"
  echo "entering gen_grib1_files..."

  # initialize
  model=$1
  rev=$2
  fcst=$3
  gribFile=$4
  GRB_NCL=$5

  # Print out value of input variables
  echo "  model           = ${model}"
  echo "  rev             = ${rev}"
  echo "  fcst            = ${fcst}"
  echo "  gribFile        = ${gribFile}"
  echo "  GRB_NCL         = ${GRB_NCL}"
  echo

  ls -l *inv.$fcst

  for invFile in *inv.$fcst              ## loop through all inventory files
  do
    echo "================================================="
    echo "Generating GRIB1 files from forecast $invFile"
    type=`echo $invFile | cut -d. -f1`
    if [ -s $invFile ]; then
      wgrib -i -grib -o ${model}.${rev}.${type}.${gribFile}.grb1 $gribFile < $invFile
    fi
  done

  # concatenate NCL and w2flds GRIB1 files together
  if [ "$GRB_NCL" == "T" ]; then
    if [[ -s ${model}.${rev}.w2flds.${gribFile}.grb1 ]] && [[ -s ${model}.${rev}.ncl.${gribFile}.grb1 ]]; then 
      cat ${model}.${rev}.w2flds.${gribFile}.grb1 >> ${model}.${rev}.ncl.${gribFile}.grb1
    else
      echo "Can't concatenate NCL and w2flds GRIB1 files together...."
      ls -l ${model}.${rev}.*.${gribFile}.grb1 
      exit 3
    fi
  fi
}

#############################################################
##   gen_grib2_files gribFile FIM_HOME fcst
##
##     this function converts multiple GRIB1 files to GRIB2
##
##      naming convention:
##        MODEL.REV.TYPE.YYYYDDDHH0000hh.grb1
##        MODEL.REV.TYPE.YYYYDDDHH0000hh.grb2
##         
##          where MODEL = FIM9, FIM8
##                REV   = revision number in repository
##                TYPE  = hyb1, hyb2, ncl, prs, sfc, w2flds
##                YYYY  = year
##                DDD   = day of year
##                HH    = runtime
##                hh    = forecast hour
##
##    inputs
##      gribFile
##      FIM_HOME 
##      fcst
##
#############################################################
gen_grib2_files()
{

  echo "*************************************************************"
  echo "entering gen_grib2_files..."

  # initialize
  gribFile=$1
  FIM_HOME=$2
  fcst=$3
  RM=/bin/rm
  grib1dir=`pwd`
  grib2dir=`echo ${grib1dir} | sed -e s/grib1/grib2/`
  echo $grib1dir
  echo $grib2dir
  if [ ! -d $grib2dir ]; then
    echo "making $grib2dir"
    mkdir $grib2dir
  fi

  # Print out value of input variables
  echo "  gribfile         = ${gribFile}"
  echo "  FIM_HOME         = ${FIM_HOME}"
  echo "  fcst             = ${fcst}"
  echo

  # convert to GRIB2 and delete GRIB1 file
  for file in *${gribFile}.grb1
  do
    echo "converting ${file} to GRIB2...."
    gfile=`echo $file | cut -d. -f5 --complement`
    cnvgrib -nv -g12 -p40 $gfile.grb1 ${grib2dir}/$gfile.grb2
    status=$?
    if [ $status -eq 0 ]; then
      if [ ${fcst} -gt 252 ]; then
        echo "correcting fcst_len in grib2 file ${grib2dir}/${gfile}.grb2"
        echo ${FIM_HOME}/util/timestamp_fix_grib2/timestamp_fix_grib2 ${grib2dir}/${gfile}.grb2
        ${FIM_HOME}/util/timestamp_fix_grib2/timestamp_fix_grib2 ${grib2dir}/${gfile}.grb2
      fi
      # link ncl file to original filename
      if [[ "$file" == *ncl* ]]; then
        #echo "removing ${grib1dir}/${gribFile}...."
        #$RM $gribFile
        echo "deleting ${file}"
        $RM $file
        #echo "linking $file to $gribFile"
        #ln -s $file $gribFile
        cd $grib2dir
        g2file=`echo $file | sed 's/grb1/grb2/'`
        ln -s $g2file $gribFile
        cd $grib1dir
      else
        echo "deleting ${file}"
        $RM $file
      fi
    else
      echo "ERROR:  Can't convert $file to GRIB2 format!"
      return $status
    fi
  done
  echo "removing ${grib1dir}/${gribFile}...."
  $RM $gribFile
}

#############################################################
##   cleanup grib1dir
##
##     this function deletes the inventory files 
##
##     inputs
##       grib1dir
##
#############################################################
cleanup()
{

  echo "*************************************************************"
  echo "entering cleanup..."

  #initialize
  RM=/bin/rm 
  grib1dir=$1

  # Print out value of input variables
  echo "  grib1dir             = ${grib1dir}"
  echo

  # delete inventory files
  cd $grib1dir
  for f in *inv*$fcst
  do
    echo "removing $f...."
    $RM -f $f
  done

}

############################
##    main
############################

# load modules
module load wgrib
module load cnvgrib/1.4.0

# Print out value of required environment variables
echo entering split_grib1_files.ksh...
echo

# check for empty strings
if [ -z "${FIM_HOME}" ]; then
    echo "FIM_HOME is unset or set to the empty string"
    exit 1
else
  echo "FIM_HOME         = ${FIM_HOME}"
  echo "MEMBER_ID        = ${MEMBER_ID}"
  echo "yyyymmddhh       = ${yyyymmddhh}"
  echo "model            = ${model}"
  echo "rev              = ${rev}"
  echo "fcst             = ${fcst}"
  echo "gribfile         = ${gribFile}"
  echo
fi

# initialize
fcst_mod=`expr $fcst % 6`            ## determine 6-hourly forecasts
grib1dir=${FIM_RUN}/${yyyymmddhh}/post_${MEMBER_ID}/174/NAT/grib1/
INV_HYB1=F
INV_HYB2=F
INV_PRS=F
INV_SFC=F
# only generate NCL and W2 grib files for 6-hourly forecasts
if [ $fcst_mod -eq 0 ] 
then
  INV_NCL=T      
  INV_W2=T
else
  INV_NCL=F
  INV_W2=F
fi
INV_HIWPP=F

cd $grib1dir
wgrib $gribFile > full_inventory.$fcst
create_grib1_inv $fcst $INV_HYB1 $INV_HYB2 $INV_NCL $INV_PRS $INV_SFC $INV_W2 $INV_HIWPP
gen_grib1_files $model $rev $fcst $gribFile $INV_NCL
gen_grib2_files $gribFile $FIM_HOME $fcst
cleanup $grib1dir
