#!/bin/ksh --login

## this script compresses the binary model output files (fim_out* ocn_out*) 
##     *.wlt   3D
##     *.szp   2D
##
##   J. Henderson    07/2016

# Print out value of required environment variables
echo entering compress_driver.ksh
echo "yyyymmddhh     = ${yyyymmddhh}"
echo "T              = ${T}"
echo "FIMRUN         = ${FIMRUN}"
echo "OUTDIR         = ${OUTDIR}"
echo "TIMEUNIT      = ${TIMEUNIT}"

DAYDIR=${FIMRUN}/fim_${yyyymmddhh}
echo "DAYDIR         = ${DAYDIR}"

# Initialize
CP=/bin/cp
MV=/bin/mv
ENC=$FIMRUN/../util/waveletcompress/enc.tcsh
export SZ_PROG=$FIMRUN/../util/waveletcompress/sz_prog
export ENCODE3D=$FIMRUN/../util/waveletcompress/encode3d
export ENCODE3D_OCN=$FIMRUN/../util/waveletcompress/encode3d_ocn

# Call compression script  (/scratch3/BMC/fim/nmme/ocn_util/waveletcompress/encall)
mkdir -p $OUTDIR
cd $OUTDIR
for f in `ls ../???_out*${T}${TIMEUNIT}`; do ln -s $f; done
if [[ $T -eq 0 ]]; then
  $CP ../../prep/icos_grid_info_level.dat .
  $CP ../../prep/output_isobaric_levs.txt .
  $CP ../FIMnamelist .
fi

for file in `ls ???_out*${T}${TIMEUNIT}`
do

  fld2D=`echo $file | awk '{print substr($0,9,2)}'`
  typ=`echo $file | awk '{print substr($0,1,3)}'`
  fld=`echo $file | awk '{print substr($0,9,4)}'`
  echo xD,fld,type: $fld2D $fld $typ
  
  if [[ $fld2D = '2D' ]]; then # 2D file for both atmosphere and ocean
    $SZ_PROG $file $file.szp 
    /bin/rm  $file
  
  elif [[ $typ = 'fim' ]] then # atmosphere data
    a1=0.04
    a2=0.1	 
    if [[ $fld = 'tmpP' ]] || [[ $fld = 'th3D' ]] || [[ $fld = 'td3D' ]] || [[ $fld = 'us3D' ]] || [[ $fld = 'vs3D' ]] || [[ $fld = 'up3P' ]] || [[ $fld = 'vp3P' ]]; then
      a1=0.05
      a2=0.1
    elif [[ $fld = 'hgtP' ]]; then
      a1=0.04
      a2=0.1
    elif [[ $fld = 'ph3D' ]]; then
      a1=0.4
      a2=1.0
    elif [[ $fld = 'rh3D' ]] || [[ $fld = 'rp3P' ]]; then
      a1=0.4
      a2=0.8
    elif [[ $fld = 'pv3D' ]] || [[ $fld = 'pr3D' ]]; then
      a1=0.2
      a2=0.5
    fi
    $ENCODE3D $file $fld $file.wlt $a1 $a2
    /bin/rm $file
  
  elif [[ $typ = 'ocn' ]] then # ocean data
    a1=0.01
    a2=0.03
    if [[ $fld = 'temp' ]] || [[ $fld = 'tave' ]] || [[ $fld = 'uave' ]] || [[ $fld = 'vave' ]] || [[ $fld = 'saln' ]] || [[ $fld = 'save' ]]; then
      a1=0.005
      a2=0.01
    fi
    $ENCODE3D_OCN $file $fld $file.wlt $a1 $a2
    /bin/rm  $file
  else
    echo warning: $file is not compressed
  fi

done
if [[ $T -eq 780 ]]; then
  $MV ${DAYDIR} ${DAYDIR}_www
fi
