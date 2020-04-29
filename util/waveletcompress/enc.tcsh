#! /bin/tcsh
#
#
set LS = /bin/ls

if ($#argv != 2) then
    echo "Usage: $0 [date filter: 00, 02, ... 32] "
    echo "Compress all fim(ocn)_out_0000[date filter]dy files."
    exit 0
endif
echo $1 $2
# create a filename list file
if ($2 == 'dy') then
	$LS -1 ???_out_????0000$1$2 > fimfiles$1
else if ($2 == 'hr') then
	$LS -1 ???_out_????000$1$2  > fimfiles$1
endif

foreach file (`cat fimfiles$1`)

  set  fld2D = `echo $file | awk '{print substr($0,9,2)}'`
  echo fld2D is $fld2D
  set  typ = `echo $file | awk '{print substr($0,1,3)}'`
  echo type is $typ
  set  fld = `echo $file | awk '{print substr($0,9,4)}'`
  echo fld is $fld
  
  # 2D file for both atmosphere and ocean
        if ($fld2D == '2D') then
          $SZ_PROG $file $file.szp 
          continue
        endif
  
  # atmosphere data
        if ($typ == 'fim') then
          set a1=0.04
          set a2=0.1	 
  	if ($fld == 'tmpP' || $fld == 'th3D' || $fld == 'td3D' || $fld == 'us3D' || $fld == 'vs3D' || $fld == 'up3P' || $fld == 'vp3P') then
            set a1=0.05
            set a2=0.1
  	else if ($fld == 'hgtP') then
            set a1=0.04
            set a2=0.1
  	else if ($fld == 'ph3D') then
            set a1=0.4
            set a2=1.0
  	else if ($fld == 'rh3D' || $fld == 'rp3P') then
            set a1=0.4
            set a2=0.8
  	else if ($fld == 'pv3D' || $fld == 'pr3D') then
            set a1=0.2
            set a2=0.5
          endif
          echo $a1
          echo $a2
          $ENCODE3D $file $fld $file.wlt $a1 $a2
  
  # ocean data
        else if ($typ == 'ocn') then
          set a1=0.01
          set a2=0.03
  	if ($fld == 'temp' || $fld == 'tave' || $fld == 'uave' || $fld == 'vave' || $fld == 'saln' || $fld == 'save') then
            set a1=0.005
            set a2=0.01
  	endif
  	$ENCODE3D_OCN $file $fld $file.wlt $a1 $a2
        else
          echo warning: $file is not compressed
        endif
end
