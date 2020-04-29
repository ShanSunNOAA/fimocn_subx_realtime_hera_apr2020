#!/bin/ksh -l

module load intel
module load nco
# Print out value of required environment variables
echo
echo "FIM_RUN = ${FIM_RUN}"
echo "YEAR    = ${YEAR}"
echo "MONTH   = ${MONTH}"
echo "DAY     = ${DAY}"
echo "HOUR    = ${HOUR}"
echo

# Set up paths to shell commands
MKDIR=/bin/mkdir
CP=/bin/cp
RM=/bin/rm
UNZIP=/usr/bin/unzip
TAR=/bin/tar

yyyy=${YEAR}
mm=${MONTH}
dd=${DAY}
hh=${HOUR}

# find the dir
wdir=${FIM_RUN}/fim${GLVL}_${NVL}_${PES}_${yyyy}${mm}${dd}${hh}
wdir_done=${FIM_RUN}/fim${GLVL}_${NVL}_${PES}_${yyyy}${mm}${dd}${hh}_done
ncfile=${NC3D_DIR}/fim${GLVL}_${yyyy}${mm}${dd}${hh}.nc
ncfile_0z=${NC3D_DIR}/fim${GLVL}_${yyyy}${mm}${dd}00.nc
icefile=${NC3D_DIR}/ice_hemi_${yyyy}${mm}${dd}${hh}.txt

end_file=$wdir/fim/timing.summary
echo "endfile= ${end_file}"
echo "nz = ${ncfile}"
echo "nc0z = ${ncfile_0z}"

if [[ -f ${end_file} && ! -f $ncfile ]]; then
  echo "qq1 do nc $ncfile ssw=$ssw"
  if [[ $hh -eq 12 ]] ; then
    st1=12
    end1=780
  elif [[ $hh -eq 18 ]] ; then
    st1=6
    end1=774
  elif [[ $hh -eq 0 ]] ; then
    st1=0
    end1=768
  elif [[ $hh -eq 6 ]] ; then
    st1=42
    end1=762
  fi

  cd ${wdir}/fim

  if [[ $hh -ne 6 ]]; then
    if [[ $ssw == 1 ]]; then
      echo "do fim2nc with ssw for ${yyyy}${mm}${dd}${hh} $wdir"
      ~Shan.Sun/bin/fim2nc_2017 -table ~Shan.Sun/bin/fim_nctable -nlfile ~Shan.Sun/myfim2nc/ssw_24hr_nl -sfile $st1 -efile $end1 >> $HOME/subx_nc.out
      ncks -Oh -d lev,34,63 -v pv3D,us3D,vs3D,mp3D,theta fim.nc pv.nc
      ncks -Oh -d intfc_index,34,64 -v pr3D              fim.nc pr.nc
      ncks -Oh -x -v pv3D,us3D,vs3D,mp3D,pr3D,theta,intfc_index,lev fim.nc nopv.nc
      ncks -h -A  pr.nc pv.nc
      ncks -h -A  pv.nc nopv.nc
      ncatted -h -a description,global,o,c,"FIMr1.1" nopv.nc
      /bin/mv nopv.nc $ncfile
      correctsize=2526080

    else
      echo "do fim2nc no ssw for ${yyyy}${mm}${dd}${hh} $wdir"
      ~Shan.Sun/bin/fim2nc_2017 -table ~Shan.Sun/bin/fim_nctable -nlfile ~Shan.Sun/myfim2nc/subx_24hr_nl -sfile $st1 -efile $end1 >> $HOME/subx_nc.out
      ncks -Oh -d lev,34,63 -v pv3D,theta fim.nc pv.nc
      ncks -Oh -x -v pv3D,theta,intfc_index,lev fim.nc final.nc
      ncks -h -A  pv.nc final.nc
      ncatted -h -a description,global,o,c,"FIMr1.1" final.nc
      /bin/mv final.nc $ncfile
      correctsize=1364064
      correctsize=1373664	# with hice
      correctsize=1192700	# value on hera
      correctsize=1185128	# value on hera
      correctsize=1179664	# value on hera
    fi
    echo "done fim2nc for ${yyyy}${mm}${dd}${hh}"
    /bin/rm pv.nc fim.nc
  elif [[ $hh -eq 6 && -f ${ncfile_0z} ]]; then
    echo "qq2 do nc $ncfile "
    if [[ $ssw == 1 ]]; then
      echo "do fim2nc with ssw for ${yyyy}${mm}${dd}${hh} $wdir"
      ~Shan.Sun/bin/fim2nc_2017 -table ~Shan.Sun/bin/fim_nctable -nlfile ~Shan.Sun/myfim2nc/ssw_24hr_nl -sfile $st1 -efile $end1 >> $HOME/subx_nc.out
      ncks -Oh -d lev,34,63 -v pv3D,us3D,vs3D,mp3D,theta fim.nc pv.nc
      ncks -Oh -d intfc_index,34,64 -v pr3D              fim.nc pr.nc
      ncks -Oh -x -v pv3D,us3D,vs3D,mp3D,pr3D,theta,intfc_index,lev fim.nc nopv.nc
      ncks -h -A  pr.nc pv.nc
      ncks -h -A  pv.nc nopv.nc
      /bin/mv nopv.nc final.nc
      correctsize=2526080
    else
      echo "do fim2nc wo ssw for ${yyyy}${mm}${dd}${hh} $wdir"
      ~Shan.Sun/bin/fim2nc_2017 -table ~Shan.Sun/bin/fim_nctable -nlfile ~Shan.Sun/myfim2nc/subx_24hr_nl -sfile $st1 -efile $end1 >> $HOME/subx_nc.out
      ncks -Oh -d lev,34,63 -v pv3D,theta fim.nc pv.nc
      ncks -Oh -x -v pv3D,theta,intfc_index,lev fim.nc final.nc
      ncks -h -A  pv.nc final.nc
      correctsize=1364064
      correctsize=1373664	# with hice
      correctsize=1192700	# value on hera
      correctsize=1185128	# value on hera
      correctsize=1179664	# value on hera
    fi # $ssw
### adding first 2 records from 0z run
    ncks -h -d time,1,2 -F $ncfile_0z part1_${yyyy}${mm}${dd}
    ncrcat -h part1_${yyyy}${mm}${dd} final.nc tmp.nc
    ncatted -h -a description,global,o,c,"FIMr1.1: first 2 records are from 0z run" tmp.nc
    /bin/mv tmp.nc $ncfile
    echo "done fim2nc for ${yyyy}${mm}${dd}${hh}"
    /bin/rm part1_${yyyy}${mm}${dd}
    echo "done fim2nc for ${yyyy}${mm}${dd}${hh}"
    /bin/rm pv.nc fim.nc
  fi # $hh -ne 6
fi

### remove run dir if nc has the right size 
if  [[ -f $ncfile && -d $wdir ]]; then
  actualsize=$(du -k $ncfile | cut -f 1)
  if  [[ $actualsize -ge $correctsize ]]; then
    ncap2 -h -s 'time[time]={0.0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5}' $ncfile tmp_$yyyy_$mm_$dd_$hh
    ncatted -Oh -a long_name,time,o,c,"time"                                      tmp_$yyyy_$mm_$dd_$hh
    ncatted -Oh -a units,time,o,c,"days since "${yyyy}"-"${mm}"-"${dd}" 00:00:00" tmp_$yyyy_$mm_$dd_$hh
    ncatted -Oh -a calendar,time,o,c,"JULIAN"                                     tmp_$yyyy_$mm_$dd_$hh
    ncatted -Oh -a ncout.F90_rev,global,d,c,                                      tmp_$yyyy_$mm_$dd_$hh
    ncatted -Oh -a case,global,d,c,                                               tmp_$yyyy_$mm_$dd_$hh
    ncatted -Oh -a interp,global,d,c,                                             tmp_$yyyy_$mm_$dd_$hh
    ncatted -Oh -a fixedgridorder,global,d,c,                                     tmp_$yyyy_$mm_$dd_$hh
    /bin/mv tmp_$yyyy_$mm_$dd_$hh $ncfile

    echo "passed nc chk ${yyyy}${mm}${dd}${hh} $actualsize $correctsize; removing run dir"
    /bin/mv ice_hemi.txt $icefile
    /bin/mv $wdir $wdir_done
    if [[ -d $wdir_done ]]; then
      /bin/rm -rf $wdir_done
    fi
  else
    echo "failed ncchk1 ${yyyy}${mm}${dd}${hh} $actualsize "
    echo "failed ncchk2 ${yyyy}${mm}${dd}${hh} $correctsize "
  fi
fi

exit $?
