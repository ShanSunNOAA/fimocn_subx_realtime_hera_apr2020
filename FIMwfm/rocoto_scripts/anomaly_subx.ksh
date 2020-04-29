#!/bin/ksh -l

module load intel
module load nco
mon=`echo $MMM | tr '[A-Z]' '[a-z]'`
date=${DAY}${mon}${YEAR}
echo "date = $date hr == $HH"

day=${YEAR}${MM}${DAY}

cd ${NC_SUBX_DIR}

if [[ $HH -eq 6 && -f zg_500_FIM_${date}*_m01.nc && -f zg_500_FIM_${date}*_m02.nc && -f zg_500_FIM_${date}*_m03.nc && -f zg_500_FIM_${date}*_m04.nc ]]; then

  wdir=${NC3D_DIR}/tmp_$day
  if [[ ! -d $wdir ]]; then
    mkdir $wdir
    echo " mkdir $wdir "
  fi
  file1=$wdir/combo1.nc
  file2=$wdir/combo2.nc
  file3=$wdir/combo3.nc
  file4=$wdir/combo4.nc
  ens=ens.nc
  wk1=w1_${day}.nc
  wk2=w2_${day}.nc
  wk3=w3_${day}.nc
  wk4=w4_${day}.nc
  wk34=w34_${day}.nc
  cl1=c1_${jday}.nc
  cl2=c2_${jday}.nc
  cl3=c3_${jday}.nc
  cl4=c4_${jday}.nc
  cl34=c34_${jday}.nc
  ano1=${NC3D_DIR}/ano_wk1_${day}.nc
  ano2=${NC3D_DIR}/ano_wk2_${day}.nc
  ano3=${NC3D_DIR}/ano_wk3_${day}.nc
  ano4=${NC3D_DIR}/ano_wk4_${day}.nc
  ano34=${NC3D_DIR}/ano_wk34_${day}.nc

  echo " grep 3 fields "
    ncks -h -A zg_500_FIM_${date}*_m01.nc $file1
    ncks -h -A zg_500_FIM_${date}*_m02.nc $file2
    ncks -h -A zg_500_FIM_${date}*_m03.nc $file3
    ncks -h -A zg_500_FIM_${date}*_m04.nc $file4

    ncks -h -A tas_2m_FIM_${date}*_m01.nc $file1
    ncks -h -A tas_2m_FIM_${date}*_m02.nc $file2
    ncks -h -A tas_2m_FIM_${date}*_m03.nc $file3
    ncks -h -A tas_2m_FIM_${date}*_m04.nc $file4

    ncks -h -A pr_sfc_FIM_${date}*_m01.nc $file1
    ncks -h -A pr_sfc_FIM_${date}*_m02.nc $file2
    ncks -h -A pr_sfc_FIM_${date}*_m03.nc $file3
    ncks -h -A pr_sfc_FIM_${date}*_m04.nc $file4

  cd $wdir
   echo " match variable names for ncdiff file1=$file1 "
   ncrename -h -v zg,h500 -v pr,r12D -v tas,t22D $file1
   ncrename -h -v zg,h500 -v pr,r12D -v tas,t22D $file2
   ncrename -h -v zg,h500 -v pr,r12D -v tas,t22D $file3
   ncrename -h -v zg,h500 -v pr,r12D -v tas,t22D $file4

   ncap2 -h -O -s 'time[time]={0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 26.5, 27.5, 28.5, 29.5, 30.5, 31.5}' $file1 tmp.nc
   echo "getting ensemble= $ens "

   nces -h -O tmp.nc $file2 $file3 $file4 tmp2.nc
   ncap -h -O -s "r12D=r12D*86400"   tmp2.nc $ens
   ncatted  -h -a long_name,r12D,o,c,daily_precipitation $ens
   ncatted  -h -a units,r12D,o,c,mm/day                  $ens
   echo "getting weekly ave forecast wk1=$wk1 ens=$ens "
# getting weekly ave; first record is first day; since 0hr data are not saved in flat files
   ncwa -h -a time -d time,4,10 -F $ens $wk1
   ncatted -h -a description,global,o,c,"1-week lead forecast (day4-10)"   $wk1
   ncwa -h -a time -d time,11,17 -F $ens $wk2
   ncatted -h -a description,global,o,c,"2-week lead forecast (day11-17)"  $wk2
   ncwa -h -a time -d time,18,24 -F $ens $wk3
   ncatted -h -a description,global,o,c,"3-week lead forecast (day18-24)"  $wk3
   ncwa -h -a time -d time,25,31 -F $ens $wk4
   ncatted -h -a description,global,o,c,"4-week lead forecast (day25-31)"  $wk4
   ncwa -h -a time -d time,18,31 -F $ens $wk34
   ncatted -h -a description,global,o,c,"week3-4 lead forecast (day18-31)" $wk34

# getting weekly climo, based on smoothed climo at /ESRL/BMC/fim/2year/subx_rt_FIMr1.1
   echo "getting weekly climo cl1=$cl1 "
   ncks -h -d  record,$jday,$jday -F ${NC3D_DIR}/smo_clim_wk1.nc $cl1
# getting anomaly 
   ncdiff -h -O $wk1 $cl1 $ano1
  
   ncks -h -d  record,$jday,$jday -F ${NC3D_DIR}/smo_clim_wk2.nc $cl2
   ncdiff -h -O $wk2 $cl2 $ano2

   ncks -h -d  record,$jday,$jday -F ${NC3D_DIR}/smo_clim_wk3.nc $cl3
   ncdiff -h -O $wk3 $cl3 $ano3

   ncks -h -d  record,$jday,$jday -F ${NC3D_DIR}/smo_clim_wk4.nc $cl4
   ncdiff -h -O $wk4 $cl4 $ano4

   ncks -h -d  record,$jday,$jday -F ${NC3D_DIR}/smo_clim_wk34.nc $cl34
   ncdiff -h -O $wk34 $cl34 $ano34

# done anomaly; clean-up
    echo " rmdir $wdir "
    /bin/rm $wdir
fi # if 6z & all 4 files exist

exit $?
