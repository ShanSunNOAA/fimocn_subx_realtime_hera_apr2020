#!/bin/csh

# csh $FIM_HOME/FIMwfm/ensplots/plot_ellipse.csh $date $tcvitalsPath $glvl_in $fcst_len_in $rundir_in $nmembers $fcst_output_int $fimensPath $FIM_HOME

setenv date $1
echo "*** in plot_ellipse.csh plotPath: $date"
setenv tcvitalsPath $2
echo "*** in plot_ellipse.csh tcvitalsPath: $tcvitalsPath"
setenv glvl_in $3
echo "*** in plot_ellipse.csh glvl_in: $glvl_in"
setenv fcst_len_in $4
echo "*** in plot_ellipse.csh fcst_len_in: $fcst_len_in"
setenv nmembers $5
echo "*** in plot_ellipse.csh nmembers: $nmembers"
setenv fcst_output_int_in $6
echo "*** in plot_ellipse.csh fcst_output_int_in: $fcst_output_int_in"
setenv fimensPath $7
echo "*** in plot_ellipse.csh fimensPath: $fimensPath"
setenv plotPath $8
echo "*** in plot_ellipse.csh plotPath: $plotPath"
setenv fim_home $9
echo "*** in plot_ellipse.csh fim_home : $fim_home"

set hr=`echo $date | cut -c9-10`
set tcvitfile=$tcvitalsPath/tcvitals.$date.txt
set storms=`cat $tcvitfile | awk '{print $2}'`
/bin/rm -rf $plotPath/ensemble-storm-data.txt
touch $plotPath/ensemble-storm-data.txt
foreach storm ($storms)
 echo storm: $storm
 set stormbasin=`echo $storm | cut -c3-3`
 if ($stormbasin == 'E' || $stormbasin == 'W' || $stormbasin == 'L') then
 if ($stormbasin == 'E') set stormbasin='EP'
 if ($stormbasin == 'W') set stormbasin='WP'
 if ($stormbasin == 'L') set stormbasin='AL'
 set stormnum=`echo $storm | cut -c1-2`
 echo "$date $stormbasin $stormnum $storm" >> $plotPath/ensemble-storm-data.txt
 endif
end
echo "in plot_ellipse.csh PWD: $PWD"
sh $fim_home/FIMwfm/ensplots/idl_ellipse.deck
