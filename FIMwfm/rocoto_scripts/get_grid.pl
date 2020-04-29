sub get_grid {
    use POSIX;
    my($file,$thisDir,$DEBUG) = @_;
    require "timelocal.pl";   #includes 'timegm' for calculations in gmt
    # get proper path for appropriate version of wgrib2
    # so do it with this terrible hack:
#    my $wgrib2;
    my $wgrib2="/apps/wgrib2/0.1.9.5.1/wgrib2";
    my($year,$julday,$hour,$fcst_proj);
   # open(S,"/home/rtfim/utilities/find_wgrib2.sh|")||
#    open(S,"find_wgrib2.sh|")||
#	die "could not do which: $!";
#    $wgrib2 = <S>;
#    close S;
#    chomp $wgrib2;
    #print "wgrib2 is $wgrib2\n";

    my($alat1,$elon1,$elonv,$alattan,$nx,$ny,$dx,$grib_type,$grid_type,
       $run_date,$fcst_proj);
    $grib_type = 0;
   
    # xue
    $grid_type = 0;

    # first, let's see if its grib(1)
    $unsafe_arg = qq{$thisDir/col_wgrib.x -V $file};
    # clean it for the taint flag
    $unsafe_arg =~ /([-\w. \/\>\|\:\']+)/;
    $arg = $1;
    open(INVENTORY, "$arg 2> /dev/null |");
    while(<INVENTORY>) {
	#print;
	if(/date (\d+) .*anl/) {
	   $run_date=$1;
	   $fcst_proj=0;
	   #print "DATE IS $run_date $fcst_proj (anx)\n";
       } elsif(/date (\d+) .* (\d+)hr fcst/) {
	   $run_date=$1;
	   $fcst_proj=$2;
	   #print "DATE IS $run_date $fcst_proj\n";
	} elsif(/Lambert Conf: Lat1 ([\-\.\d]+) Lon1 ([\-\.\d]+) Lov ([\-\.\d]+)/) {
	    #print "got lambert\n";
	    $grib_type = 1;
	    $la1 = $1;
	    $lo1 = $2;
	    $lov = $3;
	    $grid_type = 1;		# Lambert conformal conic with 1 std parallel
	    $_ = <INVENTORY>;
	    #print "Got latin1\n";
	    /Latin1 ([\.\d]+)/;
	    $latin1 = $1;
	    $_ = <INVENTORY>;
	    #print "got pole\n";
	    /North Pole \((\d+) x (\d+)\) Dx ([\.\d]+)/;
	    $nx = $1;
	    $ny = $2;
	    $dx = $3 * 1000;  # send it out in meters
	    last;
	} elsif (/latlon: lat  ([\-\.\d]+) to ([\-\.\d]+) by ([\.\d]+)/) {
	    $grib_type = 1;
	    $la1 = $1;			# lat1
	    $lo1 = $2;			# lat2
	    $dx = $3 + 0;		# make sure its a number
	    $_ = <INVENTORY>;
	    #print;
	    /long ([\-\.\d]+) to ([\-\.\d]+) by ([\.\d]+), \((\d+) x (\d+)/;
	    $lov = $1;			# lon1
	    $latin1 = $2;		# lon2
	    if($latin1 < $lov) {
		$latin1 += 360;
	    }
	    my $dy = $3 + 0;		# make sure its a number
	    $nx = $4;
	    $ny = $5;
	    if($dx == $dy) {
		if($la1 < $lo1) {
		    $grid_type = 10;	# 10 = lat-lon grid, S to N
		} else {
		    $grid_type = 11;        # 11 = lat-lon grid, N to S
		}
	    } else {
		$grid_type = 0;
	    }
	    if($DEBUG) {
		print "lat: $la1 to $lo1 by $dx\n";
		print "lon: $lov to $latin1 by $dy\n";
	    }
	    last;
	} elsif(/timerange .* nx (\d*) ny (\d*)/) {
	    $nx = $1;
	    $ny = $2;
	    if($nx == 758 && $ny == 567) {
		$grib_type = 1;
		$grid_type = 20; # define '20' as rotated lat/lon grid for RR
		last;
	    }
	}
    }
    close INVENTORY;
    print "#######################\n";
    print "Xue grid_type= $grid_type\n";
    if($grib_type != 1) {
	# its not grib(1). MustMaybe  be grib2
	$unsafe_arg = qq{$wgrib2 -grid $file};
	$unsafe_arg =~ /([\=-\w. \/]+)/;
	$arg = $1;
	if($DEBUG) {
	    print "NUMBER 2 arg is $arg\n";
	}
	open(INVENTORY,"$arg 2> /dev/null |") ||
	    print "problem with |$arg|: $!";
	$_ = <INVENTORY>;
	#print;
	$_ = <INVENTORY>;
	#print;
	if(/Lambert Conformal:/) {
	    $grib_type = 2;
	    $_ = <INVENTORY>;
	    #print;
	    /Lat1 ([\-\.\d]+) Lon1 ([\-\.\d]+) Lov ([\-\.\d]+)/i;
	    $la1 = $1;			# lat1
	    $lo1 = $2;			# lat2
	    $lov = $3;
	    #print "la1 |$la1|, lo1 $lo1\n";
	    $_ = <INVENTORY>;
	    /Latin1 ([\.\d]+)/;
	    $latin1 = $1;
	    $_ = <INVENTORY>;
	    $_ = <INVENTORY>;
	    #print;
	    /North Pole \((\d+) x (\d+)\) Dx ([\.\d]+) (.) Dy ([\.\d]+) /i;
	    $nx = $1;
	    $ny = $2;
	    $dx = $3;
	    my $mkm = $4;	# "m" or (maybe) "km"
	    $dy = $5;
	    if($mkm ne "m") {
		$dx *= 1000;  # send it out in meters
	        $dy *= 1000;
	    }
	    # a hack, because the '-grid' option USED TO NOT work on jet
	    if($dx == 0) {
		if($nx == 349 && $ny == 277) {
		    # AWIPS 221 grid
		    $dx = 32463.00;
		} else {
		    # all other grib2 grids
		    $dx = 13545.00;
		}
		$dy = $dx;
	    }
	    #print "dx |$dx|, dy |$dy|\n";
	    if($dx == $dy) {
		$grid_type = 1;	# 1 = Lamb. Conf with 1 std parallel
	    } else {
		$grid_type = 0;	# unknown grid type
	    }
        }elsif(/polar stereographic grid:/) { # start of xue
	    $grib_type = 2;
	    $grid_type = 3;
#	    $_ = <INVENTORY>;
	    #print;
#	    /lat1 ([\-\.\d]+) lon1 ([\-\.\d]+) Lov ([\-\.\d]+)/i;
	    /polar stereographic grid: \((\d+) x (\d+)\) /i;
	    $nx = $1;
	    $ny = $2;
	    print "Xue nx= $nx ny= $ny $_\n";

	    $_ = <INVENTORY>; # reading next line

#	    /North pole Lat1 ([\-\.\d]+) Lon1 ([\-\.\d]+) LatD ([\-\.\d]+) LonV ([\-\.\d]+) /i;
	    /North pole Lat1 ([\-\.\d]+) Lon1 ([\-\.\d]+) LatD ([\-\.\d]+) LonV ([\-\.\d]+) Dx ([\.\d]+) (.) Dy ([\.\d]+)/i;

#	    /North pole lat1 ([\-\.\d]+) lon1 ([\-\.\d]+) latD ([\-\.\d]+) LonV ([\-\.\d]+) dx ([\-\.\d]+) dy ([\-\.\d]+)/i;
	    $la1 = $1;			# lat1
	    $lo1 = $2;			# lat2
	    $latD = $3;			# lat2
	    $lov = $4;
            $dx = $5;
	    my $mkm = $6;	# "m" or (maybe) "km"
            $dy = $7;

	    print "Xue la1=$la1, lo1=$lo1, latD=$latD,lonv=$lov dx=$dx dy=$dy  $_\n";

	    $_ = <INVENTORY>;
	    /Latin1 ([\.\d]+)/;
	    $latin1 = $1;
	    $_ = <INVENTORY>;
	    $_ = <INVENTORY>;
	    #print;  # end of xue
	} elsif(/lat-lon grid:\((\d+) x (\d+)/) {
	    $grib_type = 2;
	    $nx = $1;
	    $ny = $2;
	    $_ = <INVENTORY>;
	    /lat (.+) to (.+) by (.+)/;
	    $la1 = $1;			# lat1
	    $lo1 = $2;			# lat2
	    $dx = $3 + 0;		# make sure its a number
	    $_ = <INVENTORY>;
	    /lon (.+) to (.+) by (.+)/;
	    $lov = $1;			# lon1
	    $latin1 = $2;		# lon2
	    my $dy = $3 + 0;		# make sure its a number
	    if($dx == $dy) {
		if($la1 < $lo1) {
		    $grid_type = 10;	# 10 = lat-lon grid, S to N
		} else {
		    $grid_type = 11;        # 11 = lat-lon grid, N to S
		}
	    } else {
		$grid_type = 0;
	    }
	} elsif(/grid_template=32769/) {
	    #print "setting grid to 20!\n";
	    $grid_type = 20;	# rotLL grid
	    $grib_type = 2;
	} else {
	    $grid_type = 0;
	}
	close INVENTORY;
	# get run date
	$unsafe_arg = qq{$wgrib2 $file};
	$unsafe_arg =~ /([\=-\w. \/]+)/;
	$arg = $1;
	print "NUMBER 3 arg is |$arg|\n";
	open(DATE,"$arg 2> /dev/null|") ||
	    print "problem with |$arg|: $!";
	my $line = <DATE>;
	($run_date) = $line =~ /d=(\d*)/;
	($fcst_proj) = $line =~ /(\d*) hour fcst/;
	$fcst_proj += 0;
	print "DATE for grib2 is $run_date, $fcst_proj\n";
	close DATE;
    }
    if($grib_type == 0) {
	#maybe a netCDF file
	my $nc_file = "NC$$.nc";
	$arg = qq{ln -s $file $nc_file; ncdump -h NC$$.nc; rm $nc_file};
	if($DEBUG) {
	    print "NUMBER 4 arg is $arg\n";
	}
	open(INVENTORY,"$arg 2> /dev/null |") ||
	    print "problem with |$arg|: $!";
	while(<INVENTORY>) {
	    if(/x = (\d+)/) {
		$nx = $1;
		$grib_type = 3;	# a netCDF file
	    }
	    if(/y = (\d+)/) {
		$ny = $1;
		last;
	    }
	}
	($year,$julday,$hour,$fcst_proj) = $file =~ m|netcdf/(..)(...)(..)(.*)|;
	$year += 2000;
	$fcst_proj += 0;	# change from string to number
	my $valid_secs = jy2secs($julday,$year)+3600*$hour;
	#print "valid_secs from file is $valid_secs\n";
	$run_date = strftime("%Y%m%d%H",gmtime($valid_secs));
	#print "run_date: $run_date\n";
	close INVENTORY;
    }
    print "grib_type2= $grib_type\n";
    # get valid_date from run_date
    ($year,$mon,$mday,$hour) = $run_date =~ /(....)(..)(..)(..)/;
    my $min=0;
    my $sec=0;
    my $valid_secs = timegm(0,0,$hour,$mday,$mon-1,$year-1900) + 3600*$fcst_proj;
    my $valid_date = sql_datetime($valid_secs);
    ($sec,$min,$hour,$maday,$mon,$year) = gmtime($valid_secs);
    if(1) {
	print "get_grid results: ".
	    "$la1,$lo1,$lov,$latin1,$nx,$ny,$dx,$grib_type,$grid_type,$valid_date,$fcst_proj\n";
    }
    return($la1,$lo1,$lov,$latin1,$nx,$ny,$dx,$grib_type,$grid_type,$valid_date,$fcst_proj);
}

sub sql_datetime {
    my $time = shift;
    my($sec,$min,$hour,$mday,$mon,$year) = gmtime($time);
    $mon++;
    $year += 1900;
    return sprintf("%4d-%2.2d-%2.2d %2.2d:%2.2d:%2.2d",
		   $year,$mon,$mday,$hour,$min,$sec);
}

sub jy2secs {
    my ($i,$julday,$leap,$timeSecs);
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst);
    ($julday,$year)=@_;
    
    #daytab holds number of days per month for regular and leap years
    my (@daytab) =(0,31,28,31,30,31,30,31,31,30,31,30,31,
	              0,31,29,31,30,31,30,31,31,30,31,30,31);
    my (@month)=(Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec);
    my (@day)=(Sun,Mon,Tue,Wed,Thu,Fri,Sat);


    #see if year was defined
    if($year == 0) {
	$timeSecs=time();
	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)
	    = gmtime($timeSecs);
	$year += 1900;
    } elsif ($year < 1000) {
	#2-digit year was (probably) input
	if($year > 70) {
	    $year += 1900;
	} else {
	    $year += 2000;
	}
    }

    #see if the year is a leap year
    $leap = ($year%4 == 0 && $year%100 != 0) || ($year%400 == 0);
    my $tt = $year%400;
    for($i=1,$mday = $julday ; $mday  > $daytab[$i + 13 * $leap]  ; $i++) {
	$mday -= $daytab[$i + 13 * $leap];
    }
    $mon=$i-1;
    my $dum;
    $timeSecs=timegm(0,0,0,$mday,$mon,$year);
    return($timeSecs);
 }
	
1;
