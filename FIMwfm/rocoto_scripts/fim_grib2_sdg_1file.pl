#!/usr/bin/perl
use strict;
my $DEBUG=0;

my $unixStart = `/bin/date +%s`;
my $startTime = `/bin/date`;

$|=1;  #force flush of buffers after each print
open (STDERR, ">&STDOUT") || die "can't dup stdout ($!)\n";
$ENV{'PATH'}="/bin";

$SIG{ALRM} = \&my_timeout;
my %month_num = (Jan => 1, Feb => 2, Mar => 3, Apr => 4, May => 5, Jun => 6,
		 Jul => 7, Aug => 8, Sep => 9, Oct =>10, Nov =>11, Dec =>12);


$ENV{'TZ'}="GMT";
my ($aname,$aid,$alat,$alon,$aelev,@description);
my ($found_airport,$lon,$lat,$lon_lat,$time);
my ($location);
my ($startSecs,$endSecs);
my ($desired_filename,$type,$fcst_len,$elev,$name,$id,$data,$bad_data);
my ($good_data,$found_sounding_data,$maps_coords,$title,$logFile);
my ($dist,$dir,$differ);
my ($loaded_soundings);
my ($dummy,$direction);
my ($sounding_file);
my ($valid_time);

#get directory and URL
use File::Basename; 
my ($basename,$thisDir,$thisURLDir,$returnAddress,$apts,$line);
($basename,$thisDir) = fileparse($0);
($basename,$thisURLDir) = fileparse($ENV{'SCRIPT_NAME'} || '.');
#untaint
$thisDir =~ /([a-zA-Z0-9\.\/\~\_\-]*)/;
$thisDir = $1;
$basename =~ /([a-zA-Z0-9\.\/\~\_\-]*)/;
$basename = $1;
$thisURLDir =~ /([a-zA-Z0-9\.\/\~\_\-]*)/;
$thisURLDir = $1;

#change to the proper directory
use Cwd 'chdir'; #use perl version so this isn't unix-dependent
chdir ("$thisDir") ||
    die "Content-type: text/html\n\nCan't cd to $thisDir: $!\n";

use lib "./";
require "timelocal.pl";   #includes 'timegm' for calculations in gmt

#get best return address
$returnAddress = "(Unknown-Requestor)";
if($ENV{REMOTE_HOST}) {
    $returnAddress = $ENV{REMOTE_HOST};
} else {
    # get the domain name if REMOTE_HOST is not set
    my $addr2 = pack('C4',split(/\./,$ENV{REMOTE_ADDR}));
    $returnAddress = gethostbyaddr($addr2,2) || $ENV{REMOTE_ADDR};
}

use DBI;
print "DBI_DSN_SOUNDINGS: $ENV{DBI_DSN_SOUNDINGS} DBI_DSN_RUC_UA: $ENV{DBI_DSN_RUC_UA} DBI_USER: $ENV{DBI_USER} DBI_PASS: $ENV{DBI_PASS}\n";

# connect to the database
# my $dbh = DBI->connect(undef,undef,undef, {RaiseError => 1});
my $dbh = DBI->connect($ENV{DBI_DSN_SOUNDINGS},$ENV{DBI_USER}, $ENV{DBI_PASS}, {RaiseError => 1});
my $model = $ENV{MODEL};
my $sth;
my $query="replace into ${model}_raob_soundings (site,time,fcst_len,s,hydro) values (?,?,?,?,?)";
my $sth_load = $dbh->prepare($query);
my $gzipped_sounding="";
my $gzipped_hydro="";
my $sql_date;
use Compress::Zlib;


# set up for fim grids
my $file = $ENV{FILENAME};
my $fim_soundings = $ENV{FIM_SOUNDINGS};
my $scripts = $ENV{SCRIPTS_DIR};
my($year,$jday,$hour,$fcst_len) = $file =~ m|(..)(...)(..)000(...)$|;
my $valid_time = jy2secs($jday,$year)+$fcst_len*3600;
# fixed by WRM 25Sep09
$valid_time += $hour*3600;
print "YJH: $year,$jday,$hour fcst_len: $fcst_len valid_secs: $valid_time\n";

#my $stationDir = "./";
my $stationDir = "$ENV{WORK_DIR}";
my $stationFile = "fim_raobs.txt";

# create raobs file
print "before create\n";

my $region = 7;
my $cmd = "${scripts}/createRaobsFile.pl $stationDir/$stationFile $region";
print "createCMD: $cmd\n";
my $exit_code = `$cmd`;

if ($exit_code != 0) {
   return ($exit_code);
}

if(-e "$file") {
    print "FILE EXISTS TO BE PROCESSED\n";
    my $runCmd = "$fim_soundings $model $file $stationDir/$stationFile";
    print "runCmd: $runCmd\n";
    open(PULL,"$runCmd |") ||
	die "couldn't run command: $!";
    $data = "";
    while(<PULL>) {
	if(/Begin sounding data for (\w*)/) {
	    $name = $1;
	    unless($bad_data) {
		$good_data=1;
		$found_sounding_data=1;
		$loaded_soundings++;
		$title = "";
		$data="";
	    }
	} elsif (/End sounding data/) {
	    $good_data=0;
	    $differ = "grid point $dist nm / $dir deg from $name:";
	    my $fcst_len2 = sprintf("%2.2d",$fcst_len);
	    $title = $model;
	    if($fcst_len == 0) {
		$title .=" analysis valid for ";
	    } else {
		$title .=" $fcst_len2 h forecast valid for ";
	    }
	    $title .= $differ;
	    print "title: $title\ndata: $data\n";
	    my $un_gzipped = "$title\n$data";
	    $gzipped_sounding = Compress::Zlib::memGzip($un_gzipped) ||
		die "cannot do memGzip for sounding data: $!\n";
	    $sql_date = sql_datetime($valid_time);
	} elsif(/Begin hydrometeor data for (\w*)/) {
	    if($1 ne $name) {
		die "big problem for hydro data: $1 ne $name\n";
	    }
	    $data="";
	    $good_data=1;
	} elsif(/End hydrometeor data/) {
	    $good_data = 0;
	    my $hydro_length = length($data);
	    if($hydro_length == 0) {
		$gzipped_hydro = undef;
	    } else {
		$gzipped_hydro = Compress::Zlib::memGzip($data) ||
		    die "cannot do memGzip for hydro data: $!\n";
	    }
	    print "$name $fcst_len $sql_date (hydro: $hydro_length)\n";
	    print "fileName: $file model: $model name: $name sql_date: $sql_date fcst_len: $fcst_len\n";
	    # $sth_load->execute($model,$name,$sql_date,$fcst_len,
	    $sth_load->execute($name,$sql_date,$fcst_len,
			       $gzipped_sounding,$gzipped_hydro);
	} elsif (/Invalid Coordinates/) {
	    $bad_data=1;
	} elsif (/Sounding data for point/) {
	    /Sounding data for point \((.*?)\).*?\(.*? (\(.*?\))/;
	    $lon_lat=$1;
	    $maps_coords=$2;
	} elsif(/delta_east= (.*) delta_north= (.*)/) {
	    my $d_east = $1;		#
	    my $d_north = $2;
	    $dist = sqrt($d_north*$d_north + $d_east*$d_east);
	    $dist = sprintf("%.1f",$dist);
	    $dir = atan2(-$d_east,-$d_north)*57.3 + 180;
	    $dir = sprintf("%.0f",$dir);
	} elsif ($good_data) {
	    $data .= $_;
	}
	if($DEBUG) {print;}
    }

} else {
    print "$file NOT AVAILABLE\n";
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
    return $timeSecs;
}
