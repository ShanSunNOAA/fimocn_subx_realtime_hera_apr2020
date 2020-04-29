#!/usr/bin/perl

##  statFIM_zeus.pl
##
##     JK Henderson
##     13 Feb 2013
##
##     this script is a modification of statFIM.pl which runs on jet. It reads 
##     from the statusAll file and updates the file to include the latest 3 
##     runs. Number of runs in file is limited to 60 (30 days). It uses
##     the timestamp of the most recent file created in the post/fim/NAT/grib1
##     directory as the completion date. The completion date will always correspond
##     to the last file created. 
##
##     Number of runs can be overridden by calling the routine with a parameter.
##        FIMstat_zeus.pl num_runs_to_process
##     
##     Assumptions:  
##       reads input file from $RTFIMHOME/statusAll
##       copies file to $RTHOME/statusAll
##       looks for GRIB1 files in 
##         FIM*/FIMrun/fim_<#lvls>_<#procs>_<#fhrs>_YYYYMMDDHHMM/post/fim/NAT/grib1
##       ???must include FIMY and FIMCO2 directories for the web to display 
##         correctly
##
##       ???Table order on web:   FIM, FIM7, FIMX, FIMZ, FIM9
##       ???Table order in file:  FIM, FIMX, [FIMY], FIM7, FIMZ, FIM9, [FIMCO2]
##
##       ???Table order in file (8/6/2013):  FIM, FIMZEUS, FIM7, FIM9ZEUS, FIMX, FIMXZEUS
##

use File::stat;
use POSIX qw(strftime);

my $rtfimhome = "/home/rtfim/";
my $rthome1 = "/scratch1/portfolios/BMC/fim/";
my $rthome2 = "/scratch2/portfolios/BMC/fim/";
my $fimDir = "$rthome1/FIM/";
my $fimzeusDir = "$rthome1/FIMZEUS/";
my $fim7Dir = "$rthome1/FIM7/";
my $fim9zeusDir = "$rthome1/FIM9ZEUS/";
my $fimxDir = "$rthome1/FIMX/";
my $fimxzeusDir = "$rthome1/FIMXZEUS/";
my $outStatus = "$rthome1/statusAll";

my @fimDirs = ("$fimDir","$fimzeusDir","$fim7Dir","$fim9zeusDir","$fimxDir","$fimxzeusDir");
my $numRuns = @fimDirs;

print "numRuns: $numRuns\n";

my $now = time;

my %months = ("Jan" => 1, "Feb" => 2, "Mar" => 3,
              "Apr" => 4,"May" => 5,"Jun" => 6,
              "Jul" => 7,"Aug" => 8,"Sep" => 9,
              "Oct" => 10,"Nov" => 11,"Dec" => 12);

if (@ARGV == 1) {
  $numruns = $ARGV[0];
}
  else {
  $numruns = 2;
}

print "number of runs to process:  $numruns\n";

my $hour = 3600;
#my $twelveHours = 12 * 3600;
my $sixHours = 6  * 3600;
my $endTime = $now - ($now % $sixHours);
my $startTime = $endTime - ($numruns * $sixHours);
my ($dir, $grepCmd, @lines, @values, $lastFileAvail, $lastFile, $completeStr, $outStr);
my ($runString,$fileDate,$fileTime,$fileName,$avail);
my ($dum0,$dum1,$dum2,$dum3);
my $runTimeFile = "$rtfimhome/statusAll";
my %data;

my $fimRunDir = "/FIMrun/fim_8_64_240_";
my $fimzeusRunDir = "/FIMrun/fim_8_64_240_";
my $fim7RunDir = "/FIMrun/fim_7_64_120_";
my $fim9zeusRunDir = "/FIMrun/fim_9_64_1200_";
my $fimxRunDir = "/FIMrun/fim_7_64_240_";
my $fimxzeusRunDir = "/FIMrun/fim_8_64_240_";

# read in current statusAll file
open (IN, $runTimeFile) || print "could not access $runTimeFile\n";
my @runs = <IN>;
close IN;

# limit to 60 runs or 30 days
my $keepruns = 60;
my $numlines = @runs;
print "$numlines lines\n";
if ($numlines > $keepruns) {
  @runs = @runs[0..$keepruns]
}

my $key;

for my $r (@runs) {
    chop $r;
    ($key,$fimRunTime,$fimRunComplete,$fimLastFcst,$fimzeusRunComplete,$fimzeusLastFcst,
     $fim7RunComplete,$fim7LastFcst,$fim9zeusRunComplete,$fim9zeusLastFcst,
     $fimxRunComplete,$fimxLastFcst,$fimxzeusRunComplete,$fimxzeusLastFcst) = split /\|/, $r;
    $data{$key}{fimRunTime} = $fimRunTime;
    $data{$key}{runComplete0} = $fimRunComplete;
    $data{$key}{lastFcst0} = $fimLastFcst;
    $data{$key}{fimzeusRunTime} = $fimzeusRunTime;
    $data{$key}{runComplete1} = $fimzeusRunComplete;
    $data{$key}{lastFcst1} = $fimzeusLastFcst;
    $data{$key}{runComplete2} = $fim7RunComplete;
    $data{$key}{lastFcst2} = $fim7LastFcst;
    $data{$key}{runComplete3} = $fim9zeusRunComplete;
    $data{$key}{lastFcst3} = $fim9zeusLastFcst;
    $data{$key}{runComplete4} = $fimxRunComplete;
    $data{$key}{lastFcst4} = $fimxLastFcst;
    $data{$key}{runComplete5} = $fimxzeusRunComplete;
    $data{$key}{lastFcst5} = $fimxzeusLastFcst;
}

print 'START:  ',  strftime("%Y-%m-%d %H:%M", localtime($startTime)), "\n";
print 'END:    ',  strftime("%Y-%m-%d %H:%M", localtime($endTime)), "\n";

for (my $time=$startTime; $time<=$endTime; $time+=$sixHours) {

    my ($sec,$min,$hour,$mday,$m,$yr,$wday,$yday,$isdst) = gmtime($time);
    $yr = ($yr + 1900);
    $year = ($yr % 2000);
    $yday++;
    $m++;
    $key = sprintf("%04d%02d%02d%02d",$yr,$m,$mday,$hour);
    my $outDateDir = sprintf("%04d%02d%02d%02d00",$yr,$m,$mday,$hour);
    my $gribFileName = sprintf("%02d%03d%02d",$year,$yday,$hour);
    my $runStr = sprintf("%02d/%02d/%04d %02dZ",$m,$mday,$yr,$hour);
    my $zeusFIMGribDir = "${fimDir}${fimRunDir}${outDateDir}/post_C/fim/NAT/grib1/";
    my $zeusFIMzeusGribDir = "${fimzeusDir}${fimzeusRunDir}${outDateDir}/post_C/fim/NAT/grib1/";
    my $zeusFIM7GribDir = "${fim7Dir}${fim7RunDir}${outDateDir}/post_C/fim/NAT/grib1/";
    my $zeusFIM9zeusGribDir = "${fim9zeusDir}${fim9zeusRunDir}${outDateDir}/post_C/fim/NAT/grib1/";
    my $zeusFIMXGribDir = "${fimxDir}${fimxRunDir}${outDateDir}/post_C/fim/NAT/grib1/";
    my $zeusFIMXzeusGribDir = "${fimxzeusDir}${fimxzeusRunDir}${outDateDir}/post_C/fim/NAT/grib1/";

    my @fimGribDirs = ("$zeusFIMGribDir","$zeusFIMzeusGribDir","$zeusFIM7GribDir","$zeusFIM9zeusGribDir","$zeusFIMXGribDir","$zeusFIMXzeusGribDir");
    my @fimRunDirs = ("$fimRunDir","$fimzeusRunDir","$fim7RunDir","$fim9zeusRunDir","$fimxRunDir","$fimxzeusRunDir");

    $outStr="";
    $data{$key}{fimRunTime} = $runStr;
    for (my $i=0; $i<$numRuns; $i++) {
        my $outDir = sprintf("%04d%02d%02d%02d",$yr,$m,$mday,$hour);
        print "==> outDir:  $outDir\n" if ($i == 0);

        # find latest grib file avail 
        my @files = <${fimGribDirs[$i]}/${gribFileName}*>;
        my $numfiles = @files;
        print " $fimGribDirs[$i]  $numfiles files \n";
        if ($numfiles > 0) {
          $lastFileAvail = @files[-1];
          print "    last file available: ",  (split '/', $lastFileAvail)[-1], "\n";
          $completeStr = strftime("%m/%d/%Y %H:%M Z", localtime(stat($lastFileAvail)->mtime));
          $lastFile = substr($lastFileAvail,-3);                       ## get forecast hour
          $data{$key}{"lastFcst${i}"} = $lastFile;
          $data{$key}{"runComplete$i"} = $completeStr;
          print "    completeStr: $completeStr\n";
        }
    } # each run
} # each time

# don't print YYYYMMDDHH to output file used by web site 
open (OUT, ">$runTimeFile") || print "could not access $runTimeFile\n";
for my $k1 ( reverse sort keys %data ) {
    my $printValues = 0;
    for (my $k2=0; $k2<$numRuns; $k2++) {
      my $kR = "runComplete${k2}";
      my $kF = "lastFcst${k2}";
      if (defined($data{ $k1 }{$kR}) || defined($data{ $k1 }{$kF})) {
        $printValues = 1; 
        last;
      }
    }
    if ($printValues) {
      print OUT "$k1|";
      print OUT "$data{$k1}{fimRunTime}|";
      for (my $k2=0; $k2<$numRuns; $k2++) {
        my $kR = "runComplete${k2}";
        my $kF = "lastFcst${k2}";
        print OUT "$data{ $k1 }{$kR}|$data{ $k1 }{$kF}|";
      }
      print OUT "\n";
   }
}
close OUT;
my $cmd = "/bin/cp -p $runTimeFile $outStatus";
`$cmd`;
