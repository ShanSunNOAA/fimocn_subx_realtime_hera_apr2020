#!/usr/bin/perl

##  statFIM.pl
##
##     JK Henderson
##     10 Feb 2013
##
##     this script is a modification of getStatusAllRunsTest.pl. It reads 
##     from the statusAll file and updates the file to include the latest 3 
##     runs. Number of runs in file is limited to 60 (30 days). It no longer
##     reads the completed timestamp from the stdout file. Instead it uses
##     the timestamp of the most recent file created in the post/fim/NAT/grib1
##     directory. The completion date will always correspond to the last file
##     created. In the other script, the date would be blank until all runs
##     were completed. 
##
##     Number of runs can be overridden by calling the routine with a parameter.
##        FIMstat.pl num_runs_to_process
##     
##     Assumptions:  
##       reads input file from /whome/rtfim/statusAll
##       copies file to /pan2/projects/fim-njet/statusAll
##       looks for GRIB1 files in 
##         FIM*/FIMrun/fim_<#lvls>_<#procs>_<#fhrs>_YYYYMMDDHHMM/post/fim/NAT/grib1
##       must include FIMY and FIMCO2 directories for the web to display 
##         correctly
##
##       Table order on web:   FIM, FIM7, FIMX, FIMZ, FIM9
##       Table order in file:  FIM, FIMX, [FIMY], FIM7, FIMZ, FIM9, [FIMCO2]
##

use File::stat;
use POSIX qw(strftime);

my $fimDir = "/whome/rtfim/FIM/";
my $fimxDir = "/whome/rtfim/FIMX/";
my $fimyDir = "/whome/rtfim/FIMY/";
my $fim7Dir = "/whome/rtfim/FIM7/";
my $fimzDir = "/whome/rtfim/FIMZ/";
my $fim9Dir = "/whome/rtfim/FIM9/";
my $fim95Dir = "/whome/rtfim/FIM95/";
my $fimCO2Dir = "/whome/rtfim/FIMCO2/";
my $outStatus = "/pan2/projects/fim-njet/statusAll";

my @fimDirs = ("$fimDir","$fimxDir","$fimyDir","$fim7Dir","$fimzDir","$fim9Dir","$fim95Dir","$fimCO2Dir");
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
my $sixHours = 6 * 3600;
my $endTime = $now - ($now % $sixHours);
my $startTime = $endTime - ($numruns * $sixHours);
my ($dir, $grepCmd, @lines, @values, $lastFileAvail, $lastFile, $completeStr, $outStr);
my ($runString,$fileDate,$fileTime,$fileName,$avail);
my ($dum0,$dum1,$dum2,$dum3);
my $runTimeFile = "/whome/rtfim/statusAll";
my %data;

#my $fimRunDir = "/FIMrun/fim_8_64_240_";
#my $fimxRunDir = "/FIMrun/fim_7_64_240_";
#my $fimyRunDir = "/FIMrun/fim_8_64_600_";
#my $fimzRunDir = "/FIMrun/fim_8_64_240_";
#my $fim7RunDir = "/FIMrun/fim_7_64_120_";
#my $fim9RunDir = "/FIMrun/fim_9_64_1600_";
#my $fimCO2RunDir = "/FIMrun/fim_7_64_240_";

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
    ($key,$fimRunTime,$fimRunComplete,$fimLastFcst,$fimxRunComplete,$fimxLastFcst,
     $fimyRunComplete,$fimyLastFcst,$fim7RunComplete,$fim7LastFcst,$fimzRunComplete,$fimzLastFcst,
     $fim9RunComplete,$fim9LastFcst,$fim95RunComplete,$fim95LastFcst,$fimCO2RunComplete,
     $fimCO2LastFcst) = split /\|/, $r;
    $data{$key}{fimRunTime} = $fimRunTime;
    $data{$key}{runComplete0} = $fimRunComplete;
    $data{$key}{lastFcst0} = $fimLastFcst;
    $data{$key}{runComplete1} = $fimxRunComplete;
    $data{$key}{lastFcst1} = $fimxLastFcst;
    $data{$key}{runComplete2} = $fimyRunComplete;
    $data{$key}{lastFcst2} = $fimyLastFcst;
    $data{$key}{runComplete3} = $fim7RunComplete;
    $data{$key}{lastFcst3} = $fim7LastFcst;
    $data{$key}{runComplete4} = $fimzRunComplete;
    $data{$key}{lastFcst4} = $fimzLastFcst;
    $data{$key}{runComplete5} = $fim9RunComplete;
    $data{$key}{lastFcst5} = $fim9LastFcst;
    $data{$key}{runComplete6} = $fim95RunComplete;
    $data{$key}{lastFcst6} = $fim95LastFcst;
    $data{$key}{runComplete7} = $fimCO2RunComplete;
    $data{$key}{lastFcst7} = $fimCO2LastFcst;
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
    my $outDateDir = sprintf("%04d%02d%02d%02d",$yr,$m,$mday,$hour);
    my $gribFileName = sprintf("%02d%03d%02d",$year,$yday,$hour);
    my $runStr = sprintf("%02d/%02d/%04d %02dZ",$m,$mday,$yr,$hour);
    my $jetFIMGribDir = "${fimDir}/FIMrun/${outDateDir}/post_C/fim/NAT/grib1/";
    my $jetFIMXGribDir = "${fimxDir}/FIMrun/${outDateDir}/post_C/fim/NAT/grib1/";
    my $jetFIMYGribDir = "${fimyDir}/FIMrun/${outDateDir}/post_C/fim/NAT/grib1/";
    my $jetFIM7GribDir = "${fim7Dir}/FIMrun/${outDateDir}/post_C/fim/NAT/grib1/";
    my $jetFIMZGribDir = "${fimzDir}/FIMrun/${outDateDir}/post_C/fim/NAT/grib1/";
    my $jetFIM9GribDir = "${fim9Dir}/FIMrun/${outDateDir}/post_C/fim/NAT/grib1/";
    my $jetFIM95GribDir = "${fim95Dir}/FIMrun/${outDateDir}/post_C/fim/NAT/grib1/";
    my $jetFIMCO2GribDir = "${fimCO2Dir}/FIMrun/${outDateDir}/post_C/fim/NAT/grib1/";

    my @fimGribDirs = ("$jetFIMGribDir","$jetFIMXGribDir","$jetFIMYGribDir","$jetFIM7GribDir","$jetFIMZGribDir","$jetFIM9GribDir","$jetFIM95GribDir","$jetFIMCO2GribDir");

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
          $completeStr = strftime("%m/%d/%Y %H:%MZ", localtime(stat($lastFileAvail)->mtime));
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
