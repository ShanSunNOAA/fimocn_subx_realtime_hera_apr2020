#!/usr/bin/perl

##
##  purgeFimout.pl 
##
##    05 Mar 2014 - J.K. Henderson
##                  started with purgeGrib1.pl 
##
##    this script removes the fimout files if they are more than 14 hours old
##
##      $fimHome/FIM*/FIMrun/fim_<GLVL>_<NLVL>_<NPROC>_YYYYMMDDHHMM/fim_C/
###############################################################################

use Time::Local;

my $fimHome = "/scratch1/portfolios/BMC/fim/";
my @rootDirs = ("$fimHome/FIMZEUS/", "$fimHome/FIM9ZEUS/", "$fimHome/FIMXZEUS/");
my @runDirs = ("fim_8","fim_9","fim_8");

my $numDirs = @rootDirs;

my $i, $cmd;

my $fhr_incr = 6;
my $now = time;
$now = $now - ($now % ($fhr_incr * 3600));
my($sec,$min,$hour,$mday,$mon,$year) = gmtime($now);

my $saveAllTime = $now - (.6 * 86400);

for ($i=0; $i<$numDirs; $i++) {

   print "***********************\n";
   print "DIR: $rootDirs[$i]/FIMrun\n";
   opendir (my $dh, "$rootDirs[$i]/FIMrun") || next;
   my @dirs = grep {/[^$runDirs[$i]][00]$/} readdir ($dh);
   closedir $dh;

   ## check for fim_C/fim_out* and delete if greater than a day old
   $cmd = "/usr/bin/find $rootDirs[$i]/FIMrun/ -maxdepth 2 -mindepth 2 -type d | grep fim_C '";
   my @dirs = `$cmd`;
   for my $dirName (@dirs) {
     chop $dirName;
     my @tmp = split /\_/, $dirName;
     my $year = substr ($tmp[4],0,4);
     $year = $year % 1900;
     my $month = substr ($tmp[4],4,2);
     $month--;
     my $day = substr ($tmp[4],6,2);
     my $hour = substr ($tmp[4],8,2);
     my $unixTime = timegm (0,0,$hour,$day,$month,$year);

     ## keep if less than 14 hours
     if ($unixTime >= $saveAllTime) {
       print "$dirName NEXT\n";
       next;
     }

     ## delete if more than 14 hours
     else {
       if ( -d $dirName ) {
         print "$dirName found....\n";
         $cmd = "/bin/rm f $dirName/fim_out*";
         print "TO BE DELETED:   $cmd\n";
##         `$cmd`;
       }
     }  
   }  # for each directory
 }  # for each model
