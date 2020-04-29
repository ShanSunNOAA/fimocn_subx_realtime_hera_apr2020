#!/usr/bin/perl

##
##  purgeFIMXhourlyFiles.pl
##     this script removes all non six-hourly forecast files
##       [01-05, 07-11, 13-17, ..., 157-161, 263-167]
##
##

use Time::Local;

my @rootDirs = ("/pan2/projects/fim-njet/FIMX/");
my @runDirs = ("fim_7");

my $numDirs = @rootDirs;

my $i, $cmd;

my $now = time;
$now = $now - ($now % (12 * 3600));
my($sec,$min,$hour,$mday,$mon,$year) = gmtime($now);

my $startTime = $now - (2 * 86400);
my $killTime = $now - 86400;

for ($i=0; $i<$numDirs; $i++) {
   print "DIR: $rootDirs[$i]/FIMrun\n";
   opendir (my $dh, "$rootDirs[$i]/FIMrun") || next;
   my @dirs = grep {/[^$runDirs[$i]][00]$/} readdir ($dh);
   closedir $dh;

   # fim_8_64_240_201006030000
   for my $dirName (@dirs) {
     print "dirName: $dirName\n";
     if ($dirName =~ /$runDirs[$i]/) {
       my @tmp = split /\_/, $dirName;
       my $year = substr ($tmp[4],0,4);
       $year = $year % 1900;
       my $month = substr ($tmp[4],4,2);
       $month--;
       my $day = substr ($tmp[4],6,2);
       my $hour = substr ($tmp[4],8,2);
       my $unixTime = timegm (0,0,$hour,$day,$month,$year);
       my $cmd;
       if ($unixTime >= $startTime && $unixTime <= $killTime) {
         my @files;
         print "directory: ${rootDirs[$i]}/FIMrun/${dirName}/fim_C\n";
         opendir (my $dir, "$rootDirs[$i]/FIMrun/${dirName}/fim_C") || next;
         @files = grep {/^fim_out/} readdir ($dir);
         close $dir;
         # fim_out_d4sP000031hr
         for my $file (@files) {
            my $hr = substr($file,15,3);
            # print "hr: $hr\n";
            if (($hr % 6) != 0) {
              my $cmd = "/bin/rm -f ${rootDirs[$i]}/FIMrun/${dirName}/fim_C/$file";
              print "cmd: $cmd\n";
              `$cmd`;
            }
         }
       }
     }
  }
}
1;
