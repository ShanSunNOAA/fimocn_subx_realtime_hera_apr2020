#!/usr/bin/perl

##
##  purgeFiles.pl
##    this script removes the FIM*/FIMrun/fim_<GLVL>_<NLVL>_<NPROC>_YYYYMMDDHHMM 
##      directories if they are more than 10 days old; 
##    for directories 2-10 days old, delete ensics/, fim_C/ prep_C/, soundings/, 
##      and verif/ sub-directories
##    for directories 2-10 days old, delete all post_C sub-directories except fim
##      (219, 236, 244, 28, europe, africa, w_pacific, floating, e_pacific)
##    if FIMX, keep fim_C directories back 4 days
##    delete logfiles older than 3 days
##

use Time::Local;

my @rootDirs = ("/scratch1/portfolios/BMC/fim/FIM/",
                "/scratch1/portfolios/BMC/fim/FIMZEUS/",
                "/scratch1/portfolios/BMC/fim/FIM7/",
                "/scratch1/portfolios/BMC/fim/FIM9ZEUS/",
                "/scratch1/portfolios/BMC/fim/FIMX/",
                "/scratch1/portfolios/BMC/fim/FIMXZEUS/");
my @runDirs = ("fim_8","fim_8","fim_7","fim_9","fim_7","fim_8");

my $numDirs = @rootDirs;

my $i, $cmd;

my $now = time;
$now = $now - ($now % (12 * 3600));
my($sec,$min,$hour,$mday,$mon,$year) = gmtime($now);

my $killTime = $now - (8 * 86400);  # previously: $now - (10 * 86400);
my $fiveDays = $now - (5 * 86400);
my $fourDays = $now - (4 * 86400);
my $saveAllTime = $now - (2 * 86400);

for ($i=0; $i<$numDirs; $i++) {
   # print "DIR: $rootDirs[$i]/FIMrun\n";
   opendir (my $dh, "$rootDirs[$i]/FIMrun") || next;
   my @dirs = grep {/[^$runDirs[$i]][00]$/} readdir ($dh);
   closedir $dh;

   # fim_8_64_240_YYYYMMDDHHMM
   for my $dirName (@dirs) {
     if ($dirName =~ /$runDirs[$i]/ && $dirName !~ /save/) {
      my @tmp = split /\_/, $dirName;
      my $year = substr ($tmp[4],0,4);
      $year = $year % 1900;
      my $month = substr ($tmp[4],4,2);
      $month--;
      my $day = substr ($tmp[4],6,2);
      my $hour = substr ($tmp[4],8,2);
      my $unixTime = timegm (0,0,$hour,$day,$month,$year);
      my $cmd;
      my $linkName = substr($tmp[4],0,10);
      print "linkName: $linkName\n";
      if ($unixTime >= $saveAllTime) {
        print "$dirName NEXT\n";
        next;
      }
      elsif ($unixTime <= $killTime) {
        $cmd = "/bin/rm -rf $rootDirs[$i]/FIMrun/$dirName";
        print "cmd: $cmd\n";
        `$cmd`;
        $cmd = "/bin/rm $rootDirs[$i]/FIMrun/$linkName";
        print "cmd: $cmd\n";
        `$cmd`;
      }
      else {
         # if ($rootDirs[$i] =~ /FIMX/ ) {
           # print "unixTime: $unixTime fourDays: $fourDays\n";
         # }
         if ($rootDirs[$i] =~ /FIMX/ && $unixTime > $fiveDays) {
            $cmd = "/usr/bin/find $rootDirs[$i]/FIMrun/$dirName  -maxdepth 1 -mindepth 1 -type d | grep -v /tracker_C | grep -v /ncldiff_C | grep -v /ncl_C | grep -v /post_C | grep -v /fim_C";
         }
         else {
            $cmd = "/usr/bin/find $rootDirs[$i]/FIMrun/$dirName  -maxdepth 1 -mindepth 1 -type d | grep -v /tracker_C | grep -v /ncldiff_C | grep -v /ncl_C | grep -v /post_C ";
          }
          print "CMD: $cmd\n";
          my @dirs = `$cmd`;
          for my $d (@dirs) {
            chop $d;
            $cmd = "/bin/rm -rf $d";
            print "TO BE DELETED CMD: $cmd\n";
             `$cmd`;
          }
          $cmd = "/usr/bin/find $rootDirs[$i]/FIMrun/$dirName/post_C  -maxdepth 1 -mindepth 1 -type d | grep -v /post_C/fim ";
          print "CMD: $cmd\n";
          my @dirs = `$cmd`;
          for my $d (@dirs) {
            chop $d;
            $cmd = "/bin/rm -rf $d";
            print "TO BE DELETED $cmd\n";
             `$cmd`;
         }
       }
      }  # good dir name
   }

   # purge log files
   $cmd = "/usr/bin/find $rootDirs[$i]/FIMwfm/log -mindepth 2 -maxdepth 2 -type f -mtime +2 | grep -v svn";
   print "CMD: $cmd\n";
   my @files = `$cmd`;
   for my $f (@files) {
     chop $f;
     $cmd = "/bin/rm -f $f";
     print "$cmd\n";
     `$cmd`;
   }

}
