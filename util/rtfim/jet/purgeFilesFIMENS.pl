#!/usr/bin/perl

use Time::Local;

  my @rootDirs = ("/pan2/projects/fim-njet/FIMENS_2013/");
  my @runDirs = ("fim_8");

my $numDirs = @rootDirs;

my $i, $cmd;

my $now = time;
$now = $now - ($now % (12 * 3600));
my($sec,$min,$hour,$mday,$mon,$year) = gmtime($now);

# my $killTime = $now - (10 * 86400);
my $killTime = $now - (6 * 86400);
my $fiveDays = $now - (5 * 86400);
my $fourDays = $now - (4 * 86400);
my $twoDays = $now - (2 * 86400);
# my $saveAllTime = $now - (6 * 86400);
my $saveAllTime = $now - (5 * 86400);

for ($i=0; $i<$numDirs; $i++) {
   # print "DIR: $rootDirs[$i]/FIMrun\n";
   opendir (my $dh, "$rootDirs[$i]/FIMrun") || next;
   my @dirs = grep {/[^$runDirs[$i]][00]$/} readdir ($dh);
   closedir $dh;

   # fim_8_64_208_YYYYMMDDHH00
   for my $dirName (@dirs) {
     if ($dirName =~ /$runDirs[$i]/) {
      my @tmp = split /\_/, $dirName;
      my $year = substr ($tmp[4],0,4);
      $year = $year % 1900;
      my $month = substr ($tmp[4],4,2);
      $month--;
      my $day = substr ($tmp[4],6,2);
      my $hour = substr ($tmp[4],8,2);
      my $unixTime = timegm (0,0,$hour,$day,$month,$year);
      my $linkName = substr($tmp[4],0,10);
      print "linkName: $linkName\n";
      my $cmd;
      if ($unixTime >= $saveAllTime) {
        print "$dirName NEXT\n";
        next;
      }
      elsif ($unixTime <= $killTime) {
        $cmd = "/bin/rm -rf $rootDirs[$i]/FIMrun/$dirName";
        print "cmd: $cmd\n";
   #     `$cmd`;
        $cmd = "/bin/rm f $rootDirs[$i]/FIMrun/$linkName";
        print "cmd: $cmd\n";
   #     `$cmd`;
      }
      else {
         if ($rootDirs[$i] =~ /FIMX/ && $unixTime > $fourDays) {
            $cmd = "/usr/bin/find $rootDirs[$i]/FIMrun/$dirName  -maxdepth 1 -mindepth 1 -type d | grep -v /tracker_ | grep -v /ncldiff_ | grep -v /ncl_ | grep -v /post_ | grep -v /fim_";
         }
         else {
         $cmd = "/usr/bin/find $rootDirs[$i]/FIMrun/$dirName  -maxdepth 1 -mindepth 1 -type d | grep -v /tracker_ | grep -v /ncldiff_ | grep -v /ncl_ | grep -v /post_ ";
          }
          print "CMD: $cmd\n";
          my @dirs = `$cmd`;
          for my $d (@dirs) {
            chop $d;
#            $cmd = "/bin/rm -rf $d";
            print "TO BE DELETED CMD: $cmd\n";
#             `$cmd`;
          }
#          $cmd = "/usr/bin/find $rootDirs[$i]/FIMrun/$dirName/post_  -maxdepth 1 -mindepth 1 -type d | grep -v /post_*/fim ";
#          print "CMD: $cmd\n";
#          my @dirs = `$cmd`;
#          for my $d (@dirs) {
#            chop $d;
#            $cmd = "/bin/rm -rf $d";
#            print "TO BE DELETED $cmd\n";
#            # `$cmd`;
         }
       }
      }  # good dir name
   }

   # purge log files
   $cmd = "/usr/bin/find $rootDirs[$i]/FIMwfm/log -mindepth 2 -maxdepth 2 -type f -mtime +10";
   print "CMD: $cmd\n";
   my @files = `$cmd`;
   for my $f (@files) {
     chop $f;
     $cmd = "/bin/rm -f $f";
     print "$cmd\n";
     `$cmd`;
   }

}
