#!/usr/bin/perl

##
##  purgeFiles.pl
##    this script removes the $fimHome/FIM*/FIMrun/fim_<GLVL>_<NLVL>_<NPROC>_YYYYMMDDHHMM 
##      directories if they are more than 11 days old; 
##    for directories 2-11 days old, delete ensics/, fim_C/ prep_C/, soundings/, 
##      and verif/ sub-directories
##    for directories 2-11 days old, delete all post_C sub-directories except fim
##      (219, 236, 244, 28, europe, africa, w_pacific, floating, e_pacific)
##    if FIMX, keep fim_C directories back 4 days
##    delete logfiles older than 6 days in $fimHome/FIM*/FIMwfm/log/
##    delete logfiles older than 7 days in $fimHome/logs/
##

use Time::Local;

my $date = `date`;
print "$date\n";

my $fimHome = "/pan2/projects/fim-njet";
my $lfs2fimHome = "/lfs2/projects/fim-njet";
my @rootDirs = ("$fimHome/FIM/", "$fimHome/FIMX/", "$fimHome/FIM9/", "$fimHome/FIM7/",
                "$fimHome/FIMY", "$fimHome/FIMZ/", "$fimHome/FIM95/", 
                "$fimHome/FIMXVolc/FIMX", "$lfs2fimHome/GFSVerification");
my @runDirs = ("fim_8","fim_7","fim_9","fim_7","fim_8","fim_8", "fim_9", "fim_7", "fim_8");

my $numDirs = @rootDirs;

my $i, $cmd;

my $now = time;
$now = $now - ($now % (12 * 3600));
my($sec,$min,$hour,$mday,$mon,$year) = gmtime($now);

my $killTime = $now - (11 * 86400);
my $fiveDays = $now - (5 * 86400);
my $fourDays = $now - (4 * 86400);
my $twoDays = $now - (2 * 86400);
my $saveAllTime = $now - (1 * 86400);

for ($i=0; $i<$numDirs; $i++) {

   print "**************************\n";
   print "DIR: $rootDirs[$i]/FIMrun\n";
   opendir (my $dh, "$rootDirs[$i]/FIMrun") || next;
   my @dirs = grep {/[^$runDirs[$i]][00]$/} readdir ($dh);
   closedir $dh;

   # fim_8_64_240_YYYYMMDDHHMM
   for my $dirName (@dirs) {
     print "dirName: $dirName\n";
     if ($dirName =~ /^$runDirs[$i]/) {     
       # assume directories begin with fim_X_XX_XXX_YYYYMMDDHH00
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
       # delete all sub-directories if greater than 11 days
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
         # find command for domains
         $cmd = "/usr/bin/find $rootDirs[$i]/FIMrun/$dirName/post_C  -maxdepth 1 -mindepth 1 -type d | grep -v /post_C/fim ";
         print "CMD: $cmd\n";
         my @dirs = `$cmd`;
         for my $d (@dirs) {
           chop $d;
           if ( ($d =~ /130/ && $unixTime < $fiveDays) || ($d !~ /130/) ) {
             $cmd = "/bin/rm -rf $d";
             print "TO BE DELETED $cmd\n";
             `$cmd`;
           }
         }  # for
       }   # else
     }  # good dir name
   }

   # purge log files
   print "***********************\n";
   print "*****  LOG FILES  *****\n";
   $cmd = "/usr/bin/find $rootDirs[$i]/FIMwfm/log -mindepth 2 -maxdepth 2 -type f -mtime +6 | grep -v svn";
   print "CMD: $cmd\n";
   my @files = `$cmd`;
   for my $f (@files) {
     chop $f;
     $cmd = "/bin/rm -f $f";
     print "$cmd\n";
     `$cmd`;
   }

}

# purge log files in $fimHome/logs
print "***********************\n";
print "Removing from $fimHome/logs:  \n";
$cmd = "find ${fimHome}/logs -mtime +7 -name \\*.log\\* -print -exec /bin/rm -f {} \\; ";
print "$cmd\n";
`$cmd`;
