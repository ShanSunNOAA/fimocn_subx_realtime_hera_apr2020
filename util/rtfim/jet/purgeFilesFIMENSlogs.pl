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
   print "DIR: $rootDirs[$i]/FIMrun\n";
   opendir (my $dh, "$rootDirs[$i]/FIMrun") || next;
   my @dirs = grep {/[^$runDirs[$i]][00]$/} readdir ($dh);
   closedir $dh;

   # purge log files
   $cmd = "/usr/bin/find $rootDirs[$i]/FIMwfm/log -mindepth 2 -maxdepth 2 -type f -mtime +7";
   print "CMD: $cmd\n";
   my @files = `$cmd`;
   for my $f (@files) {
     chop $f;
     $cmd = "/bin/rm -f $f";
     print "$cmd\n";
     `$cmd`;
   }

}
