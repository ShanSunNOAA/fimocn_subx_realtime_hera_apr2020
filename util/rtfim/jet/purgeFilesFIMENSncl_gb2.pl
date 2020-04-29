#!/usr/bin/perl

  my @rootDirs = ("/pan2/projects/fim-njet/FIMENS_2013/");
  my @runDirs = ("fim_8");

my $numDirs = @rootDirs;
my $i, $cmd;

for ($i=0; $i<$numDirs; $i++) {
   print "DIR: $rootDirs[$i]/FIMrun\n";
   opendir (my $dh, "$rootDirs[$i]/FIMrun") || next;
   my @dirs = grep {/[^$runDirs[$i]]$/} readdir ($dh);
   closedir $dh;

   # fim_8_64_208_YYYYMMDDHH00
   for my $dirName (@dirs) {
     print "processing $dirName\n";
     if ($dirName =~ /$runDirs[$i]/) {

       # delete ncl directories
       $cmd = "/usr/bin/find $rootDirs[$i]/FIMrun/$dirName -maxdepth 1 -mindepth 1 -type d | grep ncl";
       my @dirs = `$cmd`;
       for my $d (@dirs) {
         chop $d;
         $cmd = "/bin/rm -rf $d";
         print "TO BE DELETED CMD: $cmd\n";
         `$cmd`;
       }

       # delete grib2 directories
       $cmd = "/usr/bin/find $rootDirs[$i]/FIMrun/$dirName -maxdepth 4 -mindepth 4 -type d | grep grib2";
       my @dirs = `$cmd`;
       for my $d (@dirs) {
         chop $d;
         $cmd = "/bin/rm -rf $d";
         print "TO BE DELETED CMD: $cmd\n";
         `$cmd`;
       }
     }
   }
}
