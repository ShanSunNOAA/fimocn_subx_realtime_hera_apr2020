#!/usr/bin/perl

# pass in date directory and fcst
my $yyyymmddhh = $ARGV[0];
my $fcst = $ARGV[1];
die "Usage:  cp_ATCFtoTCMT_retro.pl yyyymmddhh fcst\n" unless @ARGV == 2; 

# declare variables
my $DEBUG=0;
my $BASEDIR="/pan2/projects/fim-njet/FIM9/FIMrun_sjet_p";
my $HFIPDIR="/lfs1/projects/tcmt/tier1/";
#my $HFIPDIR="/lfs2/projects/fim-njet/jhender/scratch/tier1/";

my $trackDir = "${BASEDIR}/${yyyymmddhh}/tracker_C/${fcst}/";
my $cmd;

# get names of tracker files
opendir(DIR, $trackDir);
my @filelist = grep { /\.tracker.FIM9$/ } readdir(DIR);
closedir(DIR);

# rename files to ATCF convention
#   WP07-201307090600.tracker.FIM9  ==>  awp072013_FIM9_d3162_201307090600.dat
my @newlist = ();
for my $index (0 .. $#filelist) {
  @newlist[$index] = &rename_file(@filelist[$index]);
}

if ($DEBUG) {
  print "\nelements in original and new array:  \n";
  for my $i (0 .. $#filelist) {
    print "@filelist[$i]  ==>  @newlist[$i]\n";
  }
  print "\n\n\n";
}

# copy files to TCMT project directory
print "fcst:  $fcst\n";
print "${trackDir}\n";
for my $idx (00 .. $#filelist) {
  if (-s "${trackDir}/${file}" > 1) {
    $cmd = "cp -p ${trackDir}/@filelist[$idx] ${HFIPDIR}/@newlist[$idx]";
    print "cmd: $cmd\n";
    `$cmd`;
  }
}


###==============================
### rename_file:  rename from FIM convention to ATCF convention
###
###  a<basin><cyclone number><YYYY>_<mdl>_<h,d><vvvv>_<yyyymmddhh>.dat where
###
###  <basin> = ocean basin abbreviation corresponding to location of the storm
###       al = North Atlantic basin, north of the Equator
###       sl = South Atlantic basin, south of the Equator
###       ep = North East Pacific basin, eastward of 140°W
###       cp = North Central Pacific basin, between the dateline and 140°W
###       wp = North West Pacific basin, westward of the dateline
###       io = North Indian Ocean basin, north of the Equator, between 40°E and 100°E
###       sh = South Pacific Ocean basin and South Indian Ocean basin
###  <cyclone number> = 01-49, 80-89 "training/test", 90-99 "invests" 
###  <YYYY>           = year portion of the full ATCF id for the storm - year of the hurricane
###  <mdl>            = assigned model ATCF ID   [FIM9]
###  <h,d>            = use "h" to designate historical or retrospective forecasts or "d" for demo or quasi-real-time forecasts
###  <vvvv>           = indicate the version of your model used to generate the forecast  [3162]
###  <yyyymmddhh>     = date/time of forecast initialization
###
###    example:  WP07-201307090600.tracker.FIM9  ==>  awp072013_FIM9_d3162_201307090600.dat
###

sub rename_file {
  my ($oldname) = @_;

  my $mdl = "FIM9";
  my $fcst_type = "d";
  my $version = "3162";

  my $basin = lc(substr($oldname,0,2));
  my $cyclone_number = substr($oldname,2,2);
  my $datestr = substr($oldname,5,10);
  my $year = substr($datestr,0,4);

  my $newname = "a${basin}${cyclone_number}${year}_${mdl}_${fcst_type}${version}_${datestr}.dat";
  }
###==============================
