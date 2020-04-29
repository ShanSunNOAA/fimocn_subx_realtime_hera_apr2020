#!/usr/bin/perl -I.

use DBI;
use Time::Local;

# Compute 20-wave anomaly correlation of 500 mb heights for
# GFS and FIM forecasts for 0, 24, 48, 72, 96 and 120 hour forecasts
# over both extratropics (20N - 80N and 20S-80S)
# Taken from script in diffgb.doc

my $diffgb = $ENV{DIFFGB};
my $climateFile = $ENV{CLIMATE_FILE};
my $model = $ENV{MODEL};
my $validTime = $ENV{VALID_TIME};
my $validHour = $ENV{VALID_HOUR};

# my $variables = $ENV{VARIABLES};
# my $kStrings = $ENV{K_STRINGS};
# my $levels = $ENV{LEVELS};
#
my $levels = "500:850,250:850,250";
my $variables = "HGT:UGRD:VGRD";
my $kStrings = "-t 'ano' -k '4*-1 7 100 :-t 'ano' -k '4*-1 33 100 :-t 'ano' -k '4*-1 34 100 :";

my $analysisFileName = "$ENV{ANX_DIR}/$ENV{ANX_FILE_NAME}";
my $anxModel = "$ENV{ANX_MODEL}";
my $modelFileNames = "$ENV{MODEL_FILE_NAMES}";
my @fileNames = split /:/,$modelFileNames;

print "in verif.pl modelFileNames: $modelFileNames\n";

my $dsn= $ENV{DBI_DSN};
my $user= $ENV{DBI_USER};
my $pass = $ENV{DBI_PASS};

print "dsn: $dsn\n";
print "user: $user\n";
print "pass: $pass\n";
print "analysisFileName: $analysisFileName\n";
print "anxModel: $anxModel\n";


# 255 - user defined grid
# 0 - lat/lon grid
# 144 - numLons
# 25 - numLats
# 20000 - start lat
# 0, 0 ?????
# 80000 - end lat
# -2500 - 2.5 degree
 
# numLons = 720 / 2.5 = 144
# numLats = 80 -20 = 60 / 2.5 = 24 + 1 = 25
my $ngrid= "255,0,144,25,+20000,0,0,+80000,-2500,0,0,64,0,0,0,0,0,0,0,0,255";
my $sgrid= "255,0,144,25,-80000,0,0,-20000,-2500,0,0,64,0,0,0,0,0,0,0,0,255";
# numLats = 20 - -20 = 40 / 2.5 = 16 + 1 = 117
my $tgrid= "255,0,144,17,-20000,0,0,+20000,-2500,0,0,64,0,0,0,0,0,0,0,0,255";
# full grid - .5 degree
my $ggrid= "255,0,720,361,-90000,0,0,+90000,-500,0,0,64,0,0,0,0,0,0,0,0,255";
# numLats = 90 - 70 = 20 / 2.5 = 8 + 1 = 9
# numLons = 720 / 2.5 = 144
my $npgrid = "255,0,144,9,+70000,0,0,+90000,-2500,0,0,64,0,0,0,0,0,0,0,0,255";
my $spgrid = "255,0,144,9,-70000,0,0,-90000,-2500,0,0,64,0,0,0,0,0,0,0,0,255";

my @grids = ($ngrid, $sgrid, $tgrid, $ggrid, $npgrid, $spgrid);
my @region = (10, 9, 8, 7, 11, 12);

my $gridsCt = @grids;

my ($cmd, $wacorr);

if (-s "$analysisFileName" <= 0) {
      print "$analysisFileName does not exist - program aborted\n";
      exit(1);
}

for $f (@fileNames) {
  print "fimFileName: $f\n";
}

# -h no header
# -t 'ano' - type is anomaly correlation
# 4 *-1 - for the 1st 4 kpds paramaters, use defaults
# -k '4*-1 7 100 500' - 7,100 (HGT) 500 (mb) 
# -g grid
# -C climate grid
# -c climate grid index


my @output = ();

my $numFiles = @fileNames;
print "numFiles: $numFiles\n";
for (my $j=0; $j<$numFiles; $j++) {
   print "fileName: $fileNames[$j]\n";
}

# connect to database and set up query
my $dbh = DBI->connect($dsn,$user,$pass) or die "Unable to connect $DBI::errstr\n";
print "after connected\n";
my ($sth,$query);
$query =<<"EOI";
replace into stats_retro
 (valid_date,valid_hour,model,anx_model,level,variable,fcst_len,region,wacorr)
 values(?,?,?,?,?,?,?,?,?)
EOI
;
my $sth_load = $dbh->prepare($query);

my @vars = split /:/, $variables;
my @kStrs = split /:/, $kStrings;
my @levs = split /:/, $levels;

my $numVars = @vars;

for (my $j=0; $j<$numFiles; $j++) {
  if (-e "$fileNames[$j]") { 
     print "j: $j fileName: $fileNames[$j]\n";
     # get rid of .grib1 on analysis files
     my $tmpName = $fileNames[$j];
     print "tmpName: $tmpName\n";
     $tmpName =~ s/.grib1//;
     print "tmpName: $tmpName\n";
     my $len = length ($tmpName);
     my $fc = substr($tmpName,$len-3,3);
     print "fc: $fc\n";
     $fcst = sprintf ("%d", $fc);
     print "FOUND fileName: $j: $fileNames[$j] fcst: $fcst\n";
     for (my $k=0; $k<$numVars; $k++) {
        my $variable = $vars[$k];
        my $kString = $kStrs[$k];
        my @levels = split /,/, $levs[$k];
        my $numLevels = @levels;
        for (my $l=0; $l<$numLevels; $l++) {
           $level = $levels[$l];
            for (my $i=0; $i<$gridsCt; $i++) {
               $cmd = "$diffgb -h $kString $level ' -g\"$grids[$i]\"  -C$climateFile  -z '1 20' -x $analysisFileName $fileNames[$j]";
               print "cmd: $cmd\n\n";
               my $out = `$cmd`;
               chop ($out);
               $out =~ s/^\s+//;
               print "out: $out\n";
               my @values = split /\s+/, $out;
               my $ct = @values;
               if ($ct > 11) {
                 $wacorr = $values[12] * 100;
                 $sth_load->execute("$validTime",$validHour,"$model","$anxModel",$level,"$variable",
                                   $fcst,"$region[$i]",$wacorr);
                 my $str = "$validTime:$validHour:$model:$anxModel:$level:$variable:$fcst:$region[$i]:$wacorr";
                 print "str: $str\n";

               }
            }
        }
     }
  }
  else {
   print "NOT FOUND\n";
  }
}

$dbh->disconnect();

exit $?;

1;
