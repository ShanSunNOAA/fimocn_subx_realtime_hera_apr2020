#!/usr/bin/perl

use strict;
print "in create\n";

# PREAMBLE vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

#for security, must set the PATH explicitly
$ENV{'PATH'}="/usr/bin";
    
#get directory and URL
use File::Basename;
my ($dum,$thisDir) = fileparse($0);
$thisDir =~ m|([-~.\w/]*)|;	# untaint
$thisDir = $1;

#change to the proper directory
use Cwd 'chdir'; #use perl version so this isn't unix-dependent
chdir ("$thisDir") ||
          die "Content-type: text/html\n\nCan't cd to $thisDir: $!\n";

use DBI;
#set database connection parameters


#my $dsn = $ARGV[0];
#my $user = $ARGV[1];
#my $pass = $ARGV[2];
my $out_file = $ARGV[0];
my $region = $ARGV[1];

print " out_file: $out_file region: $region\n";
# connect to the database
my $dbh = DBI->connect($ENV{DBI_DSN_RUC_UA},$ENV{DBI_USER},$ENV{DBI_PASS}, {RaiseError => 1});
print "AFTER CONNECT\n";

if (!defined ($out_file)) {
  print ("createRRraobs.pl outFileName\n");
  exit 1;
};
print "IN CREATERAOB: DIR OUTFILE: $out_file\n";
#######################!!!!!!!!!!!!!!!!!!!!!!!! PUT BACK?????
#open(OUT,">$thisDir/$out_file") ||
open(OUT,">$out_file") ||
    die "Cannot open $out_file: $!";
print "opening $out_file\n";

my $query =<<"EOQ";
select wmoid, name, lat, lon, elev 
 from metadata where
 reg like "%$region%";
EOQ

print "query is $query\n";
my $sth = $dbh->prepare($query);
$sth->execute();
my($wmoid,$name,$lat,$lon,$elev);
$sth->bind_columns(\$wmoid,\$name,\$lat,
		   \$lon,\$elev);
while($sth->fetch()) {
    printf(OUT "%d %s %6.2f %6.2f %d\n",
	   $wmoid,$name,$lat/100.0,$lon/100.0,$elev);
}
$sth->finish();
$dbh->disconnect();
