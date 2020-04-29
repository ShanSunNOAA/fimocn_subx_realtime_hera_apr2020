#!/usr/bin/perl                                                                                                                                                      
use strict;

# PREAMBLE vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv                                                                                                                       

my $DEBUG=0;  # MAKE THIS NON-ZERO TO PRODUCE EXTRA DEBUG PRINTOUT                                                                                                   

#get directory                                                                                                                                                       
use File::Basename;
my ($basename,$thisDir) = fileparse($0);
$basename =~ m|([\-\~\.\w]*)|;  # untaint                                                                                                                            
$basename = $1;
$thisDir =~ m|([\-\~\.\w\/]*)|; # untaint                                                                                                                            
$thisDir = $1;

#change to the proper directory                                                                                                                                      
use Cwd 'chdir'; #use perl version so this isn't unix-dependent                                                                                                      
chdir ("$thisDir") ||
    die "Content-type: text/html\n\nCan't cd to $thisDir: $!\n";
$thisDir = $ENV{PWD};

#useful DEBUGGING info vvvvvvvvvvvvvv                                                                                                                                
if($DEBUG) {
    foreach my $key (sort keys(%ENV)) {
        print "$key: $ENV{$key}\n";
    }
    print "thisDir is $thisDir\n";
    print "basename is $basename\n";
    print "\n";
}
# end useful DEBUGGING info ^^^^^^^^^^^^^^^^^                                                                                                                        
# END PREAMBLE^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                                                                                  


my $file = $ARGV[0];
my $out_file = $ARGV[1];
my $grib_type = $ARGV[2];
my $rh_flag = $ARGV[3];
# xue add this 2015-01-14
$rh_flag = 1;

print "weixue rh_flag = $rh_flag\n";
my %inventory;

my @var_list = qw(PRES DPT UGRD VGRD TMP RH);
if($rh_flag == 1) {
    @var_list = qw(PRES DPT UGRD VGRD TMP SPFH);
} elsif($rh_flag == 2) {
    # dummy up DPT as TMP because GFS doesn't provide DPT                                                                                                            
    @var_list = qw(PRES TMP UGRD VGRD TMP SPFH);
}
#open(I,"wgrib2 $file|") ||
open(I,"/apps/wgrib2/0.1.9.5.1/wgrib2 $file|") ||
    die "could not run wgrib2 on $file: $!";
while(<I>) {
    #print;                                                                                                                                                          
    chomp;
#    if(/(PRES):surface:/) {                                                                                                                                         
    if(/(PRES):(surface:|1 hybrid level:)/) {
        $inventory{$1} = $_;
    } else {
#       if(/((DPT|TMP|RH|SPFH)):2 m above ground:/) {                                                                                                                
        if(/((TMP|RH|SPFH)):2 m above ground:/) {
            $inventory{$1} = $_;
        } elsif(/((UGRD|VGRD)):10 m above ground:/) {
            $inventory{$1} = $_;
        } elsif(/(DPT):(2 m above ground:|1 hybrid level:)/) {
            $inventory{$1} = $_;
        }
    }
}
close I;

open(D,"|/apps/wgrib2/0.1.9.5.1/wgrib2 -i -order raw -no_header -bin $out_file $file >/dev/null 2>&1") ||
    die "could not make wgrib2 dump file: $!";
foreach my $var (@var_list) {
    print "$inventory{$var}\n";
    print D "$inventory{$var}\n";
}
close D;
