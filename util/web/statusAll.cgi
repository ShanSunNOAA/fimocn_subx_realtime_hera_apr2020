#!/usr/bin/perl -I. -I/w3/rapb/utilities

use CGI qw/:standard :html3 :form :cgi :html/;
use LWP::UserAgent;
use HTTP::Request;
use HTTP::Response;
use URI::Heuristic;
use File::stat;
use Time::localtime;

require "get_status_file_from_jet.pl";

my $statusName = "statusAll";
#my $jetDir = "/lfs0/projects/rtfim/";
my $jetDir = "/pan2/projects/fim-njet/";
my $webDir = "/w3/rapb/fim/from_jet/";


my $DEBUG = 0;
get_status_file_from_jet ("${jetDir}${statusName}","${webDir}${statusName}", $DEBUG);
my $mtime = stat("${webDir}${statusName}")->mtime;
my ($sec2,$min2,$hour2,$mday2,$mon2,$year2,$wday2,$yday2,$isdst2) = gmtime($mtime);
$year2 = ($year2 + 1900);
$mon2++;
my $lastChecked = sprintf("%02d/%02d/%04d %02d:%02dZ",$mon2,$mday2,$year2,$hour2,$min2);

print <<EOF;
Content-Type: text/html

<!--#set var="page_title" value="FIM Run Status" -->
<!--#set var="refresh" value="refresh" -->
<!--#set var="fieldTableWidth" value="85em" -->
<!--#set var="bodyWidth" value="96em" -->
<!--#include virtual="headers/FIM_header.html" -->

<h2>Last Checked: $lastChecked</h2>

EOF

my %data;
open (IN, "${webDir}${statusName}") || print "could not access ${webDir}${statusName}\n";
my @runs = <IN>;
close IN;

print "<BR>\n";

print "<table class=\"fieldTable\">\n";
print "<tr class=\"header\">";
print "<th class=\"FIMstatusHdr\">Run Time</th>";
print "<th class=\"FIMcompleteHdr\">FIM Run<br>Complete</th>";
print "<th class=\"FIMstatusHdr\">FIM<br>Last<br>Fcst</th>";
print "<th class=\"FIMcompleteHdr\">FIM7 Run<br>Complete</th>";
print "<th class=\"FIMstatusHdr\">FIM7<br>Last<br>Fcst</th>";
print "<th class=\"FIMcompleteHdr\">FIMX Run<br>Complete</th>";
print "<th class=\"FIMstatusHdr\">FIMX<br>Last<br>Fcst</th>";
print "<th class=\"FIMcompleteHdr\">FIMY Run<br>Complete</th>";
print "<th class=\"FIMstatusHdr\">FIMY<br>Last<br>Fcst</th>";
# print "<th class=\"FIMcompleteHdr\">FIMZ Run<br>Complete</th>";
# print "<th class=\"FIMstatusHdr\">FIMZ<br>Last<br>Fcst</th>";
print "<th class=\"FIMcompleteHdr\">FIM9 Run<br>Complete</th>";
print "<th class=\"FIMstatusHdr\">FIM9<br>Last<br>Fcst</th>";
# print "<th class=\"FIMcompleteHdr\">FIM95 Run<br>Complete</th>";
# print "<th class=\"FIMstatusHdr\">FIM95<br>Last<br>Fcst</th>";
print "</tr>";

my $count = 0;
for my $r (@runs) {
    # my ($fimRunString,$fimEndTime,$fimLastFcst,$fimxEndTime,$fimxLastFcst,$fimyEndTime,$fimyLastFcst,
    my ($runTime, $fimRunString,$fimEndTime,$fimLastFcst,$fimxEndTime,$fimxLastFcst,$fimyEndTime,
        $fimyLastFcst,$fim7EndTime,$fim7LastFcst,$fimzEndTime,$fimzLastFcst,$fim9EndTime,$fim9LastFcst,
        $fim95EndTime,$fim95LastFcst,$fimCO2EndTime,$fimCO2LastFcst) = split /\|/, $r;
#   print " ($runTime, $fimRunString,$fimEndTime,$fimLastFcst,$fimxEndTime,$fimxLastFcst,$fimyEndTime,$fimyLastFcst,$fim7EndTime,$fim7LastFcst,$fimzEndTime,$fimzLastFcst,$fim9EndTime,$fim9LastFcst<br>\n";
    if ($count % 2 != 0) {
               print "<tr class=\"white\">\n";
            }
            else {
               print "<tr class=\"color\">\n";
            }
   $count++;
   print "<td class=\"fieldNameFIMstatus\">$fimRunString</td>\n";
   print "<td class=\"fieldNameFIMstatus\">$fimEndTime</td>\n";
   print "<td class=\"fieldNameFIMstatus\">$fimLastFcst</td>\n";
   print "<td class=\"fieldNameFIMstatus\">$fim7EndTime</td>\n";
   print "<td class=\"fieldNameFIMstatus\">$fim7LastFcst</td>\n";
   print "<td class=\"fieldNameFIMstatus\">$fimxEndTime</td>\n";
   print "<td class=\"fieldNameFIMstatus\">$fimxLastFcst</td>\n";
   print "<td class=\"fieldNameFIMstatus\">$fimyEndTime</td>\n";
   print "<td class=\"fieldNameFIMstatus\">$fimyLastFcst</td>\n";
#  print "<td class=\"fieldNameFIMstatus\">$fimzEndTime</td>\n";
#  print "<td class=\"fieldNameFIMstatus\">$fimzLastFcst</td>\n";
   print "<td class=\"fieldNameFIMstatus\">$fim9EndTime</td>\n";
   print "<td class=\"fieldNameFIMstatus\">$fim9LastFcst</td>\n";
#  print "<td class=\"fieldNameFIMstatus\">$fim95EndTime</td>\n";
#  print "<td class=\"fieldNameFIMstatus\">$fim95LastFcst</td>\n";
#  print "<td class=\"fieldNameFIMstatus\">$fimCO2EndTime</td>\n";
#  print "<td class=\"fieldNameFIMstatus\">$fimCO2LastFcst</td>\n";
   print "</tr>\n";
}
print "</TABLE>\n";
print <<EOF;
<div class="clear"> </div>
<!--#include virtual="/footer.html" -->
<BR>
EOF
