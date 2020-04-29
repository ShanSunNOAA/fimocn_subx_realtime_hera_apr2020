#!/usr/bin/perl -I/. -I/w3/rapb/utilities

# gets a list of model files that are available on jet for
#
use strict;
use Time::Local;
use Time::HiRes "usleep";

my $DEBUG=0;

sub get_status_file_from_jet {
 

  my ($jet_filename, $local_filename, $DEBUG) = @_;
  my $from_jet_dir = "from_jet";
  my $wjetIsAvailable = 0;
  use POSIX "SIGALRM";
  my ($mask, $action, $oldaction, $alarm);
  # need to use POSIX alarms, because perl's regular alarm function doesn't
  # work with the DBI calls to mysql. (WRM 8 Oct 2009)
  $mask = POSIX::SigSet->new( SIGALRM ); # signals to mask in the handler
  $action = POSIX::SigAction->new( sub {die "die ALARM"; }, $mask);
  $oldaction = POSIX::SigAction->new();
  POSIX::sigaction(SIGALRM, $oldaction); # restore original
  eval {
    alarm 3;
    # for testing
    #sleep 6;
    if ($DEBUG) {
         print ("/usr/bin/scp -i /w3_server/.ssh/new_jetscp ".
                "'William.R.Moninger\@jetscp.rdhpcs.noaa.gov:$jet_filename' $local_filename\n");
    }
      # my $status = system("/usr/bin/scp -p -o ConnectTimeout=3 -i /w3_server/.ssh/new_jetscp ".
      #        "'moninger\@jetscp.rdhpcs.noaa.gov:/$jet_filename' $local_filename");
        my $status = `/usr/bin/scp -p -o ConnectTimeout=4 -i /w3_server/.ssh/new_jetscp 'William.R.Moninger\@jetscp.rdhpcs.noaa.gov:/$jet_filename' $local_filename 2>&1`;
      # my $status = system("/usr/bin/scp -p -o ConnectTimeout=10 -i /w3_server/.ssh/new_jetscp 'moninger\@jetscp.rdhpcs.noaa.gov:/$jet_filename' $local_filename");
      if ($DEBUG) {
         print "status: $status\n";
         print "wjetIsAvailable: $wjetIsAvailable\n";
      }
      if ($status !~ /closed/) {
         chmod(0777,$local_filename);
         $wjetIsAvailable = 1;
        if ($DEBUG) {
            print "COPIED\n";
         }
      }
  };  #eval
  if ($DEBUG) {
    warn "got alarm" if $@ and $@ =~ /Alarm/;
  }
  alarm(0);
    
    # process any interrupts
    if($@) {
        if ($DEBUG) {
           print "  'timed out' message";
        }
    }
    if ($DEBUG) {
         print "before return: wjetIsAvailable:  $wjetIsAvailable\n";
    }
  
  return ($wjetIsAvailable);
}

1;
