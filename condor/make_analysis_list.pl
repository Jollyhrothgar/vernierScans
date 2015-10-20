#!/usr/bin/perl
use strict;
use warnings;

my $NUMBER_OF_RUNS = 7;
my $run_index = 0;

for($run_index = 0; $run_index < $NUMBER_OF_RUNS; $run_index++) {
  print "root -l -b -q /direct/phenix+spin2/beaumim/vernierScans/processing/data_event_prdf_analysis/macros/Run_BeamWidthTime.C\'($run_index)\'\n";
  print "root -l -b -q /direct/phenix+spin2/beaumim/vernierScans/processing/hourglass_correction/macros/Run_HourglassData.C\'($run_index)\'\n";
  print "root -l -b -q /direct/phenix+spin2/beaumim/vernierScans/processing/data_event_prdf_analysis/macros/Run_WcmDcct.C\'($run_index)\'\n";
}
