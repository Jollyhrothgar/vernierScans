#!/usr/bin/perl
use strict;
use warnings;

my $NUMBER_OF_RUNS = 7;
my $run_index = 0;

for($run_index = 0; $run_index < $NUMBER_OF_RUNS; $run_index++) {
  system("root -l -b -q processing/data_event_prdf_analysis/macros/Run_BeamWidthTime.C\'($run_index)\'");
  system("root -l -b -q processing/hourglass_correction/macros/Run_HourglassData.C\'($run_index)\'");
  system("root -l -b -q processing/data_event_prdf_analysis/macros/Run_WcmDcct.C\'($run_index)\'");
}
