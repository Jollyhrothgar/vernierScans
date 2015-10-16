#!/usr/bin/perl -w
use strict;
use warnings;

my $dst_list_name    = $ARGV[0]; # DST list name, must contain run number
my $output_directory = $ARGV[1]; # must be of form /dir1/dir2 (no trailing "/" character)

# Hardcode these to avoid a large number of arguments fed to make_job.pl...
my $rootmacro = "/direct/phenix+spin2/beaumim/vernierScans/data_reduction/macros/Run_Reduction.C";
my $exedir    = "/direct/phenix+spin2/beaumim/vernierScans/data_reduction/macros/condor";
my $jobdir    = "/direct/phenix+spin2/beaumim/vernierScans/data_reduction/macros/condor";
my $logdir    = "/direct/phenix+spin2/beaumim/vernierScans/data_reduction/macros/condor";
my $errdir    = "/direct/phenix+spin2/beaumim/vernierScans/data_reduction/macros/condor";
my $outdir    = "$output_directory";
my $stub      = "Vernier";


my $run_number = "";
if($dst_list_name =~ /.*(\d{6}).*/m) {
  $run_number = $1;
} else {
  $run_number = "ERROR";
  die "you did not supply compatible arguments to to this script!\n";
}

my $exefile = "$exedir/$stub"."_$run_number.sh";
my $jobfile = "$exedir/$stub"."_$run_number.job";

my $jobFile =
"Universe     = vanilla
Notification = Error
Initialdir   =  
Executable   = $exefile
Arguments    =  
Log          = $jobdir/$stub"."_$run_number.log
Output       = $jobdir/$stub"."_$run_number.out
Error        = $jobdir/$stub"."_$run_number.err
Notify_user  =  \"cas\"
Queue";

my $output_file = "$output_directory/$run_number"."_reduced.root"
my $scriptFile = 
"#! /bin/bash
# THIS IS AN AUTOMATICALLY GENERATED SCRIPT
# MODIFY/RUN AT YOUR OWN RISK.
# SCRIPT: $stub.sh

root -l -b -q $rootmacro\'(\"$dst_list_name\",\"$output_file\")\'";

open EXEFILE, ">", $exefile or die $!;
open JOBFILE, ">", $jobfile or die $!;
print JOBFILE $jobFile."\n";
print EXEFILE $scriptFile."\n";

system("chmod u+x $exefile");
close EXEFILE;
close JOBFILE;
