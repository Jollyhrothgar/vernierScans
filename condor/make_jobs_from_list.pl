#!/usr/bin/perl -w 
use strict;
use warnings;

open JOBLIST, "<condor_jobs.txt" or die $!;
my $stub = "vernier_job";
my $counter = 0;

while(<JOBLIST>) {
  chomp $_;
  my $script_contents = 
"#! /bin/csh
# THIS IS AN AUTOMATICALLY GENERATED SCRIPT
# MODIFY/RUN AT YOUR OWN RISK.
# SCRIPT: $stub.csh

# Get the proper environment:
setenv HOME /phenix/u/\$LOGNAME
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
source \$i
end
source /opt/phenix/bin/phenix_setup.csh
source \$HOME/.cshrc
echo \$MYINSTALL
echo \$LD_LIBRARY_PATH

$_";
  my $out_script_name = "$stub"."_$counter.csh";
  open OUTSCRIPT, ">$out_script_name";
  print OUTSCRIPT $script_contents;
  close OUTSCRIPT;

my $jobFile =
"Universe     = vanilla
Notification = Error
Initialdir   =  
Executable   = /direct/phenix+spin2/beaumim/vernierScans/condor/$stub\_$counter.csh
Arguments    =  
Log          = /direct/phenix+spin2/beaumim/vernierScans/condor/$stub\_$counter.log
Output       = /direct/phenix+spin2/beaumim/vernierScans/condor/$stub\_$counter.out
Error        = /direct/phenix+spin2/beaumim/vernierScans/condor/$stub\_$counter.err
Notify_user  = beaumim
Queue";

  open OUTJOB,">$stub\_$counter.job";
  print OUTJOB $jobFile;
  close OUTJOB;


  $counter++;
}
