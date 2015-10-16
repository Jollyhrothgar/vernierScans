#!/usr/bin/perl -w 
use Cwd;
use strict;
use Getopt::Long;

sub make_job_file;
sub show_job_file;
sub read_config;
sub write_config;
sub make_root_call;
sub make_script;
sub dump_script;
sub make_job_file_script;

our $dir = getcwd;
# Job File Template - update with setup.conf 
our %JobFile = (
    Universe     => "vanilla"       , 
    Notification => "Error"         , 
    Initialdir   => " "              , # Full path to .job file
    Executable   => " "              , # Full path to executable (.sh) file
    Arguments    => " "              , 
    Log          => " "              , # Full path to .log file (condor log)
    Output       => " "              , # Full path to .out file (std::out)
    Error        => " "              , # Full path to .err file (std::err)
    Notify_user  => " "              , # Email address for notification
    GetEnv       => "True"          , 
    JobType      => "\"cas\""       , 
    Queue        => " "                
    );

our %config = (
    CONDOR_JOB_DIRECTORY    => " ", 
    EXECUTABLE_DIRECTORY    => " ",
    CONDOR_LOG_DIRECTORY    => " ",
    CONDOR_STDOUT_DIRECTORY => " ",
    CONDOR_STDERR_DIRECTORY => " ",
    MACRO_ARGUMENTS_LIST    => " ",
    EMAIL_NOTIFICATION      => " ",
    ROOT_MACRO              => " ",
    SCRIPT_NAME_STUB        => " ",
    INSTANCES_PER_JOB       => " "  
    );


# Default Value for Command Line Arguments
my $script = "";
my $help  = 0;
my $generate = 0;
my $make_jobs = 0;

# Please view README.txt for specific instructions, or call "make_condor_job.pl -h".
GetOptions(
    'script=s'  => \$script,
    'generate!' => \$generate,    
    'make-jobs!'=> \$make_jobs,
    'help!'     => \$help,        
    ) or die "Incorrect usage!\n";

################################
#      
#      COMMAND SWITCHES
#
################################

if( $help ){
  open HELP, "<","README.txt", or die "$! could not find README.txt. Are you in the right directory? Check it out again from CVS, if unsure. \n";
  while(<HELP>){
    print $_;
  }
  exit;
}

if( $generate ) {
  print "Generating setup.conf\n";
  write_config;
  exit;
}

# This guys will do most of the work for you.
if( $make_jobs ) {
  print "Reading in setup.conf\n";
  read_config;
  show_job_file;
  open ARGUMENTS, "<$config{'MACRO_ARGUMENTS_LIST'}", or die "$! ERROR $config{'MACRO_ARGUMENTS_LIST'} was not found!\n";
  my @arguments = <ARGUMENTS>;
  close ARGUMENTS;
  my $job_counter = 0; 
  my $file_counter = 1;
  my @root_calls = (); # this is an array of every call to send to root.
  my @root_call_list = ();
  my $root_calls_per_job = $config{'INSTANCES_PER_JOB'};
  my @job_files;
  my $call_counter = 0;

  # we can now split jobs arbitrarily. For each script we dump, we need to dump a job file as well.
  while(@arguments) {
    my $args = pop @arguments;
    @root_call_list = (); # tokens which will construct a call to root
    push @root_call_list, $config{'ROOT_MACRO'}; # the macro + path which is called after root.exe
    my @root_call_args = split(' ',$args); # the arguments given to the root macro
    push @root_call_list, @root_call_args; # an array of the macro + the macro's arguments
    my $root_call = make_root_call(@root_call_list); # concatinate the array with the right number of escape characters
    push @root_calls, $root_call;
    if( scalar(@root_calls) == $root_calls_per_job ) {
      $job_counter++;
      my @script_array = make_script($job_counter,@root_calls);
      dump_script($job_counter,@script_array);
      push @job_files, make_job_file_script($job_counter);
      @root_calls = ();
    }
    $call_counter++;
  }
  # leftovers 
  if( scalar(@root_calls) != 0 ) {
    $job_counter++;
    my @script_array = make_script($job_counter,@root_calls);
    dump_script($job_counter,@script_array);
    push @job_files, make_job_file_script($job_counter);
  }
  open CONDOR_SUBMIT, ">$dir/condor_submit.sh", or die "$! could not open condor_submit.sh for reading\n";
  print CONDOR_SUBMIT "#! /bin/bash\n";
  print CONDOR_SUBMIT "# THIS SCRIPT WILL SUBMIT ALL JOBS MADE FROM make_condor_job.pl\n";
  foreach(@job_files) {
    print CONDOR_SUBMIT "condor_submit $_\n"
  }
  close CONDOR_SUBMIT;
  system("chmod u+x $dir/condor_submit.sh");
  print "  Generated $job_counter jobs. Split over these jobs, we execute the root macro $call_counter times.\n";
  print "   Sumbit with: condor_sumbit.sh\n\n";
}


if( $script ) {
  my $stub = "";
  print "Making a condor job for: $script\n";
  if( $script =~ /(\S+)\.\S+/m)
  {
    $stub = $1;   
  }
  elsif( $script =~ /(\S+)/m)
  {
    $stub = $1;
  }
  else
  {
    die "You must supply a non-blank script name\n";
  }
  make_job_file $stub, $script, $dir; 
  exit;
}

##################################
#    
#    SUBROUTINE DEFINITIONS
# 
##################################


sub make_job_file {
  my $arg_num = scalar(@_);
  if( $arg_num != 3 ) 
  {
    die "cannot make job file";
  }
  my ($stub, $script_name, $directory) = @_;
  $JobFile{"Initialdir"  } =  "$directory/$stub.job"; # Full path to .job file
  $JobFile{"Executable"  } =  "$directory/$script_name"  ; # Full path to executable (.sh) file
  $JobFile{"Log"         } =  "$directory/$stub.log"; # Full path to .log file (condor log)
  $JobFile{"Output"      } =  "$directory/$stub.out"; # Full path to .out file (std::out)
  $JobFile{"Error"       } =  "$directory/$stub.err"; # Full path to .err file (std::err)
  $JobFile{"Notify_user" } =  ""              ; # Email address for notification

  my $out_file = $directory."/".$stub.".job";
  open JOBFILE, ">$out_file", or die "$! could not open $out_file for writing\n";
  print JOBFILE "Universe     = $JobFile{'Universe'    }\n";
  print JOBFILE "Notification = $JobFile{'Notification'}\n";
  print JOBFILE "Initialdir   = $JobFile{'Initialdir'  }\n";
  print JOBFILE "Executable   = $JobFile{'Executable'  }\n";
  print JOBFILE "Arguments    = $JobFile{'Arguments'   }\n";
  print JOBFILE "Log          = $JobFile{'Log'         }\n";
  print JOBFILE "Output       = $JobFile{'Output'      }\n";
  print JOBFILE "Error        = $JobFile{'Error'       }\n";
  print JOBFILE "Notify_user  = $JobFile{'Notify_user' }\n";
  print JOBFILE "GetEnv       = $JobFile{'GetEnv'      }\n";
  print JOBFILE "+JobType     = $JobFile{'JobType'     }\n";
  print JOBFILE "Queue\n";
  close JOBFILE;
}

sub make_root_call {
  my $root_macro = $_[0];
  my @arguments = @_[1..(scalar(@_)-1)];
  my $root_call = "root -l -b -q $root_macro\'(";
  for(my $i = 0; $i < scalar(@arguments); $i++)
  {
    if($i != scalar(@arguments) -1) {
      $root_call = $root_call."$arguments[$i], ";
    }
    else
    {
      $root_call = $root_call."$arguments[$i]";
    }
  }
  $root_call = $root_call.")\'";
  return $root_call;
}

sub make_job_file_script {
  
  my $job_number = $_[0];

  my $job_dir = $config{'CONDOR_JOB_DIRECTORY'};
  my $job_file_name = "$job_dir/$config{'SCRIPT_NAME_STUB'}_$job_number.job"; 
  open JOB_FILE, ">$job_file_name",or die "$! could not open $job_file_name for writing";
  print JOB_FILE "Universe     = $JobFile{'Universe'    }\n";
  print JOB_FILE "Notification = $JobFile{'Notification'}\n";
  print JOB_FILE "Initialdir   = $JobFile{'Initialdir'  }\n";
  print JOB_FILE "Executable   = $config{'EXECUTABLE_DIRECTORY'}/$config{'SCRIPT_NAME_STUB'}_$job_number.sh\n";
  print JOB_FILE "Arguments    = $JobFile{'Arguments'   }\n";
  print JOB_FILE "Log          = $config{'CONDOR_LOG_DIRECTORY'}/$config{'SCRIPT_NAME_STUB'}_$job_number.log\n";
  print JOB_FILE "Output       = $config{'CONDOR_STDOUT_DIRECTORY'}/$config{'SCRIPT_NAME_STUB'}_$job_number.out\n";
  print JOB_FILE "Error        = $config{'CONDOR_STDERR_DIRECTORY'}/$config{'SCRIPT_NAME_STUB'}_$job_number.err\n";
  print JOB_FILE "Notify_user  = $config{'EMAIL_NOTIFICATION'}\n";
  print JOB_FILE "GetEnv       = $JobFile{'GetEnv'      }\n";
  print JOB_FILE "+JobType     = $JobFile{'JobType'     }\n";
  print JOB_FILE "Queue"."\n";
  close JOB_FILE;
  return $job_file_name;
}

sub write_config
{
  my $check_file = "$dir/setup.conf";
  if( -e $check_file ) {
    print "Setup file already exists, delete it, then run:\n";
    print "    make_condor_job.pl --generate\n\n";
    exit;
  }
  open SETUP_FILE,">$dir/setup.conf", or die "$! could not open $dir/setup.conf for reading\n";
  foreach my $key ( keys %config ) {
    printf SETUP_FILE ("%-30s = %s\n", $key, $config{$key});
  }
  close SETUP_FILE;
}

sub read_config
{
  my $check_file = "$dir/setup.conf";
  if( -e $check_file ) {
    print "Using existing setup.conf to prepare condor jobs\n";
  }
  else {
    print "No setup file found. Please run: \n";
    print "    make_condor_job.pl --generate\n\n";
    print "to setup running multiple parallel jobs\n";
    exit;
  }
  open SETUP_FILE,"<$dir/setup.conf", or die "$! could not open $dir/setup.conf for reading\n";
  while(<SETUP_FILE>) {
    chomp $_;
    if( $_ =~ /\s*(\S+)\s*=\s*(\S+)\s*/m ) {
      $config{$1} = $2;
      print "Read in $1 as : $config{$1} \n";
      
    }
  }
}

sub make_script {
  my $job_counter = $_[0];
  my @root_calls = @_[1..(scalar(@_)-1)];
  
  my @script_array = ();
  push @script_array, "#! /bin/bash";
  push @script_array, "# THIS IS AN AUTOMATICALLY GENERATED SCRIPT from make_condor_job.pl";
  push @script_array, "# MODIFY/RUN AT YOUR OWN RISK.";
  push @script_array, "# SCRIPT: $config{'SCRIPT_NAME_STUB'}_$job_counter.sh";
  push @script_array, "";
  for(my $i = 0; $i < scalar(@root_calls); $i++) {
    push @script_array, $root_calls[$i]
  }
  return @script_array;
}

sub dump_script {
  my $job_counter = $_[0];
  my @script_array = @_[1..scalar(@_)-1];
  my $script_file = "$config{'EXECUTABLE_DIRECTORY'}/$config{'SCRIPT_NAME_STUB'}_$job_counter.sh";
  open SCRIPT_FILE, ">$script_file", or die "$! could not open $script_file for writing\n";
  foreach (@script_array) {
    print SCRIPT_FILE $_, "\n";
  }
  close SCRIPT_FILE;
  system("chmod u+x $script_file");
}

sub show_job_file {
  print "TEMPLATE JOB FILE: \n";
  print "\tUniverse     = $JobFile{'Universe'    }\n";
  print "\tNotification = $JobFile{'Notification'}\n";
  print "\tInitialdir   = $JobFile{'Initialdir'  }\n";
  print "\tExecutable   = $config{'EXECUTABLE_DIRECTORY'}/$config{'SCRIPT_NAME_STUB'}.sh\n";
  print "\tArguments    = $JobFile{'Arguments'   }\n";
  print "\tLog          = $config{'CONDOR_LOG_DIRECTORY'}/$config{'SCRIPT_NAME_STUB'}.log\n";
  print "\tOutput       = $config{'CONDOR_STDOUT_DIRECTORY'}/$config{'SCRIPT_NAME_STUB'}.out\n";
  print "\tError        = $config{'CONDOR_STDERR_DIRECTORY'}/$config{'SCRIPT_NAME_STUB'}.err\n";
  print "\tNotify_user  = $config{'EMAIL_NOTIFICATION'}\n";
  print "\tGetEnv       = $JobFile{'GetEnv'      }\n";
  print "\t+JobType     = $JobFile{'JobType'     }\n";
  print "\tQueue"."\n";
}
