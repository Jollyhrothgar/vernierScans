#!/usr/bin/perl -w
use Cwd;
use strict;
use warnings;


#
# THIS FILE RECURSIVELY CLEANS, THEN BUILDS EVERYTING IT ASSUMES YOU HAVE LEFT
# THE RELATIVE DIRECTORY STRUCTURE EXACTLY CONSISTANT WITH WHAT IS CHECKED IN TO
# GITHUB.
#

my $base_dir = getcwd;

# DIRECTORIES

my @build_dir = (
 "$base_dir/dst_analysis/build",
 "$base_dir/hourglass_correction/build",
 "$base_dir/scaler_event_prdf_analysis/build",
 "$base_dir/data_event_prdf_analysis/build",
 );

foreach my $dir (@build_dir) {
# CLEAN 
  print "CLEANING: ",$dir, "\n";
  chdir("$dir/..") or die "cannot change: $!\n";
  system("./fullClean.csh");
# BUILD
  print "BUILDING: ",$dir,"\n";
  chdir("$dir") or die "cannot change: $!\n";
  system("../src/autogen.sh --prefix=\$MYINSTALL");
  system("make -j4 install"); # hopefully your machine has at least four logical cores!
}

print "Everything is fresh and new!\n";
print "Confirm libraries are at their newest state:\n";
print "System time: ";
system("date");
print "\n";
system ("ls -ltrh \$MYINSTALL/lib | grep -i Vernier");
