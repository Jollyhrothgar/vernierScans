#! /usr/bin/perl

#for(my $i = 0; $i < 26; $i++){
for(my $i = 0; $i < 1; $i++){
  my $ii = sprintf("%02d", $i);
	my $root_call = "root -l -b -q Run_HourglassSimulation.C\'\(\"359711/seed_conf/359711_step_$ii.conf\",\"359711/step_$ii\",\"359711_step_$ii\",1,2\)\' > 359711/step_$ii/logfile.txt";
	#system($root_call);
	print $root_call, "\n";
	system($root_call);
}
