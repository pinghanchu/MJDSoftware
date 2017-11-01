#!/usr/bin/perl
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDSkim/";
my $plotdeltaT = $scriptpath."plotdeltaT";
my $enr = 67;
my $window = 5;
my $wfdata = "./data/wf_*_".$enr."_".$window.".txt";
my $wffile = "./data/wf_".$enr."_".$window.".txt";
system("cat $wfdata > $wffile"); 
system("$plotdeltaT $enr $window");
