#!/usr/bin/perl
system("./MkCookie");
sleep(5);
open(my $fin, "<", "./List/runlist/cal.list1.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], "\n";

    for($i = $array[0];$i<=$array[1];$i++){
	my $run = $i;
	my $file = $run;
	my $z   = "fill_".$run.".csh";
	my $builtfile = "./built/*/OR_run".$run.".root";
	#my $gatfile   = "./Hist/gat/mjd_run".$run.".root";
	#my $gatfile   = "./gatified/*/mjd_run".$run.".root";
	#my $histfile  = "./Hist/cal/hist_".$run.".root";
	#system("rm out.$run err.$run");
	#system("rm $gatfile");
	#system("rm $histfile");
	    
	open(my $fh, ">", $z) or die "cannot open";#
	print $fh "#!/bin/tcsh\n";
	print $fh "./process_mjd_cal $builtfile\n";
	close $fh;
	system("chmod 755 $z");
	#system("./$z");
	system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $z");
	#sleep(5);	
    }
}
