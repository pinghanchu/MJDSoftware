#!/usr/bin/perl
#system("./MkCookie");
#sleep(5);
open(my $fin, "<", "./List/runlist/cal.list4.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], "\n";

    for($i = $array[0];$i<=$array[1];$i++){
	my $run = $i;
	my $file = $run;
	my $z   = "fill_".$run.".csh";
	my $out = "err.".$run."*";
	my $builtfile = "./built/*/OR_run".$run.".root";
	
	#my $gatfile   = "./Hist/gat/mjd_run".$run.".root";
	#my $gatfile   = "./gatified/*/mjd_run".$run.".root";
	#my $histfile  = "./Hist/hist_".$run."*.root";
	#system("rm $histfile");
	#system("rm $out");
    }
}
