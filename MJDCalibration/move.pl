#!/usr/bin/perl
#system("./MkCookie");
#sleep(5);
open(my $fin, "<", "./List/runlist/cal.list5.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], "\n";
    my $startrun = $array[0];
    my $endrun = $array[1];

    my $run = $startrun."_".$endrun;
    my $file = $run;
    my $z   = "*_".$run.".csh";
    my $err = "err.".$run."*";
    my $out = "out.".$run."*";
    my $builtfile = "./built/*/OR_run".$run.".root";
    my $plotfile = "./Plot/cal/".$run."/peak_trapENF_".$run."_*.pdf";
	#my $gatfile   = "./Hist/gat/mjd_run".$run.".root";
	#my $gatfile   = "./gatified/*/mjd_run".$run.".root";
	#my $histfile  = "./Hist/hist_".$run."*.root";
	#system("rm $histfile");
    system("cp $plotfile ./temp/");
}
