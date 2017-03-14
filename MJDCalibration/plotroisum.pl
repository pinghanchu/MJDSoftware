#!/usr/bin/perl
open(my $fin, "<", "./List/runlist/cal.range.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    my $startrun = $array[0];
    my $endrun = $array[1];
    my $run = $startrun."_".$endrun;
    my $plotpath = "./Plot/cal/".$run."/";
    my $plotfiles = "reso*".$run."*.pdf";
    system("./plotroisum $array[0] $array[1] trapENFCal");
    #system("mv $plotfiles $plotpath");
}
