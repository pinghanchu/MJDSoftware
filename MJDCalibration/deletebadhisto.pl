#!/usr/bin/perl
#This script deletes bad histograms with size less than 300
#
system("ls -sl ./Hist/hist_*_*_*.root > List/list.txt");
open(my $fin, "<", "./List/list.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[5]," ", $array[9],"\n";
    if($array[5]<300){
	system("rm $array[9]");
    }
}

