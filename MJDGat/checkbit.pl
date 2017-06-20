#!/usr/bin/perl
# get run list information based on 
# https://github.com/mppmu/GAT/blob/master/Apps/DataSetInfo.hh
my $inputfile = "./List/runlist/cal.DS5.txt";
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
my $outfile = "checkbit.txt";
while(my $line = <$fin>) {
    chomp $line;
    my @array1 = split(" ",$line);
    for(my $i = $array1[0];$i<=$array1[1];$i++){
	system("./checkbit $i >> $outfile");
    #system("./checkbit $array1[0] >> $outfile");
    #system("./checkbit $array1[1] >> $outfile");
    }
}

