#!/usr/bin/perl
# get run list information based on 
# https://github.com/mppmu/GAT/blob/master/Apps/DataSetInfo.hh
my $inputfile = "test.txt";
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
my $i = 0;
my @runlist;
my $totaltime = 0;
while(my $line = <$fin>) {
    chomp $line;
    my @array1 = split(" ",$line);
    my $time = $array1[3]/1e8/60;
    $totaltime = $time + $totaltime;
}
print $totaltime,"\n";
