#!/usr/bin/perl
#This script deletes all text file in ./List/ with size equal to zero.
#
system("ls -ls ./List/*_*/*.txt > List/list.txt");
open(my $fin, "<", "./List/list.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[5]," ", $array[9],"\n";
    if($array[5] == 0){
	system("rm $array[9]");
    }
}

