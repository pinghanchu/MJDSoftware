#!/usr/bin/perl
open(my $fin, "<", "./List/runlist/cal.list.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    system("./plotcovariance $array[0] $array[1] trapENF");
}
