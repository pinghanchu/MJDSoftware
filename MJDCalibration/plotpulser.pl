#!/usr/bin/perl
open(my $fin, "<", "./List/runlist/bk.list1.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    system("./plotpulser $array[0] $array[1] trapENF");
}
