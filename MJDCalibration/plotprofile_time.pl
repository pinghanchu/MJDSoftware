#!/usr/bin/perl
open(my $fin, "<", "./List/runlist/cal.long.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    system("./profile_time $array[0] $array[1] 2614 10");
}
