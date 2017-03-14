#!/usr/bin/perl
#system("rm ./List/cal/ezfit1_*_*.txt");
open(my $fin, "<", "./List/runlist/cal.list2.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    #system("./plotezfit $array[0] $array[1] trapENF 0");
    #system("./plotezfit $array[0] $array[1] trapENF 1");
    #system("./plotezfit $array[0] $array[1] trapENF 2");
    system("./plotezfit $array[0] $array[1] trapENMCal 3");
}
