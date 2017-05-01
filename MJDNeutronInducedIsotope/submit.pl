#!/usr/bin/perl
my $inputfile = "./List/runlist/bk.DS0.txt";
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1],"\n";
    system("./search.pl $array[0] $array[1]");
}
