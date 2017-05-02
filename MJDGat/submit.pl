#!/usr/bin/perl
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDGat/";
my $inputfile = $scriptpath."List/runlist/bk.DS0.txt";
my $search = $scriptpath."searchpileupwf.pl";
my $enr = 67;
my $window = 5;
my $index = 1;
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1],"\n";
    system("$search $array[0] $array[1] $enr $window $index");
}
