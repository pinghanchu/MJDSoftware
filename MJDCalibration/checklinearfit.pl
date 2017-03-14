#!/usr/bin/perl
system("rm ./List/cal/error*.txt");
#open(my $fin, "<", "./List/runlist/cal.DS3.txt") or die "Failed to open file: $!\n";
open(my $fin, "<", "./List/error.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    
    #system("./checklinearfit $array[0] $array[1] trapENF");
    my $startrun = $array[2];
    my $endrun = $array[3];
    my @enr = split("Cal",$array[4]);
    my $energy = $enr[0];
    my $channel = $array[1];
    #my $startrun = $array[0];
    #my $endrun = $array[1];
    #my $energy = $array[4];
    #my $channel = $array[3];

    print $startrun," ",$endrun," ",$energy," ",$channel,"\n";
    system("./checklinearfit $startrun $endrun $energy $channel");
}
system("cat ./List/cal/error_*.txt >> ./List/cal/error.list.txt");
system("cat ./List/cal/errorchan_*.txt >> ./List/cal/error.chan.txt");
system("rm ./List/cal/error_*.txt ./List/cal/errorchan_*.txt");
