#!/usr/bin/perl
for(my $i=0;$i<=5;$i++){
    my $DS=$i;
    open(my $fin, "<", "./List/gettime_DS$DS.txt") or die "Failed to open file: $!\n";
    my $totaltime = 0;
    while(my $line = <$fin>) {
	chomp $line;
	my @array = split(' ',$line);
	my $time = $array[5];
	$totaltime = $totaltime + $time;
    }
    print $totaltime," " , $totaltime/60, " " , $totaltime/3600,"\n";
}
