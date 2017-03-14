#!/usr/bin/perl
system("rm gettime_*.log ./List/gettime_DS*.txt");
for(my $i=0;$i<=5;$i++){
    my $DS = $i;
    open(my $fin, "<", "./List/runlist/cal.DS".$DS.".txt") or die "Failed to open file: $!\n";
    while(my $line = <$fin>) {
	chomp $line;
	my @array = split(' ',$line);
	
	system("./gettime $DS $array[0] $array[1] >> gettime_$DS.log");
    }
}
