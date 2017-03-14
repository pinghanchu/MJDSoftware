#!/usr/bin/perl
open(my $fin, "<", "./bk.list.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    for(my $i = $array[0];$i<$array[1];$i++){
	my $run = $i;
    #my $file = "./Plot/cal/".$run."/";
	my $file1 = "*_".$run."_*.pdf";
	print $file1,"\n";
	#system("mkdir $file");
	system("mv $file1 ./Plot/pulser/");
    }
}
