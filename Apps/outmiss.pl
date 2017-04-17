#!/usr/bin/perl
system("./MkCookie");
sleep(5);
open(my $fin, "<", "./runrange.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], "\n";

    for($i = $array[0];$i<=$array[1];$i++){
	my $run = $i;
	my $file = $run;
	my $app = "fill_".$run.".csh";
	my $outfile = "out.".$run;
	if(-e $outfile){
	}else{
	    system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $app");
	}
    }
}
