#!/usr/bin/perl
system("ls -ls List/pulser/pulser_*.txt > List/miss.txt");
open(my $fin, "<", "./List/miss.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    #print $array[5]," ", $array[9],"\n";
    if($array[5]<3000){
	my @array1 = split('List/pulser/pulser_',$array[9]);
	my @array2 = split('.txt',$array1[1]);
	#my @array3 = split('_',$array2[0]);
	my $run = $array2[0];
	my $z = "run_".$run.".csh";
	#system("rm hist_$run.root\n");
	print $z,"\n";
	#system("chmod 755 *.csh");
	#system("./$z");
	system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $z")	
    }
}
