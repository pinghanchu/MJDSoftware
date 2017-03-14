#!/usr/bin/perl
system("ls -ls ./out.*_* > List/miss.txt");
open(my $fin, "<", "./List/miss.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    my @array1 = split('out.',$array[9]);
    my @array2 = split('_',$array1[1]);
    my $run1 = $array2[0];
    my $run2 = $array2[1];
    my $pos = $array2[2];
    my $run = $run1."_".$run2;
    #$print $run," ",$pos,"\n";
    #print $array[5]," ", $array[9],"\n";
    #if($array[5]<200 || $array[5]>500){	
    if($array[5]<10){
	#print $run,"\n";
	my $file = $run."_".$pos;
	my $z = "reso_".$file.".csh";
	#system("rm out.$run\n");
	print $z,"\n";
	system("chmod 755 *.csh");
	#system("./$z");
	system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $z")	
    }
}
