#!/usr/bin/perl
#system("rm ./List/cal/cal1_*_*.txt");
system("ls List/cal/cal1_*_*.txt > List/cal1.list.txt");
open(my $fin1, "<", "./List/cal1.list.txt") or die "Failed to open file: $!\n";
my @startrun;
my @endrun;
while(my $line1 = <$fin1>) {
    chomp $line1;
    my @array1 = split(' ',$line1);
    my @array2 = split('_',$array1[0]);
    my @array3 = split('.txt',$array2[2]);
#    print $array2[1]," ", $array3[0],"\n";
    push(@startrun, $array2[1]);
    push(@endrun, $array3[0]);
}
my $size = @startrun;

open(my $fin, "<", "./List/runlist/cal.list4.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    my $r1 = $array[0];
    my $r2 = $array[1];
    my $run = $r1."_".$r2;
    for(my $i = 0;$i<$size;$i++){
	if($startrun[$i]>= $r1 && $endrun[$i]<=$r2){
	    my $file = "./List/cal/cal1_".$startrun[$i]."_".$endrun[$i].".txt";
	    #print $file,"\n";
	    system("rm $file");
	}
    }
    system("./plotcalibration $r1 $r2 trapENF");
}
