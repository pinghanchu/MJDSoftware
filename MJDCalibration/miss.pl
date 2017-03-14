#!/usr/bin/perl
open(my $fin, "<", "./List/runlist/cal.DS0.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    #print $array[0]," ", $array[1], " " , $array[2], " " , $array[3], "\n";
    for(my $i = $array[0];$i<=$array[1];$i++){
	#my $file = "./Hist/hist_".$i.".root";
	#my $file = "./Hist/gat/mjd_run".$i.".root";
	#my $file = "./out.$i";
	#print $i,"\n";
	if(-e $file){
	    #print "exist!\n";
	}else{
	    my $run = $i;
	    my $z = "fill_".$run.".csh";
	    print $i,"\n";
	    #system("qsub -l projectio=1 -cwd -o out.$run -e err.$run $z");  
	}
    }
}
