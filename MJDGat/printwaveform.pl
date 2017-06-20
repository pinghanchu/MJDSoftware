#!/usr/bin/perl
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDGat/";
my $printwaveform = $scriptpath."printwaveform";
my $inputfile = "./wf.txt";
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    my $run = $array[0];
    my $entry = $array[1];
    my $chan = $array[2];
    my $enr = $array[3];
    my $nX = $array[4];
    my $ratio = $array[5];
    my $deltaT = $array[6];
    my $FFT = $array[7];
    my $aovere = $array[8];
    if($deltaT > 1000 && $ratio>1 && $ratio<10){
	print $run," ",$entry," ",$chan," " , $enr," " , $ratio," ",$deltaT,"\n";
	system("$printwaveform $run $entry $chan $enr");
    }
}



