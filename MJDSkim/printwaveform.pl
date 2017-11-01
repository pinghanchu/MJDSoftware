#!/usr/bin/perl
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/git/MJDSoftware/MJDSkim/";
my $printwaveform = $scriptpath."printwaveform";
my $inputfile = "./pileup.txt";
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(',',$line);

    my $run = $array[3];
    my $entry = $array[4];
    my $chan = $array[6];
    my $enr = $array[7];
    my $ratio = $array[8];
    my $deltaT = $array[9];
    my $avse = $array[10];
    my $dcr = $array[11];
    my $tailmin = $array[12];
    my $file = "waveform_".$run."_".$entry."_".$chan.".pdf";
    if(-e $file){
    }else{
	if($ratio>1.5 && $ratio<6 && $deltaT>1000 && $avse>-3){
	    print $run," ",$entry," ",$chan," " , $enr," " , $ratio," ",$deltaT," ",$avse, " ",$dcr," ",$tailmin,"\n";
	    system("$printwaveform $run $entry $chan $enr");
	}
    }
}




