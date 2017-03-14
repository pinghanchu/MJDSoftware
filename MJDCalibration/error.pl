#!/usr/bin/perl
open(my $fin, "<", "./List/error.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], " " , $array[2], " " , $array[3], "\n";
    my $startrun = $array[2];
    my $endrun = $array[3];
    my $channel = $array[1];
    my @enr  = split('Cal',$array[4]);
    my $energy = $enr[0];
    my $peak = $array[5];
    my $run = $startrun."_".$endrun;
    my $file = $run."_".$channel;
    my $z   = "run_".$file.".csh";
    my $histfile = "./Hist/hist_".$run.".root";
    my $plotfile = "./*_".$run."_*.pdf";
    my $plotpath = "./Plot/cal/".$run."/";
    my $z   = "run_".$file.".csh";
    if($peak>2510.){
	open(my $fh, ">", $z) or die "cannot open";#
	print $fh "#!/bin/tcsh\n";	
	print $fh "./calhistchan $startrun $endrun $startrun $endrun $energy $channel $histfile\n";	
	print $fh "mv $plotfile $plotpath\n";
	close $fh;
	system("chmod 755 $z");
	#system("./$z");
	#system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $z");    
    }
}
