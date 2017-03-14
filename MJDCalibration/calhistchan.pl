#!/usr/bin/perl
#system("./MkCookie");
#sleep(1);
#my $runrange = "2361_7635";
my $runrange = "9034_14384";
#my $runrange = "16836_18351";
#my $runrange = "60000791_60001926";
#my $runrange = "18740_22282";
#open(my $fin, "<", "./List/bad.cal_$runrange.txt") or die "Failed to open file: $!\n";
#open(my $fin, "<", "./List/offseterr.cal_$runrange.txt") or die "Failed to open file: $!\n";
#open(my $fin, "<", "./List/align.cal_$runrange.txt") or die "Failed to open file: $!\n";
open(my $fin, "<", "./List/cal/error.list.txt") or die "Failed to open file: $!\n";
#open(my $fin, "<", "./List/runlist/cal.list3.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], " " , $array[2], " " , $array[3], "\n";
    my $startrun = $array[0];
    my $endrun = $array[1];
    my $channel = $array[3];
    my $energy = $array[4];
    my $run = $startrun."_".$endrun;
    my $file = $run."_".$channel;
    my $z   = "run_".$file.".csh";
    my $histfile = "./Hist/hist_".$run.".root";
    my $plotfile = "./*_".$run."_*.pdf";
    my $plotpath = "./Plot/cal/".$run."/";
    #system("rm $plotpath*.pdf");
    my $z   = "run_".$file.".csh";
    open(my $fh, ">", $z) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";	
    print $fh "./calhistchan $startrun $endrun $startrun $endrun $energy $channel $histfile\n";	
    print $fh "mv $plotfile $plotpath\n";
    close $fh;
    system("chmod 755 $z");
    #system("./$z");
    system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $z");    
    
}
