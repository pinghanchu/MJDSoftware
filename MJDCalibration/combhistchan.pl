#!/usr/bin/perl
#system("./MkCookie");
sleep(1);
open(my $fin, "<", "./List/runlist/cal.list3.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], " " , $array[2], " " , $array[3], "\n";

    my $startrun = $array[0];
    my $endrun = $array[1];
    my $channel = $array[3];
    my $energy = $array[4];
    my $run = $startrun."_".$endrun;    
    my $z   = "run_".$run."_".$channel.".csh";
    open(my $fh, ">", $z) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";
    print $fh "./combhistchan $startrun $endrun $energy 800000 0 8000 300000 0 3000 $channel\n";
    close $fh;
    system("chmod 755 $z");
    system("./$z");
    #system("qsub -l projectio=1 -cwd -o out.$run -e err.$run $z");    
}
