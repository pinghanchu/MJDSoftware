#!/usr/bin/perl
my @pos1 = (11,12,13,14,15,16,17);
my @pos2 = (21,22,23,24,25,26,27);
my @pos3;
push(@pos3,@pos1);
push(@pos3,@pos2);
system("./MkCookie");
sleep(5);
open(my $fin, "<", "./List/runlist/cal.DS1.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], " " , $array[2], " " , $array[3], "\n";
    my $startrun = $array[0];
    my $endrun = $array[1];
    my $run = $array[0]."_".$array[1];
    my $z   = "run_".$run.".csh";
    open(my $fh, ">", $z) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";	
    print $fh "./uploadlinearfit $startrun $endrun $array[2] $array[3] trapENF\n";	
    close $fh;
    system("chmod 755 $z");
    system("./$z");
    #system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $z");    
}
