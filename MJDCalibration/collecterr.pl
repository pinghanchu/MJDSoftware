#!/usr/bin/perl
my @pos1 = (11,12,13,14,15,16,17);
my @pos2 = (21,22,23,24,25,26,27);
my @pos3;
push(@pos3,@pos1);
push(@pos3,@pos2);
#system("./MkCookie");
#sleep(1);
open(my $fin, "<", "./List/runlist/cal.DS5.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], " " , $array[2], " " , $array[3], "\n";

    my @pos;
    if($array[1]<18590){
        @pos = @pos1;
    }elsif($array[0]>=18590 && $array[1]<4500000){
        @pos = @pos3;
    }elsif($array[0]>60000000 && $array[1]<65000000){
	@pos = @pos2;
    }
    my $startrun = $array[0];
    my $endrun = $array[1];
    my $run = $startrun."_".$endrun;
    my $linearfile = "./List/cal/linearfiterr_trapENF_".$run."_*.txt";
    my $linearpath = "./List/cal/linearfiterr_".$run.".txt";
    system("cat $linearfile >> $linearpath");
}
