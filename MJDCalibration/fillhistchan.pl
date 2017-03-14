#!/usr/bin/perl
my @pos1 = (11,12,13,14,15,16,17);
my @pos2 = (21,22,23,24,25,26,27);
my @pos3;
push(@pos3,@pos1);
push(@pos3,@pos2);
system("./MkCookie");
sleep(5);
open(my $fin, "<", "./List/runlist/cal.list3.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], " " , $array[2], " " , $array[3], "\n";
    my $startrun = $array[0];
    my $endrun = $array[1];
    my $channel = $array[3];
    my $energy = $array[4];

    for($i = $startrun;$i<=$endrun;$i++){
	my $run = $i;
	my $z   = "fill_".$run."_".$channel.".csh";
	my $file = $run."_".$channel;
	open(my $fh, ">", $z) or die "cannot open";#
	print $fh "#!/bin/tcsh\n";
	print $fh "./fillhistchan $run $energy 800000 0 8000 300000 0 3000 $channel\n";
	close $fh;
	system("chmod 755 $z");
	#system("./$z");
	system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $z");	
    }
}
