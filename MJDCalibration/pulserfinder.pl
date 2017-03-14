#!/usr/bin/perl
my @pos1 = (11,12,13,14,15,16,17);
my @pos2 = (21,22,23,24,25,26,27);
my @pos3;
push(@pos3,@pos1);
push(@pos3,@pos2);
#system("./MkCookie");
#sleep(1);
open(my $fin, "<", "./List/runlist/bk.list0.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], "\n";
    my @pos;
    if($array[1]<18590){
        @pos = @pos1;
    }elsif($array[0]>=18590 && $array[1]<4500000){
        @pos = @pos3;
    }elsif($array[0]>60000000 && $array[1]<65000000){
	@pos = @pos2;
    }
    for($i = $array[0];$i<=$array[1];$i++){
	my $run = $i;
	my $z   = "run_".$run.".csh";
	my $builtfile = "./built/*/OR_run".$run.".root";
	my $gatfile   = "./gatified/*/mjd_run".$run.".root";
	#my $histfile  = "./hist_".$run.".root";
	#system("rm out.$run err.$run");
	open(my $fh, ">", $z) or die "cannot open";#
	print $fh "#!/bin/tcsh\n";
	#print $fh "rm ./List/pulser/pulser_$run.txt\n";
	foreach $ip (@pos){
	    my $file = $run."_".$ip;
	    my $z   = "pulser_".$run."_".$ip.".csh";
            open(my $fh, ">", $z) or die "cannot open";#
            print $fh "#!/bin/tcsh\n";
	    print $fh "./pulserfinder $run trapENF $ip $gatfile\n";
	    close $fh;
	    system("chmod 755 $z");
	    #system("./$z");
	    system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $z");
	}
    }
}
