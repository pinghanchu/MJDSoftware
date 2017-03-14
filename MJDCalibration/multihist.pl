#!/usr/bin/perl
my @pos1 = (11,12,13,14,15,16,17);
my @pos2 = (21,22,23,24,25,26,27);
my @pos3;
push(@pos3,@pos1);
push(@pos3,@pos2);
#system("./MkCookie");
#sleep(1);
open(my $fin, "<", "./List/runlist/cal.list3.txt") or die "Failed to open file: $!\n";
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
    my $startrun = $array[0];
    my $endrun = $array[1];
    my $run = $startrun."_".$endrun;
    my $histfile = "./Hist/hist_".$run.".root";
    my $plotfile = "./*_".$run."_*.*";
    my $plotpath = "./Plot/cal/".$run."/";
    my $listfile = "./List/cal/*_".$run.".txt";
    my $file = $run."_".$ip;

    #system("mkdir $plotpath");
    #system("rm $plotpath*.pdf");
    #system("rm $listfile");
    foreach $ip (@pos){
	my $file = $run."_".$ip;
	my $z   = "run1_".$run."_".$ip.".csh";
	open(my $fh, ">", $z) or die "cannot open";#
	print $fh "#!/bin/tcsh\n";	
	print $fh "./multihist $startrun $endrun trapE $ip\n";	
	print $fh "mv $plotfile $plotpath\n";
	close $fh;
	system("chmod 755 $z");
	#system("./$z");
	system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $z");    
    }
}
