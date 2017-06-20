#!/usr/bin/perl
system("\$MGDODIR/bin/MkCookie");
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/Apps/";
my $builtpath = "\$MJDDATADIR/surfprot/data/built/P3END/";
my $datapath = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDGat/List/runlist/";
my $datafile = $datapath."bk.DSPM.txt";
my $process = $scriptpath."process_mjd_cal";
open(my $fin, "<", $datafile) or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], "\n";

    for($i = $array[0];$i<=$array[1];$i++){
	my $run = $i;
	my $file = $run;
	my $app = "fill_".$run.".csh";
	my $builtfile = $builtpath."OR_run".$run.".root";
	my $gatfile = "mjd_run".$run.".root";
	open(my $fh, ">", $app) or die "cannot open";#
	print $fh "#!/bin/tcsh\n";
	print $fh "$process $builtfile\n";
	#print $fh "mv $gatfile GAT/\n";
	close $fh;
	system("chmod 755 $app");
	#system("./$z");
	system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $app");
	#sleep(5);	
    }
}
