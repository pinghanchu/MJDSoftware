#!/usr/bin/perl
system("./MkCookie");
sleep(5);
open(my $fin, "<", "./runrange.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], "\n";

    for($i = $array[0];$i<=$array[1];$i++){
	my $run = $i;
	my $file = $run;
	my $app = "fill_".$run.".csh";
	my $builtfile = "./built/*/OR_run".$run.".root";
	my $gatfile = "mjd_run".$run.".root";
	open(my $fh, ">", $app) or die "cannot open";#
	print $fh "#!/bin/tcsh\n";
	print $fh "./process_mjd_cal $builtfile\n";
	print $fh "mv $gatfile GAT/\n";
	close $fh;
	system("chmod 755 $app");
	#system("./$z");
	system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $app");
	#sleep(5);	
    }
}
