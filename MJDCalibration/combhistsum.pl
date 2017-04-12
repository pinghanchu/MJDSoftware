#!/usr/bin/perl
my @energy = ("trapENF","trapENM","trapE");
foreach $ienr (@energy){
    open(my $fin, "<", "./runrange.txt") or die "Failed to open file: $!\n";
    while(my $line = <$fin>) {
	chomp $line;
	my @array = split(' ',$line);
	my $startrun = $array[0];
	my $endrun = $array[1];
	my $coverstartrun = $array[2];
	my $coverendrun = $array[3];
	my $enr = $ienr;
	my $home = "./";
	my $inputpath = $home."Hist/";
	my $outputfile = $home."Hist/hist_".$startrun."_".$endrun."_".$ienr.".root";
	print $startrun," ", $endrun," ",$coverstartrun," ",$coverendrun,"\n";
	
	system("./combhistsum $startrun $endrun $ienr $inputpath $outputfile");
    }
}
