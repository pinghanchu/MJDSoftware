#!/usr/bin/perl
my $dataset = "P3LQK";
system("rm *.csh");
system("mkdir Plot");
system("mkdir List");
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
	my $inputpath = "/global/project/projectdirs/majorana/data/mjd/surfmjd/data/gatified/".$dataset."/";
	my $outputpath = $home."Hist/";
	print $startrun," ", $endrun," ",$coverstartrun," ",$coverendrun,"\n";
	
	#####Step 1: calibration###############
	system("./deletezero.pl");
	system("./deletebadhisto.pl");
	system("./calibration.pl $startrun $endrun $coverstartrun $coverendrun $enr $inputpath $outputpath");
	
	######################################################################
	
	my $calibrationfile = "calibration_".$startrun."_".$endrun."_".$enr;
	my $calibrationtex = $calibrationfile.".tex";
	my $inputpath = $home."Plot/".$startrun."_".$endrun."/";
	
	######Step 2: generate calibration.pdf#############
	#system("./calibrationdoc.pl $startrun $endrun $enr $inputpath");
	#system("cp calibration.tex $calibrationtex");
	#system("pdflatex $calibrationtex");
	
	#####Step 3: upload calibration parameters#########
	my $inputpath = $home."List/".$startrun."_".$endrun."/calibration_".$startrun."_".$endrun."_".$enr."_*.txt";
	my $outputfile = $home."List/".$startrun."_".$endrun."/calibration_".$startrun."_".$endrun.".txt";
	#system("cat $inputpath > $outputfile");
	#system("./MkCookie");
	#system("./uploadlinearfit $startrun $endrun $coverstartrun $coverendrun $enr");
    }
}
