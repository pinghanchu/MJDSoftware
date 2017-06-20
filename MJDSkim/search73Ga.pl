#!/usr/bin/perl

print "Please input [dataset] [energy] [time]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 3){
    print "You miss arguments\n";
}elsif( $numArgs >3){
    print "You have too many arguments\n";
}
my $dataset = $ARGV[0];
my $enr = $ARGV[1];
my $time = $ARGV[2];
my $timems = $time*1000;
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDSkim/";

print "dataset = ", $dataset, "; energy = ", $enr, "; time = ", $time,"\n";
my $search =$scriptpath."search73Ga";
my $datapath = "./data/";
system("mkdir ./data/");
#system("mkdir $datapath");
my $file = $dataset."_".$enr."_".$timems;
my $app  = "wf_".$file.".csh";
my $waveform = "waveform_".$file.".root";
my $wf = "wf_".$file.".txt";
my $data = "data_".$file.".txt";

open(my $fh, ">", $app) or die "cannot open";#
print $fh "#!/bin/tcsh\n";
print $fh "$search $dataset $enr $time\n";
#print $fh "mv $waveform $datapath\n";
#print $fh "mv $wf $datapath\n";
#print $fh "mv $data $datapath\n";
close $fh;
system("chmod 755 $app");
#system("./$app");
system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $app");
