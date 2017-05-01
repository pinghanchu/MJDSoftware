#!/usr/bin/perl

print "Please input [startrun] [endrun]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 2){
    print "You miss arguments\n";
}elsif( $numArgs >2){
    print "You have too many arguments\n";
}
my $startrun = $ARGV[0];
my $endrun = $ARGV[1];
#my $home = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDCalibration/";
my $home = "./";

print "startrun = ", $startrun,"; endrun = ", $endrun,"\n";

for($i = $startrun;$i<=$endrun; $i++){
    my $run = $i;
    my $file = $run;
    my $app  = "wf_".$file.".csh";
    my $waveform = "waveform_".$run.".root";
    my $wf = "wf_".$run.".txt";
    my $data = "data_".$run.".txt";

    open(my $fh, ">", $app) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";
    print $fh "./search $run\n";
    print $fh "mv $waveform ./Hist/\n";
    print $fh "mv $wf ./List/wf/\n";
    print $fh "mv $data ./List/data/\n";
    close $fh;
    system("chmod 755 $app");
    #system("./$app");
    system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $app");
}
