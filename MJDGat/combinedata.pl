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
my $datapath = "./data/";
my $outputfile = $datapath."wf_".$startrun."_".$endrun.".txt";
print "startrun = ", $startrun,"; endrun = ", $endrun,"\n";

system("rm $outputfile");
for($i = $startrun;$i<=$endrun; $i++){
    my $run = $i;
    my $file = $run;
    my $datafile  = $datapath."wf_".$file.".txt";
    system("cat $datafile >> $outputfile");
}
