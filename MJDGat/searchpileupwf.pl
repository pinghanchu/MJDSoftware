#!/usr/bin/perl

print "Please input [startrun] [endrun] [energy] [window] [index]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 5){
    print "You miss arguments\n";
}elsif( $numArgs >5){
    print "You have too many arguments\n";
}
my $startrun = $ARGV[0];
my $endrun = $ARGV[1];
my $enr = $ARGV[2];
my $window = $ARGV[3];
my $index = $ARGV[4];
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDGat/";

print "startrun = ", $startrun,"; endrun = ", $endrun,"; energy = ", $enr, "; window = ", $window,"\n";
my $search =$scriptpath."searchpileupwf";
my $datapath = "./data/".$enr."/";
#system("mkdir ./data/");
#system("mkdir $datapath");

my $entry = int(($endrun-$startrun)/$index);
for($i = 0;$i<$entry;$i++){
    my $r1 = $startrun+$index*$i;
    my $r2 = $startrun+$index*($i+1)-1;
    my $dir = $r1."_".$r2;
    my $app  = "wf_".$dir.".csh";
    open(my $fh, ">", $app) or die "cannot open";#
print $fh "#!/bin/tcsh\n";
    for($j = $r1;$j<=$r2; $j++){
	my $run = $j;
	my $waveform = "waveform_".$run.".root";
	my $wf = "wf_".$run.".txt";
	my $data = "data_".$run.".txt";
	print $fh "$search $run $enr $window\n";
	#print $fh "mv $waveform $datapath\n";
	#print $fh "mv $wf $datapath\n";
	#print $fh "mv $data $datapath\n";
    }

    close $fh;
    system("chmod 755 $app");
#system("./$app");
    system("qsub -l projectio=1 -cwd -o out.$dir -e err.$dir $app");
}
my $last = $endrun-$entry*$index;
if($last>0){
    my $r1 = $entry*$index+$startrun;
    my $r2 = $endrun;
    my $dir = $r1."_".$r2;
    my $app  = "wf_".$dir.".csh";
    open(my $fh, ">", $app) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";
    for($j = $r1;$j<=$r2; $j++){
        my $run = $j;
        my $waveform = "waveform_".$run.".root";
        my $wf = "wf_".$run.".txt";
        my $data = "data_".$run.".txt";
        print $fh "$search $run $enr $window\n";
        #print $fh "mv $waveform $datapath\n";
        #print $fh "mv $wf $datapath\n";
        #print $fh "mv $data $datapath\n";
    }

    close $fh;
    system("chmod 755 $app");
#system("./$app");
    system("qsub -l projectio=1 -cwd -o out.$dir -e err.$dir $app");
}
