#!/usr/bin/perl
#######Good detector list:##############
#DS0
#my @pos1 = (111,112,113,114,122,123,134,141,142,143,144,145,151,152,154,162,163,171,172,173);
#DS1
#my @pos1 = (112,113,114,121,122,123,132,133,134,153,161,163,164,171,172,173,174);
#DS2
#my @pos1 = (112,113,114,121,122,123,132,133,134,153,161,163,164,171,172,173,174);
#DS3
#my @pos1 = (112,113,114,121,122,123,132,133,134,141,142,143,144,145,152,153,161,163,164,171,172,173,174);
#DS4
#my @pos2 = (211,212,213,214,221,222,223,231,232,241,242,244,251,253,254,261,262,272,273,274);
#DS5
my @pos1 = (112,113,114,122,123,132,133,134,141,142,143,144,145,152,153,161,163,164,171,172,173,174);
my @pos2 = (211,212,213,214,221,222,223,231,232,241,242,244,251,253,254,261,262,273,274);
my @pos3;
push(@pos3,@pos1);
push(@pos3,@pos2);
#########################################

print "Please input [startrun] [endrun] [cover startrun] [cover endrun] [energy name] [input file path] [output file path]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 7){
    print "You miss arguments\n";
}elsif( $numArgs >7){
    print "You have too many arguments\n";
}
my $startrun = $ARGV[0];
my $endrun = $ARGV[1];
my $coverstartrun = $ARGV[2];
my $coverendrun = $ARGV[3];
my $energy = $ARGV[4];
my $gatpath = $ARGV[5];
my $histpath = $ARGV[6];
my $dir = $startrun."_".$endrun;
#my $home = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDCalibration/";
my $home = "./";

system("mkdir /List/$dir");
system("mkdir ./Plot/$dir");
print "startrun = ", $startrun,"; endrun = ", $endrun,"; cover start run = ", $coverstartrun,"; cover end run = ", $coverendrun,"\n";
print "Energy = ", $energy, "; GAT data path = ", $gatpath,"\n";

my @pos;
if($endrun<18590){
    @pos = @pos1;
}elsif($startrun>=18590 && $endrun<4500000){
    @pos = @pos3;
}elsif($startrun>60000000 && $endrun<65000000){
    @pos = @pos2;
}

for($i = $startrun;$i<=$endrun; $i++){
    my $run = $i;
    my $gatfile = $gatpath."mjd_run".$run.".root";
    foreach $ip (@pos){
	my $file = $run."_".$energy."_".$ip;
	my $histfile = $histpath."hist_".$file.".root";
	my $app  = "fill_".$file.".csh";

	open(my $fh, ">", $app) or die "cannot open";#
	print $fh "#!/bin/tcsh\n";
	print $fh "./fillhist $run $energy $ip $gatfile $histfile\n";
	close $fh;
	system("chmod 755 $app");
    }
}
system("./deletezero.pl");
my $count = 0;
foreach $ip (@pos){
    $count = 0;
    my $file = $startrun."_".$endrun."_".$energy."_".$ip;
    my $histfile = $histpath."hist_".$file.".root";
    my $app  = "run_".$file.".csh";
    open(my $fh, ">", $app) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";

    #####Step 1: fill histograms##########
    for($i = $startrun;$i<=$endrun;$i++){
	my $run = $i;
	my $runfile = $run."_".$energy."_".$ip;
	my $app = "fill_".$runfile.".csh";
        my $file1 = "./Hist/hist_".$runfile.".root";
	my $file2 = "./Hist/hist_".$startrun."_".$endrun."_".$energy."_".$ip.".root";
	if(-e $file1 || -e $file2){
	}else{
	    print $fh "./$app\n";
	    $count++;
	}
    }
    #####Step 2: combine histograms##########
    my $file = "./Hist/hist_".$startrun."_".$endrun."_".$energy."_".$ip.".root";
    if(-e $file){
    }else{
	print $fh "./combhist $startrun $endrun $energy $ip $histpath $histfile\n";
	$count++;
    }
    #####Step 3: GATMultiPeakFitter histogram#
    my $file = "./List/".$startrun."_".$endrun."/calibration_".$startrun."_".$endrun."_".$energy."_".$ip.".txt";
    if(-e $file){
    }else{
	print $fh "./multihist $startrun $endrun $coverstartrun $coverendrun $energy $ip $histfile $home\n";
	$count++;
    }

    #####Step 4: move all txt and pdf
    my $pdf = "*_".$dir."_*.pdf";
    my $txt = "*_".$dir."_*.txt";
    print $fh "mv $pdf ./Plot/$dir/\n";
    print $fh "mv $txt ./List/$dir/\n";    

    #####Step 5: delete histogram of each run
    #for($i = $startrun;$i<=$endrun;$i++){
	#my $run = $i;
        #my $file1 = "./Hist/hist_".$run."_".$ip."_".$energy.".root";
	#print $fh "rm $file1\n";	
    #}
    close $fh;

    #####Step 6: run job#############
    if($count>0){
	my $file = $startrun."_".$endrun."_".$energy."_".$ip;
	system("chmod 755 $app");
	#system("./$app");
	system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $app");
	#system("qsub -o out.$file -e err.$file $app");
	#sleep(3*60);
    }
}
