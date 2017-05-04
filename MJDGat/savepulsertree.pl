#!/usr/bin/perl
print "Please input [Data Set] [IsPulser]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 2){
    print "You miss arguments\n";
}elsif( $numArgs >2){
    print "You have too many arguments\n";
}
my $index = $ARGV[0];
my $ispulser = $ARGV[1];
my $dataset;
my $datapath;
if($index ==0){
#DS0
    $dataset = "DS0";
    $datapath = "P3JDY";
}elsif($index ==1){
#DS1
    $dataset = "DS1";
    $datapath = "P3KJR";
}elsif($index ==2){
#DS2
    $dataset = "DS2";
    $datapath = "P3KJR";
}elsif($index ==3){
#DS3
    $dataset = "DS3";
    $datapath = "P3KJR";
}elsif($index ==4){
#DS4
    $dataset = "DS4";
    $datapath = "P3LQG";
}elsif($index ==5){
#DS5
    $dataset = "DS5";
    $datapath = "P3LQK";
}

my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDGat/";
my $home = "./";
#my $inputpath = "\$MJDDATADIR/surfmjd/analysis/pulser/".$datapath."/";
#my $inputpath = "./Gat/";
#my $inputpath = $home;
my $outputpath = "./Hist/";
system("mkdir $outputpath");
my $outputname;
if($ispulser == 1){
    $outputname = "pulser";
}elsif($ispulser == 0){
    $outputname = "hist";
}

my $savepulsertree = $scriptpath."savepulsertree";
my $inputfile = $scriptpath."List/runlist/bk.".$dataset.".txt";
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1],"\n";
    
    my $startrun = $array[0]; 
    my $endrun = $array[1]; 
    for(my $i = $startrun;$i<=$endrun;$i++){
	my $outputfile = $outputpath.$outputname."_".$i;	
	system("$savepulsertree $i $outputfile $ispulser");
    }
}
