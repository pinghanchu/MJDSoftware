#!/usr/bin/perl
my $roi = 2039;
my $reso1;
my $resoerr1;
my $reso2;
my $resoerr2;
my $reso3;
my $resoerr3;
my $enr = "10";
my $nat = "00";
my $hg = "HG";
my @array;
my @array1;
my $startrun;
my $endrun;
my $dataset;
my $file;
my $channel;
my $enr;
my $reso;
my $resoerr;
open(my $fin, "<", "./List/runlist/cal.range.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    @array = split(' ',$line);
    $startrun = $array[0];
    $endrun = $array[1];
    $dataset = $array[2];
    $file = "./List/cal/roi_".$startrun."_".$endrun.".txt";
    @array="";
    open(my $fin1, "<", $file);
    while(my $line1 = <$fin1>){
	chomp $line1;
	@array1 = split(' ',$line1);
	$channel = $array1[1];
	$enr = $array1[5];
	$reso = $array1[8];
	$resoerr = $array1[9];
	$reso1 = sprintf("%.3f",$reso*2.355);
	$resoerr1 = sprintf("%.3f",$resoerr*2.355);

	print $dataset," & ", $channel," & \$", $reso1," \\pm ",$resoerr1,"\$","\n";
	print "\\hline\n";
    }
}
