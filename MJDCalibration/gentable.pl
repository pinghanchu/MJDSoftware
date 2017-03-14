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
	print $hg,",",$channel,",",$enr,",",$roi,",",$reso,",",$reso*2.355,"\n";
	if( $enr == $roi){
	    if( $channel == $nat){		
		$reso1 = sprintf("%.2f",$reso*2.355);
		$resoerr1 = $resoerr;
		#print $startrun," ",$endrun," ",$channel," ",$enr," ", $reso," " , $reso1,"\n";
            }elsif( $channel == $hg){
                $reso3 = sprintf("%.2f",$reso*2.355);
                $resoerr3 = $resoerr;
                #print $startrun," ",$endrun," ",$channel," ",$enr," ", $reso," " , $reso3,"\n";
	    }elsif( $channel == $enr || $channel == 10){
		$reso2 = sprintf("%.2f",$reso*2.355);
		$resoerr2 = $resoerr;
		#print $startrun," ",$endrun," ",$channel," ",$enr," ", $reso," " , $reso2,"\n";
	    }

	}
	@array1 = "";
    }
    print $dataset," & ", $reso1, " & ", $reso2, " & " , $reso3, "\n";
    #print $dataset," & ", $reso1," \\pm ",$resoerr1, " & ", $reso2,"\\pm ",$resoerr2, " & ",$reso3,"\\pm ",$resoerr3,"\\","\n";
    print "\\hline\n";
}
