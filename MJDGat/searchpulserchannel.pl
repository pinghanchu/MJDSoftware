#!/usr/bin/perl
my $startrun = 21533;
my $endrun = 21564;
for(my $i = $startrun;$i<=$endrun;$i++){
    system("./searchpulserchannel $i");
}
my $inputfile = "./pulser_".$startrun."_".$endrun.".txt";
system("cat data_*.txt >> $inputfile");
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
my $outfile = "pulser.txt";
my $totalcount = 0;
my $totalevent = 0;
while(my $line = <$fin>) {
    chomp $line;
    my @array1 = split(" ",$line);
    my $count = $array1[2];
    $totalcount = $totalcount + $count;
    if($count>0){
	$totalevent = $totalevent + 1;
    }
}
print $totalcount," ", $totalevent,"\n";


