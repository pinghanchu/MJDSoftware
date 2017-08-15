#!/usr/bin/perl
print "Please input [StartRun] [EndRun]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 2){
    print "You miss arguments\n";
}elsif( $numArgs >2){
    print "You have too many arguments\n";
}
my $startrun = $ARGV[0];
my $endrun =$ARGV[1];
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDGat/";
my $search = $scriptpath."searchpulserchannel";

for(my $i = $startrun;$i<=$endrun;$i++){
    system("$search $i");
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


