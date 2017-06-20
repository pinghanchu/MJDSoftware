#!/usr/bin/perl
print "Please input [Data Set]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 1){
    print "You miss arguments\n";
}elsif( $numArgs >3){
    print "You have too many arguments\n";
}
my $index = $ARGV[0];
my $file1 = "hist_".$index.".root";
my $file2 = "pulser_".$index.".root";
my $histname = "trapENFCal";
my $outputfile = "hist_".$index.".pdf";
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDGat/";
my $home = "./";
my $inputpath = $home;

my $plothisto = $scriptpath."plothisto";
system("$plothisto $file1 $file2 $histname");
system("mv hist.pdf $outputfile");
