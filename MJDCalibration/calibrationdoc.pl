#!/usr/bin/perl
#This script generates a tex file, "plot.tex", including the plots of spectrum and the deviation of calibration.
#
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
my @pos1 = (112,113,114,121,122,123,132,133,134,141,142,143,144,145,152,153,161,163,164,171,172,173,174);
my @pos2 = (211,212,213,214,221,222,223,231,232,241,242,244,251,253,254,261,262,273,274);
my @pos3;
push(@pos3,@pos1);
push(@pos3,@pos2);
print "Please input [startrun] [endrun] [energy name] [input file path]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 4){
    print "You miss arguments\n";
}elsif( $numArgs >4){
    print "You have too many arguments\n";
}
my $startrun = $ARGV[0];
my $endrun = $ARGV[1];
my $energy = $ARGV[2];
my $plotpath = $ARGV[3];
my $dir = $startrun."_".$endrun;
print "startrun = ", $startrun,"; endrun = ", $endrun, "\n";
print "Plot path = ", $plotpath,"\n";

my @pos;
if($endrun<18590){
    @pos = @pos1;
}elsif($startrun>=18590 && $endrun<4500000){
    @pos = @pos3;
}elsif($startrun>60000000 && $endrun<65000000){
    @pos = @pos2;
}
my $plotfile = $plotpath."delta_".$dir."_".$energy."_*.pdf";
system("rm plot.txt plot.tex");
system("ls $plotfile > plot.txt");
open(my $fin, "<", "./plot.txt") or die "Failed to open file: $!\n";
my $z = "plot.tex";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split('_',$line);
    my $energy = $array[4];
    my $pos = $array[5];
    my $channel = $array[6];
    my @temp = split('.pdf',$channel);
    my $channel = $temp[0];
    my $file = $dir."_".$energy."_".$pos."_".$channel;
    my $delta = $plotpath."delta_".$file.".pdf";
    my $spectrum0 = $plotpath."spectrum_".$file."_0.pdf";
    my $spectrum1 = $plotpath."spectrum_".$file."_1.pdf";
    my $spectrum2 = $plotpath."spectrum_".$file."_2.pdf";
    my $spectrum3 = $plotpath."spectrum_".$file."_3.pdf";

    open(my $fh, ">>", $z) or die "cannot open";
    print $fh "\\begin\{figure\}\[hb\]\n";
    print $fh "\\centering\n";
    print $fh "\\includegraphics\[width=0.45\\textwidth, height=0.2\\textheight\]\{$spectrum0\}\n";
    print $fh "\\includegraphics\[width=0.45\\textwidth, height=0.2\\textheight\]\{$spectrum1\}\n";
    print $fh "\\includegraphics\[width=0.45\\textwidth, height=0.2\\textheight\]\{$spectrum2\}\n";
    print $fh "\\includegraphics\[width=0.45\\textwidth, height=0.2\\textheight\]\{$spectrum3\}\n";
    print $fh "\\includegraphics\[width=0.8\\textwidth, height=0.4\\textheight\]\{$delta\}\n";
    print $fh "\\caption\{ Detector $pos, Channel $channel of run $startrun - $endrun.\}\n";
    print $fh "\\label\{fig:$file}\n";
    print $fh "\\end\{figure\}\n";
    print $fh "\\clearpage\n";
}
