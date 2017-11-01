#!/usr/bin/perl
print "Please input [dataset]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 1){
    print "You miss arguments\n";
}elsif( $numArgs >1){
    print "You have too many arguments\n";
}

my $dataset = $ARGV[0];
my $inputfile = $scriptpath."List/runlist/DS".$dataset.".callist.txt";
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
my @startsubset =();
my @endsubset=();
while(my $line = <$fin>) {
    chomp $line;
    #print $line,"\n";
    #my @array = split(' ',$line);
    print $dataset, " ", $line,"\n";
}

