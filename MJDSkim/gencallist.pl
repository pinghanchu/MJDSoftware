#!/usr/bin/perl
#print "Please input [dataset]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 1){
    print "You miss arguments\n";
}elsif( $numArgs >1){
    print "You have too many arguments\n";
}
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/git/WORK/MJDSkim/";
my $dataset = $ARGV[0];
my $inputfile = $scriptpath."List/runlist/DS".$dataset."cal.list.txt";
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
my @subset =();
while(my $line = <$fin>) {
    chomp $line;
    #print $line,"\n";
    my @array = split('skimDS',$line);
    my @array1 = split('.root',$array[1]);
    my @array2 = split('_',$array1[0]);
    my @array3 = split('run',$array2[1]);
    #print $array3[1],"\n"; 
    push(@subset, int($array3[1]));
}
my @sort_subset = sort { $a <=> $b } @subset;
my $size = @subset;
my @startrun = ();
my @endrun = ();
my $starttemp =$sort_subset[0];
push (@startrun,$sort_subset[0]);
for(my $i=1;$i<$size;$i++){
    #print $sort_subset[$i]," ",$starttemp,"\n";    
    if(($sort_subset[$i]-$starttemp)>5){
	push(@startrun,$sort_subset[$i]);
	push(@endrun,$sort_subset[$i-1]);
    }
    $starttemp = $sort_subset[$i];
}
push(@endrun,$sort_subset[$size-1]);
$size = @startrun;
for(my $i=0;$i<$size;$i++){
    print $startrun[$i]," ",$endrun[$i],"\n";
}
