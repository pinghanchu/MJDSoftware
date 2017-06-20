#!/usr/bin/perl
# get run list information based on 
# https://github.com/mppmu/GAT/blob/master/Apps/DataSetInfo.hh
my $inputfile = "list.txt";
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
my $i = 0;
my @runlist;
while(my $line = <$fin>) {
    chomp $line;
    my @array1 = split(" ",$line);
     push @runlist, $array1[0];
    $i = $i+1;
}
my @startrun;
my @endrun;
push @startrun, $runlist[0];
my $size = @runlist;
for(my $i = 1;$i<$size;$i++){
    if($runlist[$i]-$runlist[$i-1] >1){
	push @endrun, $runlist[$i-1];
	push @startrun, $runlist[$i];
    }
}
push @endrun, $runlist[$size-1];
my $size1= @startrun;
for(my $i=0;$i<$size1;$i++){
    print $startrun[$i]," ",$endrun[$i],"\n";
}
