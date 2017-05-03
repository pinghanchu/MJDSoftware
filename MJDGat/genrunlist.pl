#!/usr/bin/perl
# get run list information based on 
# https://github.com/mppmu/GAT/blob/master/Apps/DataSetInfo.hh
my $inputfile = "list.txt";
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
my $i = 0;
my @runlist;
while(my $line = <$fin>) {
    my $front = "{".$i.", {";
    my $end = "}}";
    chomp $line;
    my @array1 = split($front,$line);
    my @array2 = split($end,$array1[1]);
    my $text = $array2[0];
    $text =~ tr/,/   /;
    my @array3 = split(" ",$text);
    push @runlist, @array3;
    $i = $i+1;
}
my @list0;
my @list1;
my $i =0;
foreach $ii (@runlist){

    #print $ii,"\n";
    if($i%2==0){
	push @list0,$ii;
    }else{
	push @list1,$ii;
    }
    $i = $i+1;
}
my $size = @list0;
for(my $i=0;$i<$size;$i++){
    print $list0[$i]," ",$list1[$i],"\n";
}
