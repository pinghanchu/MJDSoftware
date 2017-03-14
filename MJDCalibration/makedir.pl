#!/usr/bin/perl
open(my $fin, "<", "./cal.list.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    my $run = $array[0]."_".$array[1];
    my $file = "./Plot/cal/".$run."/";
    my $file1 = "*_".$run."_*.pdf";
    print $file,"\n";
    system("mkdir $file");
    system("mv $file1 $file");
}
