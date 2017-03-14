#!/usr/bin/perl
open(my $fin, "<", "List/align.ezfit_18740_22282.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    my $file1 = "./Plot/cal/".$array[0]."_".$array[1]."/peak_trapENF_".$array[0]."_".$array[1]."_".$array[2]."_".$array[3]."_*.pdf";
    my $file2 = "./Plot/cal/".$array[0]."_".$array[1]."/gaus_trapENF_".$array[0]."_".$array[1]."_".$array[2]."_".$array[3]."_*.pdf";
    my $file3 = "./Plot/cal/".$array[0]."_".$array[1]."/calibrationdelta_trapENF_".$array[0]."_".$array[1]."_".$array[2]."_".$array[3].".pdf";
    system("cp $file1 ./tempezfit/");
    system("cp $file2 ./tempezfit/");
    system("cp $file3 ./tempezfit/");
}
