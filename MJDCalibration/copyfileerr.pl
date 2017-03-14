#!/usr/bin/perl
#my $runrange = "2361_7635";
my $runrange = "9034_14384";
#my $runrange = "16836_18351";
#my $runrange = "60000791_60001926";
#my $runrange = "18740_22282";
open(my $fin, "<", "List/offseterr.cal_$runrange.txt") or die "Failed to open file: $!\n";
#open(my $fin, "<", "List/cal/error.list.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    my $file1 = "./Plot/cal/".$array[0]."_".$array[1]."/peak_trapENF_".$array[0]."_".$array[1]."_".$array[2]."_".$array[3]."_*.pdf";
    my $file2 = "./Plot/cal/".$array[0]."_".$array[1]."/gaus_trapENF_".$array[0]."_".$array[1]."_".$array[2]."_".$array[3]."_*.pdf";
    my $file3 = "./Plot/cal/".$array[0]."_".$array[1]."/calibrationdelta_trapENF_".$array[0]."_".$array[1]."_".$array[2]."_".$array[3].".pdf";
    system("cp $file1 ./tempcal/");
    system("cp $file2 ./tempcal/");
    system("cp $file3 ./tempcal/");
}
