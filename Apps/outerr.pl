#!/usr/bin/perl
system("ls -ls ./out.* > miss.txt");
open(my $fin, "<", "miss.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    #print $array[0],"\n";
    my @array1 = split('./out.',$array[9]);
    my $run = $array1[1];
    if($array[5]>600 || $array[5]<200){
	my $file = $run;
	my $app = "fill_".$file.".csh";
	system("rm out.$file err.$file\n");
	print $app," ",$array[5],"\n";
	system("chmod 755 $app");
	system("qsub -o out.$file -e err.$file $app")
    }
}
