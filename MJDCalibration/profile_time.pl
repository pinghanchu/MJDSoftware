#!/usr/bin/perl
open(my $fin, "<", "./List/runlist/cal.long.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    my $startrun = $array[0];
    my $endrun = $array[1];
    my $run = $startrun."_".$endrun;
    my $z   = "profile_".$run.".csh";

    open(my $fh, ">", $z) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";
    print $fh "./profile_time $startrun $endrun 2614 10\n";
    close $fh;
    system("chmod 755 $z");
    #system("./$z");
    system("qsub -l projectio=1 -cwd -o out.$run -e err.$run $z");
}
