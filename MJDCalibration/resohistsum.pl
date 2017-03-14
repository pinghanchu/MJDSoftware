#!/usr/bin/perl
#system("./MkCookie");
#sleep(1);
open(my $fin, "<", "./List/runlist/cal.list0.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1],"\n";

    my $startrun = $array[0];
    my $endrun = $array[1];
    my $run = $startrun."_".$endrun;
    my $plotfile = "./*_".$run."_*.pdf";
    my $plotpath = "./Plot/cal/".$run."/";
    system("mkdir $plotpath");
    #system("rm $plotpath*.pdf");
    my $histfile = "./Hist/hist_".$run.".root";
    my $file = $run;
    my $z   = "resosum_".$file.".csh";
    open(my $fh, ">", $z) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";	
    print $fh "./resohistsum $startrun $endrun trapENFCal $histfile\n";	
    print $fh "mv $plotfile $plotpath\n";
    close $fh;
    system("chmod 755 $z");
    system("./$z");
    #system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $z");        
}
