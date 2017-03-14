#!/usr/bin/perl
open(my $fin, "<", "./List/runlist/cal.list0.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
    print $array[0]," ", $array[1], "\n";

    my $startrun = $array[0];
    my $endrun = $array[1];
    my $run = $startrun."_".$endrun;
    my $z   = "run_".$run.".csh";
    my $histfile = "./Hist/hist_".$run.".root";
    my $plotfile = "./*_".$run."*.pdf";
    my $plotpath = "./Plot/cal/".$run."/";    
    my $listfile = "./List/cal/*_".$run.".txt";
    #system("mkdir $plotpath");
    system("mv $plotfile $plotpath");
    #system("rm $plotpath*.pdf");
    #system("rm $listfile");
    open(my $fh, ">", $z) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";	
    print $fh "./calhistsum $startrun $endrun trapENFCal $histfile\n";	
    #print $fh "mv $plotfile $plotpath\n";
    close $fh;
    system("chmod 755 $z");
    #system("./$z");
    #system("qsub -l projectio=1 -cwd -o out.$run -e err.$run $z");    
}
