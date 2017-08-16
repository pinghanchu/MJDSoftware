#!/usr/bin/perl
my $version = "GAT-v01-06-125-gd9332b6";
my $skimpath = "/global/project/projectdirs/majorana/data/mjd/surfmjd/analysis/skim/";

for(my $i=0;$i<=5;$i++){
    my $path = $skimpath."DS".$i."cal/".$version."/";
    my $outfile = "./List/runlist/DS".$i."cal.list.txt";
    print $path,"\n";
    system("ls $path > $outfile");
}
