#!/usr/bin/perl
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDGat/";
my $home = "./";
my $totalcount = 0;
my $totaltime = 0;
for(my $i = 0;$i<6;$i++){
    my $pulserfiles = "./pulser.DS".$i.".txt";
    system("cat $pulserfiles > pulser.bk.txt");
    my $inputfile = "pulser.bk.txt";
    open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
    my $count = 0;
    
    while(my $line = <$fin>) {
	chomp $line;
	my @array = split(' ',$line);
	$count = $count+$array[1];
    }
    if($i!=2){
	$totalcount = $totalcount + $count;
	$totaltime = $totaltime + $count*20*1e-6;
	print $i," & ", $count," & ", $count*20*1e-6, "\\","\\ \n";
    }elsif($i==2){
        $totalcount = $totalcount + $count;
        $totaltime = $totaltime + $count*20*1e-6;
	print $i," & ", $count," & ", $count*20*1e-6, "\\","\\ \n";
    }
}
print " Total & ", $totalcount, " & ", $totaltime, "\\","\\ \n";
