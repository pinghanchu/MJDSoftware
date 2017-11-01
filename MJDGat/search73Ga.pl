#!/usr/bin/perl

print "Please input [dataset]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 1){
    print "You miss arguments\n";
}elsif( $numArgs >1){
    print "You have too many arguments\n";
}

my $dataset = $ARGV[0];
my $enr = 60;
my $time = 0.5;
#my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/git/WORK/MJDSkim/";
#my $inputfile = $scriptpath."List/runlist/DS".$dataset.".list.txt";
#my $inputfile = "./DS".$dataset.".list.txt";
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/git/MJDSoftware/MJDGat/";
my $search = $scriptpath."search73Ga";
my $inputfile = "run.list.txt";
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
my @runlist =();
while(my $line = <$fin>) {
    chomp $line;
    print $line,"\n";
    my @array = split(' ',$line);
    push(@runlist, $array[0]);
}


print "dataset = ", $dataset, "; energy = ", $enr, "; time = ", $time,"\n";

my $size = @runlist;
for(my $i = 0;$i<$size;$i++){
    my $run = $runlist[$i];
    my $file = $run;
    my $app  = "search_".$file.".csh";
    open(my $fh, ">", $app) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";
    print $fh "$search $run $enr $time\n";
    close $fh;
    system("chmod 755 $app");
    #system("./$app");
    system("sbatch slurm-job.sh './$app'"); 
}
