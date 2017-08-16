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
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/git/WORK/MJDSkim/";
my $inputfile = $scriptpath."List/runlist/DS".$dataset.".list.txt";
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
my @subset =();
while(my $line = <$fin>) {
    chomp $line;
    print $line,"\n";
    my @array = split('skimDS',$line);
    my @array1 = split('.root',$array[1]);
    my @array2 = split('_',$array1[0]);
    print $array2[1],"\n";
    push(@subset, $array2[1]);
}


print "dataset = ", $dataset, "; energy = ", $enr, "; time = ", $time,"\n";
my $search =$scriptpath."search73Ga";

my $size = @subset;
for(my $i = 0;$i<$size;$i++){
    my $run = $subset[$i];
    my $file = $dataset."_".$run;
    my $app  = "search_".$file.".csh";
    open(my $fh, ">", $app) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";
    print $fh "$search $dataset $run 0 $enr $time\n";
    close $fh;
    system("chmod 755 $app");
    #system("./$app");
    #system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $app");
}
