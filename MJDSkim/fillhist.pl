#!/usr/bin/perl
print "Please input [dataset]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 1){
    print "You miss arguments\n";
}elsif( $numArgs >1){
    print "You have too many arguments\n";
}
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/git/WORK/MJDSkim/";
my $search = $scriptpath."fillhist";
print $search,"\n";
my $dataset = $ARGV[0];
my $inputfile = $scriptpath."List/runlist/DS".$dataset.".callist.txt";
open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
my @startsubset =();
my @endsubset=();
while(my $line = <$fin>) {
    chomp $line;
    print $line,"\n";
    my @array = split(' ',$line);
    push(@startsubset, $array[0]);
    push(@endsubset, $array[1]);
}


#my $inputfile = $scriptpath."List/runlist/DS".$dataset."cal.list.txt";
#open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
#my @subset =();
#while(my $line = <$fin>) {
#    chomp $line;
#    print $line,"\n";
#    my @array = split('skimDS',$line);
#    my @array1 = split('.root',$array[1]);
#    my @array2 = split('_',$array1[0]);
#    my @array3 = split('run',$array2[1]);
#    print $array3[1],"\n"; 
#    push(@subset, $array3[1]);
#}

my $size = @startsubset;
for(my $i=0;$i<$size;$i++){
    my $file = $dataset."_".$startsubset[$i]."_".$endsubset[$i];
    my $app  = "fill_".$file.".csh";
    open(my $fh, ">", $app) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";
    print $fh "$search $dataset $startsubset[$i] $endsubset[$i] 1\n";   
    close $fh;
    system("chmod 755 $app");
    system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $app");
}
