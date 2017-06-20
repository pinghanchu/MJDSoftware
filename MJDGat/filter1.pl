#!/usr/bin/perl
print "Please input [run][entry][channel][energy]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 4){
    print "You miss arguments\n";
}elsif( $numArgs >4){
    print "You have too many arguments\n";
}
my $run = $ARGV[0];
my $entry = $ARGV[1];
my $channel = $ARGV[2];
my $enr = $ARGV[3];
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDGat/";
my $filter = $scriptpath."filter1";
print $filter,"\n";
system("$filter $run $entry $channel $enr");
