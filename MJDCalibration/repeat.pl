#!/usr/bin/perl
#This script re-submits "submit.pl" every 10 mins.
#
my $username = "pchu"; #Change to your username
for(my $i= 0; $i<20; $i++){
    system("qstat -u $username > state.log");
    my $myfile = "state.log";
    my $size = -s $myfile;
    if($size == 0){
	system("./submit.pl");
    }
    sleep(60*10);
}
