#!/usr/bin/perl
for(my $i= 0; $i<20; $i++){
    system("qstat -u pchu > state.log");
    my $myfile = "state.log";
    my $size = -s $myfile;
    if($size == 0){
	    system("./submit.pl");
    }
    sleep(60*5);
}

