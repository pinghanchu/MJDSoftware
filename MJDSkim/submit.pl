#!/usr/bin/perl
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/git/MJDSoftware/MJDSkim/";
my $search = $scriptpath."savehglg";
for(my $i=1;$i<=5;$i++){
    system("$search $i");
}
