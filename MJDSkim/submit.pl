#!/usr/bin/perl
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/ana/WORK/MJDSkim/";
my $search = $scriptpath."search.pl";
my $enr = 67;
my $window = 5;
for(my $i=0;$i<=5;$i++){
    system("$search $i $enr $window");
}
