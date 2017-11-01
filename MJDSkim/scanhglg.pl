#!/usr/bin/perl
my $scriptpath = "/global/projecta/projectdirs/majorana/users/pchu/git/MJDSoftware/MJDSkim/";
my $search = $scriptpath."scanhglg";
for(my $i=0;$i<40;$i++){
    my $low = 20370+$i;
    my $up = 20370+($i+1);
    my $file  = $low."_".$up;
    my $app = "scan_".$file.".csh";
    print $app,"\n";
    open(my $fh, ">", $app) or die "cannot open";#
    print $fh "#!/bin/tcsh\n";
    print $fh "$search $low $up\n";
    close $fh;
    system("chmod 755 $app");
    system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $app");
}
