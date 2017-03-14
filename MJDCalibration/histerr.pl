#!/usr/bin/perl
my @pos1 = (11,12,13,14,15,16,17);
my @pos2 = (21,22,23,24,25,26,27);
my @pos3;
push(@pos3,@pos1);
push(@pos3,@pos2);
system("ls -ls out.* > List/miss.txt");
open(my $fin, "<", "./List/miss.txt") or die "Failed to open file: $!\n";
while(my $line = <$fin>) {
    chomp $line;
    my @array = split(' ',$line);
#    print $array[5]," ", $array[9],"\n";
#    if($array[5]<1000000){
    my @array2 = split('out.',$array[9]);
    my $run = $array2[1];

    #if($array[5] < 2157 && $run<9000){
    if($array[5]<1700 && ($run<18740 || $run>60000000)){

	#my @array1 = split('hist_',$array[9]);
	#my @array2 = split('.root',$array1[1]);
	#my @array3 = split('_',$array2[0]);
	my @pos;
	if($run<18590){
	    @pos = @pos1;
	}elsif($run>=18590 && $run<4500000){
	    @pos = @pos3;
	}elsif($run>60000000 && $run<65000000){
	    @pos = @pos2;
	}
	print $run,"\n";
	my $out = "out.".$run;
	my $temp = "out.temp";
	system("grep Save $out > $temp");
	open($fin1, "<", $temp);
	my @detpos = ();
	while(my $line1 = <$fin1>){
	    chomp $line1;
	    my @array3 = split(' ',$line1);
	    #print $array3[3],"\n";
	    push(@detpos,$array3[3]);
	}
	my @posmiss=();
	foreach $n1 (@pos){
	    my $count = 0;
	    foreach $n2 (@detpos){
		#print $n1*10," ",($n1+1)*10," ",$n2,"\n";
		#$count = 0;
		if($n2>=$n1*10 && $n2<($n1+1)*10){
		    $count++;
		    print $n2," ",$count,"\n";
		}
	    }
	    if($count<=10){
		push(@posmiss,$n1);
	    }
	}
	my $z = "fill_".$run.".csh";
	system("ls -ls $z >> test.log");
	open(my $fh, ">", $z) or die "cannot open";#
        print $fh "#!/bin/tcsh\n";
        #print $fh "./process_mjd_cal $builtfile\n";
        foreach $ip (@posmiss){
	    my $z = "fill_".$run."_".$ip.".csh";
	    #system("ls -ls $z >> test.log");
	    open(my $fh, ">", $z) or die "cannot open";#
	    print $fh "#!/bin/tcsh\n";
            #print $fh "./fillhist $run trapENF 800000 0 8000 300000 0 3000 $ip\n";
            #print $fh "./fillhist $run trapENFBL 1000 -10 10 1000 -10 10 $ip\n";
            print $fh "./fillhist $run trapENFCal 300000 0 3000 300000 0 3000 $ip\n";
	    #print $fh "./fillhist $run trapECal 300000 0 3000 300000 0 3000 $ip\n";
        #}
	    close $fh;
	    system("chmod 755 $z");
	    #print $z,"\n";
	    #system("./$z");
	    system("qsub -l projectio=1 -cwd -o out.$run -e err.$run $z");
	}
	#foreach $n3 (@posmiss){
	#    print $n3,"\n";
	#}
	#if($array3[1]<18740 && $array3[1]==0){
	#    my $run = $array2[0];
	#    my $z = "fill_".$run.".csh";
	    #system("rm hist_$run.root\n");
	#    print $z,"\n";
	    #system("chmod 755 *.csh");
	    #system("./$z");
	    #system("qsub -l projectio=1 -cwd -o out.$file -e err.$file $z")	
	#}
    }
}
