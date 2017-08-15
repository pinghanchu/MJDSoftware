#!/usr/bin/perl
print "Please input [Data Set]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 1){
    print "You miss arguments\n";
}elsif( $numArgs >1){
    print "You have too many arguments\n";
}
my $index = $ARGV[0];

my @pMod1 = (111,112,113,114,121,122,123,124,131,132,133,134,141,142,143,144,145,151,152,153,154,161,162,163,164,171,172,173,174);
my @pMod2 = (211,212,213,214,221,222,223,224,225,231,232,233,241,242,243,244,245,251,252,253,254,261,262,263,264,271,272,273,274);
my @p12;
push(@p12, @pMod1);
push(@p12, @pMod2);

my @pos;
if($index == 0 || $index == 1 || $index == 2 || $index == 3){
    @pos = @pMod1;
}elsif($index >=5){
    @pos = @p12;
}elsif($index == 4){
    @pos = @pMod2;
}
my @newpos = ();
foreach $ip (@pos){
    my $file1 = $ip."0";
    my $file2 = $ip."1";
    push(@newpos,$file1);
    push(@newpos,$file2);
}
my @pulsercount = ();
my @pulserchannel = ();
foreach $ip (@newpos){
    my $file = $index."_".$ip;
    my $inputfile = "pulser_".$file.".txt";
    open(my $fin, "<", $inputfile) or die "Failed to open file: $!\n";
    my $count = 0;
    my $channel = 0;
    while(my $line = <$fin>) {
	chomp $line;
	my @array = split(' ',$line);
	$channel  = $array[2];
	$count = $count+$array[3];
    }
    push(@pulserchannel,$channel);
    push(@pulsercount,$count);
}
my $size = @pulserchannel;
my $outputfile = "pulser_DS".$index.".txt";
open(my $fout,">>", $outputfile); 
for(my $i=0;$i<$size;$i++){    
    print $fout $index," ",$newpos[$i]," ", $pulserchannel[$i]," ",$pulsercount[$i],"\n";
}
