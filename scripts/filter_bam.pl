#!/usr/bin/env perl
BEGIN{
    if (exists $ENV{MGmap_HOME}) {
        push(@INC, $ENV{MGmap_HOME} . '/modules');
    } else {
        die "Environment variable MGmap_HOME does not exists. Make it point to MGmapper's home.\n";
    }
}

use Getopt::Std;
use strict;
use MGmapper::MGmapper;

my $prog_samtools = &MGmapper::find_program('MGmap_SAMTOOLS', "samtools");
{
    my $st = `$prog_samtools 2>&1`;
    my @tmp = split(m/\n/, $st);
    @tmp = grep(m/^Version/, @tmp);
    push(@tmp, "Version: no version");
    $tmp[0] =~ m/^Version:\s+(\S+)/;
    die "Samtools $prog_samtools must be version 1.6\n" unless $1 eq '1.6';
}

getopts('hi:o:f:a:c:SN:')||Usage();

my $bamIn;
my $bamOut;
my $FMM;
my $MAS;
my $AMM;
my $NM;
my $max_edit_distance;
my $FMM_flag=0;
my $MAS_flag=0;
my $AMM_flag=0;
my $NM_flag=0;
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
    Usage();
}
sub Usage {
  print ("Usage: $0 [-h] [-i bam-file] [-o bam-file] [-f number] [-a number] [-c number] [S] [N number]\n");
  print ("Description:\n");
  print ("$0 - filter a bam file for alignment score and fraction or absolute number of matches+mismatches in cigar field\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input bam file\n");
  print ("  -o  : output bam file\n");
  print ("  -f  : FMM value i.e. fraction of matches+mismatches\n");
  print ("  -a  : minimum alignment score\n");
  print ("  -c  : absolute count of matches+mismatches\n");
  print ("  -N  : maximum edit distance\n");
  print ("  -S  : single-end bam file [off]\n");
  exit;
}

if (not defined($Getopt::Std::opt_i)){
    print STDERR "Can not read from stdin\n";
    die;
}
else{
    $bamIn = $Getopt::Std::opt_i;
}
if (not defined($Getopt::Std::opt_o)){
    print STDERR "Can not write to stdout\n";
    die;
}
else{
    $bamOut = $Getopt::Std::opt_o;
}
if (defined($Getopt::Std::opt_f)){
    $FMM=$Getopt::Std::opt_f;
    $FMM_flag=1;
}

if (defined($Getopt::Std::opt_a)){
    $MAS=$Getopt::Std::opt_a;
    $MAS_flag=1;
}

if (defined($Getopt::Std::opt_c)){
    $AMM=$Getopt::Std::opt_c;
    $AMM_flag=1;
}

if (defined($Getopt::Std::opt_N)){
    $max_edit_distance=$Getopt::Std::opt_N;
    $NM_flag=1;
}


my $cmd = "$prog_samtools view -h $bamIn |";
open(FILE,"$cmd");

open(YES,"| $prog_samtools view -Sb - > $bamOut");

my %rec=();

my $line_read;
my $line_mate;
my $reads=0;
my $headerLine=1;
my @matchRatios=();
my @matchCounts=();
my @alignmentScoreFlag=();
my @editDistanceFlag=();
my @as=();
my @xs=();
my @nm=();

my $buffer='';
my $buffer_lines=100000;
my $readsOkCount=0;

if (! defined($Getopt::Std::opt_S)){
    while (defined($_=<FILE>)){

	if (((m/^\@SQ/) || (m/^\@PG/)) && ($headerLine)){
	    # this is a header line
	    print YES "$_";
	    next;
	}
	$headerLine=0;
	$reads++;
	$line_read=$_;
	chomp($line_read);
	my @w_read=split(/\s+/,$line_read);
	
	my $len=length($w_read[9]);
	my $cigar = $w_read[5];
	
	$xs[0] = substr($w_read[14],5);
	$as[0] = substr($w_read[13],5);
	$nm[0] = substr($w_read[11],5);
	#
	# Now count number of matches in cigar string
	#
	my $matchSum=0;
	my @tmp = $cigar =~ m/(\d+)M/g;
	foreach my $s (@tmp) {
	    $matchSum += $s;
	}
	
	$matchRatios[0] = $matchSum/$len;
	$matchCounts[0] = $matchSum;
	#
	# calc the same values as above for the mate
	#
	$line_mate=<FILE>;
	chomp($line_mate);
	$reads++;
	my @w_mate=split(/\s+/,$line_mate);
	
	$len=length($w_mate[9]);
	$cigar = $w_mate[5];
	$xs[1] = substr($w_mate[14],5);
	$as[1] = substr($w_mate[13],5);
	$nm[1] = substr($w_mate[11],5);
	
	$matchSum=0;
	@tmp = $cigar =~ m/(\d+)M/g;
	foreach my $s (@tmp) {
	    $matchSum += $s;
	}
	
	$matchRatios[1] = $matchSum/$len;
	$matchCounts[1] = $matchSum;
	
	#
	# Init flags as true
	#
	my @matchRatioFlags=();
	$matchRatioFlags[0]=1;
	$matchRatioFlags[1]=1;
	
	my @matchCountFlags=();
	$matchCountFlags[0]=1;
	$matchCountFlags[1]=1;
	my @matchCountFlags=();
	$matchCountFlags[0]=1;
	$matchCountFlags[1]=1;
	
	$alignmentScoreFlag[0]=1;
	$alignmentScoreFlag[1]=1;
	
	$editDistanceFlag[0]=1;
	$editDistanceFlag[1]=1;
	for (my $i=0;$i<=1;$i++){
	    if ($FMM_flag){
		if ($matchRatios[$i] < $FMM){
		    $matchRatioFlags[$i] = 0;
		}
	    }
	    
	    if ($AMM_flag){
		if ($matchCounts[$i] < $AMM){
		    $matchCountFlags[$i] = 0;
		}
	    }
	    
	    if ($MAS_flag){
		if ($as[$i] < $MAS){
		    $alignmentScoreFlag[$i]=0;
		}
	    }

	    if ($NM_flag){
		if ($nm[$i] > $max_edit_distance){
		    $editDistanceFlag[$i]=0;
		}
	    }
	}
	
	my $readsOk = 1;
	for (my $i=0;$i<=1;$i++){
	    $readsOk *= $matchRatioFlags[$i] * $matchCountFlags[$i] * $alignmentScoreFlag[$i] * $editDistanceFlag[$i];
	}
	
	if ($readsOk){
	    $readsOkCount+=2;
	    $buffer .= "$line_read\n";
	    $buffer .= "$line_mate\n";
	    
	    if ($readsOkCount%$buffer_lines == 0){
		print YES "$buffer";
		$buffer='';
	    }
	}
    }
    close(FILE);
}
else{
    while (defined($_=<FILE>)){

	if (((m/^\@SQ/) || (m/^\@PG/)) && ($headerLine)){
	    # this is a header line
	    print YES "$_";
	    next;
	}
	$headerLine=0;
	$reads++;
	$line_read=$_;
	chomp($line_read);
	my @w_read=split(/\s+/,$line_read);
	
	my $len=length($w_read[9]);
	my $cigar = $w_read[5];
	
	$xs[0] = substr($w_read[14],5);
	$as[0] = substr($w_read[13],5);
	$nm[0] = substr($w_read[11],5);
	#
	# Now count number of matches in cigar string
	#
	my $matchSum=0;
	my @tmp = $cigar =~ m/(\d+)M/g;
	foreach my $s (@tmp) {
	    $matchSum += $s;
	}
	
	$matchRatios[0] = $matchSum/$len;
	$matchCounts[0] = $matchSum;
	#
	# Init flags as true
	#
	my @matchRatioFlags=();
	$matchRatioFlags[0]=1;
	
	my @matchCountFlags=();
	$matchCountFlags[0]=1;

	my @matchCountFlags=();
	$matchCountFlags[0]=1;

	
	$alignmentScoreFlag[0]=1;
	$editDistanceFlag[0]=1;

	
	for (my $i=0;$i<=0;$i++){
	    if ($FMM_flag){
		if ($matchRatios[$i] < $FMM){
		    $matchRatioFlags[$i] = 0;
		}
	    }
	    
	    if ($AMM_flag){
		if ($matchCounts[$i] < $AMM){
		    $matchCountFlags[$i] = 0;
		}
	    }
	    
	    if ($MAS_flag){
		if ($as[$i] < $MAS){
		    $alignmentScoreFlag[$i]=0;
		}
	    }

	    if ($NM_flag){
		if ($nm[$i] > $max_edit_distance){
		    $editDistanceFlag[$i]=0;
		}
	    }
	}
	
	my $readsOk = 1;
	for (my $i=0;$i<=0;$i++){
	    $readsOk *= $matchRatioFlags[$i] * $matchCountFlags[$i] * $alignmentScoreFlag[$i] * $editDistanceFlag[$i];
	}
	
	if ($readsOk){
	    $readsOkCount++;
	    $buffer .= "$line_read\n";

	    if ($readsOkCount%$buffer_lines == 0){
		print YES "$buffer";
		$buffer='';
	    }
	}
    }
    close(FILE);
}
if (length($buffer) >0){
    print YES "$buffer";
}
close(YES);

sub check_program{
    my ($program, @possible_path) = @_;

    foreach my $id (@possible_path){
        if (-e "$id/$program"){
            return("$id/$program");
        }
    }
    print STDERR "can not find program: '$program'\n";
    print STDERR "checked these locations: '@possible_path'\n";
    die;
}
