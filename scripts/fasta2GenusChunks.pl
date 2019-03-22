#!/usr/bin/env perl
my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
  $command .= "$ARGV[$i] ";
  $i++;
}
use Getopt::Std;
use strict;
use Cwd;
# Default parameters
*LOG=*STDERR;
my $verbose=0;
my $maxSize=10000000000;

#
# write fasta files to this temp dir - the directory will be removed af program finish
#
my $tmpDir = "/tmp/fasta2GenusChunks_$$";

#
# store sequences in memory until buffer_size is reached (1000000000 ~ 1Gb mem)
#
my $bufferSize=1000000000;
my $baseName="subset";
my $counterStart = 1;
#
# Process command line
#
getopts('hi:vl:m:o:b:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i name] [-o name] [-m number] [-b number] [-l logfile] [-v]\n");
  print ("Description:\n");
  print ("$0 - split a fasta file in chunks, where each chunk contains entries belonging to the same genus.\n");
  print ("Genus is read from fasta header i.e. the second field. If that is not present then a 'dummy' genus is made with these entries\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input file name\n");
  print ("  -o  : output base name [$baseName]\n");
  print ("  -m  : max output file size [$maxSize]\n");
  print ("  -b  : buffer size. printout buffer if length of a sequence is larger than [$bufferSize]\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -v  : Verbose [off]\n");
  print ("\n");
 exit;
} # Usage

#
# Open input
#
if (not defined($Getopt::Std::opt_i)){
  # Read from standard input
  *INP = *STDIN;
  print STDERR "Can not read from stdin\n";
  exit;
} 
else{
  # Read from file
  if (($Getopt::Std::opt_i=~/\.gz$/) || ($Getopt::Std::opt_i=~/\.Z$/)){
    open(INP,"gunzip -c $Getopt::Std::opt_i |") || die ("can't open file $Getopt::Std::opt_i: $!");
  }
  else{
    open(INP,"<$Getopt::Std::opt_i") || die ("can't open file $Getopt::Std::opt_i: $!");
  }
}

if (defined($Getopt::Std::opt_l)){
    open(LOG,">$Getopt::Std::opt_l");
}
if (defined($Getopt::Std::opt_v)){
    $verbose=1;
}
if (defined($Getopt::Std::opt_m)){
    $maxSize=$Getopt::Std::opt_m;
}
if (defined($Getopt::Std::opt_b)){
    $bufferSize=$Getopt::Std::opt_b;
}
if (defined($Getopt::Std::opt_o)){
    $baseName=$Getopt::Std::opt_o;
}
###############################################################################
# Main
#
###############################################################################
my $datestring = localtime();
print LOG "## Local date and time $datestring - start program\n" if ($verbose);
print LOG "# $command\n" if ($verbose);
my $thisDir=cwd();
print LOG "# working dir: $thisDir\n" if ($verbose);

my $line;
my $len=0;
my $fastaEntries=0;
my $fastaId;
my $genusId;
my $seq;
my $outFile;
if (! -d $tmpDir){
    system("mkdir $tmpDir");

    while (defined($_=<INP>)){
	$line=$_;
	chomp;
	$_ =~ s/\'//g;
	$_ =~ s/\[//g;
	$_ =~ s/\]//g;
	$_ =~ s/\,//g;
	$_ =~ s/\{//g;
	$_ =~ s/\}//g;
	$_ =~ s/\)//g;
	$_ =~ s/\(//g;
	$_ =~ s/\=//g;
	$_ =~ s/\;//g;
	$_ =~ s/\://g;
	
	if (/^>(\S+)/){
	    
	    if ($fastaEntries > 0){
		print GENUS "$seq";
		close(GENUS);
	    }
	    #
	    # print out info for sequence to tmp file
	    #
	    $fastaEntries++;
	    $len=0;
	    $fastaId=$1;
	    $seq='';
	    my @w=split(/\s+/);
	    if (/^>(\S+)\s(\S+)/){
		if (($2 eq 'Uncultured') && (defined($w[2]))){
		    $genusId=$w[2];
		}
		elsif (($2 eq 'PREDICTED') && (defined($w[2]))){
		    $genusId=$w[2];
		}
		else{
		    $genusId=$2;
		}
	    }
	    else{
		$genusId='Dummy';
	    }
	    $outFile="$tmpDir/$genusId.fasta";
	    open(GENUS,">>$outFile");
	}
	$seq .= $line;
	
	my $lenSeq = length($seq);
	if ($lenSeq > $bufferSize){
	    print LOG "Empty buffer as $lenSeq > $bufferSize\tfile= $outFile\n " if ($verbose); 
	    print GENUS "$seq";
	    $seq='';
	}
    }
    #
    # print the last fasta entry
    #
    print GENUS "$seq";
    close(GENUS);
}

my $genuses = `ls -1 $tmpDir | wc -l`;
chomp($genuses);
print LOG "Number of genuses: $genuses\n" if ($verbose);


my $i = $counterStart;

while ($genuses>0){
    my $out = "$baseName.$i";
    print LOG "Making file $out\n" if ($verbose);
    &makeChunks($out);
    $i++;
}

#
# Cleanup
#
system("rm -r $tmpDir");

$datestring = localtime();
print LOG "## Local date and time $datestring - end program\n" if ($verbose);

sub makeChunks{
    my ($outFile) = @_;
    open(LIST,"ls -ASl $tmpDir |");

    # Read the 'total size line'
    $_=<LIST>;
    
    my $sizeTotal=0;

    while (defined($_=<LIST>)){
	my @w=split(/\s+/);
	my $file = "$tmpDir/$w[-1]";
	my $size = $w[-5];
	if ($size >= $maxSize){
	    system("mv $file $outFile");
	    print LOG "Including $file in $outFile\n" if ($verbose);
	    $genuses--;
	    close(LIST);
	    return;
	}

	$sizeTotal += $size;
	if ($sizeTotal <= $maxSize){
	    print LOG "Including $file in $outFile\n" if ($verbose);
	    system("cat $file >> $outFile");
	    system("rm $file");
	    $genuses--;
	}
	else{
	    $sizeTotal -= $size
	}
    }
    return;
}
