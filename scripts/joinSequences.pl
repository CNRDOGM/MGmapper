#!/usr/bin/env perl
my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
  $command .= "$ARGV[$i] ";
  $i++;
}
use Getopt::Std;
use Cwd;
use strict;
# Default parameters
*LOG=*STDERR;
my $verbose=0;
my %discard_sequence_from_header_words=();
my $outFile;
#
# Process command line
#
getopts('hi:o:vl:d:D:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i name] [-o name] \n");
  print ("Description:\n");
  print ("$0 - Read and write files\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input file name [STDIN]\n");
  print ("  -o  : output file name [STDOUT]\n");
  print ("  -d  : discard sequences based on words in header\n");
  print ("  -D  : output fasta file with discarded sequences\n");
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
#
# If not file name is given, use standard output
#
# Open file to write to
if (defined($Getopt::Std::opt_o)){
    open(OUT,">$Getopt::Std::opt_o");
}
else{
    *OUT=*STDOUT;
}
if (defined($Getopt::Std::opt_l)){
    open(LOG,">$Getopt::Std::opt_l");
}
if (defined($Getopt::Std::opt_v)){
    $verbose=1;
}
if (defined($Getopt::Std::opt_d)){
    my @w = split(' ',$Getopt::Std::opt_d);
    my $i=0;
    while (defined($w[$i])){
	$discard_sequence_from_header_words{$w[$i]}=1;
	$i++;
    }
}
if (defined($Getopt::Std::opt_D)){
    open(DISCARD,">$Getopt::Std::opt_D");
}
###############################################################################
# Main
#
###############################################################################
my $datestring = localtime();
my $thisDir=cwd();
if ($verbose){
    print LOG "## Local date and time $datestring - Start program\n";
    print LOG "# $command\n";
    print LOG "# working dir: $thisDir\n\n";
}
if ($verbose){
    foreach my $key (keys %discard_sequence_from_header_words){
	print LOG "Discard sequences with this word in fasta header: $key\n"
    }
}

my $printSequence=1;
my $file;
my $header;
my $n=0;
while (defined($file=<INP>)){
    chomp($file);
    
    #
    # open the file
    #
    if ($file=~/\.gz$/){
	open(FILE,"gunzip -c $file |");
    }
    else{
	open(FILE,"<$file");
    }
    
    
    #
    # parse the fasta file
    #
    if (! defined($Getopt::Std::opt_d)){
	while (defined($_=<FILE>)){
	    print OUT "$_";
	}
	close(FILE);
    }
    else{
	while (defined($_=<FILE>)){
	    chomp;
	    if (m/^>(\S+)/){
		my $seqName=$1;
		$header=$_;
		$printSequence=1;
		foreach my $key (keys %discard_sequence_from_header_words){
		    if (($header =~ $key) && ($printSequence == 1)){
			$printSequence=0;
			print LOG "found '$key' : Discarding entry $seqName in file: $file\n" if ($verbose);
			$n++;
		    }
		}
	    }
	    if ($printSequence){
		print OUT "$_\n";
	    }
	    else{
		if (defined($Getopt::Std::opt_D)){
		    print DISCARD "$_\n";
		}
	    }
	}
	close(FILE);
    }
}
if (defined($Getopt::Std::opt_D)){
    close(DISCARD);
}
print LOG "Discarded $n sequences\n" if ($verbose);
