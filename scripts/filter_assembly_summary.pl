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
#
# Process command line
#
getopts('hi:o:vl:n:')||Usage();
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
  print ("  -n  : entries with version_status=latest, but not selected via option -o\n");
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
if (not defined($Getopt::Std::opt_o)){
  # Output goes to std output
  *OUT = *STDOUT;
} else {
  # Open file to write to
  open(OUT, ">$Getopt::Std::opt_o") || die ("can't open file $Getopt::Std::opt_o: $!");
}
if (defined($Getopt::Std::opt_l)){
    open(LOG,">$Getopt::Std::opt_l");
}
if (defined($Getopt::Std::opt_n)){
    open(NOT,">$Getopt::Std::opt_n");
}
if (defined($Getopt::Std::opt_v)){
    $verbose=1;
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
my @w=();
my %rec=();
# version_status in column 11
$rec{latest}=10;
# refseq_category in column 5
$rec{'representative genome'}=4;
$rec{'reference genome'}=4;
# assembly_level in column 12
$rec{'Complete Genome'}=11;

while (defined($_=<INP>)){
    if (/^#/){
	next;
    }
    chomp;
    @w=split(/\t/);

    # version_status
    if ($w[10] ne 'latest'){
	next;
    }

    # genome_rep
    if ($w[13] ne 'Full'){
	if (defined($Getopt::Std::opt_n)){
	    print NOT "$_\n";
	}
	next;
    }

    # assembly_level
    if ( !(($w[11] eq 'Complete Genome') || ($w[11] eq 'Chromosome')) ){
	if (defined($Getopt::Std::opt_n)){
	    print NOT "$_\n";
	}
	next;
    }

#    if (! (($w[4] eq 'representative genome') || ($w[4] eq 'reference genome')) ){
#	next;
#    }
    print OUT "$_\n"; 
}
