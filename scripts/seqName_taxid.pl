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
getopts('hi:o:vl:p:')||Usage();
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
  print ("  -p  : output file for sequences with the word plasmid in header\n");
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
if (defined($Getopt::Std::opt_p)){
    open(PLASMID,">$Getopt::Std::opt_p");
}
if (defined($Getopt::Std::opt_l)){
    open(LOG,">$Getopt::Std::opt_l");
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
my $taxid;
my $organism_name;
my $file;
my $line;
my $strain;
while (defined($_=<INP>)){
    chomp;
    @w=split(/\t/);
    $taxid=$w[0];
    $file=$w[1];
    $organism_name=$w[2];
    $strain=$w[3];
#    print LOG "$taxid\t$organism_name\t$file\n" if ($verbose);
    if (-e $file){
	open(FILE,"zgrep '^>' $file |");
	while (defined($line=<FILE>)){
	    chomp($line);
	    if ($line =~ m/^>(\S+)/){
		my $name=$1;
		print OUT "$name\t$taxid\t$strain\n";
		if (defined($Getopt::Std::opt_p)){
		    if ($line =~ m/plasmid/){
			print PLASMID "$name\t$taxid\t$strain\n";
		    }
		}
	    }
	}
	close(FILE);
    }
    else{
	print LOG "File not found: $file\n" if ($verbose);
    }
}
if (defined($Getopt::Std::opt_o)){
    close(OUT);
}
if (defined($Getopt::Std::opt_p)){
    close(PLASMID);
}
