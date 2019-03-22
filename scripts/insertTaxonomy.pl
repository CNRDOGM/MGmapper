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
use KyotoCabinet;
# Default parameters
*LOG=*STDERR;
my $verbose=0;
my $database="/home/databases/metagenomics/db/Taxonomy/update.kch";
my $safeMode=1;
#
# Process command line
#
getopts('hi:o:vl:k:s')||Usage();
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
  print ("$0 - insert sequence name and taxonomy in kch database\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input file name [STDIN]\n");
  print ("  -k  : kch database name [$database]\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -s  : Safe mode [on]\n");
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
if (defined($Getopt::Std::opt_v)){
    $verbose=1;
}
if (defined($Getopt::Std::opt_s)){
    $safeMode=0;
}
if (defined($Getopt::Std::opt_k)){
    $database=$Getopt::Std::opt_k;
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
 # create the database object
my $db = new KyotoCabinet::DB;
 
 # open the database
if (!$db->open("$database", $db->OWRITER | $db->OCREATE)) {
    printf STDERR ("open error: %s\n", $db->error);
}
my @w=();
my $value;
my $key;
my $j;
my $check;
my $n=0;
while (defined($_=<INP>)){
    if (m/^#/){next}
    chomp;
    @w=split(/\t/);
    $key=$w[0];

    $j=1;
    $value='';
    while (defined($w[$j])){
	$value .= "$w[$j]\t";
	$j++;
    }
    chomp($value);
    $n++;
    if ($safeMode){
	$check=$db->get($key);
	if (defined($check)){
	    if ($check ne $value){
		print LOG "existing\t$key\t$value\nnew\t$key\t$check\n\n" if ($verbose);
	    }
	}
	else{
	    $db->set($key,$value);
	    if ($n%1000000==0){
		$datestring = localtime();
		print LOG "Inserted $n keys\t$datestring\n" if ($verbose);
	    }
	}
    }
    else{
	$db->set($key,$value);
	if ($verbose){
	    if ($n%1000000==0){
		$datestring = localtime();
		print LOG "Inserted $n keys\t$datestring\n";
	    }
	}
    }
}
$db->close;
print LOG ("Inserted $n entries in $database\n"); 
