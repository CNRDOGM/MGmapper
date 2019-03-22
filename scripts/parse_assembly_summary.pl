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
my $mirrorPath;
my $accept;
my $downloadLog="wget.log";
#
# Process command line
#
getopts('hi:o:vl:m:w:t:')||Usage();
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
  print ("$0 - get representative genomes from assembly_summary.txt file\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input parse_assembly.txt file [STDIN]\n");
  print ("  -o  : output file name [STDOUT]\n");
  print ("  -m  : mirror path\n");
  print ("  -w  : wget logfile [$downloadLog]\n");
  print ("  -t  : output taxid and local filename\n");
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
if (defined($Getopt::Std::opt_v)){
    $verbose=1;
}
if (defined($Getopt::Std::opt_w)){
    $downloadLog=$Getopt::Std::opt_w;
}
if (defined($Getopt::Std::opt_t)){
    open(TaxidFilename,">$Getopt::Std::opt_t");
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
my @tmp=();
my $file;
my $name;
my $dir;
my $cmd;
my $url;
my $url_file;
my @w=();
my $taxid;
my $organism_name;
my $strain;
my $strain_field;
my $x;
while (defined($_=<INP>)){
    chomp;
    if (m/^#/){next}
    @w=split(/\t/);
    
    $taxid=$w[5];
    $organism_name=$w[7];
    $strain_field=$w[8];
    $x=length($strain_field);
    if ($x <1){
	$strain='unknown';
    }
    else{
	@tmp=split(/=/,$strain_field);
	$strain=$tmp[1];
    }
    $url = $w[19];
#    print STDERR "$url\n";
    $url_file='';

    @tmp=split(/\//,$url);
    $dir=$tmp[-1];
    $name = "$dir" ."_genomic.fna.gz";
    $file = "$Getopt::Std::opt_m" . "/$dir" ."_genomic.fna.gz";
    $cmd='';
    $url_file = $url . "/$dir" ."_genomic.fna.gz";
#    print STDERR "$url_file\n";
    $cmd = "wget -nH --cut-dirs=3";
    if (defined($Getopt::Std::opt_m)){
	$cmd .= " -P $Getopt::Std::opt_m";
    }
    $cmd .= " $url_file -a $downloadLog";
    if (! -e $file){
	print OUT "$cmd\n";	
    }
    else{
	my $checksum_remote = `curl $url/md5checksums.txt | grep $name | cut -f1 -d ' '`;
	chomp($checksum_remote);

	my $checksum_local = `md5sum $file | cut -f1 -d ' '`;
	chomp($checksum_local);
	if ($checksum_remote ne $checksum_local){
	    #
	    # remove old file and download the newer version
	    #
	    system("rm $file");
	    print OUT "$cmd\n";	
	    print LOG "updating file: $file\t$checksum_remote\t$checksum_local\n" if ($verbose);
	}
	else{
	    print LOG "Not updating file: $file\n" if ($verbose);
	}
    }
    #
    # make output file with taxid and filename
    #
    if (defined($Getopt::Std::opt_t)){
	print TaxidFilename "$taxid\t$file\t$organism_name\t$strain\n";
    }
}
if (defined($Getopt::Std::opt_t)){
    close(TaxidFilename);
}
