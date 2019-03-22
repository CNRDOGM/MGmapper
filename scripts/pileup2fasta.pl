#!/usr/bin/env perl
my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
    $command .= "$ARGV[$i] ";
    $i++;
}


## Program Info:
#
# Name: pileup2fasta
#
# Purpose: Intended to generate a fasta file, with minimum bias to a backbone, from a BAM file of a scaffolded
#          genome assembly.  It also generates the SNPs and INDELS from the reference -> test genome
#
# Author: John Nash
# Copyright (c) Public Health Agency of Canada, 2011,
#   all rights reserved.
#
# Licence: This script may be used freely as long as no fee is charged 
#   for use, and as long as the author/copyright attributions 
#   are not removed. It remains the property of the copyright holder.
#
# I learnt C in 1986 - I use curly brackets like an old fart. Bite me!
#
# History:
#
# 2011-02-23: v1.0 - Original...
# 2012-12-06: v1.1 - '-s' option fixed
# 2012-12-07: v1.2 - Handling mpileup files from bamfiles with multiple contigs in the reference file.
#                    (GFF output only)
# 2012-12-08: v1.3 - Fixed bug which tried to send data to a GFF file when none existed
# 2012-12-08: v1.4 -  Handling mpileup files from bamfiles with multiple contigs in the reference file.
#                    (FASTA output) 
# 2015-05-09  v2.0   Modified by Thomas Nordahl Petersen (tnp@cbs.dtu.dk). Major changes - only fasta output now.
##

## start of main function

use Getopt::Std;
use Cwd;
use strict;
use strict;
use Text::Wrap;
my $columns=70;
$Text::Wrap::columns = $columns + 1;
my $errorLines = 0;

# Init variables:
my $title = "pileup2fasta";
my $version = "2.0";
my $date = "09 December, 2012";
my $author = 'john.nash@phac-aspc.gc.ca';


# Default parameters
*LOG=*STDERR;
my $verbose=0;
# Minimum base coverage:
my $min_base_coverage = 4;

# Threshold fraction for accepting a SNP, insertion or deletion:
my $frac_limit = 0.50;

# minimum number of bases in output files
my $min_contig_len=30;

#
# Process command line
#
getopts('hi:o:vl:b:f:m:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i name] [-o name] [-b number] [-f number] [-m number][-v]\n");
  print ("Description:\n");
  print ("$0 - read a pileup file and make FASTA contig file\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input pileup file name [STDIN]\n");
  print ("  -o  : output contig file name [STDOUT]\n");
  print ("  -b  : minimum nucleotide depth to consider SNPs and INDELs [$min_base_coverage]\n");
  print ("  -f  : minimum read fraction required to support a SNP, insertion or deletion [$frac_limit]\n");
  print ("  -m  : minimum number of base pairs in output fasta file [$min_contig_len]\n");
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
  *INFILE = *STDIN;
} 
else{
  # Read from file
  if (($Getopt::Std::opt_i=~/\.gz$/) || ($Getopt::Std::opt_i=~/\.Z$/)){
    open(INFILE,"gunzip -c $Getopt::Std::opt_i |") || die ("can't open file $Getopt::Std::opt_i: $!");
  }
  else{
    open(INFILE,"<$Getopt::Std::opt_i") || die ("can't open file $Getopt::Std::opt_i: $!");
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
if (defined($Getopt::Std::opt_b)){
    $min_base_coverage = $Getopt::Std::opt_b;
}
if (defined($Getopt::Std::opt_f)){
    $frac_limit = $Getopt::Std::opt_f;
}
if (defined($Getopt::Std::opt_m)){
    $min_contig_len = $Getopt::Std::opt_m;
}

###############################################################################
# Main
#
###############################################################################
my $datestring = localtime();
my $thisDir=cwd();
if ($verbose){
    print LOG "$title version:$version working...\n" if ($verbose);
    print LOG "## Local date and time $datestring - Start program\n";
    print LOG "# $command\n";
    print LOG "# working dir: $thisDir\n\n";
}

# error messages:
my $error_msg = "Type \"$title -h\" for help.";
my $bug_msg = "I've never seen this in a samtools mpileup file, so I would love to see it in order to figure out how to deal with it. Please report to $author.";

#
# Run the main pileup parsing routing:
#
&process();
exit (0);
###########################################
#
# end of main program
#
###########################################


sub process {
# first field in an mpileup entry is the entry name:
  my $entry_name;
  my $old_entryName;

# Are we at the start of the file? 1 = YES 0 = NO
my $firstLine = 1;

# How many lines have we looked at (debug code):
  my $count;
  my $oldCount = 1;
#  my $lineCount = 0;

# final sequence:
  my $final_sequence;

# Where all of the entries are going:
  my @outgoing;

# Flag for the first line in a new entry
  my $first_entry_line=1;

# holder for the first reference position in an entry
  my $start;

# depth is the total number of bases within a contig, ave_depth is the average
  my $depth=0;
  my $ave_depth=0;

# Read the pileup file, line-by-line:
  while(<INFILE>) {

# Count the lines processed from the mpileup file:
    $count++;

# Fix up screwed up DOS/Windows CR/LF: 
   s/\r\n/\n/g;
    chomp;

# Check which entry is being processed, if the original bam file was referenced to a 
#  multiple fasta file:
    my @line = split /\t/;

# Because of the various ways that different OSes handle EOL,
#  I'm adding my own because I want to split the files downstream based on "\n":
    my $line = "$_\n";

# mpileup files (w/o consensus) have 6 fields per line:
#    die ("Not an mpileup consensus file: Error line below:\n$line") if scalar @line != 6;

    if (scalar @line != 6) {
#	print STDERR ("Error in mpileup consensus file:\nCommand: $command\nError line is shown below:\n$line");
	next;
    }
     $entry_name = $line[0];
 # Special case for first line:
    if ($firstLine) {

# This is where the first line of where the first entry_name appears:
      $old_entryName = $entry_name;
      $firstLine = 0;
    }

# Work on each entry in a multiple sequence entry, one-by-one:
    
    if ($first_entry_line){
	$start=$line[1];
	$first_entry_line=0;
    }
    if ($old_entryName eq $entry_name) {
      push @outgoing, $line;
    }
    else {
      $final_sequence =  &process_line(@outgoing);
      $start=$line[1];

# This is the first element of the new $entry_name. Zero it out.  Start afresh:
#      print STDERR "Working on: ", $entry_name, "\tLine: ", $count, "\n" if ($verbose);
      @outgoing = ();
      push @outgoing, $line;
      $oldCount = $count;
    } # end of last else
    
    $old_entryName = $entry_name;
  } # end of while(<INFILE>) 
  
# Deal with last scaffold...
  $final_sequence =  &process_line(@outgoing);

# Final status report:
#  print STDERR "\n", $count, " lines in ", $infile, " were processed.\n" if ($verbose);
  if ($verbose){
      print LOG "\n", $count, " lines in ", $Getopt::Std::opt_i, " were processed. Error lines : ",$errorLines," were found and skipped\n";
  }
  close(INFILE);

  if (defined($Getopt::Std::opt_o)){
      close(OUT);
  }
}  # end of sub process


sub process_line {

# first field in an mpileup entry is the entry name:
  my $entry_name;
  my $contig_counter=0;

# second field is the position:
  my $pos = 0;
  my $old_pos = 0;

# third field is the base:
  my $ref_sequence;

# fourth field is the number of reads:
  my $num_base_reads;

# fifth field is the encoded bases:
  my $variants;
  my $calc_var_qual = 0.0;
  my $from;
# The calculated sequnce:

  my $contig_sequence;

# sixth field is the encoded quality string. I could use it but have had no need to yet.
 
# @incoming is the array of all of the reads from each $entry_name:
  my @incoming = @_;
  my $number_of_N=0;
  my $contig_start;
  my $lines=0;
  my $depth=0;
  my $depth_min=99999;
  my $depth_max=-99999;
  my $ave_depth=0;
  my $SNPs=0;
  my $INDELs=0;
  my $insertions=0;
  my $deletions=0;
  my @chars = qw "a t c g \*";
  foreach (@incoming) {
    chomp;
    my $dumpbase_flag =0;
 
# Process each line:
    my @line = split /\t/;
    $lines++;
    if ($lines == 1){
	$contig_start = $line[1];
    }
# Assign the various mpileup entries to variables:
    $entry_name = $line[0];

# Positional checks in case the incoming mpileup file has deletions:
    $old_pos = $pos;
    $pos = $line[1];
    if ($old_pos == 0){
	$from = $pos;
    }
# To eliminate sequences which do not start at Position 1 being considered as deletions from the start, 
# i.e. for analysis of regions, we have to reset $old_pos for the first run:
    if (($old_pos == 0) && ($pos > 1)) {
      $old_pos = ($pos - 1);
    }

    $ref_sequence = $line[2];
    $num_base_reads = $line[3];

    if ($pos == ($old_pos+1)){
	$depth += $num_base_reads;
	if ($num_base_reads < $depth_min){
	    $depth_min =$num_base_reads;
	}
	if ($num_base_reads > $depth_max){
	    $depth_max =$num_base_reads;
	}
    }

    $variants = $line[4];

## Check position accuracy - part one:
    if ($pos == $old_pos) {
      die ("Error: unknown position error type 1. ", $bug_msg, "\n");
    }

## Check position accuracy - part two:
    if ($pos <= $old_pos) {
      die ("Error: unknown position error type 2. ", $bug_msg, "\n");
    }

# Find deletions and pad with * characters:
    if ($pos != ($old_pos + 1))  {
	my $len=length($contig_sequence);
	if ($len >= $min_contig_len){
	    my $end = $contig_start+$len -1;
	    $contig_counter++;
	    $ave_depth=$depth/$len;
	    my $contig_name = "$entry_name" . ".$contig_counter";
	    
	    #
	    # print fasta header
	    #
	    printf OUT ">$contig_name Start= $contig_start Len= $len depth= %.2f min_depth= %d max_depth= %d SNPs= %d INS= %d DEL= %d\n",$ave_depth,$depth_min,$depth_max,$SNPs,$insertions,$deletions;
	    
	    #
	    # print fasta sequence
	    #
	    print OUT wrap '', '', $contig_sequence, "\n";
	    if ($len % $columns == 0){
		print OUT "\n";
	    }
	}
	
	$contig_sequence='';
	$SNPs=0;
	$INDELs=0;
	

	$depth = $num_base_reads;
	$depth_min = $num_base_reads;
	$depth_max = $num_base_reads;

	$contig_start=$pos;
	my $deletion = $pos-1 - $old_pos + 1;
    } # end of if ($pos != ($old_pos + 1))
    
    my $sign;
    my $digit;
    my $match;

# If we have the desired base coverage, continue:
    if ($num_base_reads >= $min_base_coverage) {
    
# sub out the characters beginning with ^ as they are irrelevant:
      $variants =~ s/\^.{1}//g;
      $variants =~ s/\$//g;
      my $str = $variants;

# do not consider strings where there are no variants, i.e. all . or ,:
      if ($variants =~ /[.,]{$num_base_reads}/) {
	  $contig_sequence .= $ref_sequence;
	  next;
      }
      my $count = $str =~ s/([,.])/$1/gi;
      if ($count/$num_base_reads >= (1-$frac_limit)){
	  $contig_sequence .= $ref_sequence;
	  next;	  
      }
# Deal with variants:
      if ($variants =~ /([\+\-]*)(\d*)([GATC\*]+)/gi) {
	$sign = $1;
	$digit = $2;
	$match = $3;

# Process INSERTIONS and DELETIONS:
	if ($sign) {
	    my $var_count = 0;
	    $var_count++ while $variants =~ /[\+\-]\d[GATC\*]+/gi;
	    $calc_var_qual = $var_count / $num_base_reads;
	    #
	    # Threshold level for insertions:
	    #
	    if ($calc_var_qual > $frac_limit) {
		
		if ($sign eq '+'){
		    $contig_sequence .= lc($match);
		    $insertions++;
		    $dumpbase_flag=1;
		    if ($verbose){
			print LOG "Insertion accepted: insertion fraction: $calc_var_qual >= $frac_limit\tbase depth: $num_base_reads >= $min_base_coverage\n$_\n\n";
		    }
		}
	    }
	    if (($dumpbase_flag) || ($sign eq '-')){
		next;
	    }
	}
# Process SNPs and deletions:
	my $count;
	my $winner;
	my $most_counts=0;
	foreach my $id (@chars){
	    $count = $str =~ s/($id)/$1/gi;
	    if ($count > $most_counts){
		$most_counts=$count;
		$winner=uc $id;
	    }
	}
	if ($most_counts >= $min_base_coverage){
	    if ($winner eq '\*'){
		$calc_var_qual = $most_counts / $num_base_reads;
		if ($calc_var_qual > $frac_limit) {
		    $deletions++;
		    if ($verbose){
			print LOG "Deletion: DEL frac: $calc_var_qual > $frac_limit\tbase depth: $most_counts > $min_base_coverage\n$_\n\n";
		    }
		    $dumpbase_flag=1;
		}
	    }
	    else{
		$calc_var_qual = $most_counts / $num_base_reads;
		if ($calc_var_qual > $frac_limit) {
		    $SNPs++;			
		    $contig_sequence .= lc($winner);
		    if ($verbose){
			print LOG "SNP '$winner': SNP frac: $calc_var_qual > $frac_limit\tbase depth: $most_counts > $min_base_coverage\n$_\n\n";
		    }
		    $dumpbase_flag = 1;
		}
	    }
	}
	if ($dumpbase_flag == 0) {
	    $contig_sequence .= $ref_sequence;
	    $dumpbase_flag = 1;
	}
      }
    } 
    
    # Handle anything processed which failed... call it "reference":
    if ($dumpbase_flag == 0) {
      $contig_sequence .= $ref_sequence;
      $dumpbase_flag = 1;
    }
  } # end of foreach (@incoming)

# Return the differential sequence from the contig:
  
  my $len=length($contig_sequence);
  if ($len >= $min_contig_len){
      my $end = $contig_start+$len -1;
      $contig_counter++;
      $ave_depth=$depth/$len;
      my $contig_name = "$entry_name" . ".$contig_counter";
      printf OUT ">$contig_name Start= $contig_start Len= $len depth= %.2f min_depth= %d max_depth= %d SNPs= %d INS= %d DEL= %d\n",$ave_depth,$depth_min,$depth_max,$SNPs,$insertions,$deletions;
      print OUT wrap '', '', $contig_sequence, "\n";
      if ($len % $columns == 0){
	  print OUT "\n";
      }
      $depth=0;
  }
  return($contig_sequence);
} # end of sub process_lines
