#!/usr/bin/env perl
BEGIN{
    if (exists $ENV{MGmap_HOME}) {
        push(@INC, $ENV{MGmap_HOME} . '/modules');
    } else {
        die "Environment variable MGmap_HOME does not exists. Make it point to MGmapper's home.\n";
    }
}



my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
  $command .= "$ARGV[$i] ";
  $i++;
}
use Getopt::Std;
use MGmapper::MGmapper;
use strict;


my $prog_samtools = &MGmapper::find_program('MGmap_SAMTOOLS', "samtools");
{
    my $st = `$prog_samtools 2>&1`;
    my @tmp = split(m/\n/, $st);
    @tmp = grep(m/^Version/, @tmp);
    push(@tmp, "Version: no version");
    $tmp[0] =~ m/^Version:\s+(\S+)/;
    die "Samtools $prog_samtools must be version 1.6\n" unless $1 eq '1.6';
}

*LOG=*STDERR;
my $bam;
#my $prog_samtools = "samtools";
#
# Process command line
#
getopts('hi:o:vl:a:b:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i name] [-o name] [-a indexFile] [-l logfile] [-v]\n");
  print ("Description:\n");
  print ("$0 - Read and write files\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input mpileup file name [STDIN]\n");
  print ("  -a  : bwa index file\n");
  print ("  -b  : bam file\n");
  print ("  -o  : output file name [STDOUT]\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -v  : Verbose [off]\n");
  print ("\n");
 exit;
} # Usage

#
# Open input
#
if (! defined($Getopt::Std::opt_a)){
    print "Option -a has not been defined\n";
    exit;
}
else{
    open(IDX,"<$Getopt::Std::opt_a") || die ("File not found: $Getopt::Std::opt_a");
}

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
if (defined($Getopt::Std::opt_b)){
    $bam = $Getopt::Std::opt_b;
    if (! -e $bam){
	print LOG "File not found: '$bam'\n";
    }
}
else{
    print "bam file not defined with option -b\n";
    die;
}
###############################################################################
# Main
#
###############################################################################
if (defined($Getopt::Std::opt_v)){
    print LOG ("# $command\n");
}
my %rec=();
#
# Read info in Index file
#
$_=<IDX>;
while (! eof (IDX)){
    $_=<IDX>;
    chomp;
    my $name;
    my $desc;
    my $size;
    if (m/^\d+ (\S+) (.+)/){
	$name=$1;
	$desc=$2;
	$rec{$name}{name}=$name;
	$rec{$name}{desc}=$desc;
    }
    elsif (m/^\d+ (\S+)/){
	$name=$1;
	$rec{$name}{name}=$name;
	$rec{$name}{desc}='NotDefined';
    }
    else{
	print LOG "Error parsing file: $Getopt::Std::opt_a\n";
    }
    $_=<IDX>;
    chomp;
    if (m/^\d+ (\d+) \d+/){
	$size=$1;
    }
    else{
	print LOG  "Error parsing File: $Getopt::Std::opt_a\n";
    }
    $rec{$name}{size}=$size;
}
close(IDX);

#
# Read the mpileup file
#
while (defined($_=<INP>)){
    chomp;
    if (m/^(\S+)\s+(\d+)\s+\w+\s+(\d+)\s+(\S+)\s+\S+/){
	my $name = $1;
	my $pos = $2;
	my $count = $3;
	my $line_ALT = $4;
	if ($count >= 1){
	    $rec{$name}{covered_positions}++;
	}
	$rec{$name}{nucleotides} += $count;
    }
}
#
# read the bam file and get the readcount and readCount where AS>XS i.e. uniquely mapped reads
#
open(BAM,"$prog_samtools view $bam |");
my $xs;
my $as;
my $nm;
my $name;
while (defined($_=<BAM>)){
    my @w=split(/\s+/);
    
    # target name
    $name=$w[2];

    # nedit distance i.e. number of changes to a read to make it identical to the reference sequence
    $nm = substr($w[11],5);

    # alignment score
    $as = substr($w[13],5);

    # secondary alignment score i.e. second best alignment score
    $xs = substr($w[14],5);
    $rec{$name}{reads}++;
    if ($as > $xs){
	$rec{$name}{uniqReadCount}++;
    }
    $rec{$name}{nm} += $nm;
}
close(BAM);

foreach my $key (keys %rec){
    if (exists ($rec{$key}{covered_positions})){
	my $reads = $rec{$key}{reads};
	my $uniqReadCount=$rec{$key}{uniqReadCount};
	my $covered_positions = $rec{$key}{covered_positions};

	# number of bp in reference sequence
	my $size=$rec{$key}{size};
	my $desc=$rec{$key}{desc};
	printf OUT "%s\t", $rec{$key}{name};
	printf OUT "%d\t", $size;
	printf OUT "%d\t", $rec{$key}{nucleotides};
	printf OUT "%d\t", $covered_positions;
	printf OUT "%d\t", $reads;
	printf OUT "%d\t", $uniqReadCount;
	printf OUT "%d\t", $rec{$key}{nm};
	printf OUT "%s\n",$desc;
    }
}
if (defined($Getopt::Std::opt_o)){
    close(OUT);
}
