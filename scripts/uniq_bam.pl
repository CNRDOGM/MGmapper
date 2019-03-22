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



getopts('hi:o:b:S')||Usage();


my $bamIn;
my $bamOut;
my $headerLine=1;
my @as=();
my @xs=();
my $buffer='';
my $buffer_lines=100000;
my $reads=0;

if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
    Usage();
}
sub Usage {
  print ("Usage: $0 [-h] [-i bam-file] [-o bam-file] [-b number] [-S]\n");
  print ("Description:\n");
  print ("$0 - extract uniqly mapped reads from a bam file, based on bam fields AS and XS\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input bam file\n");
  print ("  -o  : output bam file\n");
  print ("  -b  : lines in buffer before print to file [$buffer_lines]\n");
  print ("  -S  : Single-end reads [off]\n");
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
if (defined($Getopt::Std::opt_b)){
    $buffer_lines = $Getopt::Std::opt_b;
}

open(BamIn,"$prog_samtools view -h $bamIn |");
open(BamOut,"| $prog_samtools view -h -Sb - > $bamOut");

if (! defined($Getopt::Std::opt_S)){
    while (defined($_=<BamIn>)){
	if (((m/^\@SQ/) || (m/^\@PG/)) && ($headerLine)){
	    # this is a header line                                                                                                                          
	    print BamOut "$_";
	    next;
	}
	my @w=split(/\s+/);
	$xs[0] = substr($w[14],5);
	$as[0] = substr($w[13],5);
	
	my $line_mate = <BamIn>;
	my @w_mate=split(/\s+/,$line_mate);
	$xs[1] = substr($w_mate[14],5);
	$as[1] = substr($w_mate[13],5);
	
	my $uniqPair=0;
	if ($as[0] > $xs[0]){
	    $uniqPair=1;
	}
	if ($as[1] > $xs[1]){
	    $uniqPair=1;
	}
	
	if ($uniqPair){
	    $reads +=2;
	    if ($reads%$buffer_lines == 0){
		print BamOut "$buffer";
		$buffer='';
	    }
	    else{
		$buffer .= $_;
		$buffer .= $line_mate;
	    }
	}
    }
    close(BamIn);
    if (length($buffer) >0){
	print BamOut "$buffer";
    }
    close(BamOut);
}
else{
    while (defined($_=<BamIn>)){
	if (((m/^\@SQ/) || (m/^\@PG/)) && ($headerLine)){
	    # this is a header line                                                                                                                          
	    print BamOut "$_";
	    next;
	}
	my @w=split(/\s+/);
	$xs[0] = substr($w[14],5);
	$as[0] = substr($w[13],5);
	
	
	my $uniqPair=0;
	if ($as[0] > $xs[0]){
	    $uniqPair=1;
	}
	
	if ($uniqPair){
	    $reads++;
	    if ($reads%$buffer_lines == 0){
		print BamOut "$buffer";
		$buffer='';
	    }
	    else{
		$buffer .= $_;
	    }
	}
    }
    close(BamIn);
    if (length($buffer) >0){
	print BamOut "$buffer";
    }
    close(BamOut);
}

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
