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
my $taxDir = "/home/databases/metagenomics/mirror/taxonomy";
my @types = qw(fasta gi taxid scientific_name);
my $type='fasta';
#
# Process command line
#
getopts('hi:o:vl:p:t:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i name] [-o name] [-p name] [-t name]\n");
  print ("Description:\n");
  print ("$0 - Read and write files\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input file name [STDIN]\n");
  print ("  -o  : output file name [STDOUT]\n");
  print ("  -p  : taxonomy path [$taxDir]\n");
  print ("  -t  : input data type: 'fasta' or 'gi' or 'taxid' or 'scientific_name' \(default: $type\)\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -v  : Verbose [off]\n");
  print ("\n");
  print ("      type='fasta'\t'>gi\|213123\|' (a file with fasta entries)\n");
  print ("      type='gi'\t\t213123 (a file with 1 column of gi-numbers)\n");
  print ("      type='taxid'\tsomeIdentifies\ttaxid\n");
  print ("      type='scientific_name'\tsomeIdentifies\tscientific name\n");
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
if (defined($Getopt::Std::opt_p)){
    $taxDir=$Getopt::Std::opt_p;
}
if (defined($Getopt::Std::opt_t)){
    $type= $Getopt::Std::opt_t;
    my $found=0;
    foreach my $id (@types) {
	if ($id eq $type){
	    $found = 1;
	}
    }
    if ($found == 0){
	print STDERR "Unknown input type: '$Getopt::Std::opt_t'\n";
	print STDERR "Acceptable options: @types\n";
	exit;
    }
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
my $names="$taxDir/names.dmp";
my $nodes="$taxDir/nodes.dmp";
my $gitaxid="$taxDir/gi_taxid_nucl.dmp";

my %input=();
my %input_tax=();
my %input_strain=();
my %input_scientific=();
my $n=0;
my $Qname;
if ($type eq 'fasta'){
    #
    # Fasta format
    #    
    if (defined($Getopt::Std::opt_i)){
	print LOG "Reading fasta file: $Getopt::Std::opt_i\n" if ($verbose);
    }
    else{
	print LOG "Reading fasta file: stdin\n" if ($verbose);
    }
   while (defined($_=<INP>)){
	chomp;
	if (m/>(\S+)/){
	    $Qname=$1;
	    if (m/>gi\|(\d+)\|/){
		my @w=split(/\s+/);
		my $gi=$1;
		$w[0] =~ s/\>//g;
		$input{$gi}=$w[0];
#		print LOG "$Qname\t$gi\n";
		$n++;
	    }
	}
    }
}
elsif ($type eq 'gi'){
    #
    # gi numbers
    #
    while (defined($_=<INP>)){
	chomp;
	$input{$_}=$_;
	$n++;
    }
}
elsif ($type eq 'scientific_name'){
    #
    # scientific or synonym names
    #
    while (defined($_=<INP>)){
	chomp;
	my @w=split(/\t/);
	$w[0] =~ s/\s+$//;
	$w[1] =~ s/\s+$//;
	my $uc=uc($w[1]);
	$input_scientific{$w[0]}=$uc;
	$n++;
    }
}
else{
    #
    # taxid numbers
    #
    while (defined($_=<INP>)){
	chomp;
	my @w=split(/\t/);
	$input_tax{$w[0]}=$w[1];
	$input{$w[0]}=$w[1];
	$input_strain{$w[0]}=$w[2];
	$n++;
    }
}
print LOG "Input identifiers: $n\n" if ($verbose);

#
# reading taxid and parent
#
print LOG "Reading file: $nodes\n" if ($verbose);
open(NODES,"<$nodes") || die ("File not found: '$nodes'\n");
my %parent=();
my %rank=();
while (defined($_=<NODES>)){
    chomp;
    my @w=split(/\t/,$_);
    my $tax=$w[0];
    my $parentID=$w[2];
    my $rankID=$w[4];
    $parent{$tax} = $parentID;
    $rank{$tax}=$rankID;
}
close(NODES);
#
# reading taxonomy name given a taxid
#
print LOG "Reading file: $names\n" if ($verbose);
open(NAMES,"<$names") || die ("File not found: '$names'\n");
my %names=();
my %scientific=();
my %synonym=();
while (defined($_=<NAMES>)){
    chomp;
    my @w=split(/\t/,$_);
    my $name_class=$w[6];
    my $name_txt;
    my $tax_id;
    if ($type eq 'scientific_name'){
	if ($name_class eq 'scientific name' ||
	    $name_class eq 'equivalent name'){
	    $tax_id=$w[0];
	    $name_txt=$w[2];
	    $names{$tax_id}=$name_txt;
	    my $uc = uc($name_txt);
	    $scientific{$uc}=$tax_id;
	}
	elsif ($name_class =~ m/synonym/ || $name_class eq 'misspelling'){
	    $tax_id=$w[0];
	    $name_txt=$w[2];
	    my $uc = uc($name_txt);
	    $synonym{$uc}=$tax_id;
	}
    }
    else{
	if ($name_class ne 'scientific name'){
	    next;
	}
	my $tax_id=$w[0];
	my $name_txt=$w[2];
	$names{$tax_id}=$name_txt;
    }
}
close(NAMES);


my $continue=$n;
my $giFound=0;
if (($type eq 'fasta') || ($type eq 'gi')){
    #
    # reading gi_tax file 
    #
    print LOG "Reading file: $gitaxid\n" if ($verbose);
    open(GITAXID,"<$gitaxid") || die ("File not found: '$gitaxid'\n");
    my $n=0;
    while ((defined($_=<GITAXID>)) && ($continue)){
	$n++;
	if (m/(\d+)\t(\d+)/){
	    if (exists($input{$1})){
		$Qname=$input{$1};
		$input_tax{$Qname}=$2;
		$continue--;
		$giFound++;
	    }
	}
	if ($n%10000000 == 0){
	    print LOG "read $n lines, found $giFound gi-taxid relations\n" if ($verbose);
	}
    }
    close(GITAXID);
    print LOG "gi-taxid relations found: $giFound\n" if ($verbose);
}
elsif ($type eq 'scientific_name'){
    foreach my $key (keys %input_scientific){
	my $taxon = $input_scientific{$key};
	my @sp = split(/\s+/,$taxon);
	if (exists($scientific{$taxon})){
	    $input_tax{$key} = $scientific{$taxon};
	}
	elsif (exists($synonym{$taxon})){
	    $input_tax{$key} = $synonym{$taxon};
	}
	elsif ($#sp > 0 && exists($scientific{$sp[0]})){
	    $input_tax{$key} = $scientific{$sp[0]};
	}
	elsif ($#sp > 0 && exists($synonym{$sp[0]})){
	    $input_tax{$key} = $synonym{$sp[0]};
	}
	else{
	    print LOG "Unknown Scientific/synonym name '$taxon' (has been upper cased) for seqId '$key'\n" if ($verbose);
	    $input_tax{$key}=0;
	}
    }
}
else{
    #
    # storing relation between input identifier and taxid
    #
    foreach my $key (keys %input){
	$input_tax{$key}=$input{$key};
    }    
}

print LOG "Establish taxonomy path\n" if ($verbose);
print OUT "# Qname\t";
print OUT "superkingdom\tsuperkingdomTax\t";
print OUT "phylum\tphylumTax\t";
print OUT "class\tclassTax\t";
print OUT "order\torderTax\t";
print OUT "family\tfamilyTax\t";
print OUT "genus\tgenusTax\t";
print OUT "species\tspeciesTax\t";
print OUT "strain\tstrainTax\n";

foreach my $key (keys %input_tax){
    my $Qname = $key;
    my $tax = $input_tax{$key};
    my $taxID = $tax;

    my $superkingdom = "Unknown";
    my $superkingdomTax = "Unknown";

    my $phylum = "Unknown";
    my $phylumTax = "Unknown";

    my $class = "Unknown";
    my $classTax = "Unknown";

    my $order = "Unknown";
    my $orderTax = "Unknown";

    my $family = "Unknown";
    my $familyTax = "Unknown";

    my $genus = "Unknown";
    my $genusTax = "Unknown";

    my $species = "Unknown";
    my $speciesTax = "Unknown";

    my $no_rank = "Unknown";
    my $no_rankTax = "Unknown";

    my $first=1;
    while ($taxID > 1){	
	my $parentID= $parent{$taxID};
	my $rankID=$rank{$taxID};

	if ($rankID eq 'superkingdom'){
	    $superkingdom = $names{$taxID};
	    $superkingdomTax = $taxID;
	}
	elsif ($rankID eq 'phylum'){
	    $phylum = $names{$taxID};
	    $phylumTax = $taxID;
	}
	elsif ($rankID eq 'class'){
	    $class = $names{$taxID};
	    $classTax = $taxID;
	}
	elsif ($rankID eq 'order'){
	    $order = $names{$taxID};
	    $orderTax = $taxID;
	}
	elsif ($rankID eq 'family'){
	    $family = $names{$taxID};
	    $familyTax = $taxID;
	}
	elsif ($rankID eq 'genus'){
	    $genus = $names{$taxID};
	    $genusTax = $taxID;
	}
	elsif ($rankID eq 'species'){
	    $species = $names{$taxID};
	    $speciesTax = $taxID;
	}
	elsif (($rankID eq 'no rank') && ($first)){
	    #
	    # try to get the strain name
	    #
	    $no_rank = $names{$taxID};
	    $no_rankTax = $taxID;
	    if ($type eq 'taxid'){
		if (exists($input_strain{$key})){
		    $no_rank=$input_strain{$key};
		}
	    }
	}
	$first=0;
	$taxID=$parentID;
    }
    print OUT "$Qname\t";
    print OUT "$superkingdom\t$superkingdomTax\t";
    print OUT "$phylum\t$phylumTax\t";
    print OUT "$class\t$classTax\t";
    print OUT "$order\t$orderTax\t";
    print OUT "$family\t$familyTax\t";
    print OUT "$genus\t$genusTax\t";
    print OUT "$species\t$speciesTax\t";
    print OUT "$no_rank\t$no_rankTax\n";
}
$datestring = localtime();
print LOG "## Local date and time $datestring - End program\n" if ($verbose);
if (defined($Getopt::Std::opt_l)){
    close(LOG);
}
if (defined($Getopt::Std::opt_o)){
    close(OUT);
}
