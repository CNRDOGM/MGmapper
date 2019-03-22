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
my $fh;
my $verbose=0;
my $minimum_readCount= 10;
my $clade = 'strain';
my $minimum_relative_abundance = "0.01";
my $collapse_databases=0;

# minimum ratio of uniqReadCount/ReadCount
my $ratio_cutoff="0.005";

# maximum ratio of nm/nucleotides i.e. fraction of errors or mismatches compared to number of basepairs
my $max_nm_nuc_ratio="0.01";

my $min_coverage=0;
my $single_end_reads=0;
my $Nreads=0;
#
# A flag to determine if failed hits are printed to file defined by option -f
#
my $failed_flag=0;
#
# Process command line
#
getopts('hi:o:vl:n:a:c:Hr:dsm:g:f:N:F')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i name] [-o name] [-s] [-n number] [-a number]  [-r number] [-m number] [-F] [-g number] [-c name] [-f name] [-d] [-H] [-N number] [-l name]  [-v]\n");
  print ("Description:\n");
  print ("$0 - filter and merge .annot files\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input file name [STDIN]\n");
  print ("  -o  : output file name [STDOUT]\n");
  print ("  -s  : read count are treated as single-end reads [off]\n");
  print ("  -n  : minimum number of mapped reads [$minimum_readCount]\n");
  print ("  -a  : minimum abundance % abundance [$minimum_relative_abundance]\n");
  print ("  -r  : minimum uniqueReadCount/ReadCount ratio [$ratio_cutoff]\n");
  print ("  -m  : maximum nucleotide error fraction [$max_nm_nuc_ratio]\n");
  print ("  -F  : Filters -n, -a, -r and -m are all ([on]|off)\n");
  print ("  -g  : minimum coverage [$min_coverage]\n");
  print ("  -c  : collapse at clade level (strain, species, genus, family, order, class, phylum, superfamily) - default [$clade]\n");
  print ("  -f  : output File with failed hits\n");
  print ("  -d  : collapse clades from different databases [off]\n");
  print ("  -H  : print Header [off]\n");
  print ("  -N  : total number of reads in sample to calculate R_Abundance [$Nreads]\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -v  : Verbose [off]\n\n");
  print ("       Examples:\n");
  print ("       cat stat.Bacteria.annot stat.Plasmid.annot | MGmapper_classify.pl -m 0.01 -a 0.01 -r 0.005 -n 10 -c strain| less -S\n");
  print ("       Same result is obtained below as those above -m, -a, -r and -n values are default\n");
  print ("       cat stat.Bacteria.annot stat.Plasmid.annot | MGmapper_classify.pl -m 0.01 -a 0.01 -r 0.005 -n 10 -c strain| less -S\n");
  print ("       An easy way to turn off the default values for -m, -a, -r and -n is to specify -F (equal to setting -m 1 -a 0 -r 0 -n 0)\n");
  print ("       cat stat.Bacteria.annot stat.Plasmid.annot | MGmapper_classify.pl -F -c strain| less -S\n");
  print ("       #\n");
  print ("       # At strain -r 0.005 should be used\n");
  print ("       # At species -r 0.005 should be used and option -v prints lowest accepted relative abundance which is used in a second round with -r 0\n");
  print ("       cat stat.Bacteria.annot stat.Plasmid.annot | MGmapper_classify.pl -m 0.01 -a 0.01 -r 0.005 -n 10 -c species -v -l logfile| less -S\n");
  print ("       read new abundance in logfile\n");
  print ("       cat stat.Bacteria.annot stat.Plasmid.annot | MGmapper_classify.pl -m 0.01 -a \$lowest_from_round1 -r 0 -n 10 -c species -v -l logfile -f negative_hits.txt| less -S\n");
  print ("       # At all other levels i.e. genus phylum etc -r 0 should be used i.e. uniq read count ratio is not used\n");
  print ("       #\n");
  print ("       # option -n 10 sets a lower readCount cutoff. Two strains belonging to the same species with each say 5 reads mapped to each,\n");
  print ("       # results in none of them being selected at strain level, but at species level there are 10 reads and the species will be accepted as a true hit.\n");
  print ("       MGmapper_classify.pl -i file.annot -c genus -r 0 | less -S\n");
  print ("       Collapse clades from different databases:\n");
  print ("       Ex. cat *.annot | MGmapper_classify.pl -c genus -d -r 0| less -S\n");
  print ("       Collapse clades from different databases and add header:\n");
  print ("       Ex. cat *.annot | MGmapper_classify.pl -c genus -d -H -r 0| less -S\n");
  print ("       coverage= 'covered_positions'/'size'\n");
  print ("       paired-end reads: abundance= (100*read_count)/(size*2)\n");
  print ("       single-end reads: abundance= (100*read_count)/(size), specify -s for single-end reads\n");
  print ("       uniqReadCount ratio (option r) = uniqReadCount/readCount\n");
  print ("       misMatch ratio (option m) = misMatches/nucleotides\n");
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
if (defined($Getopt::Std::opt_d)){
    $collapse_databases=1;
}
if (defined($Getopt::Std::opt_n)){
    $minimum_readCount=$Getopt::Std::opt_n;
}
if (defined($Getopt::Std::opt_c)){
    $clade=$Getopt::Std::opt_c;
}
if (defined($Getopt::Std::opt_a)){
    $minimum_relative_abundance=$Getopt::Std::opt_a;
}
if (defined($Getopt::Std::opt_r)){
    $ratio_cutoff=$Getopt::Std::opt_r;
}
if (defined($Getopt::Std::opt_g)){
    $min_coverage=$Getopt::Std::opt_g;
}
if (defined($Getopt::Std::opt_m)){
    $max_nm_nuc_ratio=$Getopt::Std::opt_m;
}
if (defined($Getopt::Std::opt_s)){
    $single_end_reads=1;
}
if (defined($Getopt::Std::opt_N)){
    $Nreads=$Getopt::Std::opt_N;
}
if (defined($Getopt::Std::opt_f)){
    open(FAILED,">$Getopt::Std::opt_f");
    $failed_flag=1;
}
if (defined($Getopt::Std::opt_F)){
    # change default option -m (max edit distance) 
    $max_nm_nuc_ratio=1;

    # change default option -r (min uniq readCount ratio)
    $ratio_cutoff=0;

    # change default option -a (min size normalized abundance)
    $minimum_relative_abundance=0;

    # change default option -n (min readCount)
    $minimum_readCount=0;
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
my $db;
my $strain;
my $size;
my $nucleotides;
my $covered_positions;
my $readCount;
my $readCountUniq;
my $nm;

my $superfamily;
my $superfamilyTax;
my $phylum;
my $phylumTax;
my $class;
my $classTax;
my $order;
my $orderTax;
my $family;
my $familyTax;
my $genus;
my $genusTax;
my $species;
my $speciesTax;
my $strain;
my $strainTax;

my $coverage;
my $depth;
my $relative_abundance;
my $desc;

my %r=();
my %collapse_db=();
# store database names
my %d=();
my @databases;

my $taxid;
my $col;
my $lastCol;
my $taxid_column;
my @columns= qw(superfamily superfamilyTax phylum phylumTax class classTax order orderTax family familyTax genus genusTax species speciesTax strain strainTax);
if ($clade eq 'superfamily'){
    $col=9;
    $taxid_column=10;
    $lastCol=12;
}
elsif ($clade eq 'phylum'){
    $col=11;
    $taxid_column=12;
    $lastCol=14;
}
elsif ($clade eq 'class'){
    $col=13;
    $taxid_column=14;
    $lastCol=16;
}
elsif ($clade eq 'order'){
    $col=15;
    $taxid_column=16;
    $lastCol=18;
}
elsif ($clade eq 'family'){
    $col=17;
    $taxid_column=18;
    $lastCol=20;
}
elsif ($clade eq 'genus'){
    $col=19;
    $taxid_column=20;
    $lastCol=22;
}
elsif ($clade eq 'species'){
    $col=21;
    $taxid_column=22;
    $lastCol=24;
}
elsif ($clade eq 'strain'){
    $col=1;
    $taxid_column=24;
    $lastCol=26;
}
else{
    print LOG "unknown value for option -c\n";
    exit;
}
my $name;
if (defined($Getopt::Std::opt_H)){
    if (!$collapse_databases){
	print OUT "# Database\tName\tS_Abundance (%)\tR_Abundance (%)\tSize (bp)\tSeq_count\tNucleotides\tCovered_positions\tCoverage\tDepth\tReads\tReads_uniq\tEdit_distance\tDescription";
	if ($failed_flag){
	    print FAILED "# Database\tName\tS_Abundance (%)\tR_Abundance (%)\tSize (bp)\tSeq_count\tNucleotides\tCovered_positions\tCoverage\tDepth\tReads\tReads_uniq\tEdit_distance\tDescription";
	}
    }
    else{
	print OUT "# Name\tS_Abundance (%)\tR_Abundance (%)\tSize (bp)\tSeq_count\tNucleotides\tCovered_positions\tCoverage\tDepth\tReads\tReads_uniq\tEdit_distance\tDescription";
	if ($failed_flag){
	    print FAILED "# Name\tS_Abundance (%)\tR_Abundance (%)\tSize (bp)\tSeq_count\tNucleotides\tCovered_positions\tCoverage\tDepth\tReads\tReads_uniq\tEdit_distance\tDescription";
	}
    }

    my $end = $lastCol - 11;
    for (my $i=0;$i<=$end;$i++){
	print OUT "\t$columns[$i]";
	if ($failed_flag){
	    print FAILED "\t$columns[$i]";
	}
    }
    print OUT "\n";
    if ($failed_flag){
	print FAILED "\n";	
    }

}
my $aveReadLen;
my $total_reads=0;
my $rejected_reads=0;
while (defined($_=<INP>)){
    chomp;
    @w=split(/\t/);
    
    $db=$w[0];
    if (! exists($d{$db})){
	print LOG "# databases: $db\n" if ($verbose);
	push(@databases,$db);
	$d{$db}=$db;
    }
    $strain=$w[1];
    $size=$w[2];
    $nucleotides=$w[3];
    $covered_positions=$w[4];
    $readCount = $w[5];
    $total_reads += $readCount;
    $readCountUniq = $w[6];
    $nm=$w[7];
    $desc=$w[8];
    for (my $i=9;$i<=24;$i++){
	if ($w[$i] eq 'Unknown'){
	    $w[$i] = '-';
	}
    }
    
    #
    # Taxonomy
    #
    
    $superfamily= $w[9];
    $superfamilyTax= $w[10];

    $phylum=$w[11];
    $phylumTax=$w[12];

    $class=$w[13];
    $classTax=$w[14];

    $order=$w[15];
    $orderTax=$w[16];

    $family=$w[17];
    $familyTax=$w[18];
    
    $genus=$w[19];
    $genusTax=$w[20];

    $species=$w[21];
    $speciesTax=$w[22];

    $strain=$w[23];
    $strainTax=$w[24];

    $coverage=$covered_positions/$size;
    $depth=$nucleotides/$size;

    
    $name=$w[$col];
    $taxid=$w[$taxid_column];

    if ($clade eq 'strain'){
	$taxid=$name;
	$r{$taxid}{$db}{desc}=$desc;
	if ($collapse_databases){
	    $collapse_db{$taxid}{desc}=$desc;
	}
    }
    if ($taxid eq '-'){
	$name="Unclassified";
    }

    if ($clade ne 'strain'){
	$r{$taxid}{$db}{desc}="collapsed at $clade level";
	if ($collapse_databases){
	    $collapse_db{$taxid}{desc}="collapsed at $clade level";
	}
    }

    if (!$collapse_databases){
	$r{$taxid}{$db}{seqCount}++;
	$r{$taxid}{$db}{taxid} = $taxid;
	$r{$taxid}{$db}{name} = $name;
	$r{$taxid}{$db}{size} += $size;
	$r{$taxid}{$db}{nucleotides} += $nucleotides;
	$r{$taxid}{$db}{covered_positions} += $covered_positions;
	$r{$taxid}{$db}{readCount} += $readCount;
	$r{$taxid}{$db}{readCountUniq} += $readCountUniq;
	$r{$taxid}{$db}{nm} += $nm;
	if ($single_end_reads){
	    $r{$taxid}{$db}{relative_abundance} += ($readCount*100)/$size;
	}
	else{
	    $r{$taxid}{$db}{relative_abundance} += ($readCount*100)/($size*2);
	}
	
	$r{$taxid}{$db}{superfamily} = $superfamily;
	$r{$taxid}{$db}{superfamilyTax} = $superfamilyTax;
	
	$r{$taxid}{$db}{phylum} = $phylum;
	$r{$taxid}{$db}{phylumTax} = $phylumTax;
	
	$r{$taxid}{$db}{class} = $class;
	$r{$taxid}{$db}{classTax} = $classTax;
	
	$r{$taxid}{$db}{order} = $order;
	$r{$taxid}{$db}{orderTax} = $orderTax;
	
	$r{$taxid}{$db}{family} = $family;
	$r{$taxid}{$db}{familyTax} = $familyTax;
	
	$r{$taxid}{$db}{genus} = $genus;
	$r{$taxid}{$db}{genusTax} = $genusTax;
	
	$r{$taxid}{$db}{species} = $species;
	$r{$taxid}{$db}{speciesTax} = $speciesTax;
	
	$r{$taxid}{$db}{strain} = $strain;
	$r{$taxid}{$db}{strainTax} = $strainTax;	
    }
    else{
	$collapse_db{$taxid}{seqCount}++;
	$collapse_db{$taxid}{taxid} = $taxid;
	$collapse_db{$taxid}{name} = $name;
	$collapse_db{$taxid}{size} += $size;
	$collapse_db{$taxid}{nucleotides} += $nucleotides;
	$collapse_db{$taxid}{covered_positions} += $covered_positions;
	$collapse_db{$taxid}{readCount} += $readCount;
	$collapse_db{$taxid}{readCountUniq} += $readCountUniq;
	$collapse_db{$taxid}{nm} += $nm;

	if ($single_end_reads){
	    $collapse_db{$taxid}{relative_abundance} += ($readCount*100)/$size;
	}
	else{
	    $collapse_db{$taxid}{relative_abundance} += ($readCount*100)/($size*2);
	}
	
	$collapse_db{$taxid}{superfamily} = $superfamily;
	$collapse_db{$taxid}{superfamilyTax} = $superfamilyTax;
	
	$collapse_db{$taxid}{phylum} = $phylum;
	$collapse_db{$taxid}{phylumTax} = $phylumTax;
	
	$collapse_db{$taxid}{class} = $class;
	$collapse_db{$taxid}{classTax} = $classTax;
	
	$collapse_db{$taxid}{order} = $order;
	$collapse_db{$taxid}{orderTax} = $orderTax;
	
	$collapse_db{$taxid}{family} = $family;
	$collapse_db{$taxid}{familyTax} = $familyTax;
	
	$collapse_db{$taxid}{genus} = $genus;
	$collapse_db{$taxid}{genusTax} = $genusTax;
	
	$collapse_db{$taxid}{species} = $species;
	$collapse_db{$taxid}{speciesTax} = $speciesTax;
	
	$collapse_db{$taxid}{strain} = $strain;
	$collapse_db{$taxid}{strainTax} = $strainTax;	
    }
}

my $uniq_ratio;
my $nm_over_nuc_ratio;
my $depth;
my $expected_coverage;

#
# Databases are not collapsed
#

my $observed_min_rel_abundance=0;
my $observed_max_rel_abundance=0;
my $lines=0;
my $hits=0;
if (!$collapse_databases){
    foreach my $database (@databases){
	foreach my $key (sort {$r{$b}{$database}{relative_abundance} <=> $r{$a}{$database}{relative_abundance}} keys %r){
	    if (! exists($r{$key}{$database}{size})){
		next;
	    }

	    $lines++;
	    if ($r{$key}{$database}{relative_abundance} < $minimum_relative_abundance){
		$rejected_reads += $r{$key}{$database}{readCount};
		if ($failed_flag){
		    $fh = *FAILED;
		    &printout($fh,$key,$database);
		}
		next;
	    }

	    #
	    # store some min and max observed relative abundances to be written to logfile only
	    #
	    if ($lines == 1){
		$observed_min_rel_abundance = $r{$key}{$database}{relative_abundance};
	    }
	    	    


	    if ($r{$key}{$database}{readCount} < $minimum_readCount){
		$rejected_reads += $r{$key}{$database}{readCount};
		if ($failed_flag){
		    $fh = *FAILED;
		    &printout($fh,$key,$database);
		}
		next;
	    }

	    $uniq_ratio=$r{$key}{$database}{readCountUniq}/$r{$key}{$database}{readCount};
	    if ($uniq_ratio < $ratio_cutoff){
		$rejected_reads += $r{$key}{$database}{readCount};
		if ($failed_flag){
		    $fh = *FAILED;
		    &printout($fh,$key,$database);
		}
		next;
	    }

	    $nm_over_nuc_ratio=$r{$key}{$database}{nm}/$r{$key}{$database}{nucleotides};
	    if ($nm_over_nuc_ratio > $max_nm_nuc_ratio){
		$rejected_reads += $r{$key}{$database}{readCount};
		if ($failed_flag){
		    $fh = *FAILED;
		    &printout($fh,$key,$database);
		}
		next;
	    }

	    $coverage=$r{$key}{$database}{covered_positions}/$r{$key}{$database}{size};
	    if ($coverage < $min_coverage){
		$rejected_reads += $r{$key}{$database}{readCount};
		if ($failed_flag){
		    $fh=*FAILED;
		    &printout($fh,$key,$database);
		}
		next;
	    }

	    if ($r{$key}{$database}{relative_abundance} > $observed_max_rel_abundance){
		$observed_max_rel_abundance = $r{$key}{$database}{relative_abundance};
	    }
	    if ($r{$key}{$database}{relative_abundance} < $observed_min_rel_abundance){
		$observed_min_rel_abundance = $r{$key}{$database}{relative_abundance};
	    }
	    $hits++;
	    $fh=*OUT;
	    &printout($fh,$key,$database);

	}
    }
}
else{
    foreach my $key (sort {$collapse_db{$b}{relative_abundance} <=> $collapse_db{$a}{relative_abundance}} keys %collapse_db){
	if (! exists($collapse_db{$key}{size})){
	    next;
	}
	$lines++;

	if ($collapse_db{$key}{relative_abundance} < $minimum_relative_abundance){
	    $rejected_reads += $collapse_db{$key}{readCount};
	    if ($failed_flag){
		$fh=*FAILED;
		&printout_collapse($fh,$key);
	    }
	    next;
	}

	#
	# store some min and max observed relative abundances to be written to logfile only
	#
	if ($lines == 1){
	    $observed_min_rel_abundance = $collapse_db{$key}{relative_abundance};
	}
	
	
	if ($collapse_db{$key}{readCount} < $minimum_readCount){
	    $rejected_reads += $collapse_db{$key}{readCount};
	    if ($failed_flag){
		$fh=*FAILED;
		&printout_collapse($fh,$key);
	    }
	    next;
	}

	$uniq_ratio=$collapse_db{$key}{readCountUniq}/$collapse_db{$key}{readCount};
	if ($uniq_ratio < $ratio_cutoff){
	    $rejected_reads += $collapse_db{$key}{readCount};
	    if ($failed_flag){
		$fh=*FAILED;
		&printout_collapse($fh,$key);
	    }
	    next;
	}
	
	$nm_over_nuc_ratio=$collapse_db{$key}{nm}/$collapse_db{$key}{nucleotides};
	if ($nm_over_nuc_ratio > $max_nm_nuc_ratio){
	    $rejected_reads += $collapse_db{$key}{readCount};
	    if ($failed_flag){
		$fh=*FAILED;
		&printout_collapse($fh,$key);
	    }
	    next;
	}


	$depth=$collapse_db{$key}{nucleotides}/$collapse_db{$key}{size};
	$expected_coverage=1-exp(-$depth);
	my $coverage = $collapse_db{$key}{covered_positions}/$collapse_db{$key}{size};
	
	if ($coverage < $min_coverage){
	    $rejected_reads += $collapse_db{$key}{readCount};
	    if ($failed_flag){
		$fh=*FAILED;
		&printout_collapse($fh,$key);
	    }
	    next;
	}
	$fh=*OUT;

	if ($collapse_db{$key}{relative_abundance} > $observed_max_rel_abundance){
	    $observed_max_rel_abundance = $collapse_db{$key}{relative_abundance};
	}
	if ($collapse_db{$key}{relative_abundance} < $observed_min_rel_abundance){
	    $observed_min_rel_abundance = $collapse_db{$key}{relative_abundance};
	}
	$hits++;
	&printout_collapse($fh,$key);
    }
}

if ($verbose){
    my $perc=0;
    if ($total_reads > 0){
	$perc=$rejected_reads*100/$total_reads;
    }
    print LOG "# Number of hits processed:\t$lines\n";
    print LOG "# Number of hits accepted:\t$hits\n";
    printf LOG "# Total_reads= $total_reads\trejected_reads= $rejected_reads\t % rejected= %.3f\n",$perc;
    print LOG ("# Maximum accepted abundance value=\t$observed_max_rel_abundance\n");
    print LOG ("# Minimum accepted abundance value=\t$observed_min_rel_abundance\n");
}
if (defined($Getopt::Std::opt_i)){
    close(INP);
}
if (defined($Getopt::Std::opt_o)){
    close(OUT);
}
if (defined($Getopt::Std::opt_l)){
    close(LOG);
}

sub printout_collapse{
    my ($fh, $key)=@_;

    printf $fh "%s",$collapse_db{$key}{name};
    printf $fh "\t%.3f",$collapse_db{$key}{relative_abundance};
    if (defined($Getopt::Std::opt_N)){
	my $value=0;
	if ($Nreads > 0){
	    $value = $collapse_db{$key}{readCount}*100/$Nreads;
	}
	printf $fh "\t%.3f",$value;
    }
    else{
	print $fh "\t-";
    }
    printf $fh "\t%d",$collapse_db{$key}{size};
    printf $fh "\t%d",$collapse_db{$key}{seqCount};
    printf $fh "\t%d",$collapse_db{$key}{nucleotides};
    
    
    printf $fh "\t%d",$collapse_db{$key}{covered_positions};

    $coverage=$collapse_db{$key}{covered_positions}/$collapse_db{$key}{size};
    printf $fh "\t%.3f",$coverage;
    
    my $depth=$collapse_db{$key}{nucleotides}/$collapse_db{$key}{size};
    printf $fh "\t%.3f",$depth;
    
    printf $fh "\t%d",$collapse_db{$key}{readCount};
    printf $fh "\t%d",$collapse_db{$key}{readCountUniq};
    printf $fh "\t%d",$collapse_db{$key}{nm};
    printf $fh "\t%s",$collapse_db{$key}{desc};
    
    my $column=12;
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{superfamily};
	printf $fh "\t%s", $collapse_db{$key}{superfamilyTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{phylum};
	printf $fh "\t%s", $collapse_db{$key}{phylumTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{class};
	printf $fh "\t%s", $collapse_db{$key}{classTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{order};
	printf $fh "\t%s", $collapse_db{$key}{orderTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{family};
	printf $fh "\t%s", $collapse_db{$key}{familyTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{genus};
	printf $fh "\t%s", $collapse_db{$key}{genusTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{species};
	printf $fh "\t%s", $collapse_db{$key}{speciesTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $collapse_db{$key}{strain};
	printf $fh "\t%s", $collapse_db{$key}{strainTax};
    }
    print $fh "\n"
}

sub printout{
    my ($fh,$key,$database) = @_;
    printf $fh "%s",$database;
    printf $fh "\t%s",$r{$key}{$database}{name};
    
    printf $fh "\t%.3f",$r{$key}{$database}{relative_abundance};

    if (defined($Getopt::Std::opt_N)){
	my $value=0;
	if ($Nreads > 0){
	    $value = $r{$key}{$database}{readCount}*100/$Nreads;
	}
	printf $fh "\t%.3f",$value;
    }
    else{
	print $fh "\t-";
    }
    
    printf $fh "\t%d",$r{$key}{$database}{size};
    printf $fh "\t%d",$r{$key}{$database}{seqCount};
    printf $fh "\t%d",$r{$key}{$database}{nucleotides};
    
    
    printf $fh "\t%d",$r{$key}{$database}{covered_positions};
    

    my $depth=$r{$key}{$database}{nucleotides}/$r{$key}{$database}{size};
    my $coverage = $r{$key}{$database}{covered_positions}/$r{$key}{$database}{size};

    printf $fh "\t%.3f",$coverage;
    printf $fh "\t%.3f",$depth;    
    printf $fh "\t%d",$r{$key}{$database}{readCount};
    printf $fh "\t%d",$r{$key}{$database}{readCountUniq};
    printf $fh "\t%d",$r{$key}{$database}{nm};
    printf $fh "\t%s",$r{$key}{$database}{desc};
    
    my $column=12;
    if ($column <= $lastCol){
	printf $fh "\t%s", $r{$key}{$database}{superfamily};
	printf $fh "\t%s", $r{$key}{$database}{superfamilyTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $r{$key}{$database}{phylum};
	printf $fh "\t%s", $r{$key}{$database}{phylumTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $r{$key}{$database}{class};
	printf $fh "\t%s", $r{$key}{$database}{classTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $r{$key}{$database}{order};
	printf $fh "\t%s", $r{$key}{$database}{orderTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $r{$key}{$database}{family};
	printf $fh "\t%s", $r{$key}{$database}{familyTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $r{$key}{$database}{genus};
	printf $fh "\t%s", $r{$key}{$database}{genusTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $r{$key}{$database}{species};
	printf $fh "\t%s", $r{$key}{$database}{speciesTax};
	$column +=2;
    }
    if ($column <= $lastCol){
	printf $fh "\t%s", $r{$key}{$database}{strain};
	printf $fh "\t%s", $r{$key}{$database}{strainTax};
    }
    print $fh "\n"
}
