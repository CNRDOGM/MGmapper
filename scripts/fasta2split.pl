#!/usr/bin/env perl
my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
  $command .= "$ARGV[$i] ";
  $i++;
}
use Getopt::Std;
use strict;
use Cwd;
# Default parameters
*LOG=*STDERR;
my $verbose=0;
my $Verbose=0;
my $bp_max=10000000000;
#
# store sequences in memory until buffer_size is reached (1000000000 ~ 1Gb mem)
#
my $buffer_size=5000000000;
my $baseName="subset";
my $dummy_collapsed_class = 'collapsed_class_' . "$$";
my $dummy_collapsed_class_flag = 1;
#
# Process command line
#
getopts('hi:vVl:m:o:cb:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i name] [-o name] [-m number] [-b number] [-c] [-l logfile] [-s filename][-v]\n");
  print ("Description:\n");
  print ("$0 - split a fasta file in chunks of x basepairs in each subset - Sequences with same second field (genus) end up in same chunk until max number of basepairs is reached.\n");
  print ("Collapsed_name (often genus) is read as the second field in the fasta header, if no second field is present, the sequence name is used.\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input file name\n");
  print ("  -o  : output base name [$baseName]\n");
  print ("  -m  : max number basepairs in each chunk [$bp_max]\n");
  print ("  -b  : buffer size. Store sequences in memory until buffer_size is reached (5.000.000.000 ~ 5Gb mem) default [$buffer_size]\n");
  print ("  -c  : if no second field is present, use a dummy class name instead of sequence name [on]\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -v  : Verbose [off]\n");
  print ("  -V  : Additional verbose [off]\n");
  print ("\n");
 exit;
} # Usage

#
# Open input
#
if (not defined($Getopt::Std::opt_i)){
  # Read from standard input
  *INP = *STDIN;
  print STDERR "Can not read from stdin\n";
  exit;
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

if (defined($Getopt::Std::opt_l)){
    open(LOG,">$Getopt::Std::opt_l");
}
if (defined($Getopt::Std::opt_v)){
    $verbose=1;
}
if (defined($Getopt::Std::opt_V)){
    $Verbose=1;
}
if (defined($Getopt::Std::opt_m)){
    $bp_max=$Getopt::Std::opt_m;
}
if (defined($Getopt::Std::opt_b)){
    $buffer_size=$Getopt::Std::opt_b;
}
if (defined($Getopt::Std::opt_o)){
    $baseName=$Getopt::Std::opt_o;
}
if (defined($Getopt::Std::opt_c)){
    $dummy_collapsed_class_flag=0;
}
###############################################################################
# Main
#
###############################################################################
my $datestring = localtime();
print LOG "## Local date and time $datestring - start program\n" if ($verbose);
print LOG "# $command\n" if ($verbose);
my $thisDir=cwd();
print LOG "# working dir: $thisDir\n" if ($verbose);
if ($dummy_collapsed_class_flag){
    print LOG "# Entries with an absent second field being a word string, are given the dummy 'Collapsed_name': $dummy_collapsed_class\n" if ($verbose);
}
print LOG "# Genuses are sorted by size(bp) and selected from size sorted list. '+' = reading from top of list, '-' = reading from reverse end of list\n" if ($verbose);
print LOG "# Genus\tSeq_id\tSeq_id_len\tsum basepairs\toutput file\tgenus selected from '+' or '-' of sorted list\n" if ($verbose);
my $seq='';
my $line;
my $genus;
my $name;
my %s=();
my $len;
my %genuses=();
my $genusCount=0;
my %entry=();
my $fasta_entries=0;
my $tmp_line;
while (defined($_=<INP>)){
    chomp;
    $line=$_;
    $tmp_line=$line;
    $tmp_line =~ s/\'//g;
    
    if ($tmp_line =~ m/^>(\S+)/){
	$name=$1;
	$fasta_entries++;
	if ($Verbose){
	    if ($fasta_entries%10000==0){
		print LOG "# Processed $fasta_entries entries. Number of genuses= $genusCount\n";
	    }
	}
	if ($tmp_line =~ m/>(\S+)\s+(\S+)/){	    
	    $genus=$2;

	    if ($genus =~ m/\d+/){
		if ($dummy_collapsed_class_flag){
		    $genus=$dummy_collapsed_class;
		}
		else{
		    $genus=$name;
		}
	    }
	}
	else{
	    if ($dummy_collapsed_class_flag){
		$genus=$dummy_collapsed_class;
	    }
	    else{
		$genus=$name;
	    }
	}

	if (! exists($genuses{$genus})){
	    $genuses{$genus}=1;
	    $genusCount++;
	}
	$entry{$genus}{len}{$name} = 0;
    }
    else{
	$len = length($line);
	$s{$genus}{len} +=$len;
	$entry{$genus}{len}{$name} += $len;
    }
}
close(INP);
print LOG "# Done Reading. Processed $fasta_entries entries. Number of genuses= $genusCount\n" if ($Verbose);
my $bp=0;
#
# Store genus name
#
my %selected=();
my $tmp;
my $outFile;
my $unselected_genuses=$genusCount;
my %sequences=();
my %fastaChunks=();

my $bp_in_fasta_chunk=0;
my $i=1;
while ($unselected_genuses !=0){
    ($bp_in_fasta_chunk, $i) = &assign_sequence_to_fasta_chunk($bp_in_fasta_chunk, $i);
}

for (my $j=1;$j<=$i;$j++){
    $sequences{$j}='';
    my $size=length($sequences{$j});
    print LOG "Initializing variable for fasta chunk $j. size=$size\n" if ($Verbose);
    my $outFile= "$baseName" . ".$j";
    if (-e "$outFile"){
	print LOG "removing file: $outFile\n" if ($Verbose);
	system("rm $outFile");
    }
}

print LOG "# Genus count: $genusCount\n" if ($verbose);
#
# Read input file again and now print fasta entries in pre-defined output files
#
if (($Getopt::Std::opt_i=~/\.gz$/) || ($Getopt::Std::opt_i=~/\.Z$/)){
    open(INP,"gunzip -c $Getopt::Std::opt_i |") || die ("can't open file $Getopt::Std::opt_i: $!");
}
else{
    open(INP,"<$Getopt::Std::opt_i") || die ("can't open file $Getopt::Std::opt_i: $!");
}

my $size;
my $fastaChunk;
while (defined($_=<INP>)){
    $line=$_;
    $tmp_line=$line;
    $tmp_line =~ s/\'//g;

    if ($tmp_line =~ m/^>(\S+)/){
	$name=$1;
	if ($tmp_line =~ m/>(\S+)\s+(\S+)/){
	    $genus=$2;
	    if ($genus =~ m/\d+/){
		if ($dummy_collapsed_class_flag){
		    $genus=$dummy_collapsed_class;
		}
		else{
		    $genus=$name;
		}
	    }
	}
	else{
	    if ($dummy_collapsed_class_flag){
		$genus=$dummy_collapsed_class;
	    }
	    else{
		$genus=$name;
	    }
	}

	if (! exists($entry{$genus}{out}{$name})){
	    print LOG "A problem occured. output file for genus='$genus' id='$name' has not been defined\n";
	    die;
	}
	if ($entry{$genus}{out}{$name} == 0){
	    print LOG "A problem occured. len=0 for genus='$genus' id='$name'\n";
	    die;
	}
	$fastaChunk=$entry{$genus}{out}{$name};

    }

    $size += length($line);
    $sequences{$fastaChunk} .= $line;
    
    if ($size >= $buffer_size){
	foreach my $key (keys %sequences){
	    my $file = "$baseName" . '.' . "$key";
	    
	    my $variable_size = length($sequences{$key});
	    if ($variable_size > 0){
		if (-e "$file"){
		    open(OUT,">>$file");
		}
		else{
		    open(OUT,">$file");
		}
		print LOG "size= '$size, size variable $key= '$variable_size'. Printing to file: $file\n" if ($Verbose);
		print OUT ("$sequences{$key}");
		$sequences{$key}='';
		close(OUT);
	    }
	}
	$size=0;
    }
}
#
# flush sequence variables
#
foreach my $key (keys %sequences){
    my $file = "$baseName" . '.' . "$key";
    
    my $variable_size = length($sequences{$key});
    if ($variable_size > 0){
	if (-e "$file"){
	    open(OUT,">>$file");
	}
	else{
	    open(OUT,">$file");
	}
	print LOG "printing to file: $file\n" if ($Verbose);
	print OUT ("$sequences{$key}");
	$sequences{$key}='';
	close(OUT);
    }
}

print LOG ("Done\n") if ($verbose);
#
# End main program
#
sub assign_sequence_to_fasta_chunk{
    my ($bp_in_fasta_chunk, $i) = @_;

    foreach my $key (sort {$s{$b}{len} <=> $s{$a}{len}} keys %s){
	if (exists($selected{$key})){
	    next;
	}
	$bp = $s{$key}{len} + $bp_in_fasta_chunk;
	
	
	#
	# Does a single genus contain more bases than the bp_max limit?
	#
	if ($s{$key}{len} > $bp_max){
	    #
	    # insert all sequences belonging to one genus into 2 or more fasta chunks
	    #
	    ($bp_in_fasta_chunk, $i) = &insert_genus_in_multiple_chunks($key, $i);
	    print LOG "Unselected genuses= $unselected_genuses\n" if ($Verbose);

	    #
	    # Now fill up the last fasta chunk with a smaller genuses
	    #
	    ($bp_in_fasta_chunk, $i) = &fill_up_with_smaller_genuses($bp_in_fasta_chunk, $i);
	    print LOG "Unselected genuses= $unselected_genuses\n" if ($Verbose);
	    next;
	}
	else{
	    if ($bp <= $bp_max){
		($bp_in_fasta_chunk, $i) = &insert_genus($key);
		print LOG "Unselected genuses= $unselected_genuses\n" if ($Verbose);
		next;	
	    }
	    else{
		($bp_in_fasta_chunk, $i) = &fill_up_with_smaller_genuses($bp_in_fasta_chunk, $i);
		print LOG "Unselected genuses= $unselected_genuses\n" if ($Verbose);
	    }
	}
    }
    return($bp_in_fasta_chunk, $i);
}

sub insert_genus_in_multiple_chunks{
    my ($key, $i) = @_;
    my $outFile;
    
    foreach my $id (keys %{$entry{$key}{len}}){

	#
	# is a single fasta entry bigger that the bp limit?
	#
	if ($entry{$key}{len}{$id} > $bp_max){
	    if (exists($fastaChunks{$i})){
		$i++;
	    }
	    $bp_in_fasta_chunk=$entry{$key}{len}{$id};
	    $entry{$key}{out}{$id} = $i;
	    $fastaChunks{$i}=1;
	    $outFile= "$baseName" . ".$i";
	    print LOG "$key\t$id\t$entry{$key}{len}{$id}\t$bp_in_fasta_chunk\t$outFile\t'+'\n" if ($verbose);
	    $bp_in_fasta_chunk=0;
	    $i++;
	    next;
	}

	$bp_in_fasta_chunk += $entry{$key}{len}{$id};
	if ($bp_in_fasta_chunk <= $bp_max){
	    $fastaChunks{$i}=1;
	    $entry{$key}{out}{$id} = $i;
	    $outFile= "$baseName" . ".$i";
	    print LOG "$key\t$id\t$entry{$key}{len}{$id}\t$bp_in_fasta_chunk\t$outFile\t'+'\n" if ($verbose);
	}
	else{
	    if (exists($fastaChunks{$i})){
		$i++;
	    }
	    $fastaChunks{$i}=1;
	    $bp_in_fasta_chunk = $entry{$key}{len}{$id};
	    $entry{$key}{out}{$id} = $i;
	    $outFile= "$baseName" . ".$i";
	    print LOG "$key\t$id\t$entry{$key}{len}{$id}\t$bp_in_fasta_chunk\t$outFile\t'+'\n" if ($verbose);
	}
    }
    #
    # indicate that all species belonging to same genus have been selected
    #
    $selected{$key}=1;
    $unselected_genuses--;
    return($bp_in_fasta_chunk, $i);
}

sub insert_genus{
    my ($key) = @_;
    #
    # one entire set of sequences belonging to the same genus can be placed in one fasta chunk
    #
    my $outFile;
    foreach my $id (keys %{$entry{$key}{len}}){
	$bp_in_fasta_chunk += $entry{$key}{len}{$id};
	$fastaChunks{$i}=1;
	$entry{$key}{out}{$id} = $i;
	$outFile= "$baseName" . ".$i";
	print LOG "$key\t$id\t$entry{$key}{len}{$id}\t$bp_in_fasta_chunk\t$outFile\t'+'\n" if ($verbose);
    }
    $selected{$key}=1;
    $unselected_genuses--;
    return($bp_in_fasta_chunk, $i);
}

sub fill_up_with_smaller_genuses{
    my ($bp_in_fasta_chunk, $i) = @_;
    my $outFile;

    foreach my $low (sort {$s{$a}{len} <=> $s{$b}{len}} keys %s){
	if (exists($selected{$low})){
	    next;
	}
	my $bp=$s{$low}{len};
	my $tmp = $bp_in_fasta_chunk + $bp;
	if ($tmp <= $bp_max){
	    foreach my $id (keys %{$entry{$low}{len}}){
		$fastaChunks{$i}=1;
		$bp_in_fasta_chunk += $entry{$low}{len}{$id};
		$entry{$low}{out}{$id} = $i;
		$outFile= "$baseName" . ".$i";
		print LOG "$low\t$id\t$entry{$low}{len}{$id}\t$bp_in_fasta_chunk\t$outFile\t'-'\n" if ($verbose);
	    }
	    $selected{$low}=1;
	    $unselected_genuses--;
	}
	else{
	    if (exists($fastaChunks{$i})){
		$i++;
		$bp_in_fasta_chunk=0;
	    }
	}
    }
    return($bp_in_fasta_chunk, $i);
}
