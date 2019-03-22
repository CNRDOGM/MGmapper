#!/usr/bin/env perl
my $version = "SE_2.24";

BEGIN{
    if (exists $ENV{MGmap_HOME}) {
        push(@INC, $ENV{MGmap_HOME} . '/modules');
    } else {
        die "Environment variable MGmap_HOME does not exists. Make it point to MGmapper's home.\n";
    }
    *LOG=*STDERR;
}

umask 0022;
my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
  $command .= "$ARGV[$i] ";
  $i++;
}
#use warnings;
use Getopt::Std;
use Time::HiRes;
use Cwd;
use strict;
use Parallel::ChildManager;
use MGmapper::MGmapper;

my $origDir = cwd();
#
# Process command line
#
getopts('hi:f:d:kKq:m:F:C:c:t:Te:EpSn:ws:a:vP:Q:B:RA:VMN:r:U:D:L:1:')||Usage();
if (defined($Getopt::Std::opt_V)){
    print "$version\n";
    exit;
}

#
# program calls
#
my $prog_xlsWrite = "$ENV{MGmap_HOME}/scripts/MGcsv2xls.py";
my $prog_pileup2fasta = "$ENV{MGmap_HOME}/scripts/pileup2fasta.pl";
my $prog_mpileup2stat = "$ENV{MGmap_HOME}/scripts/mpileup2stat.pl";
my $prog_filterBam = "$ENV{MGmap_HOME}/scripts/filter_bam.pl";
my $prog_MGmapper_classify = "$ENV{MGmap_HOME}/MGmapper_classify.pl";

if (! -e $prog_xlsWrite){
    die "Program not found: $prog_xlsWrite\n";
}
if (! -e $prog_pileup2fasta){
    die "Program not found: $prog_pileup2fasta\n";
}
if (! -e $prog_mpileup2stat){
    die "Program not found: $prog_mpileup2stat\n";
}
if (! -e $prog_filterBam){
    die "Program not found: $prog_filterBam\n";
}
#
# parameter files
#
my $parameterFile = "$ENV{MGmap_HOME}/databases.txt";
my $parameterFile_adapter = "$ENV{MGmap_HOME}/adapter.txt";

#
# third-party programs
#
my $prog_python = &MGmapper::find_program('MGmap_PYTHON', "python");
my $prog_bwa = &MGmapper::find_program('MGmap_BWA', "bwa");
my $prog_samtools = &MGmapper::find_program('MGmap_SAMTOOLS', "samtools");
#my $prog_samtools = "/services/tools/samtools/0.1.18/bin/samtools";
my $prog_bedtools = &MGmapper::find_program('MGmap_BEDTOOLS', "bedtools");
my $prog_pigz = &MGmapper::find_program('MGmap_PIGZ', "pigz");
my $prog_cutadapt = &MGmapper::find_program('MGmap_CUTADAPT', "cutadapt");
{
    my $st = `$prog_samtools 2>&1`;
    my @tmp = split(m/\n/, $st);
    @tmp = grep(m/^Version/, @tmp);
    push(@tmp, "Version: no version");
    $tmp[0] =~ m/^Version:\s+(\S+)/;
    die "Samtools $prog_samtools must be version 1.6\n" unless $1 eq '1.6';
}

my $datestring;
my $databaseCounter=0;
my %databaseNames=();
my %db_version=();

my @dbNames=();
my %idx=();
my %databaseCount=();
my %countEntries=();
my %runningMode=();
my %info=();

&readParameterFile;

#
# Default parameters
#
my $cm;
$i=0;
my $wwwMode=0;
my $workDir = "mapper_$$";
my $cleanFile;
my $redo=0;
my $commonName;
my $qualityCut = 30;
my $minLen = 30;
my $cores = 1;
my $matchRatio = 0.8;
my $matchRatioFlag = 1;
my $absoluteMatchCount = $minLen;
my $absoluteMatchCountFlag = 0;
my $mapNotPhiXReads=0;
my $sampleName = 'sample';
my $SNP_threshold=4;
my $QUALITY_BASE = 33;
my $alignment_score_min=30;
my $max_edit_distance=5;

#
# parameters for MGmapper_classify.pl
#
my $min_relative_abundance=0.01;
my $min_readCount = 20;
my $min_coverage=0;
my $min_readCountUniq_ratio=0.005;
my $max_mismatch_ratio=0.15;

my $unix_cmd1='';
#
# Total number of reads that are ok either from the start or after cutadapt
#

my $biological_relevant_reads=0;
my @databases_fullmode=();
my @databases_chainmode=();
my @databases_all=();
my @markdup_all=();
my $markdup_count=0;
my %abundance=();
my %fastaOut=();
my @rmFiles;
my $verbose=1;
my $cleanup=1;
my $cleanupCutadapt=1;
my $noclean="$workDir/noclean";
my $slimMode=0;
my $remove_PCR_duplicates=0;
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}
sub Usage {
    print ("Usage: $0 [-h] [-P parameter file] [-Q parameter file] [-i fastq file] [-f fastq file list] [-d outDir] [-q number] [-m number] [-F num serarate with comma] [-C numbers separate with comma] [-c cores] [-t number] [-T] [-e integer] [-E] [-p] [-s number] [-k] [-K] [-S] [-x] [-n name] [-A number] [-N number] [-M] [-r number] [-U number] [-D number] [-L number] [-V] [-1 unix_command]\n");
    print ("Description:\n");
    print ("$0 - Map SE reads against genomic sequence databases\n");
    print ("\n");
    print ("Options:\n");
    print ("  -h  : display this message\n");
    print ("  -V  : show version \($version\) and exit\n");
    print ("  -P  : parameter file with reference to database [$parameterFile]\n");
    print ("  -Q  : parameter file with adapters removed by cutadapt [$parameterFile_adapter]\n");
    print ("  -c  : cores to use [$cores]\n");
    print ("  Input files:\n");
    print ("  -i  : input fastq file\n");
    print ("  -f  : input list with fastq file names\n");
    print ("  Parameters to filter bwa mem hits for a single read:\n");
    print ("  -t  : true hits have a Match/readLen >= fraction [$matchRatio]\n");
    print ("  -T  : option -t is active [on]\n");
    print ("  -e  : minimum number of Matches+Mismatches for a valid read [$absoluteMatchCount] - only active if -E is defined\n");
    print ("  -E  : option -e is active [off]\n");
    print ("  -A  : minumum alignment score for a read [$alignment_score_min]\n");
    print ("  -N  : maximum edit distance\n");
    print ("  -p  : remove PCR duplicates [off]\n");
    print ("  Parameters for post-processing of reads assigned to a specific strain or species:\n");
    print ("  -r  : minimum size normalized abundance [$min_relative_abundance]\n");
    print ("  -U  : minimum read count [$min_readCount]\n");
    print ("  -D  : minimum Reads_uniq/Reads ratio [$min_readCountUniq_ratio]\n");
    print ("  -L  : maximum Edit_distance/Nucleotides ratio [$max_mismatch_ratio] i.e. max nucleotide mis-match ratio\n");
    print ("  Cutadapt parameters:\n");
    print ("  -S  : skip cutadapt [off]\n");
    print ("  -m  : minimum read length in cutadapt [$minLen]\n");
    print ("  -q  : Quality q cutoff in cutadapt [$qualityCut]\n");
    print ("  -B  : QUALITY_BASE, quality values are encoded as ascii(quality + QUALITY_BASE) [$QUALITY_BASE]\n");
    print ("  Mapping mode and database selections:\n");
    print ("  -F  : map reads in Full mode against these databases - comma separated numbers\n");
    print ("  -C  : map reads in Best mode against these databases - comma separated numbers - the order matters\n");
    print ("  Output files\/directory\n");
    print ("  -d  : output directory [$workDir]\n");
    print ("  -n  : Sample name [$sampleName]\n");
    print ("  Assembled fasta files and parameters:\n");
    print ("  -a  : make assembled fasta files for reads mapping to these databases - comma separated numbers - order dont matter ex.'1,2'\n");
    print ("  -s  : Min number of occurences to accept a variation in contig fasta file [$SNP_threshold]\n");
    print ("  -R  : make Read count matrices [off]\n");
    print ("  -k  : keep all files [off]\n");
    print ("  -K  : keep all cutadapt files [off]\n");
    print ("  -M  : slim mode i.e. remove all bam and unmapped.fq [off]\n");
    print ("  Command to be executed after everything else has finished:\n");
    print ("  -1  : unix command\n");

    print ("---Databases---\n");
    for (my $i=1;$i<=$databaseCounter;$i++){
	printf ("   $i\t%-20s$idx{desc}{$i}\n",$databaseNames{$i});
    }
    exit;
} # Usage
&parseOptions;
###############################################################################
#
# Main program start
#
###############################################################################
my $thisDir=cwd();
$datestring = localtime();

#
# save running time for most system calls
#
my $fh_runtime=\*RUNTIME; 
MGmapper::init_runtime_file($fh_runtime, "$workDir/log/runtime.log", $datestring);


my $fh_total_time= \*TIME;
my $fh_filestat=\*FILESTAT;
if ($wwwMode){
    #
    # save a simple file with start and end times for running MGmapper - only used in www-mode
    #
    MGmapper::init_time_file($fh_total_time, "$workDir/misc/time.tab", $workDir, $datestring);


    #
    # A file with readcount statistics - only used in www-mode
    #
    MGmapper::init_filestat($fh_filestat, "$workDir/misc/filestat.tab");
}


my $start_time1 = [Time::HiRes::gettimeofday()];
my ($user1, $system1, $child_user1, $child_system1) = times;

MGmapper::make_run_info_file($workDir, $version, $command, $datestring);

my $fh_stat = \*STAT;
MGmapper::init_MGmapper_summary($fh_stat, "$workDir/log/MGmapper.summary", $thisDir, $command);

my $cmd;
my %readCount=();
my $readCountFile = "$workDir/misc/readCount.txt";


my $sumReadsBefore=0;
my $sumReadsAfter=0;
my $fastqFiles=0;


if (! defined($Getopt::Std::opt_S)){
    if ($wwwMode){
	print $fh_filestat "#Counter\tFileName\tReads in\tReads after trimming\n";
    }
}
else{
    if ($wwwMode){
	print $fh_filestat "#Counter\tFileName\tReads in\n";
    }
    print $fh_stat "Counter\tFileName\tReads\n";
}

my $cleanFileFlag=0;
if (-e "$workDir/bam/cleaned.nophiX.bam"){
    $cleanFileFlag=1;
    if (! exists($readCount{$cleanFile})){
	$biological_relevant_reads = &bamCountReads($cleanFile);
    }
}


if (-e $readCountFile){
    open(COUNT,"<$readCountFile");
    while (defined($_=<COUNT>)){
	chomp;
	my @w=split(/\t+/);
	my $file = $w[0];
	my $reads = $w[1];
	$readCount{$file}=$reads;
	print LOG "# Reading from $readCountFile: $_\n" if ($verbose);
    }
    close(COUNT);
}

my %readsBefore=();
my %readsAfter=();
my %cutadaptFile=();


my @fileHolder=();
my $cutFile = "$workDir/misc/cutadapt.all.fq";
if (defined($Getopt::Std::opt_i)){
    push(@fileHolder,$Getopt::Std::opt_i);
}
elsif (defined($Getopt::Std::opt_f)){
    @fileHolder = MGmapper::readFastqList($Getopt::Std::opt_f);
}

foreach my $file (@fileHolder){
    if (exists($readCount{$file})){
	$readsBefore{$file}=$readCount{$file};
    }
    else{
	open(COUNT,">>$readCountFile");
	$readsBefore{$file}=MGmapper::countReads($file,\*LOG, $verbose);
	print COUNT "$file\t$readsBefore{$file}\n";
	close(COUNT);
    }
    $sumReadsBefore += $readsBefore{$file};
}


if (! defined($Getopt::Std::opt_S)){
    &run_cutadapt();
}
else{
    foreach my $file (@fileHolder){
	$fastqFiles++;
	if ($wwwMode){
	    my @FS=split(/\//,$file);
	    print $fh_filestat "$fastqFiles\t$FS[-1]\t$readsBefore{$file}\n";
	}
	print $fh_filestat "$fastqFiles\t$file\t$readsBefore{$file}\n";
    }
}

if (-e $cleanFile){
    &mapping($cleanFile);
}
elsif ($#fileHolder >= 0){
    #
    # many input fastq files
    #
    if (! -e $cutFile){
	foreach my $file (@fileHolder){
	    if (-e $cutadaptFile{$file}){
		if (($cutadaptFile{$file} =~/\.gz$/) || ($cutadaptFile{$file}=~/\.Z$/)){
		    $cmd="gunzip -c $cutadaptFile{$file} >> $cutFile";
		}
		else{
		    $cmd="cat $cutadaptFile{$file} >> $cutFile";
		}
	    }
	    elsif (($file =~/\.gz$/) || ($file=~/\.Z$/)){
		$cmd = "gunzip -c $file >> $cutFile";
	    }
	    else{
		$cmd = "cat $file >> $cutFile";
	    }
	    print LOG "# Doing: $cmd\n" if ($verbose);
	    system("$cmd");
	}
    }
    system"gzip $cutFile";
    $cutFile="$cutFile.gz";
    if (($cleanupCutadapt) || ($slimMode)){
	push(@rmFiles, $cutFile);
    }
    if (! defined($Getopt::Std::opt_S)){
	print $fh_stat "Total\tFastq files\tReads read\tReads after trimming\tAfter cleaning\n";
	print $fh_stat "Total\t$fastqFiles\t$sumReadsBefore\t$sumReadsAfter\n";
    }
    else{
	print $fh_stat "Total\tFastq files\tReads read\n";
	print $fh_stat "Total\t$fastqFiles\t$sumReadsBefore\n";    
    }
    if ($wwwMode){
	open(FILESTAT2,">$workDir/misc/filestat2.tab");
	if (! defined($Getopt::Std::opt_S)){
	    print FILESTAT2 "#Fastq files\tReads read\tReads after trimming\tAfter cleaning\n";
	    print FILESTAT2 "$fastqFiles\t$sumReadsBefore\t$sumReadsAfter\t";
	}
	else{
	    print FILESTAT2 "#Fastq files\tReads read\tAfter cleaning\n";
	    print FILESTAT2 "$fastqFiles\t$sumReadsBefore\t";
	}
    }
    &mapping($cutFile);
}
else{
    print LOG "Empty fileholders - I die now\n";
    print LOG "Done!\n";
    die;
}

############################################################################
#
# Make overall abundance statistics
#
############################################################################
open(DBCount,">$workDir/misc/dbCount.tab");
print DBCount "#Database\tVersion\tNo of sequence\tNo of nucleotides\n";


print $fh_stat "# All mapped reads\n";
print LOG "mapNotPhiXReads = $biological_relevant_reads used as reference to max possible reads that can map to a database\n" if ($verbose);

my $version_databases="$workDir/misc/version.databases";;
open(VERSION_DATABASES,">$version_databases");
my $file;
my $reads;
foreach my $key (@databases_all){
    if (exists($db_version{$key})){
        print VERSION_DATABASES "$key\t$db_version{$key}\n";
    }
    print DBCount "$key\t$db_version{$key}\t$countEntries{$key}{sequences}\t$countEntries{$key}{nuc}\n";
}
close(VERSION_DATABASES);
close(DBCount);


#
# add annotations to stat files
#

#
# KyotoCabinet not installed on cge, so use grep instead
#
#&add_taxonomy;
&add_taxonomy;

#
# Join nucleotide statistics file with read abundance statistics file
#
my $xlsWrite_cmd = "$prog_python $prog_xlsWrite $workDir/$sampleName.xlsx $workDir/abundance.databases.txt";
my $xlsWrite_negative_cmd = "$prog_python $prog_xlsWrite $workDir/$sampleName.negative.xlsx ";

&excel_and_post_processing(); 

my $unmappedBam = "$workDir/bam/unmapped.bam";
&unmapped_SE();
my $unmapped=&bamCountReads($unmappedBam);
push(@rmFiles, $unmappedBam);

#
# printout unmapped bam file
#
&printout_unmapped();


#
# make abundance.databases.txt
#
my $overall_abundance="$workDir/abundance.databases.txt";
open(AB,">$overall_abundance");
&print_overall_abundance(\*AB, $unmapped);
close(AB);

#
# Also print abundance databases info to STAT (log/MGmapper.summary);
# the STAT filehandle is already open
&print_overall_abundance($fh_stat, $unmapped);


if ($wwwMode){
    open(FILESTAT3,">$workDir/misc/filestat3.tab");
    print FILESTAT3 "#Mapping mode\tDatabase\tPercent mapped\tNo reads\n";
    &print_overall_abundance(\*FILESTAT3, $unmapped);
    close(FILESTAT3);

    #
    # Make snippet files with only 5 best abundances for all databases - only at strain level
    # For historical reasons outputFile Chainmode.tab is kept allthougt it is Bestmode
    #
    if ($#databases_chainmode >=0){
	MGmapper::make_www_tables($workDir, "$workDir/misc/Chainmode.tab", @databases_chainmode);
    }
    if ($#databases_fullmode >=0){
	MGmapper::make_www_tables($workDir, "$workDir/misc/Fullmode.tab", @databases_fullmode);
    }
}


print LOG "\n# $xlsWrite_cmd\n" if ($verbose);
print LOG "\n# $xlsWrite_negative_cmd\n" if ($verbose);

if (-e "$workDir/$sampleName.xlsx"){
    system("rm $workDir/$sampleName.xlsx");
}
if (-e "$workDir/$sampleName.negative_hits.xlsx"){
    system("rm $workDir/$sampleName.negative_hits.xlsx");
}
$datestring = localtime();
print $fh_runtime "$datestring\n";
print $fh_runtime "# $xlsWrite_cmd\n";
print $fh_runtime "# $xlsWrite_negative_cmd\n";
system("$xlsWrite_cmd");
system("$xlsWrite_negative_cmd");
#
# Make a header file that explains columns in $database.summary files
#
MGmapper::print_README("$workDir/Readme", $Getopt::Std::opt_R);

&cleanup_files();

#
# gzip all fastq files
#
$cmd="$prog_pigz -p $cores $workDir/misc/*.fq $workDir/cutadapt/*.cutadapt";
print LOG "Doing: $cmd\n" if ($verbose);
system("$cmd");

$datestring = localtime();
if ($wwwMode){
    print TIME "$datestring\n";
    close(TIME);
}
print $fh_stat "## Local date and time $datestring - done\n";
if ($remove_PCR_duplicates){
    print $fh_stat ("PCR duplicates - excluded from unmapped fastq files(s): $markdup_count\n");
}
print $fh_stat "$datestring\n";
print $fh_stat "End main programe\n";
close($fh_stat);

print $fh_runtime "$datestring\n";
print $fh_runtime "End main programe\n";
close($fh_runtime);
if ($wwwMode){
    print STDERR "Done!\n";
}
if (defined($Getopt::Std::opt_1)){
    print LOG "Executing unix command: '$unix_cmd1'\n" if ($verbose);
    system("$unix_cmd1");
}
if ($verbose){
    close(LOG);
}

exit(0);
###############################################################################
#
# Main program end
#
###############################################################################
sub run_cutadapt{
    my $first=1;
    my @commands=();
    
    my $command;    
    my $command_part1 = &make_cutadapt_command();
    my $command_part2;
    
    foreach my $file (@fileHolder){
	if ($first){
	    $datestring = localtime();	    
	    print $fh_stat "## Local date and time $datestring - cutadapt start\n";
	    $first=0;
	}
	
	$command_part2 = &cutadapt_SE($file);
	$command = "$command_part1" . "$command_part2";
	
	if ($command =~ m/-o (\S+) /){
	    $cutadaptFile{$file}=$1;
	}
	else{
	    print LOG "problem extracting output name from cutadapt command:\n$command\n";
	    exit;
	}
	
	if (-e "${cutadaptFile{$file}}.gz"){
	    $cutadaptFile{$file}= "${cutadaptFile{$file}}.gz";
	}
	
	if ($command ne ''){
	    if ((! -e $cutadaptFile{$file}) || (-z $cutadaptFile{$file})){
		push(@commands,$command);
	    }
	}
    }
    

    if (! -e $cleanFile){
	foreach my $command (@commands){
	    print LOG "$command\n" if ($verbose);
	    $cm->start("$command");		
	}
	#
	# Wait for all jobs to finish
	#
	$cm->wait_all_children;
	print $fh_runtime "Commands run in parallel:\n";
	foreach my $id (@commands){
	    print $fh_runtime "$id\n";
	}
	$datestring = localtime();
	print $fh_runtime "$datestring\n";
    }
    foreach my $file (@fileHolder){
	$fastqFiles++;
	if (exists($readCount{$cutadaptFile{$file}})){
	    $readsAfter{$file}=$readCount{$cutadaptFile{$file}};
	}
	else{
	    open(COUNT,">>$readCountFile");
	    $readsAfter{$file}=MGmapper::countReads($cutadaptFile{$file}, \*LOG, $verbose);
	    print COUNT "$cutadaptFile{$file}\t$readsAfter{$file}\n";
	    close(COUNT);
	}
	
	$sumReadsAfter += $readsAfter{$file};
	
	if ($wwwMode){
	    my @FS=split(/\//,$file);
	    print $fh_filestat "$fastqFiles\t$FS[-1]\t$readsBefore{$file}\t$readsAfter{$file}\n";	    
	}
	print $fh_stat "$fastqFiles\t$file\t$readsBefore{$file}\t$readsAfter{$file}\n";
    }
    return(0);
}

sub make_cutadapt_command{
    my $cutadapt_cmd = "$prog_cutadapt --quality-base=$QUALITY_BASE -f fastq -q $qualityCut -m $minLen";
    if (! -e "$parameterFile_adapter"){
	print LOG "can't open file $parameterFile_adapter\n";
	print LOG "Done!\n";
	die;
    }
    open(ADAPTER,"<$parameterFile_adapter") || die ("can't open file $parameterFile_adapter: $!");
    print LOG "Adapters removed by cutadapt:\n" if ($verbose);
    my $i=0;
    my $adapter;
    while (defined($_=<ADAPTER>)){
	if (/^\#/){
	    next;
	}
	chomp;
	if (m/(\S+)/){
	    $adapter = $1;
	}
	else{
	    next;
	}
	my $len=length($adapter);
	if ($len >1){
	    $i++;
	    print LOG "$i\t$adapter\n" if ($verbose);
	    $cutadapt_cmd .= " -b $adapter";
	}
    }
    $cutadapt_cmd .= ' ';
    print LOG "\n" if ($verbose);
    close(ADAPTER);
    return($cutadapt_cmd);
}

sub printout_unmapped{
    if ((! $slimMode) && ($unmapped >0)){
	my $readsUnmapped = "$workDir/misc/unmapped.fq";
	my $cmd = "$prog_bedtools bamtofastq -i $unmappedBam -fq $readsUnmapped";
	print LOG "# Doing: $cmd\n";
	system($cmd);
	
	print $fh_stat ("File with unmapped reads : $readsUnmapped\n");
	system("gzip $readsUnmapped");
	print $fh_stat ("Number of unmapped reads : $unmapped\n\n");
    }
    return(0);
}

sub excel_and_post_processing{
    foreach my $db (@databases_all){
	my $annot = "$workDir/misc/stat.$db.annot";
	my $strain_positive = "$workDir/stat/positive.strain.$db.txt";
	my $strain_negative = "$workDir/stat/negative.strain.$db.txt";
	
	my $strain_positive_noHeader = "$workDir/misc/positive.strain.$db.txt";
	my $strain_negative_noHeader = "$workDir/misc/negative.strain.$db.txt";
	
	my $strain_log = "$workDir/stat/classify.strain.$db.log";
	
	my $species_positive = "$workDir/stat/positive.species.$db.txt";
	my $species_negative = "$workDir/stat/negative.species.$db.txt";
	
	my $species_positive_noHeader = "$workDir/misc/positive.species.$db.txt";
	my $species_negative_noHeader = "$workDir/misc/negative.species.$db.txt";
	
	my $species_log = "$workDir/stat/classify.species.$db.log";
	
	
	$cmd = "$prog_MGmapper_classify -i $annot -c strain -a $min_relative_abundance -r $min_readCountUniq_ratio -m $max_mismatch_ratio -n $min_readCount -g $min_coverage -o $strain_positive -f $strain_negative -v -l $strain_log -s -H -N $biological_relevant_reads";
	
	print LOG "Doing: $cmd\n" if ($verbose);
	print $fh_runtime "# $cmd\n";
	system("$cmd");
	
	#
	# Remove header for the excel sheet
	#
	system("grep -v '^#' $strain_positive > $strain_positive_noHeader");
	system("grep -v '^#' $strain_negative > $strain_negative_noHeader");
	
	if (-z $strain_positive_noHeader){
	    print LOG "removing empty file: $strain_positive\n" if ($verbose);
	    push(@rmFiles,$strain_positive);
	}
	push(@rmFiles,$strain_positive_noHeader);
	
	if (-z $strain_negative_noHeader){
	    print LOG "removing empty file: $strain_negative\n" if ($verbose);
	    push(@rmFiles,$strain_negative);
	}
	push(@rmFiles,$strain_negative_noHeader);
	
	$cmd = "$prog_MGmapper_classify -i $annot -c species -a $min_relative_abundance -r $min_readCountUniq_ratio -m $max_mismatch_ratio -n $min_readCount -g $min_coverage -o $species_positive -v -l $species_log -s -H -f $species_negative -N $biological_relevant_reads";
	print LOG "Doing: $cmd\n" if ($verbose);
	print $fh_runtime "# $cmd\n";
	system("$cmd");
	
        #
        # now find the lowest accepted abundance and rerun
        #
        my $new_lowest_abundance = $min_relative_abundance;
        open(TMP,"<$species_log");
	my $species_count=0;
	while (defined ($_=<TMP>)){
	    chomp;
	    if (/^# Minimum accepted abundance value=/){
		my @w=split(/\t/);
		$new_lowest_abundance = $w[1];
	    }
	    if (/# Number of hits accepted:/){
		my @w=split(/\t/);
		$species_count= $w[1];
	    }
	}
        close(TMP);


	#
	# if any species was found with default parameters as above, then rerun with -a = lowest abundance among those hits and no uniq_read_ratio i.e. -r 0
	# This will find more species with abundance higher than what was previoiusly found, but no uniq_read_ratio criteria
	#
	if ($species_count >0){
	    #
	    # not pretty to subtract the very small number, but needed due to rounding off problems 
	    #
	    $new_lowest_abundance -= 0.0000000001;
	    $cmd = "$prog_MGmapper_classify -i $annot -c species -a $new_lowest_abundance -r 0 -m $max_mismatch_ratio -n $min_readCount -g $min_coverage -o $species_positive -f $species_negative -v -l $species_log -s -H -N $biological_relevant_reads";

	    print LOG "Doing: $cmd\n" if ($verbose);
	    print $fh_runtime "# $cmd\n";
	    system("$cmd");
	}
	    
	#
	# Remove header for the excel sheet
	#
	system("grep -v '^#' $species_positive > $species_positive_noHeader");
	system("grep -v '^#' $species_negative > $species_negative_noHeader");
	
	if (-z $species_positive_noHeader){
	    print LOG "removing empty file: $species_positive\n" if ($verbose);
	    push(@rmFiles,$species_positive);
	}
	push(@rmFiles,$species_positive_noHeader);
	
	if (-z $species_negative_noHeader){
	    print LOG "removing empty file: $species_negative\n" if ($verbose);
	    push(@rmFiles,$species_negative);
	}
	push(@rmFiles,$species_negative_noHeader);
	
	
	if ((! -z $strain_positive_noHeader) && (-e $strain_positive_noHeader)){
	    $xlsWrite_cmd .= " $strain_positive_noHeader";
	}
	if ((! -z $species_positive_noHeader) && (-e $species_positive_noHeader)){
	    $xlsWrite_cmd .= " $species_positive_noHeader";
	}
	
	if ((! -z $strain_negative_noHeader) && (-e $strain_negative_noHeader)){
	    $xlsWrite_negative_cmd .= " $strain_negative_noHeader";
	}
	if ((! -z $species_negative_noHeader) && (-e $species_negative_noHeader)){
	    $xlsWrite_negative_cmd .= " $species_negative_noHeader";
	}
    }
    return(0);
}

sub print_overall_abundance{
    my ($fh, $unmapped) = @_;

    printf $fh "Fullmode\tnotPhiX\t100.00\t%d\n",$biological_relevant_reads;

    my $file;
    my $reads;
    my $perc;
    foreach my $key (@databases_all){
	$file = "$workDir/bam/mapTo.$key.bam";
	$reads = &bamCountReads($file);
	$perc = $reads*100/$biological_relevant_reads;
	printf $fh "%s\t%s\t%7.3f\t%d\n",$runningMode{$key},$key,$perc,$reads;
    }
    my $str="Unmapped";
    $perc = $unmapped*100/$biological_relevant_reads;
    printf $fh "-\t%s\t%7.3f\t%d\n",$str,$perc,$unmapped;
    return(0);
}



sub cleanup_files{
    if ($cleanupCutadapt){
        system("rm -r $workDir/cutadapt");
    }
    if ($cleanup){
        my $cmd;
        if (! -e "$noclean"){
            foreach my $id (@rmFiles){
                if (-e $id){
                    $cmd = "rm $id";
                    print LOG "Doing: $cmd\n" if ($verbose);
                    system("$cmd");
                }
            }
        }
    }
    return(0);
}


sub unmapped_SE{
    my %rec=();
    my $str;
    my $revcomp;
    my $hash_str_F;
    my $hash_str_R;
    my @w=();
    my $characterLen=30;
    foreach my $db (@databases_all){
        my $file = "$workDir/bam/mapTo.$db.bam";
        open (FILE,"$prog_samtools view $file | cut -f1,10 |");
        while (defined($_=<FILE>)){
            chomp;
            @w=split(/\s+/);
            $str=$w[1];
            $revcomp = reverse($str);
            $revcomp =~ tr/ACGTacgt/TGCAtgca/;

            $hash_str_F = $w[0] . $str;
            $hash_str_R = $w[0] . $revcomp;
            $rec{$hash_str_F}=1;
            $rec{$hash_str_R}=1;
        }
        close(FILE);
    }

    if ($remove_PCR_duplicates){
	foreach my $file (@markdup_all){
	    open (FILE,"$prog_samtools view $file | cut -f1,10 |");
	    while (defined($_=<FILE>)){
		chomp;
		@w=split(/\s+/);
		$str=$w[1];
		$revcomp = reverse($str);
		$revcomp =~ tr/ACGTacgt/TGCAtgca/;
		$hash_str_F = $w[0] . $str;
		$hash_str_R = $w[0] . $revcomp;
		$rec{$hash_str_F}=2;
		$rec{$hash_str_R}=2;
		$markdup_count++;
	    }
	    close(FILE);
	}
	print LOG "Total number of PCR duplicates: $markdup_count\n";
    }

    open(BamIn,"$prog_samtools view $cleanFile |");
    open(BamOut,"| $prog_samtools view -h -Sb - > $unmappedBam");


    while (defined($_=<BamIn>)){
        @w=split(/\s+/);
        $str=$w[9];
        $revcomp = reverse($str);
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;

        $hash_str_F = $w[0] . $str;
        $hash_str_R = $w[0] . $revcomp;
        if ((! exists($rec{$hash_str_F})) && (! exists($rec{$hash_str_R}))){
            print BamOut "$_";
        }
    }

    close(BamIn);
    close(BamOut);
    return(0);
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

sub parseOptions{
    if (defined($Getopt::Std::opt_Q)){
	$parameterFile_adapter = $Getopt::Std::opt_Q;
    }

    if (defined($Getopt::Std::opt_p)){
	$remove_PCR_duplicates = 1;
    }

    #
    # execute a unix command after everything else has finished
    #
    if (defined($Getopt::Std::opt_1)){
	$unix_cmd1=$Getopt::Std::opt_1;
    }
    if (defined($Getopt::Std::opt_D)){
	$min_readCountUniq_ratio = $Getopt::Std::opt_D;
    }
    if (defined($Getopt::Std::opt_L)){
	$max_mismatch_ratio = $Getopt::Std::opt_L;
    }
#
# A working directory for all output files
#
    if (defined($Getopt::Std::opt_d)){
	$workDir = $Getopt::Std::opt_d;
    }
    
    if (! -d $workDir){
	system ("mkdir -p $workDir; chmod 0755 $workDir");
    }
    if (! -d "$workDir/bam"){
	system("mkdir -p $workDir/bam; chmod 0755 $workDir/bam");
    }
    if (! -d "$workDir/log"){
	system("mkdir -p $workDir/log; chmod 0755 $workDir/log");
    }
    if (! -d "$workDir/misc"){
	system("mkdir -p $workDir/misc; chmod 0755 $workDir/misc");
    }
    if (! -d "$workDir/stat"){
	system("mkdir -p $workDir/stat; chmod 0755 $workDir/stat");
    }
    
    if (defined($Getopt::Std::opt_A)){
	$alignment_score_min = $Getopt::Std::opt_A;
    }
    if (defined($Getopt::Std::opt_N)){
	$max_edit_distance = $Getopt::Std::opt_N;
    }
    if (defined($Getopt::Std::opt_r)){
	$min_relative_abundance = $Getopt::Std::opt_r;
    }
    if (defined($Getopt::Std::opt_U)){
	$min_readCount = $Getopt::Std::opt_U;
    }

    
    if ($verbose){
	open(LOG,">$workDir/log/MGmapper.log");
    }
    if (defined($Getopt::Std::opt_w)){
	$wwwMode=1;
    }
    if (defined($Getopt::Std::opt_s)){
	$SNP_threshold=$Getopt::Std::opt_s;
    }
    
    if (defined($Getopt::Std::opt_n)){
	$sampleName=$Getopt::Std::opt_n;
    }
    
    if (defined($Getopt::Std::opt_c)){
	$cores = $Getopt::Std::opt_c;
    }
    if (defined($Getopt::Std::opt_M)){
	$slimMode = 1;
    }
    $cm = new ChildManager($cores);
    
    if (defined($Getopt::Std::opt_i)){
	if (! -e $Getopt::Std::opt_i){
	    print LOG "File not found: $Getopt::Std::opt_i\n";
	    print STDERR "File not found: $iGetopt::Std::opt_i\nDone!\n";
	    exit;
	}
    }
    elsif (defined($Getopt::Std::opt_f)){
	open(LIST,"<$Getopt::Std::opt_f");
    }
    else{
	print STDERR "No single input file -i or list of files -f was defined\n";
	exit;
    }
    
#
# Quality score q in cutadapt 
#
    if (defined($Getopt::Std::opt_q)){
	$qualityCut = $Getopt::Std::opt_q;
    }
#
# minimum length of read after cutadapt
#
    if (defined($Getopt::Std::opt_m)){
	$minLen = $Getopt::Std::opt_m;
    }
#
# option to activate option -t, default is on
#
    if (defined($Getopt::Std::opt_T)){
	$matchRatioFlag=0;
    }
    if (defined($Getopt::Std::opt_e)){
	$absoluteMatchCount=$Getopt::Std::opt_e;
    }
    
#
# activate option -e, default is off
#
    if (defined($Getopt::Std::opt_E)){
	$absoluteMatchCountFlag=1;
    }
#
# Databases to map against phiX
#
    if (defined($Getopt::Std::opt_F)){
	@databases_fullmode=();
	my @choices = split(/\,/,$Getopt::Std::opt_F);
	my $i=0;
	while (defined $choices[$i]){
	    if ($choices[$i] eq '0'){
		@databases_fullmode=();
		$i++;
		next;
	    }
	    if (! exists $databaseNames{$choices[$i]}){	    
		print "Error with options to -F\n";
		print STDERR "Error with options to -F\nDone!";
		exit;
	    }
	    $runningMode{$databaseNames{$choices[$i]}}='Fullmode';
	    push(@databases_fullmode,$databaseNames{$choices[$i]});
	    my $j=0;
	    while (defined($databases_all[$j])){
		if ($databases_all[$j] eq $databaseNames{$choices[$i]}){
		    print LOG "Database $databaseNames{$choices[$i]} has been selected more than 1 time - i exit now\n";
		    print STDERR "Done!\n";
		    exit;
		}
		$j++;
	    }
	    push(@databases_all,$databaseNames{$choices[$i]});
	    $i++;
	}
    }
    
#
# Databases to map against in best mode
#
    if (defined($Getopt::Std::opt_C)){
	@databases_chainmode=();
	my @choices = split(/\,/,$Getopt::Std::opt_C);
	my $i=0;
	while (defined $choices[$i]){
	    if ($choices[$i] eq '0'){
		@databases_chainmode=();
		$i++;
		next;
	    }
	    if (! exists $databaseNames{$choices[$i]}){
		print LOG "Error with options to -C\n";
		print STDERR "Error with options to -C\nDone!\n";
		exit;
	    }
	    $runningMode{$databaseNames{$choices[$i]}}='Bestmode';
	    push(@databases_chainmode,$databaseNames{$choices[$i]});
	    my $j=0;
	    while (defined($databases_all[$j])){
		if ($databases_all[$j] eq $databaseNames{$choices[$i]}){
		    print LOG "Database $databaseNames{$choices[$i]} has been selected more than 1 time - i exit now\n";
		    print STDERR "Done!\n";
		    exit;
		}
		$j++;
	    }
	    push(@databases_all,$databaseNames{$choices[$i]});
	    $i++;
	}
    }
    
    if (defined($Getopt::Std::opt_t)){
	$matchRatio=$Getopt::Std::opt_t;
    }
    
    if (defined($Getopt::Std::opt_k)){
	$cleanup=0;
    }
    if (defined($Getopt::Std::opt_K)){
	$cleanupCutadapt=0;
    }
    if (defined($Getopt::Std::opt_B)){
	$QUALITY_BASE=$Getopt::Std::opt_B;
    }
    
    print LOG ("# command: $command\n\n");
    print LOG ("# Current directory: $origDir\n");
    if (defined($Getopt::Std::opt_i)){
	print LOG "# -i: file with reads: $Getopt::Std::opt_i\n";
    }
    
    if ($#databases_fullmode >= 0){
	print LOG "# -F: Full mode: map reads against these databases: @databases_fullmode\n";
    }
    if ($#databases_chainmode >= 0){
	print LOG "# -C: Best mode: map reads against these databases: @databases_chainmode\n";
    }
    if ($verbose){
	print LOG "# -c: Number of cores: $cores\n";
	print LOG "# -t: Matches\/readLength: $matchRatio\n";
	print LOG "# -v: logfile: $workDir/MGmapper.log\n";
    }
    if (defined($Getopt::Std::opt_S)){
	print LOG "# -S: Skip running cutadapt\n";
    }
    
    print LOG "# -d: Output directory: $workDir\n\n";
    
    if (defined($Getopt::Std::opt_a)){
	if (! -d "$workDir/fasta"){
	    system("mkdir -p $workDir/fasta; chmod 755 $workDir/fasta");
	}
	my @w=split(/\,/,$Getopt::Std::opt_a);
	my $i=0;
	while (defined($w[$i])){
	    if (! exists($databaseNames{$w[$i]})){
		print STDERR "Assembled fasta files requested for an unexisting database referenced with: '$w[$i]'\n";
		die;
	    }
	    $fastaOut{$databaseNames{$w[$i]}}=1;
	    print LOG "fasta: $w[$i] $databaseNames{$w[$i]}\n" if ($verbose);
	    $i++;
	}
    }
    return(0);
}

sub readParameterFile{

    if (defined($Getopt::Std::opt_P)){
	$parameterFile = $Getopt::Std::opt_P;
    }
    if (! -e "$parameterFile"){
	print LOG "can't open file $parameterFile\n";
	print LOG "Done!\n";
	die;
    }

    open(PARAMETER,"<$parameterFile");
    
    while (defined($_=<PARAMETER>)){
	if (/^\#/){
	    next;
	}
	chomp;
	my @w=split(/\s+/);
	my $file = $w[0];
	my $id = $w[1];
	$databaseNames{$databaseCounter} = $id;
	$databaseCount{$id} = $databaseCounter;
	push(@dbNames,$id);
	
	my $i=1;
	my $splitFile = $file . ".$i";
	my $flag=1;
	#
	# Big fasta files can be split into smaller chunks and called with a suffix being 1 2 3 ... N i.e. Bacteria.1 Bacteria.2  
	# Each fasta file must of-course be indexed with both bwa index and samtools faidx
	#
	if ((! -e $file) && (! -e $splitFile)){
	    die "\nError in file 'parameter file with reference to database' $parameterFile\nFile not found '$file\n";
	}

	if (-e $splitFile){
	    while ($flag){
		my $split = "$file" . ".$i";
		if (-e $split){
		    $idx{index}{$id}[$i] = $split;
		    $idx{ann}{$id}[$i] = "$split" . '.ann';
		}
		else{
		    $flag=0;
		}
		$i++;
	    }
	}
	else{
	    $idx{index}{$id}[$i] = $file;
	    $idx{ann}{$id}[$i] = "$file" . '.ann';	
	}
	my $dbNum = $databaseCount{$id};
	my $dbName = $databaseNames{$dbNum};
	$idx{dbNum}{$dbName}=$dbNum;
	$idx{dbName}{$dbNum}=$dbName;
	
	#
	# info in file like Bacteria.info
	#
	my $infoFile = "$file.info";
	$info{date}='undefined';
	if (-e $infoFile){
	    open(INFO,"<$infoFile");
	    while (defined($_=<INFO>)){
		chomp;
		if (m/^Date:/){
		    $info{date}=$_;
		}
#		if (m/^Sequences:/){
#		    $info{sequences}=$_;
#		}
		if (m/^Remark:/){
		    $info{remark}=$_;
		}
#		if (m/^Filesize:/){
#		    $info{filesize}=$_;
#		}
		if (m/^Species:/){
		    $info{species}=$_;
		}
	    }
	    close(INFO);
	    if (exists ($info{date})){
		$idx{desc}{$dbNum} .= " $info{date};";
                my @t=split(/\s+/,$info{date});
                $db_version{$id} = $t[1];
            }
            else{
                $db_version{$id} = 'undefined';
            }
	    
	    if (exists ($info{species})){
		$idx{desc}{$dbNum} .= " $info{species};";
	    }
	    if (exists ($info{sequences})){
		$idx{desc}{$dbNum} .= " $info{sequences};";
	    }
	    if (exists ($info{filesize})){
		$idx{desc}{$dbNum} .= " $info{filesize};";
	    }
	    if (exists ($info{remark})){
		$idx{desc}{$dbNum} .= " $info{remark};";
	    }
	    undef %info;
	}
	else{
	    my $i=2;
	    $idx{desc}{$dbNum} = ' ';
	    while (defined($w[$i])){
		$idx{desc}{$dbNum} .= "$w[$i] ";
		$i++;
	    }
	}
	
	#
	# name of kyotocabinet file
	#
	my $kch = "$file.kch";
	if (-e $kch){
	    $idx{kch}{$dbName} = $kch;
	}
	
	my $tax = "$file.tax";
        if (-e $tax){
            $idx{tax}{$dbName} = $tax;
        }

	#
	# A logfile showing redundant entries that were removed from the database
	#
	my $fastauniq = $file . '.fastauniq.log';
	if (-e $fastauniq){
	    $idx{fastauniq}{$dbName} = $fastauniq;
	}
	$databaseCounter++;
    }
    $databaseCounter -= 1;
    close(PARAMETER);
    return(0);
}
sub add_taxonomy{
    use KyotoCabinet;
    my $value;
    my $seqId;
    my $kch;
    my $db;

    my $file;
    my $stat;
    print LOG "Adding annotations to all misc/stat.db.raw files\n" if ($verbose);

    my $tax_unknown = '';
    for (my $i=1;$i<=16;$i++){
        $tax_unknown .= "\tUnknown";
    }

    my $kch_flag=1;
    foreach my $database (@databases_all){
        $file = "$workDir/misc/stat.$database.raw";
        $stat = "$workDir/misc/stat.$database.annot";

        if (exists($idx{kch}{$database})){
            $db = new KyotoCabinet::DB;
            $kch="$idx{kch}{$database}";
            #
            # open the database
            #
            if (!$db->open("$kch", $db->OREADER)) {
                printf STDERR ("open error: %s\n", $db->error);
            }
            $kch_flag=1;
            print LOG "Using annotaiotion database: $kch\n" if ($verbose);
        }
        else{
            $kch_flag=0;
        }
        open(FILE,"<$file");
        open(TMP,">$stat");

        my @w=();
        while (defined($_=<FILE>)){
            chomp;
            if ($kch_flag){
                @w=split(/\s+/);
                $seqId=$w[0];
                $value = $db->get("$seqId");
                if (defined($value)){
                    print TMP "$database\t$_\t$value\n";
                }
            }
            else{
                print TMP "$database\t$_";
                print TMP "$tax_unknown\n";
                print LOG "No annotation found for: '$seqId'\n" if ($verbose);
            }
        }
        close(FILE);
        close(TMP);

        if ($kch_flag){
            if (!$db->close){
                printf LOG ("close error: %s\n", $db->error);
            }
        }
    }
}


sub cutadapt_SE{
    my ($fastq) = @_;
    #
    # Identify input file for forward reads 
    #
    my $pos1 = rindex($fastq,'/') +1;
    my $name = substr $fastq , $pos1;

    my $cutadaptDir = "$workDir/cutadapt";
    if (! -d $cutadaptDir){
        system("mkdir $cutadaptDir");
    }

    #
    # remove adapter, trimm for forward reads
    #
    my $cutAdaptFile = "$cutadaptDir/$name" . '.cutadapt';
    my $cutAdaptLog = "$cutadaptDir/$name" . '.cutadapt.log';
    my $gz = $cutAdaptFile . '.gz';
    if (-e $gz){
	$cutAdaptFile = $gz;
    }
    
    if ($cleanupCutadapt){
	push(@rmFiles,$cutAdaptFile);
    }
    my $cmd = " -o $cutAdaptFile $fastq  > $cutAdaptLog";

    return("$cmd");
}


sub mapping{
    my ($readF) = @_;
    #
    # bwa against phiX
    #
    $cleanFile = "$workDir/bam/cleaned.nophiX.bam";

    my $sortBamFile;
    my $summaryFile;

    if (! -e $cleanFile){
	print LOG "Calling cleanData with parameters: $readF , $cleanFile\n" if ($verbose);
	$datestring = localtime();
	print $fh_stat "## Local date and time $datestring - start mapping against phiX\n";
	cleanData($readF, $cleanFile);
    }
    else{
	print LOG "Skipping 'cleanData' as file exists: $cleanFile\n" if ($verbose);
	$biological_relevant_reads = &bamCountReads($cleanFile);
    }
    
    if ($slimMode){
        push(@rmFiles, $cleanFile);
    }

    print $fh_stat "Reads not mapped to db= phiX:\t$biological_relevant_reads\n\n";

    if ($wwwMode){
	print FILESTAT2 "$reads\n";
	close(FILESTAT2);
    }
    #
    # map reads in Full mode against all databases
    #
    my $size_databases_all = $#databases_all;
    
    foreach my $db (@databases_all){
	#
	# Start mapping
	#
	&mapToDatabases($cleanFile,$db);
    }

    my @filter_hit_commands=();

    foreach my $db (@databases_all){
	my @list=();
	@list=&post_process_single_database_cmd($db);
	push(@filter_hit_commands, @list);
    }

    #
    # clean bam files for low alignment scores and low match+mismatches
    # output are mapFilter bam files
    #
    foreach my $command (@filter_hit_commands){
	print LOG "$command\n";
	$cm->start("$command");		
    }
    $cm->wait_all_children;
    if (@filter_hit_commands){
	print $fh_runtime "Commands run in parallel:\n";
	foreach my $id (@filter_hit_commands){
	    print $fh_runtime "# $id\n";
	}
	$datestring = localtime();
	print $fh_runtime "$datestring\n";
    }

    #
    # make mapTo bam files from mapFilter files
    # where best hits based on alignment scores are used to
    # determine in which database specific bam file a read and its pair should be 
    #

    &best_hit;
    print $fh_runtime "Command: best_hit\n";
    $datestring = localtime();
    print $fh_runtime "$datestring\n";

    my @sortBam_commands=();
    my @pileup_commands=();
    my @pileupToDepth_commands=();
    my @pileupToFasta_commands=();
    my @sortByName_commands=();
    my @fixmate_commands=();
    my @markdup_commands=();
    my @markdupRemove_commands=();

    foreach my $db (@databases_all){
	my $i=1;
	my $annot= "$workDir/misc/stat.$db.annot";

	if ((! -e $annot) || (-z $annot)){
	    while (-e "$workDir/bam/mapTo.$db.$i.bam"){
		my $mapTo = "$workDir/bam/mapTo.$db.$i.bam";
		my $mapSort = "$workDir/bam/mapTo.$db.$i.sort.bam";

		#
		# filenames only used in case of removal of PCR duplicates
		#
		my $mapSortByName = "$workDir/bam/mapTo.$db.$i.SortByName.bam";
		my $mapFixmate = "$workDir/bam/mapTo.$db.$i.fixmate.bam";
		my $mapFixmateSort = "$workDir/bam/mapTo.$db.$i.fixmate.sort.bam";
		my $mapMarkdup = "$workDir/bam/mapTo.$db.$i.markdup.bam";

		if ($slimMode){
		    push(@rmFiles, $mapTo);
		}

		my $pileup ="$workDir/misc/pileup.$db.$i";
		my $summaryFile ="$workDir/misc/stat.$db.raw.$i";

		if ($remove_PCR_duplicates){
		    #
		    # procedure is to do samtools fixmate, samtools markdup, remove duplicates and over-write previous mapSort filename
		    # to produce a sorted bam file without the PCR duplicates
		    #

		    # Sort file by readName
		    my $cmd = &samtoolsSortByName_cmd($mapTo, $mapSortByName);
		    push(@sortByName_commands, $cmd);
		    push(@rmFiles, $mapSortByName);

		    # fixmate i.e. add ms TAG to bamfile
		    $cmd = &samtoolsFixmate_cmd($mapSortByName, $mapFixmate);
		    push(@fixmate_commands, $cmd);
		    push(@rmFiles, $mapFixmate);

		    # sort fixmate bam file by coordinate
		    $cmd = &sortBam_cmd($mapFixmate, $mapFixmateSort);
		    push(@sortBam_commands, $cmd);
		    push(@rmFiles, $mapFixmateSort);

		    # mark PCR duplicates
		    $cmd = &samtoolsMarkdup_cmd($mapFixmateSort, $mapMarkdup);
		    push(@markdup_commands, $cmd);
		    push(@rmFiles, $mapMarkdup);

		    # mark PCR duplicates
		    $cmd = &markdupRemove_cmd($mapMarkdup, $mapSort);
		    push(@markdupRemove_commands, $cmd);

		    #
		    # save the markdup files such that unmapped reads do not include PCR duplicates
		    #
		    push(@markdup_all, $mapMarkdup);
		}
		else{
		    my $cmd = &sortBam_cmd($mapTo, $mapSort);
		    push(@sortBam_commands, $cmd);
		}
		$cmd = &bamToPileup_cmd($mapSort, $i, $db, $pileup);
		push(@pileup_commands, $cmd);

		$cmd = &pileupToDepth_cmd($pileup, $i, $db, $mapSort, $summaryFile);
		if (! -z $pileup){
		    push(@pileupToDepth_commands, $cmd);
		}
		
		if (exists($fastaOut{$db})){
		    my $fsa ="$workDir/fasta/$sampleName.$db.$i.fasta";
		    my $cmd = &pileupToFasta_cmd($pileup, $fsa);
		    push(@pileupToFasta_commands, $cmd);
		}
		$i++;
	    }
	}
    }

    if ($remove_PCR_duplicates){
	#
	# sort by name
	#
	if (@sortByName_commands){
	    foreach my $command (@sortByName_commands){
		print LOG "$command\n";
		$cm->start("$command");
	    }
	    $cm->wait_all_children;
	    print $fh_runtime "Commands run in parallel: samtools sort by name:\n";
	    foreach my $id (@sortByName_commands){
		print $fh_runtime "$id\n";
	    }
	    $datestring = localtime();
	    print $fh_runtime "$datestring\n";
	}

	#
	# Fixmate
	#
	if (@fixmate_commands){
	    foreach my $command (@fixmate_commands){
		print LOG "$command\n";
		$cm->start("$command");
	    }
	    $cm->wait_all_children;
	    print $fh_runtime "Commands run in parallel: samtools fixmate:\n";
	    foreach my $id (@fixmate_commands){
		print $fh_runtime "$id\n";
	    }
	    $datestring = localtime();
	    print $fh_runtime "$datestring\n";
	}

	#
	# samtools sort by coordinates
	#
	if (@sortBam_commands){
	    foreach my $command (@sortBam_commands){
		print LOG "$command\n";
		$cm->start("$command");
	    }
	    $cm->wait_all_children;
	    print $fh_runtime "Commands run in parallel: samtools sort:\n";
	    foreach my $id (@sortBam_commands){
		print $fh_runtime "$id\n";
	    }
	    $datestring = localtime();
	    print $fh_runtime "$datestring\n";
	}

	#
	# samtools markdup
	#
	if (@markdup_commands){
	    foreach my $command (@markdup_commands){
		print LOG "$command\n";
		$cm->start("$command");
	    }
	    $cm->wait_all_children;
	    print $fh_runtime "Commands run in parallel: samtools markdup:\n";
	    foreach my $id (@markdup_commands){
		print $fh_runtime "$id\n";
	    }
	    $datestring = localtime();
	    print $fh_runtime "$datestring\n";
	}

	#
	# remove duplicates
	#
	if (@markdupRemove_commands){
	    foreach my $command (@markdupRemove_commands){
		print LOG "$command\n";
		$cm->start("$command");
	    }
	    $cm->wait_all_children;
	    print $fh_runtime "Commands run in parallel: remove duplicates:\n";
	    foreach my $id (@markdupRemove_commands){
		print $fh_runtime "$id\n";
	    }
	    $datestring = localtime();
	    print $fh_runtime "$datestring\n";
	}
    }
    else{
	#
	# sort the mapTo bam files
	#
	if (@sortBam_commands){
	    foreach my $command (@sortBam_commands){
		print LOG "$command\n";
		$cm->start("$command");		
	    }
	    $cm->wait_all_children;
	    print $fh_runtime "Commands run in parallel:\n";
	    foreach my $id (@sortBam_commands){
		print $fh_runtime "# $id\n";
	    }
	    $datestring = localtime();
	    print $fh_runtime "$datestring\n";
	}
    }

    #
    # make pileup files via samtools mpileup
    #
    if (@pileup_commands){
	foreach my $command (@pileup_commands){
	    print LOG "$command\n";
	    $cm->start("$command");
	}
	$cm->wait_all_children;
	print $fh_runtime "Commands run in parallel:\n";
	foreach my $id (@pileup_commands){
	    print $fh_runtime "# $id\n";
	}
	$datestring = localtime();
	print $fh_runtime "$datestring\n";
    }

    #
    # make summary files i.e stat.$db.$i.raw
    #
    if (@pileupToDepth_commands){
	foreach my $command (@pileupToDepth_commands){
	    print LOG "$command\n";
	    $cm->start("$command");		
	}
	$cm->wait_all_children;
	print $fh_runtime "Commands run in parallel:\n";
	foreach my $id (@pileupToDepth_commands){
	    print $fh_runtime "# $id\n";
	}
	$datestring = localtime();
	print $fh_runtime "$datestring\n";
    }

    if (@pileupToFasta_commands){
	#
	# make fasta contig files
	#
	foreach my $command (@pileupToFasta_commands){
	    print LOG "$command\n";
	    $cm->start("$command");		
	}
	$cm->wait_all_children;
	print $fh_runtime "Commands run in parallel:\n";
	foreach my $id (@pileupToFasta_commands){
	    print $fh_runtime "# $id\n";
	}
	$datestring = localtime();
	print $fh_runtime "$datestring\n";
    }

    foreach my $db (@databases_all){
	my $mapToFinal = "$workDir/bam/mapTo.$db.bam";

        if ($slimMode){
            push(@rmFiles, $mapToFinal);
	}

	my $i=1;
	my @mapToList=();
	while (-e "$workDir/bam/mapTo.$db.$i.bam"){
	    my $mapTo = "$workDir/bam/mapTo.$db.$i.bam";
	    my $mapSort = "$workDir/bam/mapTo.$db.$i.sort.bam";

	    push(@mapToList,$mapSort);
	    push(@rmFiles,$mapTo);

            if ($slimMode){
                push(@rmFiles, $mapSort);
            }
	    $i++;
	}


	my $merge = "$workDir/bam/merge.$db.sam";
	&mergeBam(\@mapToList,$merge,$mapToFinal);

	my $summaryFileMain ="$workDir/misc/stat.$db" . ".raw";
	push(@rmFiles, $summaryFileMain);
	my $fsaMain ="$workDir/fasta/$sampleName.$db.fasta";
	&make_files($db, $summaryFileMain,$fsaMain,\@mapToList);
	
        #
	# Make readCount matrices
        #
        if (defined($Getopt::Std::opt_R)){
	    if (! -d "$workDir/matrix"){
		system("mkdir -p $workDir/matrix; chmod 755 $workDir/matrix");
	    }
            my $readCountMatrix ="$workDir/matrix/readCount.$db.matrix";
	    &make_readCount_matrix($readCountMatrix, $db, \@mapToList);
	}
    }
}

sub make_readCount_matrix{
    my ($matrix,$db, $ra_list)=@_;
    
    if (! -e "$matrix"){
	my $unsorted = $matrix . ".Unsorted";
	my $i=0;
	while (defined($ra_list->[$i])){
	    my $j=$i+1;
	    &readCountAnalysis($ra_list->[$i], $j, $db, $unsorted);
	    $i++;
	}
	#
	# sort the readCount matrix according to fastaEntry name
	#
	
	$cmd = "sort -k 1 $unsorted > $matrix";
	print LOG "# Doing: $cmd\n" if ($verbose);
	
	system($cmd);
	system("rm $unsorted");
    }
    return(0);
}



sub make_files{
    my ($db, $summaryFileMain, $fsaMain, $ra_mapToList) = @_;

    if ((! -e $summaryFileMain) || ((! -e $fsaMain) && (exists($fastaOut{$db}))) ){
	my @summaryFiles=();
	my @Fsa=();
	my $i=0;
	while (defined($ra_mapToList->[$i])){
	    my $j=$i+1;
	    my $summaryFile ="$workDir/misc/stat.$db.raw.$j";
	    push(@summaryFiles,$summaryFile);
	    
	    my $pileup ="$workDir/misc/pileup.$db.$j";
	    push(@rmFiles,$pileup);

	    if (! -z $pileup){
		#                                                                                                                                         
		# make assembled fasta files if requested - option -a                                                                                     
		#                                                                                                                                         
		if (exists($fastaOut{$db})){
		    my $fsa ="$workDir/fasta/$sampleName.$db.$j.fasta";
		    push(@Fsa,$fsa);
		}
	    }
	    $i++;
	}

	#
	# make assembled fasta files if requested - option -a
	#
	if (exists($fastaOut{$db})){
	    my $cmd='';
	    if ($#Fsa == 0){
		$cmd = "mv $Fsa[0] $fsaMain";
	    }
	    elsif ($#Fsa >0){
		my $i=0;
		$cmd = "cat";
		while (defined($Fsa[$i])){
		    $cmd .= " $Fsa[$i]";
		    $i++;
		}
		$cmd .= " > $fsaMain";
	    }
	    print LOG "Doing: $cmd\n" if ($verbose);
	    system("$cmd");
	    if (-z $fsaMain){
		print LOG "Empty file: $fsaMain - I remove it\n";
		push(@rmFiles, $fsaMain);
	    }
	    push(@rmFiles,@Fsa);
	}
	#
	# cat all summaryfiles for database db
	#
        my $cmd = "cat ";
        foreach my $id (@summaryFiles){
            chomp($id);
            $cmd .= "$id ";
        }
        $cmd .= " > $summaryFileMain";
        print LOG "# Doing: $cmd\n" if ($verbose);

        system($cmd);

	if (-z $summaryFileMain){
	    print LOG "$summaryFileMain is empty - I remove it\n";
	    system("rm $summaryFileMain");
	}
	if ($cleanup){
	    foreach my $id (@summaryFiles){
		if (-e $id){
		    system("rm $id");
		}
	    }
	}
    }
    return(0);
}



sub mergeBam{
    my ($ra, $sam, $bam, $symbolicLink)=@_;

    if (-e $sam){
        system("rm $sam");
    }
    my $i=0;
    my $cmd;
    while (defined($ra->[$i])){
        $cmd = "$prog_samtools view -H $ra->[$i] >> $sam";
        print LOG "# In mergeBam: $cmd\n" if ($verbose);

        system("$cmd");
        $i++;
    }
    $i=0;
    while (defined($ra->[$i])){
        $cmd = "$prog_samtools view $ra->[$i] >> $sam";
        print LOG "# In mergeBam: $cmd\n" if ($verbose);

        system("$cmd");
        $i++;
    }

    $cmd = "$prog_samtools view $sam -Sb > $bam";
    print LOG "# In mergeBam: $cmd\n" if ($verbose);

    system("$cmd");
    
    if (-e $sam){
	$cmd="rm $sam";
	print LOG "# In mergeBam: $cmd\n" if ($verbose);
	system("$cmd");
    }
    return(0);
}

sub mergeSummaryFiles{
    my ($in,$out,$db) = @_;
    open(TMP,"<$in");
    open(TMPOUT,"> $out");
    while (defined($_=<TMP>)){
        chomp;
        printf TMPOUT "%s\t$_\n",$db;
    }
    close(TMPOUT);
    close(TMP);
    $cmd = "rm $in";
    print LOG "# Doing: $cmd\n" if ($verbose);
    system("$cmd");

    return(0);
}

sub mapToDatabases{
    my ($bamIn, $db) = @_;

    my $i=1;
    my $cmd;


    while (defined($idx{index}{$db}[$i])){
	my $map = "$workDir/bam/map.$db.$i.bam";

	#
	# map against databases
	#
	my $index = $idx{index}{$db}[$i];
	my $annFile = $idx{ann}{$db}[$i];

	if (! -e $map){
	    &mapToIndex($bamIn,$index,$map);

	    my $annotationFile = $index . '.ann';
	    open(ANN,"head -1 $annotationFile |");
	    $_=<ANN>;
	    chomp;
	    my @w=split(/\s+/);
	    close(ANN);
	    $countEntries{$db}{sequences} += $w[1];
	    $countEntries{$db}{nuc} += $w[0];	    
	}
	$i++;
    }
}


sub mapToIndex{
    my ($inFile, $index, $outFile) = @_;

    #
    # Only accept properly pair
    #
    my $cmd;
    
    $cmd = "$prog_bedtools bamtofastq -i $inFile -fq /dev/stdout | $prog_bwa mem -t $cores -M $index";
    $cmd .= " - | $prog_samtools view -F 260 -Sb - > $outFile";
    print LOG "# Doing: $cmd\n" if ($verbose);

    system("$cmd");
    print $fh_runtime "# $cmd\n";
    $datestring = localtime();
    print $fh_runtime "$datestring\n";
    return(0);
}


sub post_process_single_database_cmd {
    my ($db) = @_;
    
    my $i=1;
    my $map;
    my $mapFilter;
    my $mapTo;
    my @cmdList=();

    while (-e "$workDir/bam/map.$db.$i.bam"){
	$map = "$workDir/bam/map.$db.$i.bam";
	$mapFilter = "$workDir/bam/mapFilter.$db.$i.bam";
	$mapTo = "$workDir/bam/mapTo.$db.$i.bam";
	if ((! -e $mapFilter) && (! -e $mapTo)){
            my $cmd ="$prog_filterBam -i $map -o $mapFilter -S";
            if (defined($Getopt::Std::opt_E) && defined($Getopt::Std::opt_T) && defined($Getopt::Std::opt_e)){
                $cmd .= " -c $absoluteMatchCount";
            }
            if (! defined($Getopt::Std::opt_T)){
                $cmd .= " -f $matchRatio";
            }
	    $cmd .= " -a $alignment_score_min";
	    if (defined($Getopt::Std::opt_N)){
		$cmd .= " -N $max_edit_distance";
	    }
	    push(@cmdList,$cmd);
	}
	$i++;
	push(@rmFiles, $map);
    }
    return(@cmdList);
}

sub best_hit{
    my $cmd;

    my $mapFilter;
    my $mapTo;
    my @bamFiles=();
    my @bamFilesRm=();
    my $db;

    foreach $db (@databases_fullmode){ 
	my $i=1;
	@bamFiles=();
	while (-e "$workDir/bam/mapFilter.$db.$i.bam"){
	    $mapFilter = "$workDir/bam/mapFilter.$db.$i.bam";
	    if (-e $mapFilter){
		push(@bamFiles, $mapFilter);
	    }
	    $i++;
	}
	my $fileCount= $#bamFiles +1;
	print LOG "Bamfiles in Fullmode, db= $db: $fileCount\n";
	if ($fileCount == 1){
	    $cmd="cp $workDir/bam/mapFilter.$db.1.bam $workDir/bam/mapTo.$db.1.bam";
	    print LOG "# Doing $cmd\n";
	    system("$cmd");
	}
	else{
	    #
	    # Only keep the best alignment score for duplicate hits
	    #
	    my %rec=();
	    my $i=0;
	    my $j=0;
	    foreach my $file (@bamFiles){
		$i++;
		open(FILE,"$prog_samtools view $file |");
		while (defined($_=<FILE>)){
		    my @w=split(/\s+/);
		    my $read = $w[0];
		    my $as = substr($w[13],5);
		    
		    if (! exists ($rec{$read}{sumAS})){
			$rec{$read}{sumAS}=$as;
			$rec{$read}{partition}=$i;
		    }
		    else{
			if ($as > $rec{$read}{sumAS}){
#			    if ($j <=10){
#				print LOG "# as sum before:$rec{$read}{partition}\t$read\t$rec{$read}{sumAS}\n";
#				print LOG "# as sum better: $i\t$read\t$as\n\n";
#				$j++;
#			    }
			    $rec{$read}{sumAS}=$as;
			    $rec{$read}{partition}=$i;
			}
		    }
		}
		close(FILE);
	    }
	    $i=0;
	    foreach my $file (@bamFiles){
		my $outfile = $file;
		$outfile =~ s/mapFilter/mapTo/g;

		$i++;
		open(FILE,"$prog_samtools view -h $file |");
		open(BAM,"| $prog_samtools view -h -Sb - > $outfile");
		my $headerLine=1;
		while (defined($_=<FILE>)){
		    if (((m/^\@SQ/) || (m/^\@PG/)) && ($headerLine)){
			print BAM "$_";
			next
		    }
		    else{
			$headerLine=0;
		    }
		    my @w=split(/\s+/);
		    my $read = $w[0];
		    
		    if ($rec{$read}{partition}==$i){
			print BAM "$_";
		    }
		}
		close(FILE);
		close(BAM);
		my $count1 = &bamCountReads($file);
		my $count2 = &bamCountReads($outfile);
		print LOG "read count $file: $count1\n";
		print LOG "read count $outfile: $count2\n";
	    }
	    %rec=();
	}
	push(@bamFilesRm, @bamFiles);
    }

    @bamFiles=();
    foreach $db (@databases_chainmode){
	my $i=1;
	while (-e "$workDir/bam/mapFilter.$db.$i.bam"){
	    $mapFilter = "$workDir/bam/mapFilter.$db.$i.bam";
	    push(@bamFiles,$mapFilter);
	    $i++;
	}
    }
    my $fileCount= $#bamFiles +1;
    print LOG "Bamfiles in bestmode: $fileCount\n";
    print LOG "@bamFiles\n";
    if ($fileCount == 1){
	foreach $db (@databases_chainmode){
	    $cmd="cp $workDir/bam/mapFilter.$db.1.bam $workDir/bam/mapTo.$db.1.bam";
	    print LOG "# Doing $cmd\n";
	    system("$cmd");
	}
    }
    else{
	#
	# Only keep the best alignment score for duplicate hits
	#
	my %rec=();
	my $i=0;
	my $j=0;
	foreach my $file (@bamFiles){
	    $i++;
	    print LOG "$i\t$file\n";

	    open(FILE,"$prog_samtools view $file |");
	    while (defined($_=<FILE>)){
		my @w=split(/\s+/);
		my $read = $w[0];
		my $as = substr($w[13],5);
		
		if (! exists ($rec{$read}{sumAS})){
		    $rec{$read}{sumAS}=$as;
		    $rec{$read}{partition}=$i;
		}
		else{
		    if ($as > $rec{$read}{sumAS}){
			if ($j%100000 == 0){
			    print LOG "# $j times a better alignment score is found\n";
			    print LOG "# as sum before:$rec{$read}{partition}\t$read\t$rec{$read}{sumAS}\n";
			    print LOG "# as sum better: $i\t$read\t$as\n\n";
			}
			$j++;
			$rec{$read}{sumAS}=$as;
			$rec{$read}{partition}=$i;
		    }
		}
	    }
	    close(FILE);
	}

	$i=0;
	foreach my $file (@bamFiles){
	    my $outfile = $file;
	    $outfile =~ s/mapFilter/mapTo/g;
	    $i++;
	    open(FILE,"$prog_samtools view -h $file |");
	    open(BAM,"| $prog_samtools view -h -Sb - > $outfile");
	    my $headerLine=1;
	    while (defined($_=<FILE>)){
		if (((m/^\@SQ/) || (m/^\@PG/)) && ($headerLine)){
		    print BAM "$_";
		    next
		}
		else{
		    $headerLine=0;
		}
		my @w=split(/\s+/);
		my $read = $w[0];
		
		if ($rec{$read}{partition} == $i){
		    print BAM "$_";
		}
	    }
	    close(FILE);
	    close(BAM);
	    my $count1 = &bamCountReads($file);
	    my $count2 = &bamCountReads($outfile);
	    print LOG "read count $file: $count1\n";
	    print LOG "read count $outfile: $count2\n";
	}
	push(@bamFilesRm, @bamFiles);
	%rec=();
    }
    foreach my $id (@bamFilesRm){
	if (-e $id){
	    my $cmd = "rm $id";
	    print LOG "# Doing: $cmd\n";
	    system("$cmd");
	}
    }
    return(0);
}


sub cleanData{
    my ($fastq, $outFile) = @_;

    my $dbName = $idx{index}{phiX174}[1];

    $cmd = "$prog_bwa mem -t $cores $dbName $fastq | $prog_samtools view -f4 -Sb - > $outFile";
    print LOG "# Doing: $cmd\n\n" if ($verbose);

    system("$cmd");
    print $fh_runtime "# $cmd\n";
    $datestring = localtime();
    print $fh_runtime "$datestring\n";

    $biological_relevant_reads = &bamCountReads($outFile);
    print LOG "Number of biological relevant reads: $biological_relevant_reads\n" if ($verbose);
    return(0);
}

sub cutadapt{
    my ($inFile, $outFile) = @_;

    my $log = "$outFile.log";
    my $cmd = "$prog_cutadapt -o $outFile $inFile > $log";
    
    return($cmd);
}


sub sortBam_cmd{
    my ($inFile, $outFile) = @_;
    my $cmd;

    if (! -e $outFile){
	$cmd = "$prog_samtools sort $inFile -o $outFile";
    }
    return($cmd);
}

sub samtoolsFixmate_cmd{
    my ($inFile, $outFile) = @_;

    if (! -e $outFile){
	$cmd = "$prog_samtools fixmate -m $inFile $outFile";
	return("$cmd");
    }
    return;
}

sub samtoolsSortByName_cmd{
    my ($inFile, $outFile) = @_;

    if (! -e $outFile){
	$cmd = "$prog_samtools sort -n $inFile -o $outFile";
	return("$cmd");
    }
    return;
}

sub samtoolsMarkdup_cmd{
    my ($inFile, $outFile) = @_;

    #
    # mark PCR duplicates
    #
    if (! -e $outFile){
	$cmd = "$prog_samtools markdup $inFile $outFile";
	return("$cmd");
    }
    return;
}

sub markdupRemove_cmd{
    my ($inFile, $outFile) = @_;

    #
    # remove PCR duplicates
    #
    if (! -e $outFile){
	$cmd = "$prog_samtools view -F 1024 $inFile -o $outFile -b";
	return("$cmd");
    }
    return;
}

sub bamToPileup_cmd{
    my ($bam, $i, $db, $outFile, $indexFile) = @_;
    my $index;
    if (-e $indexFile){
	$index=$indexFile;
    }
    else{
	$index = $idx{index}{$db}[$i];
    }
    my $cmd = '';

    # the -A option in samtools mpileup has the effect that also not properly paired reads are used
    # that can be usefull if bamtools convert actually make fastq files where reads are reversed if stated in bam file
    # and when its reversed it can no longer be properly paired as that only mean reads like this: ----->   <------
    if ((! -e $outFile) || (-z $outFile)){
	$cmd = "$prog_samtools mpileup -Q0 -f $index $bam > $outFile";
    }

    return($cmd);
}

sub pileupToDepth{
    my ($inFile, $i, $db, $summaryFile, $annFile) = @_;
    my $ann;
    if (-e $annFile){
	$ann=$annFile;
    }
    else{
	$ann = $idx{ann}{$db}[$i];
    }

    my $cmd = "$prog_mpileup2stat -i $inFile -a $ann > $summaryFile";
    
    if (! -e $summaryFile){
	print LOG "# Doing: $cmd\n" if ($verbose);
	system("$cmd");
        print $fh_runtime "# $cmd\n";
	$datestring = localtime();
	print $fh_runtime "$datestring\n";
    }
    return(0);
}

sub pileupToDepth_cmd{
    my ($inFile, $i, $db, $bam, $summaryFile, $annFile) = @_;
    my $ann;
    if (-e $annFile){
	$ann=$annFile;
    }
    else{
	$ann = $idx{ann}{$db}[$i];
    }

    my $cmd = "$prog_mpileup2stat -i $inFile -a $ann -b $bam -o  $summaryFile";
    
    return($cmd);
}


sub pileupToFasta_cmd{
    my ($inFile, $fasta) = @_;

    $cmd = "$prog_pileup2fasta -i $inFile -o $fasta -b $SNP_threshold";
    return($cmd);
}

sub readCountAnalysis{
    my ($inFile, $i, $db, $matrix) = @_;
    my $annFile = $idx{ann}{$db}[$i];
    my %rec=();


    open(ANN,"<$annFile");
    if (! eof ANN){
        # Read header line
	$_=<ANN>;
    }
    my $fastaName;
    my $desc;
    my @w=();
    while (defined($_=<ANN>)){
	if (m/^(\d+) (\S+) (.*)/){
	    $fastaName=$2;
	    $desc=$3;
	}

        if (! exists($rec{$fastaName}{count})){
            $rec{$fastaName}{count}=0;
            $rec{$fastaName}{desc}=$desc;
        }

        # read nest line which include length of that fasta entry
	$_=<ANN>;
        @w=split(/\s+/);
        $rec{$fastaName}{len}=$w[1];
    }
    close(ANN);

    my $cmd="$prog_samtools view $inFile | cut -f3 | sort | uniq -c |";
    print LOG "# Doing: $cmd\n" if ($verbose);

    open(READCOUNT,"$cmd ");
    while (defined($_=<READCOUNT>)){
        my @w=split(' ',$_);
        my $count=$w[0];
        my $fastaName=$w[1];
        if (! exists($rec{$fastaName}{count})){
            print LOG "Error In readCountAnalysis: $inFile,$i,$db\n";
            print LOG "A new fasta entry name is found '$fastaName' which was not in the file: $annFile\n";
            print LOG "$_";
        }
        else{
            $rec{$fastaName}{count}=$count;
        }
    }
    close(READCOUNT);
    if ($i == 1){
	print LOG "# writing to $matrix\n" if ($verbose);
        open(MATRIX,"> $matrix");
    }
    else{
	print LOG "# adding to $matrix\n" if ($verbose);
        open(MATRIX,">> $matrix");
    }

    foreach my $id (keys %rec){
        print MATRIX "$id\t$rec{$id}{count}\t$rec{$id}{len}\t$rec{$id}{desc}\n";
    }
    close(MATRIX);
    return(0);
}


sub bamCountReads{
    my ($bam) = @_;

    my $count;
    chomp($bam);
    if (exists($readCount{$bam})){
        $count=$readCount{$bam};
    }
    else{
	print LOG "# Doing: $prog_samtools flagstat $bam | head -1 | cut -f1 -d' '\n" if ($verbose);

        $count=`$prog_samtools flagstat $bam | head -1 | cut -f1 -d' '`;
        chomp($count);
        open(TMP,">>$readCountFile");
        print TMP "$bam\t$count\n";
        close(TMP);
	$readCount{$bam}=$count;
    }
    return($count);
}
