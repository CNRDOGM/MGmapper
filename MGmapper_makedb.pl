#!/usr/bin/env perl

# environment variables for bwa and samtools (MGmap_BWA, MGmap_SAMTOOLS)
# set path in MGmapper.init and do this: source MGmapper.init
# 

BEGIN{
    if (exists $ENV{MGmap_HOME}) {
        push(@INC, $ENV{MGmap_HOME} . '/modules');
    } else {
        die "Environment variable MGmap_HOME does not exists. Make it point to MGmapper's home.\n";
    }
    *LOG=*STDERR;
}

my $prog_bwa = &MGmapper::find_program('MGmap_BWA', "bwa");
my $prog_samtools = &MGmapper::find_program('MGmap_SAMTOOLS', "samtools");

####################
#
# Nothing needs to change below
#
###################

my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
  $command .= "$ARGV[$i] ";
  $i++;
}
use Getopt::Std;
use Cwd;
use strict;
use MGmapper::MGmapper;
# Default parameters
*LOG=*STDERR;
my $verbose=0;


#
# main dir where databases are installed
#
#my $mainDir = $ENV{MGmap_MAIN_DB};
my $main_db_dir = $ENV{MGmap_MAIN_DB};

my @databases=qw(taxonomy bacteria virus fungi human plant vertebrate_mammals vertebrate_other invertebrate nt archaea protozoa);

my %runningTime=();
my %runningMem=();
#
# running time on computerome
#
$runningTime{taxonomy}='1:00:00';
$runningMem{taxonomy}='1gb';

$runningTime{bacteria}='2:00:00:00';
$runningMem{bacteria}='30gb';

$runningTime{archaea}='2:00:00:00';
$runningMem{archaea}='30gb';

$runningTime{protozoa}='2:00:00:00';
$runningMem{protozoa}='30gb';

$runningTime{virus}='2:00:00:00';
$runningMem{virus}='30gb';

$runningTime{fungi}='2:00:00:00';
$runningMem{fungi}='30gb';

$runningTime{human}='2:00:00:00';
$runningMem{human}='30gb';

$runningTime{plant}='7:00:00:00';
$runningMem{plant}='50gb';

$runningTime{vertebrate_mammals}='7:00:00:00';
$runningMem{vertebrate_mammals}='50gb';

$runningTime{vertebrate_other}='7:00:00:00';
$runningMem{vertebrate_other}='50gb';

$runningTime{invertebrate}='7:00:00:00';
$runningMem{invertebrate}='50gb';

$runningTime{nt}='7:00:00:00';
$runningMem{nt}='50gb';

################################################
#
# Nothing needs to be changed below this line
#
################################################
my $DATE=`date +%Y%m%d`;
chomp($DATE);
#my $runningSystem='unix';
#
# HPC system at CBS 
#
my $runningSystem='computerome';
#
# Process command line
#
getopts('hvl:d:p:cr:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-d numbers] [-p path] [-c] [-l logfile] [-v]\n");
  print ("Description:\n");
  print ("$0 - download fasta files and update MGmapper databases\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -r  : run on system (computerome|unix) [$runningSystem]\n");
  print ("  -d  : comma separated list of database id\'s like '0,2,3'. See numbers below\n");
  print ("  -p  : main path were databases are installed [$main_db_dir]\n");
  print ("  -c  : check only - will find programs used and write commands to STDOUT but not execute them[off]\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -v  : Verbose [off]\n");
  print ("\n");
  print ("  Ex. to update taxonomy, bacteria and fungi datbases.\n");
  print ("  NB! when updating the bacteria db, 3 databases will be made: bacteria, bacteria_draft and plasmid\n");
  print ("  update.pl -v d 0,1,3 -c (only do a check (-c), shows the commands to be run without doing it)\n");
  print ("  update.pl -v d 0,1,3 (will update taxonomy, bacteria and fungi databases)\n");
  print ("\n");
  my $i=0;
  print ("== Database id\'s ==\n");
  while (defined($databases[$i])){
      print ("    $i:\t($databases[$i])\n");
      $i++;
  }
  
 exit;
} # Usage

if (defined($Getopt::Std::opt_l)){
    open(LOG,">$Getopt::Std::opt_l");
}
if (defined($Getopt::Std::opt_v)){
    $verbose=1;
}
if (defined($Getopt::Std::opt_p)){
    $main_db_dir=$Getopt::Std::opt_p;
}
my $mirrorDir = "$main_db_dir/mirror";
my $taxonomyDir = "$mirrorDir/taxonomy";

if (defined($Getopt::Std::opt_r)){
    $runningSystem=$Getopt::Std::opt_r;
    if (($runningSystem ne 'computerome') && ($runningSystem ne 'unix')){
	print "unknown system: '$runningSystem'\n";
	die;
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
my $cmd;
my @updateDatabases=();
my $runTaxonomy=0;

if (defined($Getopt::Std::opt_d)){
    my @w=split(/\,/,$Getopt::Std::opt_d);
    my $i=0;
    while (defined($w[$i])){
	my $db=$databases[$w[$i]];
	if ($db eq 'taxonomy'){
	    $runTaxonomy=1;
	    $i++;
	    next;
	}
	push(@updateDatabases, $db);
	$i++;
    }
}
else{
    print "No databases have been selected\n";
    die;
}

if (! -e "$taxonomyDir/names.dmp"){
    $runTaxonomy=1;
}


my $cmd;
my $command;
my $logDir;
my $updateLogDir;
my $updateLogFile;
if ($runTaxonomy){
    &update_taxonomy();
}

foreach my $type (@updateDatabases){
    $command='';

    if ($runningSystem eq 'computerome'){
	$logDir = "$main_db_dir/logs/${type}/${DATE}";

	if (! defined($Getopt::Std::opt_c)){
	    if (! -d $logDir){
		$cmd ="mkdir -p $logDir";
		print LOG "$cmd\n";
		system($cmd);
	    }
	}
	my $re = "${type}.e";
	my $ro = "${type}.o";
	my $jobName = "update." . "$type"; 
	$command = "xqsub -V -d $logDir -N $jobName -e $re -o $ro -l nodes=1:ppn=1,mem=$runningMem{$type},walltime=$runningTime{$type}  -de ";
	
    }
    $command .= "$ENV{MGmap_HOME}/MGmapper_makedb/$runningSystem/update.$type.sh $main_db_dir $runningSystem";

    if (! defined($Getopt::Std::opt_c)){
	$updateLogDir = "$main_db_dir/logs/update";
	if (! -d $updateLogDir){
	    $cmd = "mkdir -p $updateLogDir";
	    print LOG "$cmd\n" if ($verbose);
	    system("$cmd");
	}
        #
        # print update command and database info to file
        #
	$updateLogFile = "$updateLogDir/update.$DATE.log";
	if (-e $updateLogFile){
	    open(UPDATE_LOG,">>$updateLogFile");
	}
	else{
	    open(UPDATE_LOG,">$updateLogFile");
	}
	
	print UPDATE_LOG "command\t$command\tdate\t$DATE\n";
	close(UPDATE_LOG);
	print LOG "$command\n" if ($verbose);
	system("$command");
    }
    else{
	print "$command\n";
    }
}

if (! defined($Getopt::Std::opt_c)){
    if ($runningSystem eq 'computerome'){
	print LOG "# All started\n" if ($verbose);
    }
    else{
	print LOG "# Finished\n" if ($verbose);
    }
}

#############################################
#
# End main program
#
#############################################
sub update_taxonomy{
    $command='';
    my $type= 'taxonomy';
    if ($runningSystem eq 'computerome'){
	$logDir = "$main_db_dir/logs/${type}/${DATE}";

	if (! defined($Getopt::Std::opt_c)){
	    if (! -d $logDir){
		$cmd ="mkdir -p $logDir";
		print LOG "$cmd\n";
		system($cmd);
	    }
	}
	my $re = "$logDir/${type}.e";
	my $ro = "$logDir/${type}.o";
	my $jobName = "update." . "$type"; 
#	$command = "xqsub -V -d $mainDir/scripts -N $jobName -e $re -o $ro  -l nodes=1:ppn=1,mem=$runningMem{$type},walltime=$runningTime{$type}  -de ";
	$command = "xqsub -V -d $logDir -N $jobName -e $re -o $ro  -l nodes=1:ppn=1,mem=$runningMem{$type},walltime=$runningTime{$type}  -de ";
	
    }

    my $checkFile = "$taxonomyDir/done";
    if (! defined($Getopt::Std::opt_c)){
	if (-e "$checkFile"){
	    my $cmd = "rm $checkFile";
	    print LOG "$cmd\n" if ($verbose);
	    system("$cmd");
	}
    }

#    $command .="$ENV{MGmap_HOME}/scripts/$runningSystem/update.taxonomy.sh $mainDir";
    $command .="$ENV{MGmap_HOME}/MGmapper_makedb/$runningSystem/update.taxonomy.sh $main_db_dir";
    print LOG "$command\n" if ($verbose);
    if (defined($Getopt::Std::opt_c)){
	print "$command\n";
    }
    else{
	system($command);
    }

    if (! defined($Getopt::Std::opt_c)){
	while (! -e "$checkFile"){
	    my $sec=30;	
	    print LOG "# Sleeping $sec seconds until file $checkFile is created\n" if ($verbose);
	    sleep($sec);
	}
	print LOG "$checkFile created. Now remove it and continue\n" if ($verbose);
	my $cmd = "rm $checkFile";
	print LOG "$cmd\n" if ($verbose);
	system("$cmd");
    }
    return(0);
}

