package MGmapper;
use strict;

sub find_program {
    my ($envvar, $prog) = @_;
    if (exists $ENV{$envvar}) {
        $prog = $ENV{$envvar};
        unless (-e $prog) {
            die "Program given by environment $envvar = $prog, does not exists\n";
        }
    } else {
        my $st = `which $prog 2> /dev/null`;
        if ($st) {
            chomp $st;
            $prog = $st;
        } else {
            die "Necessary program $prog is not in the path nor given in environment by $envvar\n";
        }
    }
    return $prog;
}


sub make_run_info_file{
    my ($workDir, $version, $command, $string) = @_;
    open(RUN_INFO,">$workDir/misc/run.info");
    print RUN_INFO "Version\t$version\n";
    print RUN_INFO "Command\t$command\n";
    print RUN_INFO "WorkDir\t$workDir\n";
    print RUN_INFO "Date\t$string\n";
    close(RUN_INFO);
    close(0);
    return(0);
}

sub init_MGmapper_summary{
    my ($fh, $file, $string, $startDir, $cmd) = @_;
    if (-e "$file"){
	open($fh,">>$file");	
    }
    else{
	open($fh,">$file");
    }
    print $fh "## Local date and time $string - start main program\n";
    print $fh "# command: $cmd\n";
    print $fh "# working dir: $startDir\n\n";
    return(0);
}

sub init_runtime_file{
    my ($fh, $file, $string) = @_;
    if (-e "$file"){
	open($fh,">>$file");
    }
    else{
	open($fh,">$file");
    }
    print $fh "# Start main program\n$string\n";
    return(0);
}

sub init_time_file{
    my ($fh, $file, $workDir, $string) = @_;
    open($fh,">$file");
    print $fh "#Job start\tJob end\n";
    print $fh "$string\t";
    return(0);
}

sub init_filestat{
    my ($fh, $file) = @_;
    open($fh,"> $file");
}

sub print_README{
    my ($file, $readCountMatrix) = @_;

    my $fh=\*README;
    open($fh,">$file");

    print $fh "Files in 'stat/*.txt are post-processing results made via the program MGmapper_classify.pl\n";
    print $fh "\n";
    print $fh "The files positive.* are those passing the MGmapper_classify.pl filtering criteria):\n";
    print $fh "MGmapper_classify.pl -i ../misc/stat.Bacteria.annot -c strain (input are paired-end reads)\n\n";
    print $fh "MGmapper_classify.pl -i ../misc/stat.Bacteria.annot -s -c strain (input are single-end reads)\n\n";
    print $fh "Default parameters for MGmapper_classify.pl are:\n\n";
    print $fh "-a 0.01 (abundance of at least 0.01%)\n";
    print $fh "-r 0.005 (readCountUniq/readCount)\n";
    print $fh "-m 0.01 (maximum mismatch ratio nucleotides/mismatched_nucleotides, where mismatch_nucleotides are number of changes need to make a read match the reference sequence, also known as the edit_distance)\n";
    print $fh "-n 10 (minimum number of mapped reads)#\n";
    print $fh "# At strain and species level -r 0.05 should be used\n";
    print $fh "# At all other levels i.e. genus phylum etc -r 0 should be used i.e. uniq read count ratio is not used\n";
    print $fh "#\n";
    print $fh "-g 0 (minimum coverage)\n";
    print $fh "-c strain (collapse at clade levels = strain\n\n";
    print $fh "-s (single-end reads, effecting the abundance calculation. if -s is not defined then abundance is divided by 2 due to pair-end reads)\n";
    print $fh "Possible clades are: superfamily, phylum, class, order, family, genus, species and strain\n";
    print $fh "Also, the data from different databases can be collapsed via option -d\n";
    print $fh "\n";
    print $fh "A special file is 'abundance.databases.txt' which shows read aboundance for all databases\n";
    print $fh "Percentages are calculated based on reads available after removing those that map to PhiX \(notPhiX\). notPhiX count is set\n";
    print $fh "Unmapped reads is calculated only on runs in Bestmode.\n";
    print $fh "\n";
    if (defined($readCountMatrix)){
        print $fh "A readcount matrix file: matrix/readCount.database.matrix\n";
        print $fh "A matrix showing ALL entries that are present in a database, also entries where not reads mapped to that entry\n";
        print $fh "Column 1: Name of reference sequence\n";
        print $fh "Column 2: Number of reads that mapped to that reference sequence\n";
        print $fh "Column 3: Length of reference sequence\n";
        print $fh "Column 4: Description of reference sequence\n";
        print $fh "\n";
    }
    print $fh "All files are combined into two excel workbook named 'Sample.xlsx' and 'Sample.negative.xlsx', where the latter includes \n";
    print $fh "hits at strain and species level that did NOT pass the filters when running MGmapper_classify.pl\n";
    close($fh);
    return(0);
}

sub readFastqList{
    my ($fileList) = @_;
    open(FileList,"<$fileList") || die "File not found: '$fileList'\n";
    my @LIST=();

    while (defined($_=<FileList>)){
	if (m/^\#/){
	    next;
	}
	chomp;
	my $len = length($_);
	if ($len < 1){next}
	my $file = $_;
	if (! -e $file){
	    print LOG "File not found; $file\n";
	    print STDERR "File not found; $file\nDone!\n";
	    die;
	}
	push(@LIST,$file);
    }
    close(FileList);
    return(@LIST);
}

sub make_www_tables{
    my ($workDir, $outFile, @databases) = @_;
    my $maxShowLines=5;

    open(TAB,">$outFile");

    my $delete_flag=1;

    my $statFile;
    my $cmd;
    foreach my $key (@databases){
        print LOG "Making snippet files for $key to output file= $outFile\n";
        $statFile = "$workDir/stat/positive.strain.$key.txt";
	if (-z $statFile){
	    $cmd="rm $statFile";
	    system("$cmd");
	    next;
	}
        if ((-e $statFile) && (! -z $statFile)){
            print TAB "#$key\n";
            print TAB "#Database\tRef Seq\tS_Abundance (%)\tR_Abundance (%)\tSize (bp)\tNucleotides\tCovered positions\tCoverage\tDepth\tReadCount\tReadCount uniq\tMismatches\tDescription\t";
	    print TAB "SuperFamily\tSuperFamily taxid\tPhylum\tPhylum taxid\tClass\tClass taxid\tOrder\tOrder taxid\tFamily\tFamily taxid\tGenus\tGenus taxid\tSpecies\tSpecies taxid\tStrain\tStrain taxid\n";
            open(A,"<$statFile");
            my $i=0;
            while (defined($_=<A>) && ($i<$maxShowLines)){
		if (/^#/){
		    next;
		}		
                print TAB "$_";
                $i++;
		$delete_flag=0;
            }
            close(A);
            print TAB "\n";
        }
    }
    close(TAB);
    if ($delete_flag){
	system("rm $outFile");
    }
    return(0);
}

sub countReads{
    my ($fastq, $fh, $verbose) = @_;

    my $cmd;
    chomp($fastq);
    if ($fastq =~ '.gz$'){
	$cmd = "gunzip -c $fastq | wc -l";
    }
    else{
	$cmd = "wc -l $fastq";
    }
    print $fh "# Doing: $cmd\n" if ($verbose);

    my $reads = `$cmd`;
    chomp($reads);

    $reads /= 4;

    return($reads);
}

1;
