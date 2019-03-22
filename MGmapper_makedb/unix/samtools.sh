#!/usr/bin/bash

db=$1
workDir=$2

# The environment variable MGmap_SAMTOOLS is set by sourcing init file in main install dir for MGmapper: source MGmapper.init
samtools=$MGmap_SAMTOOLS


file="${workDir}/${db}.1"
if [ -e "$file" ]; then
    for (( i=1;i<=100;i++ ))
    do
	fasta=${workDir}/${db}.$i
	if [ -e "$fasta" ]; then
	    log="$workDir/log/samtools.$i.e"
	    $samtools faidx $fasta >& $log
	fi
    done
else
    fasta=${workDir}/${db}
    
    if [ -e "$fasta" ]; then
	log="$workDir/log/samtools.e"
	$samtools faidx $fasta >& $log
    fi
fi

