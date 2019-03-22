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
	re=${workDir}/log/${db}.${i}.samtools.e
	ro=${workDir}/log/${db}.${i}.samtools.o
	jobName=samtools.${db}.${i}
	
	if [ -e "$fasta" ]; then
	    xqsub -V -d $workDir -e $re -o $ro -N $jobName -l nodes=1:ppn=1,mem=40gb,walltime=24:00:00 -de $samtools faidx $fasta
	    sleep 2
	fi
    done
else
    fasta=${workDir}/${db}
    re=${workDir}/log/${db}.samtools.e
    ro=${workDir}/log/${db}.samtools.o
    jobName=samtools.${db}
    
    if [ -e "$fasta" ]; then
	xqsub -V -d $workDir -e $re -o $ro -N $jobName -l nodes=1:ppn=1,mem=40gb,walltime=24:00:00 -de $samtools faidx $fasta
    fi
fi

