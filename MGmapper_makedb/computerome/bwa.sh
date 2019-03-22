#!/usr/bin/bash

db=$1
workDir=$2

# The environment variable MGmap_BWA is set by sourcing init file in main install dir for MGmapper: source MGmapper.init
bwa=$MGmap_BWA

file="${workDir}/${db}.1"
if [ -e "$file" ]; then
    for (( i=1;i<=100;i++ ))
    do
	fasta=${workDir}/${db}.$i
	re=${workDir}/log/${db}.${i}.bwa.e
	ro=${workDir}/log/${db}.${i}.bwa.o
	jobName=bwa.${db}.${i}
	
	if [ -e "$fasta" ]; then
	    xqsub -V -d $workDir -e $re -o $ro -N $jobName -l nodes=1:ppn=1,mem=40gb,walltime=48:00:00 -de $bwa index $fasta
	    sleep 2
	fi
    done
else
    fasta=${workDir}/${db}
    re=${workDir}/log/${db}.bwa.e
    ro=${workDir}/log/${db}.bwa.o
    jobName=bwa.${db}
    
    if [ -e "$fasta" ]; then
	xqsub -V -d $workDir -e $re -o $ro -N $jobName -l nodes=1:ppn=1,mem=40gb,walltime=48:00:00 -de $bwa index $fasta
    fi
fi
