#!/bin/bash
fasta=$1
workDir=$2
name=$3
log="$workDir/log/blast.$name.log"

if [ ! -e "$fasta" ]; then
  echo "File not found: $fasta" >> $log
  exit
fi

re="$workDir/log/${name}.makeblastdb.e"
ro="$workDir/log/${name}.makeblastdb.o"
jobName="blast.$name"
xqsub -V -d $workDir -e $re -o $ro -N $jobName -l nodes=1:ppn=1,mem=50gb,walltime=1:00:00:00 -de $MGmap_MAKEBLASTDB -in $fasta -out $fasta -dbtype nucl -title $name -logfile $log -max_file_sz 10GB
