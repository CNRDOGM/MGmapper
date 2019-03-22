#!/bin/bash

# mainDir is location where the following directories are located:
# mirror (location for mirror of organism)
# logs (logfiles for download process and parsing of files before making files located in db (see db below)
# db location of MGmapper formatted files for an organism
mainDir=$1
system=$2
DATE=`date +%Y%m%d`

# update is a speciel directory to see commands by running $MGmap_HOME/MGmapper_makedb.pl
updateDir=$mainDir/logs/update

mirrorPath=$mainDir/mirror
taxonomyDir=$mainDir/mirror/taxonomy

# MGmap_HOME is the environment variable for main dir for MGmapper scripts
perl_scripts=$MGmap_HOME/scripts
sh_scripts=$MGmap_HOME/MGmapper_makedb/$system

organism=vertebrate_mammals

metagenomicsDB=$mainDir/db/${organism}_${DATE}
metagenomicsLOG=$mainDir/db/${organism}_${DATE}/log
downloadLOG=$mainDir/logs/$organism/$DATE
#####################################
#
# check if perl scripts exists
#
#####################################
if [ ! -e ${perl_scripts}/parse_assembly_summary.pl ]; then
  echo "program not found: ${perl_scripts}/parse_assembly_summary.pl"
  exit
fi

if [ ! -e ${perl_scripts}/seqName_taxid.pl ]; then
  echo "program not found: ${perl_scripts}/seqName_taxid.pl"
  exit
fi

if [ ! -e ${perl_scripts}/tax.pl ]; then
  echo "program not found: ${perl_scripts}/tax.pl"
  exit
fi

if [ ! -e ${perl_scripts}/insertTaxonomy.pl ]; then
  echo "program not found: ${perl_scripts}/insertTaxonomy.pl"
  exit
fi

if [ ! -e ${perl_scripts}/joinSequences.pl ]; then
  echo "program not found: ${perl_scripts}/joinSequences.pl"
  exit
fi

if [ ! -e ${perl_scripts}/fastauniq.pl ]; then
  echo "program not found: ${perl_scripts}/fastauniq.pl"
  exit
fi

if [ ! -e ${perl_scripts}/fasta2split.pl ]; then
  echo "program not found: ${perl_scripts}/fasta2split.pl"
  exit
fi

#####################################################################
#   CREATE DIRS IF NECESSARY
#####################################################################

if [ ! -d "${metagenomicsLOG}" ]; then
    mkdir -p ${metagenomicsLOG}
fi

if [ ! -d "$mirrorPath/$organism/genbank" ]; then
  mkdir -p $mirrorPath/$organism/genbank
fi

#
# path for assembly summary files
#
if [ ! -d "$downloadLOG" ]; then
  mkdir -p $downloadLOG
fi
#
# path for logging the update commands
#
if [ ! -d "$updateDir" ]; then
  mkdir -p $updateDir
fi
updateLog="$updateDir/update.$DATE.log"
#####################################################################
#   MIRROR DATABASES
#####################################################################
#
# get assembly summary file
#
subset=${downloadLOG}/assembly_summary.filtered
assembly=${downloadLOG}/assembly_summary.txt

wget -nH --cut-dirs=3 -P $downloadLOG ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/assembly_summary.txt -o $downloadLOG/assembly_summary.log
cat  $assembly | perl -ane '{@w=split(/\t/); if ( (($w[10] eq "latest") && ($w[13] eq "Full")) && ( ($w[11] eq "Complete Genome") || ($w[4] eq "representative genome") || ($w[4] eq "reference genome")) ) {print "$_"}}' > $subset


echo "Date: $DATE" > ${metagenomicsDB}/${organism}.info
echo "Remark: genome_rep=Full; assembly_lvl=Complete Genome or refseq_cat=ref genome,rep genome" >> ${metagenomicsDB}/${organism}.info


${perl_scripts}/parse_assembly_summary.pl -i $subset -o $downloadLOG/genbank_wget.sh -w $downloadLOG/wget.log -m $mirrorPath/$organism/genbank -l ${downloadLOG}/parse_assembly_summary.log -v -t $downloadLOG/taxid_filename.txt

chmod a+rx $downloadLOG/genbank_wget.sh
/bin/tcsh $downloadLOG/genbank_wget.sh

${perl_scripts}/seqName_taxid.pl -i $downloadLOG/taxid_filename.txt -o $downloadLOG/seqName_taxid.txt -v -l $downloadLOG/seqName_taxid.log

taxFileNew=${metagenomicsDB}/${organism}.tax

${perl_scripts}/tax.pl -p $taxonomyDir -i $downloadLOG/seqName_taxid.txt -t taxid -o $taxFileNew -v -l ${metagenomicsLOG}/tax.log

${perl_scripts}/insertTaxonomy.pl -i $taxFileNew -k ${metagenomicsDB}/${organism}.kch -v -l ${metagenomicsLOG}/insertTaxonomy.log -s


cut -f 2 $downloadLOG/taxid_filename.txt | ${perl_scripts}/joinSequences.pl -v -o ${metagenomicsDB}/${organism}.fsa -l ${metagenomicsLOG}/joinSequences.log

#
# remove identical fasta entries
#
${perl_scripts}/fastauniq.pl -i ${metagenomicsDB}/$organism.fsa -o ${metagenomicsDB}/$organism -v -l ${metagenomicsLOG}/fastauniq.log
rm ${metagenomicsDB}/$organism.fsa
species=`grep '^>' ${metagenomicsDB}/${organism} | cut -f2,3 -d ' ' | sort -u | wc -l`
echo "Species: $species" >> ${metagenomicsDB}/${organism}.info

#
# split big fasta files into chunks containing approx 10 billion nucleotides ~10gb file size
#
fileSize=`du -b "${metagenomicsDB}/$organism" | cut -f1`
sizeLimit=10000000000
if [ "$fileSize" -gt "$sizeLimit" ]; then
    ${perl_scripts}/fasta2split.pl -i ${metagenomicsDB}/$organism -m $sizeLimit -v -l ${metagenomicsLOG}/fasta2split.log -o ${metagenomicsDB}/$organism
fi
fileSize=`du -h "${metagenomicsDB}/${organism}" | cut -f1`
sequences=`grep -c '^>' ${metagenomicsDB}/${organism} | cut -f1`
echo "Filesize: $fileSize" >> ${metagenomicsDB}/${organism}.info
echo "Sequences: $sequences" >> ${metagenomicsDB}/${organism}.info

echo -e "database\t$organism\tdate\t$DATE\tsequences\t$sequences" >> $updateLog
#
# make bwa and samtools index files via xqsub
#

${sh_scripts}/bwa.sh $organism ${metagenomicsDB} >& /dev/null
${sh_scripts}/samtools.sh $organism ${metagenomicsDB} >& /dev/null
${sh_scripts}/makeblastdb.sh ${metagenomicsDB}/$organism ${metagenomicsDB} $organism >& /dev/null

#
# change permissions 
#
chmod 755 ${metagenomicsDB}
chmod 644 ${metagenomicsDB}/*

chmod 755 ${metagenomicsLOG}
chmod 644 ${metagenomicsLOG}/*
