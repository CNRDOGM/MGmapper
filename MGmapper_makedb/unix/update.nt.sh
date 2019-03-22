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


# nt database
organism=nt

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

if [ ! -e $MGmap_BLASTDBCMD ]; then
  echo "Environment variable MGmap_BLASTDBCMD for blastdbcmd do not exists: $MGmap_BLASTDBCMD"
  exit
fi

#####################################################################
#   CREATE DIRS IF NECESSARY
#####################################################################
mirror=$mirrorPath/$organism
if [ ! -d "$mirror" ]; then
  mkdir -p $mirror
fi

# nt
if [ ! -d "${metagenomicsLOG}" ]; then
    mkdir -p ${metagenomicsLOG}
fi

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
#   DOWNLOAD DATABASE FILES
#####################################################################
#
# mirror nt tar files
#
wget -nH --cut-dirs=3 -m ftp://ftp.ncbi.nih.gov/blast/db/nt.'*'.tar.gz -P $mirror -o $downloadLOG/nt.log

echo "Date: $DATE" > ${metagenomicsDB}/${organism}.info

# unpack in dir $metagenomicsDB
for file in `ls -1 $mirror/nt.*.tar.gz`
do
gunzip -c $file | tar xvf - -C $metagenomicsDB
done


db=$metagenomicsDB/nt
fasta=$metagenomicsDB/nt
tax=$downloadLOG/seqName_taxid.txt

# get one big fasta file
if [ ! -f "$fasta" ]; then
    $MGmap_BLASTDBCMD -db $db -entry all > $fasta
fi

# get accession, taxid
if [ ! -f "$tax" ]; then
    $MGmap_BLASTDBCMD -db $db -entry all -outfmt "%a %T" | tr " " "\t" > $tax
fi

# Make a text file with taxid's at all levels up to superfamile
taxFileNew=${metagenomicsDB}/${organism}.tax
${perl_scripts}/tax.pl -p $taxonomyDir -i $downloadLOG/seqName_taxid.txt -t taxid -o $taxFileNew -v -l ${metagenomicsLOG}/tax.log

# Make a Kyoto Cabinet database file with taxonomies
${perl_scripts}/insertTaxonomy.pl -i $taxFileNew -k ${metagenomicsDB}/${organism}.kch -v -l ${metagenomicsLOG}/insertTaxonomy.log -s

#species=`grep '^>' ${metagenomicsDB}/${organism} | cut -f2,3 -d ' ' | sort -u| wc -l`
#echo "Species: $species" >> ${metagenomicsDB}/${organism}.info

# split big fasta files into chunks containing approx 10 billion nucleotides ~10gb file size
fileSize=`du -b "${metagenomicsDB}/${organism}" | cut -f1`
sizeLimit=10000000000
if [ "$fileSize" -gt "$sizeLimit" ]; then
    ${perl_scripts}/fasta2split.pl -i ${metagenomicsDB}/${organism} -m $sizeLimit -v -l ${metagenomicsLOG}/fasta2split.log -o ${metagenomicsDB}/$organism -b 2000000000
fi

fileSize=`du -h "${metagenomicsDB}/${organism}" | cut -f1`
sequences=`grep -c '^>' ${metagenomicsDB}/${organism} | cut -f1`
echo "Filesize: $fileSize" >> ${metagenomicsDB}/${organism}.info
echo "Sequences: $sequences" >> ${metagenomicsDB}/${organism}.info

echo -e "database\t$organism\tdate\t$DATE\tsequences\t$sequences" >> $updateLog

# make bwa and samtools index files via xqsub
${sh_scripts}/bwa.sh $organism ${metagenomicsDB} >& /dev/null
${sh_scripts}/samtools.sh $organism ${metagenomicsDB} >& /dev/null
#${sh_scripts}/makeblastdb.sh ${metagenomicsDB}/$organism ${metagenomicsDB} $organism>& /dev/null
