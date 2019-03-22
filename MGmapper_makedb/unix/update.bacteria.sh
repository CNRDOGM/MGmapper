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


# Bacteria
organism=bacteria

metagenomicsDB=$mainDir/db/${organism}_${DATE}
metagenomicsLOG=$mainDir/db/${organism}_${DATE}/log
downloadLOG=$mainDir/logs/$organism/$DATE

#
# plasmids
#
organismPlasmid=plasmid

metagenomicsDBPlasmid=$mainDir/db/${organismPlasmid}_${DATE}
metagenomicsLOGPlasmid=$mainDir/db/${organismPlasmid}_${DATE}/log
downloadLOGPlasmid=$mainDir/logs/$organismPlasmid/$DATE

# Bacteria_draft
organismDraft=bacteria_draft

metagenomicsDBDraft=$mainDir/db/${organismDraft}_${DATE}
metagenomicsLOGDraft=$mainDir/db/${organismDraft}_${DATE}/log
downloadLOGDraft=$mainDir/logs/$organismDraft/$DATE
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
if [ ! -e ${perl_scripts}/filter_assembly_summary.pl ]; then
  echo "program not found: ${perl_scripts}/fasta2split.pl"
  exit
fi

#####################################################################
#   CREATE DIRS IF NECESSARY
#####################################################################
if [ ! -d "$mirrorPath/$organism/genbank" ]; then
  mkdir -p $mirrorPath/$organism/genbank
fi

# bacteria
if [ ! -d "${metagenomicsLOG}" ]; then
    mkdir -p ${metagenomicsLOG}
fi


if [ ! -d "$downloadLOG" ]; then
  mkdir -p $downloadLOG
fi

# Plasmid
if [ ! -d "${metagenomicsLOGPlasmid}" ]; then
    mkdir -p ${metagenomicsLOGPlasmid}
fi


if [ ! -d "$downloadLOGPlasmid" ]; then
  mkdir -p $downloadLOGPlasmid
fi

# Bacteria_draft
if [ ! -d "${metagenomicsLOGDraft}" ]; then
    mkdir -p ${metagenomicsLOGDraft}
fi

if [ ! -d "$downloadLOGDraft" ]; then
  mkdir -p $downloadLOGDraft
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
# get assembly summary file
#
subset=${downloadLOG}/assembly_summary.filtered
subset2=${downloadLOGDraft}/assembly_summary.Bacteria.not_selected
draft=${downloadLOGDraft}/assembly_summary.draft
wget -nH --cut-dirs=3 -P $downloadLOG ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt -o $downloadLOG/assembly_summary.log
${perl_scripts}/filter_assembly_summary.pl -i $downloadLOG/assembly_summary.txt -o $subset -v -l $downloadLOG/filter_assembly_summary.log -n $subset2

cat  $subset2 | perl -ane '{@w=split(/\t/); if ( ($w[10] eq "latest") && ($w[13] eq "Full") && ($w[4] eq "representative genome") ) {print "$_"}}' > $draft

# Bacteria
echo "Date: $DATE" > ${metagenomicsDB}/${organism}.info
echo "Remark: genome_rep=Full; assembly_lvl=Complete Genome,Chr" >> ${metagenomicsDB}/${organism}.info

# Plasmid
echo "Date: $DATE" > ${metagenomicsDBPlasmid}/${organismPlasmid}.info
echo "Remark: genome_rep=Full; assembly_lvl=Complete Genome,Chr or refseq_cat=ref genome,rep genome" >> ${metagenomicsDBPlasmid}/${organismPlasmid}.info

# Bacteria_draft
echo "Date: $DATE" > ${metagenomicsDBDraft}/${organismDraft}.info
echo "Remark: genome_rep=Full; refseq_cat=rep genome" >> ${metagenomicsDBDraft}/${organismDraft}.info


# Bacteria
${perl_scripts}/parse_assembly_summary.pl -i $subset -o $downloadLOG/genbank_wget.sh -w $downloadLOG/wget.log -m $mirrorPath/$organism/genbank -l ${downloadLOG}/parse_assembly_summary.log -v -t $downloadLOG/taxid_filename.txt

# Bacteria_draft
${perl_scripts}/parse_assembly_summary.pl -i $draft -o $downloadLOGDraft/genbank_wget.sh -w $downloadLOGDraft/wget.log -m $mirrorPath/$organism/genbank -l ${downloadLOGDraft}/parse_assembly_summary.log -v -t $downloadLOGDraft/taxid_filename.txt

#
# Download new files or files where local and remote md5sum's differ
#

# Bacteria 
chmod a+rx $downloadLOG/genbank_wget.sh
/bin/tcsh $downloadLOG/genbank_wget.sh

# Baceria_draft
chmod a+rx $downloadLOGDraft/genbank_wget.sh
/bin/tcsh $downloadLOGDraft/genbank_wget.sh


# Bacteria
${perl_scripts}/seqName_taxid.pl -i $downloadLOG/taxid_filename.txt -o $downloadLOG/seqName_taxid.txt -v -l $downloadLOG/seqName_taxid.log -p $downloadLOGPlasmid/seqName_taxid.Bacteria.txt

# Bacteria_draft
${perl_scripts}/seqName_taxid.pl -i $downloadLOGDraft/taxid_filename.txt -o $downloadLOGDraft/seqName_taxid.txt -v -l $downloadLOGDraft/seqName_taxid.log -p $downloadLOGPlasmid/seqName_taxid.Bacteria_draft.txt

cat $downloadLOGPlasmid/seqName_taxid.Bacteria.txt $downloadLOGPlasmid/seqName_taxid.Bacteria_draft.txt > $downloadLOGPlasmid/seqName_taxid.txt

# Bacteria
taxFileNew=${metagenomicsDB}/${organism}.tax
${perl_scripts}/tax.pl -p $taxonomyDir -i $downloadLOG/seqName_taxid.txt -t taxid -o $taxFileNew -v -l ${metagenomicsLOG}/tax.log
${perl_scripts}/insertTaxonomy.pl -i $taxFileNew -k ${metagenomicsDB}/${organism}.kch -v -l ${metagenomicsLOG}/insertTaxonomy.log -s

# Bacteria_draft
taxFileNewDraft=${metagenomicsDBDraft}/${organismDraft}.tax
${perl_scripts}/tax.pl -p $taxonomyDir -i $downloadLOGDraft/seqName_taxid.txt -t taxid -o $taxFileNewDraft -v -l ${metagenomicsLOGDraft}/tax.log
${perl_scripts}/insertTaxonomy.pl -i $taxFileNewDraft -k ${metagenomicsDBDraft}/${organismDraft}.kch -v -l ${metagenomicsLOGDraft}/insertTaxonomy.log -s

# plasmid
taxFileNewPlasmid=${metagenomicsDBPlasmid}/${organismPlasmid}.tax
${perl_scripts}/tax.pl -p $taxonomyDir -i $downloadLOGPlasmid/seqName_taxid.txt -t taxid -o $taxFileNewPlasmid -v -l ${metagenomicsLOGPlasmid}/tax.log
${perl_scripts}/insertTaxonomy.pl -i $taxFileNewPlasmid -k ${metagenomicsDBPlasmid}/${organismPlasmid}.kch -v -l ${metagenomicsLOGPlasmid}/insertTaxonomy.log -s

# Bacteria
cut -f 2 $downloadLOG/taxid_filename.txt | ${perl_scripts}/joinSequences.pl -d plasmid -v -o ${metagenomicsDB}/${organism}.fsa -l ${metagenomicsLOG}/joinSequences.log -D ${metagenomicsDBPlasmid}/${organismPlasmid}.bacteria.fsa

# Bacteria_draft
cut -f 2 $downloadLOGDraft/taxid_filename.txt | ${perl_scripts}/joinSequences.pl -d plasmid -v -o ${metagenomicsDBDraft}/${organismDraft}.fsa -l ${metagenomicsLOGDraft}/joinSequences.log -D ${metagenomicsDBPlasmid}/${organismPlasmid}.bacteria_draft.fsa

# Plasmid
cat ${metagenomicsDBPlasmid}/${organismPlasmid}.bacteria.fsa ${metagenomicsDBPlasmid}/${organismPlasmid}.bacteria_draft.fsa > ${metagenomicsDBPlasmid}/${organismPlasmid}.fsa
rm ${metagenomicsDBPlasmid}/${organismPlasmid}.bacteria.fsa ${metagenomicsDBPlasmid}/${organismPlasmid}.bacteria_draft.fsa

#
# remove identical fasta entries
#

# Bacteria
${perl_scripts}/fastauniq.pl -i ${metagenomicsDB}/$organism.fsa -o ${metagenomicsDB}/$organism -v -l ${metagenomicsLOG}/fastauniq.log
rm ${metagenomicsDB}/${organism}.fsa

# Bacteria_draft
${perl_scripts}/fastauniq.pl -i ${metagenomicsDBDraft}/$organismDraft.fsa -o ${metagenomicsDBDraft}/$organismDraft -v -l ${metagenomicsLOGDraft}/fastauniq.log
rm ${metagenomicsDBDraft}/$organismDraft.fsa

# Plasmid
${perl_scripts}/fastauniq.pl -i ${metagenomicsDBPlasmid}/$organismPlasmid.fsa -o ${metagenomicsDBPlasmid}/$organismPlasmid -v -l ${metagenomicsLOGPlasmid}/fastauniq.log
rm ${metagenomicsDBPlasmid}/$organismPlasmid.fsa

# Bacteria
species=`grep '^>' ${metagenomicsDB}/${organism} | cut -f2,3 -d ' ' | sort -u| wc -l`
echo "Species: $species" >> ${metagenomicsDB}/${organism}.info

# Bacteria_draft
speciesDraft=`grep '^>' ${metagenomicsDBDraft}/${organismDraft} | cut -f2,3 -d ' ' | sort -u| wc -l`
echo "Species: $speciesDraft" >> ${metagenomicsDBDraft}/${organismDraft}.info

# plasmid
speciesPlasmid=`grep '^>' ${metagenomicsDBPlasmid}/${organismPlasmid} | cut -f2,3 -d ' ' | sort -u| wc -l`
echo "Species: $speciesPlasmid" >> ${metagenomicsDBPlasmid}/${organismPlasmid}.info

#
# split big fasta files into chunks containing approx 10 billion nucleotides ~10gb file size
#

# Bacteria
fileSize=`du -b "${metagenomicsDB}/${organism}" | cut -f1`
sizeLimit=10000000000
if [ "$fileSize" -gt "$sizeLimit" ]; then
    ${perl_scripts}/fasta2split.pl -i ${metagenomicsDB}/${organism} -m $sizeLimit -v -l ${metagenomicsLOG}/fasta2split.log -o ${metagenomicsDB}/$organism
fi

# Bacteria_draft
fileSize=`du -b "${metagenomicsDBDraft}/${organismDraft}" | cut -f1`
sizeLimit=10000000000
if [ "$fileSize" -gt "$sizeLimit" ]; then
    ${perl_scripts}/fasta2split.pl -i ${metagenomicsDBDraft}/${organismDraft} -m $sizeLimit -v -l ${metagenomicsLOGDraft}/fasta2split.log -o ${metagenomicsDBDraft}/$organismDraft
fi

# plasmid
fileSize=`du -b "${metagenomicsDBPlasmid}/${organismPlasmid}" | cut -f1`
sizeLimit=10000000000
if [ "$fileSize" -gt "$sizeLimit" ]; then
    ${perl_scripts}/fasta2split.pl -i ${metagenomicsDBPlasmid}/${organismPlasmid} -m $sizeLimit -v -l ${metagenomicsLOGPlasmid}/fasta2split.log -o ${metagenomicsDBPlasmid}/$organismPlasmid
fi


# Bacteria
fileSize=`du -h "${metagenomicsDB}/${organism}" | cut -f1`
sequences=`grep -c '^>' ${metagenomicsDB}/${organism} | cut -f1`
echo "Filesize: $fileSize" >> ${metagenomicsDB}/${organism}.info
echo "Sequences: $sequences" >> ${metagenomicsDB}/${organism}.info
echo -e "database\t$organism\tdate\t$DATE\tsequences\t$sequences" >> $updateLog

# Bacteria_draft
fileSize=`du -h "${metagenomicsDBDraft}/${organismDraft}" | cut -f1`
sequences=`grep -c '^>' ${metagenomicsDBDraft}/${organismDraft} | cut -f1`
echo "Filesize: $fileSize" >> ${metagenomicsDBDraft}/${organismDraft}.info
echo "Sequences: $sequences" >> ${metagenomicsDBDraft}/${organismDraft}.info
echo -e "database\t$organismDraft\tdate\t$DATE\tsequences\t$sequences" >> $updateLog

# plasmid
fileSizePlasmid=`du -h "${metagenomicsDBPlasmid}/${organismPlasmid}" | cut -f1`
sequencesPlasmid=`grep -c '^>' ${metagenomicsDBPlasmid}/${organismPlasmid} | cut -f1`
echo "Filesize: $fileSizePlasmid" >> ${metagenomicsDBPlasmid}/${organismPlasmid}.info
echo "Sequences: $sequencesPlasmid" >> ${metagenomicsDBPlasmid}/${organismPlasmid}.info
echo -e "database\t$organismPlasmid\tdate\t$DATE\tsequences\t$sequences" >> $updateLog

#
# make bwa and samtools index files
#

# Bacteria
sh ${sh_scripts}/bwa.sh $organism ${metagenomicsDB} $bwa >& /dev/null
sh ${sh_scripts}/samtools.sh $organism ${metagenomicsDB} $samtools >& /dev/null

# Bacteria_draft
sh ${sh_scripts}/bwa.sh $organismDraft ${metagenomicsDBDraft} $bwa >& /dev/null
sh ${sh_scripts}/samtools.sh $organismDraft ${metagenomicsDBDraft} $samtools >& /dev/null

# plasmid
sh ${sh_scripts}/bwa.sh $organismPlasmid ${metagenomicsDBPlasmid}  $bwa >& /dev/null
sh ${sh_scripts}/samtools.sh $organismPlasmid ${metagenomicsDBPlasmid}  $samtools >& /dev/null

#
# change permissions 
#

# Bacteria
chmod 755 ${metagenomicsDB}
chmod 644 ${metagenomicsDB}/*

chmod 755 ${metagenomicsLOG}
chmod 644 ${metagenomicsLOG}/*

# Bacteria_draft
chmod 755 ${metagenomicsDBDraft}
chmod 644 ${metagenomicsDBDraft}/*

chmod 755 ${metagenomicsLOGDraft}
chmod 644 ${metagenomicsLOGDraft}/*

# plasmid
chmod 755 ${metagenomicsDBPlasmid}
chmod 644 ${metagenomicsDBPlasmid}/*

chmod 755 ${metagenomicsLOGPlasmid}
chmod 644 ${metagenomicsLOGPlasmid}/*

