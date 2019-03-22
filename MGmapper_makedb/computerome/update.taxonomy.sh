#!/bin/bash

mainDir=$1
DATE=`date +%Y%m%d`

dbPath=$mainDir/mirror/taxonomy
logDir=$mainDir/logs/taxonomy/$DATE

# Download NCBI's taxonomic data and GenBank ID taxonomic
# assignments.

if [ ! -d "${dbPath}" ]; then
    mkdir -p ${dbPath}
fi

if [ ! -d "${logDir}" ]; then
    mkdir -p ${logDir}
fi


########################################################################
# DO DOWNLOAD
########################################################################

## Download taxdump
TAXDUMP="taxdump.tar.gz"
find $dbPath -name "*.dmp" -print0 | xargs -0 rm --force
#rm --force $dbPath/gc.prt $dbPath/readme.txt $dbPath/${TAXDUMP}
wget -nH --cut-dirs=2 -P $dbPath ftp://ftp.ncbi.nih.gov/pub/taxonomy/${TAXDUMP} -o $logDir/download.log
tar -zxvf $dbPath/${TAXDUMP} -C $dbPath
rm --force $dbPath/${TAXDUMP} $dbPath/{citations,division,gencode,merged,delnodes}.dmp $dbPath/gc.prt $dbPath/readme.txt >& /dev/null

## Download taxid
TAXID="gi_taxid_nucl.dmp.gz"
rm --force $dbPath/${TAXID/.gz/}* >& /dev/null
wget -nH --cut-dirs=2 -P $dbPath ftp://ftp.ncbi.nih.gov/pub/taxonomy/${TAXID} >& /dev/null
gunzip $dbPath/${TAXID}

touch $dbPath/done;
chmod 755 ${dbPath}
chmod 644 ${dbPath}/*
exit 0

trap '' SIGPIPE
