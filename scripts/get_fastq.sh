#!/bin/bash
SRA=($(cut -f1 -d '	' $1 |grep -v "sra_accession"))
for (( i=0; i<${#SRA[@]}; i++ ))
 do 
    echo downloading ${SRA[i]}
     wget -r -np -nd -k -P $2 ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${SRA[i]:0:6}/${SRA[i]}/${SRA[i]}.sra
    
    sra_file=$2${SRA[i]}.sra
    echo converting $sra_file to fastq
    fastq-dump  $sra_file --split-files -O $2
    echo removing $sra_file
    rm $sra_file
done    
