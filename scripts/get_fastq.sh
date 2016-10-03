#!/bin/bash
SRA=($(cut -f1 -d '	' $1 |grep -v "sra_accession"))
    for (( i=0; i<${#SRA[@]}; i++ )); do 
      echo downloading ${SRA[i]}
      prefetch ${SRA[i]}
      sra_file=$(ls ~/ncbi/public/sra/${SRA[i]}.sra) 
      echo converting $sra_file to fastq
      fastq-dump  $sra_file --split-files -O $2
      echo removing $sra_file
      rm $sra_file
      done    
