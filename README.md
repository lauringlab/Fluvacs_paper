# Fluvacs Paper

This repository holds the analysis used in *citation,2016* it relies heavily on our other repository [variant_pipeline](https://github.com/lauringlab/variant_pipeline).

#Overview
--------

    project
    |- README          # the top level description of content
    |
    |- data            # raw and primary data, are not changed once created
    |  |- reference/  # reference fasta files to be used in in sequence alignment
    |  |- raw/         # raw data, will not be altered. In practice this is where the raw fastq files go. They will need to be downloaded from the SRA. 
    |  |- process/     # cleaned data, will not be altered once created. During anlysis this includes intermediate files used in variant calling as well as consensus sequence processing
    |
    |- scripts/           # any programmatic code
    |- results         # all output from workflows and analyses
    |  |- figures/     # graphs, likely designated for manuscript figures
    |  |- figures.Rmd  # exicutable R markdown to make the figures 
    +- Makefile        # executable Makefile for this study
    
  --------
# Dependencies    
The analysis expects the variant_pipeline repository to be your home directory. The published analysis was done using the variant_pipeline as it was on commit 3b7dd2b45e0ce392cacc60c94daf681709d7d83b.
```
    cd ~/
    git clone https://github.com/lauringlab/variant_pipeline.git
    git checkout 3b7dd2b45e0ce392cacc60c94daf681709d7d83b
```
This pipeline requires [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and the R package [DeepSNV](https://www.bioconductor.org/packages/release/bioc/html/deepSNV.html).

More information on the dependencies needed to run the varaint\_pipline can be found [here](https://github.com/lauringlab/variant_pipeline)

The python code used in this analysis has been tested in version 2.7.11 (although it will probably work in other versions) and relies on the following modules

```
argparse
Bio.Seq 
Bio 
Bio.Alphabet 
glob
os
os.path
pandas 
shutil
sys
subprocess
tempfile
```

The session information on the Rmarkdown file that reproduces the figures is below.
```
R version 3.2.3 (2015-12-10)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.11.6 (El Capitan)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggdendro_0.1-17 reshape2_1.4.1  ggplot2_2.1.0   plyr_1.8.3      knitr_1.13     

loaded via a namespace (and not attached):
 [1] MASS_7.3-45      colorspace_1.2-6 scales_0.4.0     magrittr_1.5     tools_3.2.3     
 [6] gtable_0.2.0     Rcpp_0.12.5      stringi_1.0-1    stringr_1.0.0    munsell_0.4.3 


```

The consensus sequence analysis and antigenic analysis relies on [muscle](http://www.drive5.com/muscle/downloads.htm). It expects the exicutable to be named "muscle". Please change the path to this exicutable in the Makefile.

# Reproducing the analysis

We can use the commands in the Makefile to do everyhing from downloading the fastq files from the SRA all the way through making the figures. There are 3 stages to the analysis. 1) Downloading the fastq files from the SRA, 2) Primary analsysis - Calling iSNV in each sample using the variant_pipeline repository referenced above, 3) Secondary analysis - Filtering the iSNV for quality, parsing the consensus files (deepSNV concatenates the 8 segments) and making the iSNV figures and tables. All the intermediate files needed to run the secondary analysis are included in this repository. 

Please note that in many places we refer to the genomic segment "NA" as "NR". This is to avoid complications in R as "NA" is a special term. 

## Downloading raw data

Becasue we use a plasmid control to estimate the lane specific error rates used to identify iSNV it is important that the fastq files from the SRA are split into separate directories for each Hiseq lane. This stage makes a file list for each run (if needed) and then downloads the .sra files using wget. The sra files are then converted to fastq files and are renamed to match the sample names used in the rest of the analysis. All of this is achieved my running the "download" phony target.

```
make download
```


## Processing the raw data

The following command launches the analyis pipeline from https://github.com/lauringlab/variant_pipeline. This pipeline is based in [bpipe](http://bpipe-test-documentation.readthedocs.io/en/latest/) and is smart enough to remember what commands have been run in the event of failure. It also logs the commands that were run("./data/processed/\*/commandlog.txt"). These steps are time and memory intensive. These commands output a number of large intermediate files. The important files that are used in down stream analysis are the concatenated variant calls "./data/processed/\*/Variants/all.sum.csv". The concatenated coverage files "./data/processed/\*/deepSNV/all.coverage.csv", and the deepSNV concatenated consensus sequences "./data/processed/\*/deepSNV/\*.fa". This stage of the analysis can be run using the "primary" phony target, and requires that the fastq files are downloaded and named appropriately.
```
make primary
```


## Secondary analysis

The secondary analysis is broken up into 3 separate stages. Each can be run separately or all three can be run at the same time. 

All stages are run by the phony target "secondary"

```
make secondary
```

### Consensus analysis

DeepSNV concatenates the consenus file of each sample. This stage uses a bpipe pipeline to deconcatenate the fasta files according the positions provided by deepSNV in the coverage file. These segments are then trimmed to match the coding regions provided by a separate fasta file. Finally we extract the HA and NA segements from each sample and output an HA fasta and an NA fasta file for the 2007-2008 season and the 2004-2005 & 2005-2006 seasons. In this work as with all the code here that involves alignments draw heavily (in many cases verbatim) from the HA_number script deleveped by Jesse Bloom at https://github.com/jbloomlab/HA_numbering/blob/master/HA_numbering.py. Many thanks to him for that great work.

```
make consensus
```

### iSNV analysis

The raw iSNV calls produced in the primary stage have not been filtered for quality. In this stage we filter the variants by p-value, frequency, mapping quality, phred score, and read position. At this point we narrow our focus to the coding regions of the genome. We also require that each variant be found in duplicate sequencing runs when the genome titer of the sample is between 10^3 and 10^5 genomes/ul. This is done by ./results/processing_var.Rmd. This file outputs a csv containing all iSNV for each season, a csv with iSNV found in HA alone (this is used to find putative antigenic variants below), and an updated meta data file that is helpful in naming the consensus files above. For these reasons the ./results/processing_var.Rmd script is referenced by all secondary targets. Once the iSNV are filtered the figures are made using ./results/Fluvacs_figures.Rmd. Both of these scripts are r markdown files and contain commentary on the analysis. The output from these analysis can be viewed here online in "results/Fluvacs_figures.md".

```
make isnv
```

### Antigenic analysis

This analysis combines output from the previous two targets. Each minor iSNV found in HA is placed in the contex of that sample's consensus seqeunce and each sequence is translated to determine where or not the mutation is synonymous or nonsynonymous. Nonsynonymous variants are then passed through an updated version of the Bloom labs [HA_number.py](https://github.com/jbloomlab/HA_numbering/blob/master/HA_numbering.py.) script. This aligns the varaints to common antigenic regions identified in H1N1 and H3N2 viruses.

```
make antigenic
```
