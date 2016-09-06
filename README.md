# Fluvacs Paper

This repository holds the analysis used in *citation,2016* it relies heavily on our other repo variant_pipeline.

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
The analysis expects you to have the variant_pipeline repository cloned and saved in your home directory.
```
    cd ~/
    git clone https://github.com/lauringlab/variant_pipeline.git
```
This pipeline requires [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and the R package [DeepSNV](https://www.bioconductor.org/packages/release/bioc/html/deepSNV.html).

Information on the dependencies needed to run the varaint\_pipline can be found [here](https://github.com/lauringlab/variant_pipeline)


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

#Reproducing the analysis

We can use the commands in the Makefile to do everyhing from download the fastq files from the SRA all the way through making the figures. The Makefile contains comments that may be useful if you run into problems.

## Downloading raw data
The following commands download the necessary .sra files into the default directory ~/ncbi/public/sra/,  move and convert them to fastq files in the appropriate data/raw directory and remove the .sra file. These targets also give the fastq files meaningful names based on the meta data avaialable from the SRA.

```
   Under development 
```

## Processing the raw data

The following commands launch the analyis pipeline from https://github.com/lauringlab/variant_pipeline. This pipeline is based in [bpipe](http://bpipe-test-documentation.readthedocs.io/en/latest/) and is smart enough to remember what commands have been run in the event of failure. It also logs the commands that were run. These steps are time and memory intensive. These commands output a number of intermediate files as well as putative variants and consensus sequences for each sample. 

```
make ./data/processed/Run_1293/Variants/all.sum.csv
make ./data/processed/Run_1304/Variants/all.sum.csv
make ./data/processed/2007-2008/Variants/all.sum.csv
make ./data/processed/Run_1412/Variants/all.sum.csv
make ./data/processed/2004-2005/Variants/all.sum.csv
make ./data/processed/2005-2006/Variants/all.sum.csv
```


## Making the SNV figures

The following command calls the an intermediate vairant filtering stage and then "knits" an RMD file to create the figures based on the single nucleotide variant analysis.

I usually knit the RMD in R studio; however it should be possible using the following commands. Note : *I'm not sure how these commands work in the context of making a github_document output. You can change this line in figures.Rmd as needed.* 

```
make ./results/FluVacs_Figures.md
```

## Preliminary Consensus work and antigenic variants.

The following command deconcatinates the consensus sequence of each sample (the default DeepSNV settings produce a concatenated sequence), trims to just yield the coding sequences, identifies nonsynomomous variants and highlights putuatative antigenic variants.

```
make ./results/2007-2008.HA.aa.csv
make ./resutls/2007-2008.putative.antigenic.csv
```



## Running the analysis in it's entirety

The analysis can be run in it's entirety (from fastq to final figures) using the following command.

```
make write.paper
```

A pictoral summary of the analyis is provided below.

<img src="write.paper.png"/>
