#all : ./results.pdf

# utility function to print various variables. For example, running the
# following at the command line:
#
#	make print-BAM
#
# will generate:
#	BAM=data/raw_june/V1V3_0001.bam data/raw_june/V1V3_0002.bam ...
print-%:
	@echo '$*=$($*)'


####################################################################################
#
#
#
#	Part 1  : Get the raw data
#
#
#
###################################################################################


#Download the fastq files for

#####################################################################################
#
#
#
#	Part 2  : Primary analysis - alignment, processing, and variant calling
#		This section relies heavily on a variant calling pipeline developed
#		in Bpipe. It requires Bowtie2 and DeepSNV. The pipeline can be
#		downloaded from https://github.com/lauringlab/variant_pipeline.git.
#		The following commonands expect the variant_pipeline directory to be
#		in your home directory.
#
#
#
#######################################################################################


./data/processed/Run_1293/Variants/all.sum.csv: #data/raw/Run_1293/*.fastq
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1293/ -o ./data/processed/Run_1293/ -r ./data/reference/Brisbane_seq_untranslated -p bris -d two.sided -m fisher -a 0.9

./data/processed/Run_1304/Variants/all.sum.csv: #data/raw/Run_1304/*.fastq
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1304/ -o ./data/processed/Run_1304/ -r ./data/reference/Brisbane_seq_untranslated -p Bris -d two.sided -m fisher -a 0.9

./data/processed/2007-2008/Variants/all.sum.csv: #data/raw/2007-2008/*.fastq
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/2007-2008/ -o ./data/processed/2007-2008/ -r ./data/reference/Brisbane_seq_untranslated -p Brisbane -d two.sided -m fisher -a 0.9

#####################################################################################
#
#
#
#	Part 3 : Secondary analysis make the figures!
#
#
#
#####################################################################################


# Get the difference between the whole genome and coding regions

./data/processed/bis_difference.csv:
	./scripts/trim_to_coding.py ~/muscle3.8.31/ ./data/reference/Brisbane_seq_untranslated.fa ./data/Brisbane_H3N2_plasmids.fa -csv ./data/processed/bis_difference.csv

# Make figures and run R scripts
./results/FluVacs_Figures.md ./data/processed/2007_2008.ha.csv: ./data/processed/bis_difference.csv ./data/processed/Run_1293/Variants/all.sum.csv ./data/processed/Run_1304/Variants/all.sum.csv ./data/processed/2007-2008/Variants/all.sum.csv
	cd ./results/ ; Rscript -e "knitr::knit('./FluVacs_Figures.Rmd')"


# get the start and stop of each segment in the concatenated genome.
./data/concat_pos_bris.csv: data/processed/Run_1293/deepSNV/all.sum.csv
	Rscript --vanilla scripts/get_concat_pos.R data/processed/Run_1293/deepSNV/all.coverage.csv data/processed/Run_1304/deepSNV/all.coverage.csv ./data/concat_pos_bris.csv

./data/processed/Run_1293/HA.fa ./data/processed/Run_1293/NR.fa: ./data/processed/Run_1293/Variants/all.sum.csv
	./scripts/consensus.pipe.py data/processed/Run_1293/deepSNV/ ./data/reference/Brisbane_H3N2_plasmids.fa ./data/concat_pos_bris.csv
./data/processed/Run_1304/HA.fa ./data/processed/Run_1304/NR.fa: ./data/processed/Run_1304/Variants/all.sum.csv
	./scripts/consensus.pipe.py data/processed/Run_1304/deepSNV/ ./data/reference/Brisbane_H3N2_plasmids.fa ./data/concat_pos_bris.csv

# Combine the concatenated coding regions
./data/processed/2007-2008.HA.fa : ./data/processed/Run_1293/HA.fa ./data/processed/Run_1304/HA.fa
	cat ./data/processed/Run_1293/HA.fa ./data/processed/Run_1304/HA.fa > ./data/processed/2007-2008.HA.fa
./data/processed/2007-2008.NR.fa : ./data/processed/Run_1293/NR.fa ./data/processed/Run_1304/NR.fa
	cat ./data/processed/Run_1293/NR.fa ./data/processed/Run_1304/NR.fa > ./data/processed/2007-2008.NR.fa

## Antigenic analysis

./data/processed/2007-2008.putative.antigenic.csv ./data/processed/2007-2008.HA.csv: ./data/processed/2007-2008.HA.fa ./data/processed/2007_2008.ha.csv
	./scripts/antigenic_sites.py 
