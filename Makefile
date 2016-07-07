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





write.paper :	./results/FluVacs_Figures.md\
							./results/2007-2008.putative.antigenic.csv\
							./results/2007-2008.HA.fa\
							./results/2007-2008.NR.fa
	echo 'Done!'




# Make figures and run R scripts
./results/FluVacs_Figures.md ./results/2007-2008.HA.csv: ./data/processed/bis_difference.csv\
 																																./data/processed/CalH3N2_difference.csv\
																																./data/processed/Run_1293/Variants/all.sum.csv\
																																./data/processed/Run_1304/Variants/all.sum.csv\
																																./data/processed/2007-2008/Variants/all.sum.csv\
																																./data/processed/Run_1293/deepSNV/all.coverage.csv\
																																./data/processed/Run_1304/deepSNV/all.coverage.csv\
																																./data/processed/2007-2008/deepSNV/all.coverage.csv\
																																./data/processed/Run_1412/Variants/all.sum.csv\
																																./data/processed/Run_1412/deepSNV/all.coverage.csv\
																																./data/processed/2004-2005/Variants/all.sum.csv\
																																./data/processed/2004-2005/deepSNV/all.coverage.csv\
																																./data/processed/2005-2006/Variants/all.sum.csv\
																																./data/processed/2005-2006/deepSNV/all.coverage.csv

	cd ./results/ ; Rscript -e "knitr::knit('./FluVacs_Figures.Rmd')"

# Get the difference between the whole genome and coding regions

./data/processed/bis_difference.csv:
	python ./scripts/trim_to_coding.py ~/muscle3.8.31/ ./data/reference/Brisbane_seq_untranslated.fa ./data/reference/Brisbane_H3N2_plasmids.fa -csv ./data/processed/bis_difference.csv

./data/processed/CalH3N2_difference.csv:
	python ./scripts/trim_to_coding.py ~/muscle3.8.31/ ./data/reference/CalH3N2_untranslated.fa ./data/reference/CalH3N2SeqCont.fa -csv ./data/processed/CalH3N2_difference.csv


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
#The if statements don't run if the files already exist.
#This is a hack incase you want to rerun the secondary analysis only and
# already have the results of the primary analysis

./data/processed/Run_1293/Variants/all.sum.csv ./data/processed/Run_1293/deepSNV/all.coverage.csv: #data/raw/Run_1293/*.fastq
	if[ ! -f ./data/processed/Run_1293/Variants/all.sum.csv ];then python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1293/ -o ./data/processed/Run_1293/ -r ./data/reference/Brisbane_seq_untranslated -p bris -d two.sided -m fisher -a 0.9;fi

./data/processed/Run_1304/Variants/all.sum.csv ./data/processed/Run_1304/deepSNV/all.coverage.csv: #data/raw/Run_1304/*.fastq
	if[ ! -f ./data/processed/Run_1304/Variants/all.sum.csv ];then python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1304/ -o ./data/processed/Run_1304/ -r ./data/reference/Brisbane_seq_untranslated -p Bris -d two.sided -m fisher -a 0.9;fi

./data/processed/2007-2008/Variants/all.sum.csv ./data/processed/2007-2008/deepSNV/all.coverage.csv: #data/raw/2007-2008/*.fastq
	if[ ! -f ./data/processed/2007-2008/Variants/all.sum.csv ];then python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/2007-2008/ -o ./data/processed/2007-2008/ -r ./data/reference/Brisbane_seq_untranslated -p Brisbane -d two.sided -m fisher -a 0.9;fi

./data/processed/Run_1412/Variants/all.sum.csv ./data/processed/Run_1412/deepSNV/all.coverage.csv: #data/raw/Run_1412/*.fastq
	if[ ! -f ./data/processed/Run_1412/Variants/all.sum.csv ];then python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1412/ -o ./data/processed/Run_1412/ -r ./data/reference/CalH3N2_untranslated -p CAcontrol -d two.sided -m fisher -a 0.9;fi

./data/processed/2004-2005/Variants/all.sum.csv ./data/processed/2004-2005/deepSNV/all.coverage.csv: #data/raw/2004-2005/*.fastq
	if [ ! -f ./data/processed/2004-2005/Variants/all.sum.csv ];then python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/2004-2005/ -o ./data/processed/2004-2005/ -r ./data/reference/CalH3N2_untranslated -p Cal_H3N2 -d two.sided -m fisher -a 0.9;fi
./data/processed/2005-2006/Variants/all.sum.csv ./data/processed/2005-2006/deepSNV/all.coverage.csv: #data/raw/2005-2006/*.fastq
	if[ ! -f ./data/processed/2005-2006/Variants/all.sum.csv ];then python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/2005-2006/ -o ./data/processed/2005-2006/ -r ./data/reference/CalH3N2_untranslated -p Cal_H3N2 -d two.sided -m fisher -a 0.9;fi








./results/2007-2008.putative.antigenic.csv ./results/2007-2008.HA.aa.csv: ./results/2007-2008.HA.fa #| ./results/2007-2008.HA.csv
		python ./scripts/antigenic_sites.py

# Combine the concatenated coding regions
./results/2007-2008.HA.fa : ./data/processed/Run_1293/HA.fa ./data/processed/Run_1304/HA.fa
	cat ./data/processed/Run_1293/HA.fa ./data/processed/Run_1304/HA.fa > ./results/2007-2008.HA.fa
./results/2007-2008.NR.fa : ./data/processed/Run_1293/NR.fa ./data/processed/Run_1304/NR.fa
	cat ./data/processed/Run_1293/NR.fa ./data/processed/Run_1304/NR.fa > ./results/2007-2008.NR.fa

## Get consensus sequences
./data/processed/Run_1293/HA.fa ./data/processed/Run_1293/NR.fa: ./data/processed/Run_1293/Variants/all.sum.csv ./data/processed/concat_pos_bris.csv
	python ./scripts/consensus.pipe.py data/processed/Run_1293/deepSNV/ ./data/reference/Brisbane_H3N2_plasmids.fa ./data/processed/concat_pos_bris.csv
./data/processed/Run_1304/HA.fa ./data/processed/Run_1304/NR.fa: ./data/processed/Run_1304/Variants/all.sum.csv ./data/processed/concat_pos_bris.csv
	python ./scripts/consensus.pipe.py data/processed/Run_1304/deepSNV/ ./data/reference/Brisbane_H3N2_plasmids.fa ./data/processed/concat_pos_bris.csv

./data/processed/Run_1412/HA.fa ./data/processed/Run_1412/NR.fa : ./data/processed/Run_1412/deepSNV/all.sum.csv ./data/processed/concat_pos_CalH3N2.csv
	python ./scripts/consensus.pipe.py data/processed/Run_1412/deepSNV/ ./data/reference/CalH3N2SeqCont.fa ./data/processed/concat_pos_CalH3N2.csv

./data/processed/2004-2005/HA.fa ./data/processed/2004-2005/NR.fa : ./data/processed/2004-2005/deepSNV/all.sum.csv ./data/processed/concat_pos_CalH3N2.csv
	python ./scripts/consensus.pipe.py data/processed/2004-2005/deepSNV/ ./data/reference/CalH3N2SeqCont.fa ./data/processed/concat_pos_CalH3N2.csv

./data/processed/2005-2006/HA.fa ./data/processed/2005-2006/NR.fa : ./data/processed/2005-2006/deepSNV/all.sum.csv ./data/processed/concat_pos_CalH3N2.csv
	python ./scripts/consensus.pipe.py data/processed/2005-2006/deepSNV/ ./data/reference/CalH3N2SeqCont.fa ./data/processed/concat_pos_CalH3N2.csv

# get the start and stop of each segment in the concatenated genome.
./data/processed/concat_pos_bris.csv: data/processed/Run_1293/deepSNV/all.coverage.csv data/processed/Run_1304/deepSNV/all.coverage.csv
	Rscript --vanilla scripts/get_concat_pos.R data/processed/Run_1293/deepSNV/all.coverage.csv data/processed/Run_1304/deepSNV/all.coverage.csv ./data/processed/concat_pos_bris.csv

./data/processed/concat_pos_CalH3N2.csv: data/processed/Run_1412/deepSNV/all.coverage.csv data/processed/2004-2005/deepSNV/all.coverage.csv  data/processed/2005-2006/deepSNV/all.coverage.csv
	Rscript --vanilla scripts/get_concat_pos.R data/processed/Run_1412/deepSNV/all.coverage.csv data/processed/2004-2005/deepSNV/all.coverage.csv  data/processed/2005-2006/deepSNV/all.coverage.csv ./data/processed/concat_pos_CalH3N2.csv
