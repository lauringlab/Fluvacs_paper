#all : ./results.pdf

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


./data/processed/Run_1293/Variants/all.sum.csv: data/raw/Run_1293/*.fastq
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1293/ -o ./data/processed/Run_1293/ -r ./data/reference/Brisbane_seq_untranslated -p bris -d two.sided -m fisher -a 0.9

./data/processed/Run_1304/Variants/all.sum.csv: data/raw/Run_1304/*.fastq
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1304/ -o ./data/processed/Run_1304/ -r ./data/reference/Brisbane_seq_untranslated -p Bris -d two.sided -m fisher -a 0.9

./data/processed/2007-2008/Variants/all.sum.csv: data/raw/2007-2008/*.fastq
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/2007-2008/ -o ./data/processed/2007-2008/ -r ./data/reference/Brisbane_seq_untranslated -p Brisbane -d two.sided -m fisher -a 0.9

./data/processed/Run_1412/Variants/all.sum.csv: data/raw/Run_1412/*.fastq
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1412/ -o ./data/processed/Run_1412/ -r ./data/reference/CalH3N2_untranslated -p CAcontrol -d two.sided -m fisher -a 0.9


./data/processed/2004-2005/Variants/all.sum.csv: data/raw/2004-2005/*.fastq
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/2004-2005/ -o ./data/processed/2004-2005/ -r ./data/reference/CalH3N2_untranslated -p Cal_H3N2 -d two.sided -m fisher -a 0.9


./data/processed/2005-2006/Variants/all.sum.csv: data/raw/2005-2006/*.fastq
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/2005-2006/ -o ./data/processed/2005-2006/ -r ./data/reference/CalH3N2_untranslated -p Cal_H3N2 -d two.sided -m fisher -a 0.9






#####################################################################################
#
#
#
#	Part 3 : Secondary analysis make the figures!
#
#
#
#####################################################################################




# get the start and stop of each segment in the concatenated genome.
data/concat_pos_bris.csv: data/processed/Run_1293/deepSNV/all.sum.csv
	Rscript --vanilla scripts/get_concat_pos.R data/processed/Run_1293/deepSNV/all.coverage.csv data/processed/Run_1304/deepSNV/all.coverage.csv ./data/concat_pos_bris.csv



# Parse the consensus files
#
#./scripts/parse_consensus.py data/processed/Run_1293/deepSNV/ ./data/concat_pos_bris.csv data/processed/Run_1293/parsed_fa
#
#mkdir ./data/processed/Run_1293/coding_fa
#for FILE in ./data/concat_pos_bris.csv data/processed/Run_1293/parsed_fa/*.fa; do
#	NAME=$(basename $FILE)
#	./scripts/trim_to_coding.py ~/muscle3.8.31/ $FILE ./data/reference/Brisbane_H3N2_plasmids.fa ./data/processed/Run_1293/coding_fa/$NAME
#done
#
#./scripts/concat_seg.py ./data/processed/Run_1293/coding_fa/ HA ./data/raw/2007_2008.meta.HAgm.csv ./data/processed/Run_1293/HA.fa
#./scripts/concat_seg.py ./data/processed/Run_1293/coding_fa/ NR ./data/raw/2007_2008.meta.HAgm.csv ./data/processed/Run_1293/NR.fa
#
#
#
#
#
#./scripts/parse_consensus.py data/processed/Run_1304/deepSNV/ ./data/concat_pos_bris.csv data/processed/Run_1304/parsed_fa
#
#mkdir ./data/processed/Run_1304/coding_fa
#for FILE in ./data/concat_pos_bris.csv data/processed/Run_1304/parsed_fa/*.fa; do
#	NAME=$(basename $FILE)
#	./scripts/trim_to_coding.py ~/muscle3.8.31/ $FILE ./data/reference/Brisbane_H3N2_plasmids.fa ./data/processed/Run_1304/coding_fa/$NAME
#done
#
#./scripts/concat_seg.py ./data/processed/Run_1304/coding_fa/ HA ./data/raw/2007_2008.meta.HAgm.csv ./data/processed/Run_1304/HA.fa
#./scripts/concat_seg.py ./data/processed/Run_1304/coding_fa/ NR ./data/raw/2007_2008.meta.HAgm.csv ./data/processed/Run_1304/NR.fa
#
#cat ./data/processed/Run_1293/HA.fa ./data/processed/Run_1304/HA.fa > ./data/processed/2007-2008.HA.fa
#cat ./data/processed/Run_1293/NR.fa ./data/processed/Run_1304/NR.fa > ./data/processed/2007-2008.NR.fa
#


./results/figures.md:
	cd ./results/ ; Rscript -e "knitr::knit('./figures.Rmd')"
