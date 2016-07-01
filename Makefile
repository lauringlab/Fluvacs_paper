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


data/raw/2014-5-30/%.fastq: data/raw/2014-5-30/SraRunInfo.csv
	./scripts/get_fastq.sh data/raw/2014-5-30/SraRunInfo.csv data/raw/2014-5-30/
	python scripts/change_names_sra.py -s data/raw/2014-5-30/ -f data/raw/2014-5-30/ -k data/raw/2014-5-30/SraRunInfo.csv -run
data/raw/2015-6-23/%.fastq: data/raw/2015-6-23/SraRunInfo.csv
	./scripts/get_fastq.sh data/raw/2015-6-23/SraRunInfo.csv data/raw/2015-6-23/
	python scripts/change_names_sra.py -s data/raw/2015-6-23/ -f data/raw/2015-6-23/ -k data/raw/2015-6-23/SraRunInfo.csv -run
data/raw/2015-11-14/%.fastq: data/raw/2015-11-14/SraRunInfo.csv
	./scripts/get_fastq.sh data/raw/2015-11-14/SraRunInfo.csv data/raw/2015-11-14/
	python scripts/change_names_sra.py -s data/raw/2015-11-14/ -f data/raw/2015-11-14/ -k data/raw/2015-11-14/SraRunInfo.csv -run


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

#####################################################################################
#
#
#
#	Part 3 : Secondary analysis make the figures!
#
#
#
#####################################################################################	



./results/figures.md: 
	cd ./results/ ; Rscript -e "knitr::knit('./figures.Rmd')"	
