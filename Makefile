#all : ./results.pdf

# fastq files


finalResults= ./results/FluVacs_Figures.md ./results/2007-2008.putative.antigenic.csv ./results/2007-2008.HA.fa ./results/2007-2008.NR.fa ./results/2007-2008.HA.aa.csv

musclePath= /sw/lsa/centos7/muscle/3.8.31/bin/

#musclePath= ~/muscle3.8.31/
write.paper :	$(finalResults)
	echo 'Done!'




# Make figures and run R scripts
./results/FluVacs_Figures.md: ./results/2007-2008.wg.csv ./results/2004-2005.wg.csv ./results/2005-2006.wg.csv
	cd ./results/ ; Rscript -e "knitr::knit('./FluVacs_Figures.Rmd')"

pipelineAndCodingOutput= ./data/processed/bis_difference.csv ./data/processed/CalH3N2_difference.csv ./data/processed/Run_1293/Variants/all.sum.csv ./data/processed/Run_1304/Variants/all.sum.csv ./data/processed/2007-2008/Variants/all.sum.csv ./data/processed/2004-2005/Variants/all.sum.csv ./data/processed/2005-2006/Variants/all.sum.csv ./data/processed/Run_1412/Variants/all.sum.csv
./results/2007-2008.HA.csv ./results/2007-2008.wg.csv ./results/2004-2005.wg.csv ./results/2005-2006.wg.csv : $(pipelineAndCodingOutput)
	cd ./results/ ; Rscript -e "knitr::knit('./processing_var.Rmd')"



# Get the difference between the whole genome and coding regions

./data/processed/bis_difference.csv:./data/reference/Brisbane_seq_untranslated.fa ./data/reference/Brisbane_H3N2_plasmids.fa
	python ./scripts/trim_to_coding.py $(musclePath) ./data/reference/Brisbane_seq_untranslated.fa ./data/reference/Brisbane_H3N2_plasmids.fa -csv ./data/processed/bis_difference.csv

./data/processed/CalH3N2_difference.csv: ./data/reference/CalH3N2_untranslated.fa ./data/reference/CalH3N2SeqCont.fa
	python ./scripts/trim_to_coding.py $(musclePath) ./data/reference/CalH3N2_untranslated.fa ./data/reference/CalH3N2SeqCont.fa -csv ./data/processed/CalH3N2_difference.csv



####################################################################################
#
#
#	Part 1 : Downloading files. The following steps use a metadata file and 
#  		the sratoolkit. The files are downloaded automatically to 
#  		~/ncbi/public/sra/ which is the default location for the prefetch
#  		command. The script then converts this file to fastq files in the raw
#  		directory and removes the intermidiate .sra file from ~/ncbi/public/sra/.
#  		Finally we rename these fastq file to match the sampe names used in 
#  		downstream analysis. 
#
#
#
#
#####################################################################################

# Make reference file for each run

./data/raw/Run_1293/Run-1293.info.csv:./data/reference/SRR2fastq.csv
	head -n 1 data/reference/SRR2fastq.csv > ./data/raw/Run_1293/Run-1293.info.csv
	 more data/reference/SRR2fastq.csv | grep 'Run-1293'>>./data/raw/Run_1293/Run-1293.info.csv

./data/raw/Run_1304/Run-1304.info.csv:./data/reference/SRR2fastq.csv
	head -n 1 data/reference/SRR2fastq.csv > ./data/raw/Run_1304/Run-1304.info.csv
	more data/reference/SRR2fastq.csv | grep 'Run-1304'>>./data/raw/Run_1304/Run-1304.info.csv

./data/raw/Run_1412/Run-1412.info.csv:./data/reference/SRR2fastq.csv
	head -n 1 data/reference/SRR2fastq.csv > ./data/raw/Run_1412/Run-1412.info.csv
	more data/reference/SRR2fastq.csv | grep 'Run-1412'>>./data/raw/Run_1412/Run-1412.info.csv

./data/raw/2004-2005/2004-2005.info.csv:./data/reference/SRR2fastq.csv
	head -n 1 data/reference/SRR2fastq.csv > ./data/raw/2004-2005/2004-2005.info.csv
	more data/reference/SRR2fastq.csv | grep '2004-2005'>>./data/raw/2004-2005/2004-2005.info.csv


./data/raw/2005-2006/2005-2006.info.csv:./data/reference/SRR2fastq.csv
	head -n 1 data/reference/SRR2fastq.csv > ./data/raw/2005-2006/2005-2006.info.csv
	more data/reference/SRR2fastq.csv | grep '2005-2006'>>./data/raw/2005-2006/2005-2006.info.csv


./data/raw/2007-2008/2007-2008.info.csv:./data/reference/SRR2fastq.csv
	head -n 1 data/reference/SRR2fastq.csv > ./data/raw/2007-2008/2007-2008.info.csv
	more data/reference/SRR2fastq.csv | grep '2007-2008'>> ./data/raw/2007-2008/2007-2008.info.csv




# Download files and make SRA_files_are_Downloaded.txt 


./data/raw/Run_1293/SRA_files_are_Downloaded.txt: ./data/raw/Run_1293/Run-1293.info.csv
	./scripts/get_fastq.sh ./data/raw/Run_1293/Run-1293.info.csv ./data/raw/Run_1293/
	touch ./data/raw/Run_1293/SRA_files_are_Downloaded.txt
	
./data/raw/Run_1304/SRA_files_are_Downloaded.txt: ./data/raw/Run_1304/Run-1304.info.csv
	./scripts/get_fastq.sh ./data/raw/Run_1304/Run-1304.info.csv ./data/raw/Run_1304/
	touch ./data/raw/Run_1304/SRA_files_are_Downloaded.txt


./data/raw/Run_1412/SRA_files_are_Downloaded.txt: ./data/raw/Run_1412/Run-1412.info.csv
	./scripts/get_fastq.sh ./data/raw/Run_1412/Run-1412.info.csv ./data/raw/Run_1412/
	touch ./data/raw/Run_1412/SRA_files_are_Downloaded.txt


./data/raw/2004-2005/SRA_files_are_Downloaded.txt: ./data/raw/2004-2005/2004-2005.info.csv
	./scripts/get_fastq.sh ./data/raw/2004-2005/2004-2005.info.csv ./data/raw/2004-2005/
	touch ./data/raw/2004-2005/SRA_files_are_Downloaded.txt


./data/raw/2005-2006/SRA_files_are_Downloaded.txt: ./data/raw/2005-2006/2005-2006.info.csv
	./scripts/get_fastq.sh ./data/raw/2005-2006/2005-2006.info.csv ./data/raw/2005-2006/
	touch ./data/raw/2005-2006/SRA_files_are_Downloaded.txt


./data/raw/2007-2008/SRA_files_are_Downloaded.txt: ./data/raw/2007-2008/2007-2008.info.csv
	./scripts/get_fastq.sh ./data/raw/2007-2008/2007-2008.info.csv ./data/raw/2007-2008/
	touch ./data/raw/2007-2008/SRA_files_are_Downloaded.txt


# Rename files and move them to the proper directories
allFastq = $(fastq20042005) $(fastq20052006) $(fastq20072008) $(fastqRun1293) $(fastqRun1304) $(fastqRun1412)

./data/raw/Run_1293/renaming_log.txt : ./data/raw/Run_1293/SRA_files_are_Downloaded.txt
	python scripts/change_names_sra.py -s ./data/raw/Run_1293/ -f ./data/raw/Run_1293/ -k ./data/raw/Run_1293/Run-1293.info.csv  -run	


./data/raw/Run_1304/renaming_log.txt : ./data/raw/Run_1304/SRA_files_are_Downloaded.txt
	python scripts/change_names_sra.py -s ./data/raw/Run_1304/ -f ./data/raw/Run_1304/ -k ./data/raw/Run_1304/Run-1304.info.csv  -run	

./data/raw/Run_1412/renaming_log.txt : ./data/raw/Run_1412/SRA_files_are_Downloaded.txt
	python scripts/change_names_sra.py -s ./data/raw/Run_1412/ -f ./data/raw/Run_1412/ -k ./data/raw/Run_1412/Run-1412.info.csv  -run	


./data/raw/2004-2005/renaming_log.txt: ./data/raw/2004-2005/SRA_files_are_Downloaded.txt
	python scripts/change_names_sra.py -s data/raw/2004-2005/ -f data/raw/2004-2005/ -k data/raw/2004-2005/2004-2005.info.csv  -run	

./data/raw/2005-2006/renaming_log.txt: ./data/raw/2005-2006/SRA_files_are_Downloaded.txt
	python scripts/change_names_sra.py -s data/raw/2005-2006/ -f data/raw/2005-2006/ -k data/raw/2005-2006/2005-2006.info.csv  -run	

./data/raw/2007-2008/renaming_log.txt: ./data/raw/2007-2008/SRA_files_are_Downloaded.txt
	python scripts/change_names_sra.py -s data/raw/2007-2008/ -f data/raw/2007-2008/ -k data/raw/2007-2008/2007-2008.info.csv  -run	


#####################################################################################
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

./data/processed/Run_1293/Variants/all.sum.csv : ./data/raw/Run_1293/renaming_log.txt
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1293/ -o ./data/processed/Run_1293/ -r ./data/reference/Brisbane_seq_untranslated -p Brisbane -d two.sided -m fisher -a 0.9

./data/processed/Run_1304/Variants/all.sum.csv : ./data/raw/Run_1304/renaming_log.txt
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1304/ -o ./data/processed/Run_1304/ -r ./data/reference/Brisbane_seq_untranslated -p Brisbane -d two.sided -m fisher -a 0.9

./data/processed/2007-2008/Variants/all.sum.csv : ./data/raw/2007-2008/renaming_log.txt
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/2007-2008/ -o ./data/processed/2007-2008/ -r ./data/reference/Brisbane_seq_untranslated -p Brisbane -d two.sided -m fisher -a 0.9

./data/processed/Run_1412/Variants/all.sum.csv : ./data/raw/Run_1412/renaming_log.txt
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1412/ -o ./data/processed/Run_1412/ -r ./data/reference/CalH3N2_untranslated -p Cal-H3N2 -d two.sided -m fisher -a 0.9

./data/processed/2004-2005/Variants/all.sum.csv:  ./data/raw/2004-2005/renaming_log.txt
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/2004-2005/ -o ./data/processed/2004-2005/ -r ./data/reference/CalH3N2_untranslated -p Cal-H3N2 -d two.sided -m fisher -a 0.9

./data/processed/2005-2006/Variants/all.sum.csv : ./data/raw/2005-2006/renaming_log.txt
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/2005-2006/ -o ./data/processed/2005-2006/ -r ./data/reference/CalH3N2_untranslated -p Cal-H3N2 -d two.sided -m fisher -a 0.9


bowtie_alignments:
	bowtie2-build ./data/reference/Brisbane_seq_untranslated.fa ./data/reference/Brisbane_seq_untranslated
	bowtie2-build ./data/reference/CalH3N2_untranslated.fa ./data/reference/CalH3N2_untranslated




./results/2007-2008.putative.antigenic.csv ./results/2007-2008.HA.aa.csv : ./results/2007-2008.HA.fa ./results/2007-2008.HA.csv
		python ./scripts/antigenic_sites.py $(musclePath)

# Combine the concatenated coding regions
./results/2007-2008.HA.fa : ./data/processed/Run_1293/HA.fa ./data/processed/Run_1304/HA.fa
	cat ./data/processed/Run_1293/HA.fa ./data/processed/Run_1304/HA.fa > ./results/2007-2008.HA.fa
./results/2007-2008.NR.fa : ./data/processed/Run_1293/NR.fa ./data/processed/Run_1304/NR.fa
	cat ./data/processed/Run_1293/NR.fa ./data/processed/Run_1304/NR.fa > ./results/2007-2008.NR.fa
./results/2004-2006.HA.fa : ./data/processed/Run_1412/HA.fa
	cp ./data/processed/Run_1412/HA.fa results/2004-2006.HA.fa
./results/2204-2006.NR.fa : ./data/processed/Run_1412/NR.fa
	cp ./data/processed/Run_1412/NR.fa results/2004-2006.NR.fa

## Get consensus sequences
./data/processed/Run_1293/HA.fa ./data/processed/Run_1293/NR.fa: ./data/processed/Run_1293/Variants/all.sum.csv ./data/processed/concat_pos_bris.csv
	python ./scripts/consensus.pipe.py ./data/processed/Run_1293/deepSNV/ ./data/reference/Brisbane_H3N2_plasmids.fa ./data/processed/concat_pos_bris.csv $(musclePath)

./data/processed/Run_1304/HA.fa ./data/processed/Run_1304/NR.fa: ./data/processed/Run_1304/Variants/all.sum.csv ./data/processed/concat_pos_bris.csv
	python ./scripts/consensus.pipe.py ./data/processed/Run_1304/deepSNV/ ./data/reference/Brisbane_H3N2_plasmids.fa ./data/processed/concat_pos_bris.csv $(musclePath)

./data/processed/Run_1412/HA.fa ./data/processed/Run_1412/NR.fa : ./data/processed/Run_1412/Variants/all.sum.csv ./data/processed/concat_pos_CalH3N2.csv
	cp -v ./data/processed/2005-2006/deepSNV/536.removed.fa ./data/processed/Run_1412/deepSNV/536.removed.fa
	python ./scripts/consensus.pipe.py ./data/processed/Run_1412/deepSNV/ ./data/reference/CalH3N2SeqCont.fa ./data/processed/concat_pos_CalH3N2.csv $(musclePath)

./data/processed/2004-2005/HA.fa ./data/processed/2004-2005/NR.fa : ./data/processed/2004-2005/Variants/all.sum.csv ./data/processed/concat_pos_CalH3N2.csv
	python ./scripts/consensus.pipe.py data/processed/2004-2005/deepSNV/ ./data/reference/CalH3N2SeqCont.fa ./data/processed/concat_pos_CalH3N2.csv $(musclePath)

./data/processed/2005-2006/HA.fa ./data/processed/2005-2006/NR.fa : ./data/processed/2005-2006/Variants/all.sum.csv ./data/processed/concat_pos_CalH3N2.csv
	python ./scripts/consensus.pipe.py data/processed/2005-2006/deepSNV/ ./data/reference/CalH3N2SeqCont.fa ./data/processed/concat_pos_CalH3N2.csv $(musclePath)

# get the start and stop of each segment in the concatenated genome.
./data/processed/concat_pos_bris.csv: ./data/processed/Run_1293/Variants/all.sum.csv ./data/processed/Run_1304/Variants/all.sum.csv
	Rscript --vanilla scripts/get_concat_pos.R ./data/processed/Run_1293/deepSNV/all.coverage.csv ./data/processed/Run_1304/deepSNV/all.coverage.csv ./data/processed/concat_pos_bris.csv

./data/processed/concat_pos_CalH3N2.csv: ./data/processed/Run_1412/Variants/all.sum.csv data/processed/2004-2005/Variants/all.sum.csv  data/processed/2005-2006/Variants/all.sum.csv
	Rscript --vanilla scripts/get_concat_pos.R ./data/processed/Run_1412/deepSNV/all.coverage.csv data/processed/2004-2005/deepSNV/all.coverage.csv  data/processed/2005-2006/deepSNV/all.coverage.csv ./data/processed/concat_pos_CalH3N2.csv
