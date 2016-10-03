#all : ./results.pdf

# fastq files

fastq20042005 =./data/raw/2004-2005/400.1.1.fastq ./data/raw/2004-2005/400.2.1.fastq ./data/raw/2004-2005/402.1.1.fastq ./data/raw/2004-2005/402.2.1.fastq ./data/raw/2004-2005/403.1.1.fastq ./data/raw/2004-2005/403.2.1.fastq ./data/raw/2004-2005/406.1.1.fastq ./data/raw/2004-2005/406.2.1.fastq ./data/raw/2004-2005/409.1.1.fastq ./data/raw/2004-2005/409.2.1.fastq ./data/raw/2004-2005/414.1.1.fastq ./data/raw/2004-2005/414.2.1.fastq ./data/raw/2004-2005/426.1.1.fastq ./data/raw/2004-2005/426.2.1.fastq ./data/raw/2004-2005/427.1.1.fastq ./data/raw/2004-2005/427.2.1.fastq ./data/raw/2004-2005/Cal-H3N2.1.1.fastq ./data/raw/2004-2005/Cal-H3N2.2.1.fastq

 fastq20052006 = ./data/raw/2005-2006/504.1.1.fastq  ./data/raw/2005-2006/504.2.1.fastq  ./data/raw/2005-2006/508.1.1.fastq  ./data/raw/2005-2006/508.2.1.fastq  ./data/raw/2005-2006/529.1.1.fastq  ./data/raw/2005-2006/529.2.1.fastq  ./data/raw/2005-2006/530.1.1.fastq  ./data/raw/2005-2006/530.2.1.fastq  ./data/raw/2005-2006/Cal-H3N2.1.1.fastq  ./data/raw/2005-2006/Cal-H3N2.2.1.fastq

fastq20072008 = ./data/raw/2007-2008/10.1.1.fastq ./data/raw/2007-2008/10.2.1.fastq ./data/raw/2007-2008/100.1.1.fastq ./data/raw/2007-2008/100.2.1.fastq ./data/raw/2007-2008/101.1.1.fastq ./data/raw/2007-2008/101.2.1.fastq ./data/raw/2007-2008/103.1.1.fastq ./data/raw/2007-2008/103.2.1.fastq ./data/raw/2007-2008/104.1.1.fastq ./data/raw/2007-2008/104.2.1.fastq ./data/raw/2007-2008/105.1.1.fastq ./data/raw/2007-2008/105.2.1.fastq ./data/raw/2007-2008/106.1.1.fastq ./data/raw/2007-2008/106.2.1.fastq ./data/raw/2007-2008/11.1.1.fastq ./data/raw/2007-2008/11.2.1.fastq ./data/raw/2007-2008/12.1.1.fastq ./data/raw/2007-2008/12.2.1.fastq ./data/raw/2007-2008/14.1.1.fastq ./data/raw/2007-2008/14.2.1.fastq ./data/raw/2007-2008/15.1.1.fastq ./data/raw/2007-2008/15.2.1.fastq ./data/raw/2007-2008/16.1.1.fastq ./data/raw/2007-2008/16.2.1.fastq ./data/raw/2007-2008/17.1.1.fastq ./data/raw/2007-2008/17.2.1.fastq ./data/raw/2007-2008/18.1.1.fastq ./data/raw/2007-2008/18.2.1.fastq ./data/raw/2007-2008/19.1.1.fastq ./data/raw/2007-2008/19.2.1.fastq ./data/raw/2007-2008/2.1.1.fastq ./data/raw/2007-2008/2.2.1.fastq ./data/raw/2007-2008/22.1.1.fastq ./data/raw/2007-2008/22.2.1.fastq ./data/raw/2007-2008/25.1.1.fastq ./data/raw/2007-2008/25.2.1.fastq ./data/raw/2007-2008/26.1.1.fastq ./data/raw/2007-2008/26.2.1.fastq ./data/raw/2007-2008/27.1.1.fastq ./data/raw/2007-2008/27.2.1.fastq ./data/raw/2007-2008/28.1.1.fastq ./data/raw/2007-2008/28.2.1.fastq ./data/raw/2007-2008/29.1.1.fastq ./data/raw/2007-2008/29.2.1.fastq ./data/raw/2007-2008/3.1.1.fastq ./data/raw/2007-2008/3.2.1.fastq ./data/raw/2007-2008/31.1.1.fastq ./data/raw/2007-2008/31.2.1.fastq ./data/raw/2007-2008/33.1.1.fastq ./data/raw/2007-2008/33.2.1.fastq ./data/raw/2007-2008/34.1.1.fastq ./data/raw/2007-2008/34.2.1.fastq ./data/raw/2007-2008/39.1.1.fastq ./data/raw/2007-2008/39.2.1.fastq ./data/raw/2007-2008/4.1.1.fastq ./data/raw/2007-2008/4.2.1.fastq ./data/raw/2007-2008/41.1.1.fastq ./data/raw/2007-2008/41.2.1.fastq ./data/raw/2007-2008/42.1.1.fastq ./data/raw/2007-2008/42.2.1.fastq ./data/raw/2007-2008/43.1.1.fastq ./data/raw/2007-2008/43.2.1.fastq ./data/raw/2007-2008/46.1.1.fastq ./data/raw/2007-2008/46.2.1.fastq ./data/raw/2007-2008/47.1.1.fastq ./data/raw/2007-2008/47.2.1.fastq ./data/raw/2007-2008/49.1.1.fastq ./data/raw/2007-2008/49.2.1.fastq ./data/raw/2007-2008/5.1.1.fastq ./data/raw/2007-2008/5.2.1.fastq ./data/raw/2007-2008/50.1.1.fastq ./data/raw/2007-2008/50.2.1.fastq ./data/raw/2007-2008/51.1.1.fastq ./data/raw/2007-2008/51.2.1.fastq ./data/raw/2007-2008/54.1.1.fastq ./data/raw/2007-2008/54.2.1.fastq ./data/raw/2007-2008/55.1.1.fastq ./data/raw/2007-2008/55.2.1.fastq ./data/raw/2007-2008/56.1.1.fastq ./data/raw/2007-2008/56.2.1.fastq ./data/raw/2007-2008/57.1.1.fastq ./data/raw/2007-2008/57.2.1.fastq ./data/raw/2007-2008/6.1.1.fastq ./data/raw/2007-2008/6.2.1.fastq ./data/raw/2007-2008/61.1.1.fastq ./data/raw/2007-2008/61.2.1.fastq ./data/raw/2007-2008/64.1.1.fastq ./data/raw/2007-2008/64.2.1.fastq ./data/raw/2007-2008/66.1.1.fastq ./data/raw/2007-2008/66.2.1.fastq ./data/raw/2007-2008/67.1.1.fastq ./data/raw/2007-2008/67.2.1.fastq ./data/raw/2007-2008/7.1.1.fastq ./data/raw/2007-2008/7.2.1.fastq ./data/raw/2007-2008/71.1.1.fastq ./data/raw/2007-2008/71.2.1.fastq ./data/raw/2007-2008/72.1.1.fastq ./data/raw/2007-2008/72.2.1.fastq ./data/raw/2007-2008/73.1.1.fastq ./data/raw/2007-2008/73.2.1.fastq ./data/raw/2007-2008/74.1.1.fastq ./data/raw/2007-2008/74.2.1.fastq ./data/raw/2007-2008/76.1.1.fastq ./data/raw/2007-2008/76.2.1.fastq ./data/raw/2007-2008/78.1.1.fastq ./data/raw/2007-2008/78.2.1.fastq ./data/raw/2007-2008/79.1.1.fastq ./data/raw/2007-2008/79.2.1.fastq ./data/raw/2007-2008/8.1.1.fastq ./data/raw/2007-2008/8.2.1.fastq ./data/raw/2007-2008/80.1.1.fastq ./data/raw/2007-2008/80.2.1.fastq ./data/raw/2007-2008/81.1.1.fastq ./data/raw/2007-2008/81.2.1.fastq ./data/raw/2007-2008/82.1.1.fastq ./data/raw/2007-2008/82.2.1.fastq ./data/raw/2007-2008/84.1.1.fastq ./data/raw/2007-2008/84.2.1.fastq ./data/raw/2007-2008/85.1.1.fastq ./data/raw/2007-2008/85.2.1.fastq ./data/raw/2007-2008/86.1.1.fastq ./data/raw/2007-2008/86.2.1.fastq ./data/raw/2007-2008/88.1.1.fastq ./data/raw/2007-2008/88.2.1.fastq ./data/raw/2007-2008/89.1.1.fastq ./data/raw/2007-2008/89.2.1.fastq ./data/raw/2007-2008/94.1.1.fastq ./data/raw/2007-2008/94.2.1.fastq ./data/raw/2007-2008/95.1.1.fastq ./data/raw/2007-2008/95.2.1.fastq ./data/raw/2007-2008/98.1.1.fastq ./data/raw/2007-2008/98.2.1.fastq ./data/raw/2007-2008/Brisbane.1.1.fastq ./data/raw/2007-2008/Brisbane.2.1.fastq

fastqRun1293 = ./data/raw/Run_1293/10.1.1.fastq ./data/raw/Run_1293/10.2.1.fastq ./data/raw/Run_1293/15.1.1.fastq ./data/raw/Run_1293/15.2.1.fastq ./data/raw/Run_1293/16.1.1.fastq ./data/raw/Run_1293/16.2.1.fastq ./data/raw/Run_1293/17.1.1.fastq ./data/raw/Run_1293/17.2.1.fastq ./data/raw/Run_1293/18.1.1.fastq ./data/raw/Run_1293/18.2.1.fastq ./data/raw/Run_1293/19.1.1.fastq ./data/raw/Run_1293/19.2.1.fastq ./data/raw/Run_1293/2.1.1.fastq ./data/raw/Run_1293/2.2.1.fastq ./data/raw/Run_1293/20.1.1.fastq ./data/raw/Run_1293/20.2.1.fastq ./data/raw/Run_1293/22.1.1.fastq ./data/raw/Run_1293/22.2.1.fastq ./data/raw/Run_1293/23.1.1.fastq ./data/raw/Run_1293/23.2.1.fastq ./data/raw/Run_1293/24.1.1.fastq ./data/raw/Run_1293/24.2.1.fastq ./data/raw/Run_1293/25.1.1.fastq ./data/raw/Run_1293/25.2.1.fastq ./data/raw/Run_1293/26.1.1.fastq ./data/raw/Run_1293/26.2.1.fastq ./data/raw/Run_1293/27.1.1.fastq ./data/raw/Run_1293/27.2.1.fastq ./data/raw/Run_1293/28.1.1.fastq ./data/raw/Run_1293/28.2.1.fastq ./data/raw/Run_1293/29.1.1.fastq ./data/raw/Run_1293/29.2.1.fastq ./data/raw/Run_1293/3.1.1.fastq ./data/raw/Run_1293/3.2.1.fastq ./data/raw/Run_1293/30.1.1.fastq ./data/raw/Run_1293/30.2.1.fastq ./data/raw/Run_1293/31.1.1.fastq ./data/raw/Run_1293/31.2.1.fastq ./data/raw/Run_1293/33.1.1.fastq ./data/raw/Run_1293/33.2.1.fastq ./data/raw/Run_1293/34.1.1.fastq ./data/raw/Run_1293/34.2.1.fastq ./data/raw/Run_1293/39.1.1.fastq ./data/raw/Run_1293/39.2.1.fastq ./data/raw/Run_1293/4.1.1.fastq ./data/raw/Run_1293/4.2.1.fastq ./data/raw/Run_1293/40.1.1.fastq ./data/raw/Run_1293/40.2.1.fastq ./data/raw/Run_1293/41.1.1.fastq ./data/raw/Run_1293/41.2.1.fastq ./data/raw/Run_1293/42.1.1.fastq ./data/raw/Run_1293/42.2.1.fastq ./data/raw/Run_1293/43.1.1.fastq ./data/raw/Run_1293/43.2.1.fastq ./data/raw/Run_1293/46.1.1.fastq ./data/raw/Run_1293/46.2.1.fastq ./data/raw/Run_1293/47.1.1.fastq ./data/raw/Run_1293/47.2.1.fastq ./data/raw/Run_1293/49.1.1.fastq ./data/raw/Run_1293/49.2.1.fastq ./data/raw/Run_1293/5.1.1.fastq ./data/raw/Run_1293/5.2.1.fastq ./data/raw/Run_1293/50.1.1.fastq ./data/raw/Run_1293/50.2.1.fastq ./data/raw/Run_1293/51.1.1.fastq ./data/raw/Run_1293/51.2.1.fastq ./data/raw/Run_1293/52.1.1.fastq ./data/raw/Run_1293/52.2.1.fastq ./data/raw/Run_1293/54.1.1.fastq ./data/raw/Run_1293/54.2.1.fastq ./data/raw/Run_1293/55.1.1.fastq ./data/raw/Run_1293/55.2.1.fastq ./data/raw/Run_1293/56.1.1.fastq ./data/raw/Run_1293/56.2.1.fastq ./data/raw/Run_1293/57.1.1.fastq ./data/raw/Run_1293/57.2.1.fastq ./data/raw/Run_1293/58.1.1.fastq ./data/raw/Run_1293/58.2.1.fastq ./data/raw/Run_1293/59.1.1.fastq ./data/raw/Run_1293/59.2.1.fastq ./data/raw/Run_1293/6.1.1.fastq ./data/raw/Run_1293/6.2.1.fastq ./data/raw/Run_1293/61.1.1.fastq ./data/raw/Run_1293/61.2.1.fastq ./data/raw/Run_1293/64.1.1.fastq ./data/raw/Run_1293/64.2.1.fastq ./data/raw/Run_1293/66.1.1.fastq ./data/raw/Run_1293/66.2.1.fastq ./data/raw/Run_1293/67.1.1.fastq ./data/raw/Run_1293/67.2.1.fastq ./data/raw/Run_1293/68.1.1.fastq ./data/raw/Run_1293/68.2.1.fastq ./data/raw/Run_1293/7.1.1.fastq ./data/raw/Run_1293/7.2.1.fastq ./data/raw/Run_1293/70.1.1.fastq ./data/raw/Run_1293/70.2.1.fastq ./data/raw/Run_1293/71.1.1.fastq ./data/raw/Run_1293/71.2.1.fastq ./data/raw/Run_1293/72.1.1.fastq ./data/raw/Run_1293/72.2.1.fastq ./data/raw/Run_1293/73.1.1.fastq ./data/raw/Run_1293/73.2.1.fastq ./data/raw/Run_1293/74.1.1.fastq ./data/raw/Run_1293/74.2.1.fastq ./data/raw/Run_1293/75.1.1.fastq ./data/raw/Run_1293/75.2.1.fastq ./data/raw/Run_1293/76.1.1.fastq ./data/raw/Run_1293/76.2.1.fastq ./data/raw/Run_1293/77.1.1.fastq ./data/raw/Run_1293/77.2.1.fastq ./data/raw/Run_1293/78.1.1.fastq ./data/raw/Run_1293/78.2.1.fastq ./data/raw/Run_1293/79.1.1.fastq ./data/raw/Run_1293/79.2.1.fastq ./data/raw/Run_1293/8.1.1.fastq ./data/raw/Run_1293/8.2.1.fastq ./data/raw/Run_1293/bris.1.1.fastq ./data/raw/Run_1293/bris.2.1.fastq

fastqRun1304 = ./data/raw/Run_1304/100.1.1.fastq ./data/raw/Run_1304/100.2.1.fastq ./data/raw/Run_1304/101.1.1.fastq ./data/raw/Run_1304/101.2.1.fastq ./data/raw/Run_1304/103.1.1.fastq ./data/raw/Run_1304/103.2.1.fastq ./data/raw/Run_1304/104.1.1.fastq ./data/raw/Run_1304/104.2.1.fastq ./data/raw/Run_1304/105.1.1.fastq ./data/raw/Run_1304/105.2.1.fastq ./data/raw/Run_1304/106.1.1.fastq ./data/raw/Run_1304/106.2.1.fastq ./data/raw/Run_1304/11.1.1.fastq ./data/raw/Run_1304/11.2.1.fastq ./data/raw/Run_1304/12.1.1.fastq ./data/raw/Run_1304/12.2.1.fastq ./data/raw/Run_1304/13.1.1.fastq ./data/raw/Run_1304/13.2.1.fastq ./data/raw/Run_1304/14.1.1.fastq ./data/raw/Run_1304/14.2.1.fastq ./data/raw/Run_1304/36.1.1.fastq ./data/raw/Run_1304/36.2.1.fastq ./data/raw/Run_1304/37.1.1.fastq ./data/raw/Run_1304/37.2.1.fastq ./data/raw/Run_1304/38.1.1.fastq ./data/raw/Run_1304/38.2.1.fastq ./data/raw/Run_1304/44.1.1.fastq ./data/raw/Run_1304/44.2.1.fastq ./data/raw/Run_1304/80.1.1.fastq ./data/raw/Run_1304/80.2.1.fastq ./data/raw/Run_1304/81.1.1.fastq ./data/raw/Run_1304/81.2.1.fastq ./data/raw/Run_1304/82.1.1.fastq ./data/raw/Run_1304/82.2.1.fastq ./data/raw/Run_1304/83.1.1.fastq ./data/raw/Run_1304/83.2.1.fastq ./data/raw/Run_1304/84.1.1.fastq ./data/raw/Run_1304/84.2.1.fastq ./data/raw/Run_1304/85.1.1.fastq ./data/raw/Run_1304/85.2.1.fastq ./data/raw/Run_1304/86.1.1.fastq ./data/raw/Run_1304/86.2.1.fastq ./data/raw/Run_1304/88.1.1.fastq ./data/raw/Run_1304/88.2.1.fastq ./data/raw/Run_1304/89.1.1.fastq ./data/raw/Run_1304/89.2.1.fastq ./data/raw/Run_1304/90.1.1.fastq ./data/raw/Run_1304/90.2.1.fastq ./data/raw/Run_1304/91.1.1.fastq ./data/raw/Run_1304/91.2.1.fastq ./data/raw/Run_1304/92.1.1.fastq ./data/raw/Run_1304/92.2.1.fastq ./data/raw/Run_1304/93.1.1.fastq ./data/raw/Run_1304/93.2.1.fastq ./data/raw/Run_1304/96.1.1.fastq ./data/raw/Run_1304/96.2.1.fastq ./data/raw/Run_1304/98.1.1.fastq ./data/raw/Run_1304/98.2.1.fastq ./data/raw/Run_1304/Bris.1.1.fastq ./data/raw/Run_1304/Bris.2.1.fastq

fastqRun1412 = ./data/raw/Run_1412/400.1.1.fastq ./data/raw/Run_1412/400.2.1.fastq ./data/raw/Run_1412/401.1.1.fastq ./data/raw/Run_1412/401.2.1.fastq ./data/raw/Run_1412/402.1.1.fastq ./data/raw/Run_1412/402.2.1.fastq ./data/raw/Run_1412/403.1.1.fastq ./data/raw/Run_1412/403.2.1.fastq ./data/raw/Run_1412/404.1.1.fastq ./data/raw/Run_1412/404.2.1.fastq ./data/raw/Run_1412/406.1.1.fastq ./data/raw/Run_1412/406.2.1.fastq ./data/raw/Run_1412/407.1.1.fastq ./data/raw/Run_1412/407.2.1.fastq ./data/raw/Run_1412/409.1.1.fastq ./data/raw/Run_1412/409.2.1.fastq ./data/raw/Run_1412/411.1.1.fastq ./data/raw/Run_1412/411.2.1.fastq ./data/raw/Run_1412/412.1.1.fastq ./data/raw/Run_1412/412.2.1.fastq ./data/raw/Run_1412/414.1.1.fastq ./data/raw/Run_1412/414.2.1.fastq ./data/raw/Run_1412/415.1.1.fastq ./data/raw/Run_1412/415.2.1.fastq ./data/raw/Run_1412/426.1.1.fastq ./data/raw/Run_1412/426.2.1.fastq ./data/raw/Run_1412/427.1.1.fastq ./data/raw/Run_1412/427.2.1.fastq ./data/raw/Run_1412/501.1.1.fastq ./data/raw/Run_1412/501.2.1.fastq ./data/raw/Run_1412/503.1.1.fastq ./data/raw/Run_1412/503.2.1.fastq ./data/raw/Run_1412/504.1.1.fastq ./data/raw/Run_1412/504.2.1.fastq ./data/raw/Run_1412/506.1.1.fastq ./data/raw/Run_1412/506.2.1.fastq ./data/raw/Run_1412/507.1.1.fastq ./data/raw/Run_1412/507.2.1.fastq ./data/raw/Run_1412/508.1.1.fastq ./data/raw/Run_1412/508.2.1.fastq ./data/raw/Run_1412/509.1.1.fastq ./data/raw/Run_1412/509.2.1.fastq ./data/raw/Run_1412/510.1.1.fastq ./data/raw/Run_1412/510.2.1.fastq ./data/raw/Run_1412/511.1.1.fastq ./data/raw/Run_1412/511.2.1.fastq ./data/raw/Run_1412/512.1.1.fastq ./data/raw/Run_1412/512.2.1.fastq ./data/raw/Run_1412/515.1.1.fastq ./data/raw/Run_1412/515.2.1.fastq ./data/raw/Run_1412/517.1.1.fastq ./data/raw/Run_1412/517.2.1.fastq ./data/raw/Run_1412/520.1.1.fastq ./data/raw/Run_1412/520.2.1.fastq ./data/raw/Run_1412/524.1.1.fastq ./data/raw/Run_1412/524.2.1.fastq ./data/raw/Run_1412/525.1.1.fastq ./data/raw/Run_1412/525.2.1.fastq ./data/raw/Run_1412/528.1.1.fastq ./data/raw/Run_1412/528.2.1.fastq ./data/raw/Run_1412/529.1.1.fastq ./data/raw/Run_1412/529.2.1.fastq ./data/raw/Run_1412/530.1.1.fastq ./data/raw/Run_1412/530.2.1.fastq ./data/raw/Run_1412/CAcontrol.1.1.fastq ./data/raw/Run_1412/CAcontrol.2.1.fastq


finalResults= ./results/FluVacs_Figures.md ./results/2007-2008.putative.antigenic.csv ./results/2007-2008.HA.fa ./results/2007-2008.NR.fa ./results/2007-2008.HA.aa.csv

musclePath= /sw/lsa/centos7/muscle/3.8.31/bin/

#musclePath= ~/muscle3.8.31/
write.paper :	$(finalResults)
	echo 'Done!'




# Make figures and run R scripts
./results/FluVacs_Figures.md: ./results/2007-2008.wg.csv
	cd ./results/ ; Rscript -e "knitr::knit('./FluVacs_Figures.Rmd')"

pipelineAndCodingOutput= ./data/processed/bis_difference.csv ./data/processed/CalH3N2_difference.csv ./data/processed/Run_1293/Variants/all.sum.csv ./data/processed/Run_1304/Variants/all.sum.csv ./data/processed/2007-2008/Variants/all.sum.csv ./data/processed/Run_1293/deepSNV/all.coverage.csv ./data/processed/Run_1304/deepSNV/all.coverage.csv ./data/processed/2007-2008/deepSNV/all.coverage.csv ./data/processed/Run_1412/Variants/all.sum.csv ./data/processed/Run_1412/deepSNV/all.coverage.csv ./data/processed/2004-2005/Variants/all.sum.csv ./data/processed/2004-2005/deepSNV/all.coverage.csv ./data/processed/2005-2006/Variants/all.sum.csv ./data/processed/2005-2006/deepSNV/all.coverage.csv
./results/2007-2008.HA.csv 	./results/2007-2008.wg.csv: $(pipelineAndCodingOutput)
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

./data/raw/Run_1293/Run-1293.info.tsv:./data/reference/SRA_metadata.txt.tsv
	head -n 1 data/reference/SRA_metadata.txt.tsv > ./data/raw/Run_1293/Run-1293.info.tsv
	 more data/reference/SRA_metadata.txt.tsv | grep 'Run-1293'>>./data/raw/Run_1293/Run-1293.info.tsv

./data/raw/Run_1304/Run-1304.info.tsv:./data/reference/SRA_metadata.txt.tsv
	head -n 1 data/reference/SRA_metadata.txt.tsv > ./data/raw/Run_1304/Run-1304.info.tsv
	more data/reference/SRA_metadata.txt.tsv | grep 'Run-1304'>>./data/raw/Run_1304/Run-1304.info.tsv

./data/raw/Run_1412/Run-1412.info.tsv:./data/reference/SRA_metadata.txt.tsv
	head -n 1 data/reference/SRA_metadata.txt.tsv > ./data/raw/Run_1412/Run-1412.info.tsv
	more data/reference/SRA_metadata.txt.tsv | grep 'Run-1412'>>./data/raw/Run_1412/Run-1412.info.tsv

./data/raw/2004-2005/2004-2005.info.tsv:./data/reference/SRA_metadata.txt.tsv
	head -n 1 data/reference/SRA_metadata.txt.tsv > ./data/raw/2004-2005/2004-2005.info.tsv
	more data/reference/SRA_metadata.txt.tsv | grep '2004-2005'>>./data/raw/2004-2005/2004-2005.info.tsv


./data/raw/2005-2006/2005-2006.info.tsv:./data/reference/SRA_metadata.txt.tsv
	head -n 1 data/reference/SRA_metadata.txt.tsv > ./data/raw/2005-2006/2005-2006.info.tsv
	more data/reference/SRA_metadata.txt.tsv | grep '2005-2006'>>./data/raw/2005-2006/2005-2006.info.tsv


./data/raw/2007-2008/2007-2008.info.tsv:./data/reference/SRA_metadata.txt.tsv
	head -n 1 data/reference/SRA_metadata.txt.tsv > ./data/raw/2007-2008/2007-2008.info.tsv
	more data/reference/SRA_metadata.txt.tsv | grep '2007-2008'>> ./data/raw/2007-2008/2007-2008.info.tsv




# Download files and make SRA_files_are_Downloaded.txt 


./data/raw/Run_1293/SRA_files_are_Downloaded.txt: ./data/raw/Run_1293/Run-1293.info.tsv
	./scripts/get_fastq.sh ./data/raw/Run_1293/Run-1293.info.tsv ./data/raw/Run_1293/
	touch ./data/raw/Run_1293/SRA_files_are_Downloaded.txt
	
./data/raw/Run_1304/SRA_files_are_Downloaded.txt: ./data/raw/Run_1304/Run-1304.info.tsv
	./scripts/get_fastq.sh ./data/raw/Run_1304/Run-1304.info.tsv ./data/raw/Run_1304/
	touch ./data/raw/Run_1304/SRA_files_are_Downloaded.txt


./data/raw/Run_1412/SRA_files_are_Downloaded.txt: ./data/raw/Run_1412/Run-1412.info.tsv
	./scripts/get_fastq.sh ./data/raw/Run_1412/Run-1412.info.tsv ./data/raw/Run_1412/
	touch ./data/raw/Run_1412/SRA_files_are_Downloaded.txt


./data/raw/2004-2005/SRA_files_are_Downloaded.txt: ./data/raw/2004-2005/2004-2005.info.tsv
	./scripts/get_fastq.sh ./data/raw/2004-2005/2004-2005.info.tsv ./data/raw/2004-2005/
	touch ./data/raw/2004-2005/SRA_files_are_Downloaded.txt


./data/raw/2005-2006/SRA_files_are_Downloaded.txt: ./data/raw/2005-2006/2005-2006.info.tsv
	./scripts/get_fastq.sh ./data/raw/2005-2006/2005-2006.info.tsv ./data/raw/2005-2006/
	touch ./data/raw/2005-2006/SRA_files_are_Downloaded.txt


./data/raw/2007-2008/SRA_files_are_Downloaded.txt: ./data/raw/2007-2008/2007-2008.info.tsv
	./scripts/get_fastq.sh ./data/raw/2007-2008/2007-2008.info.tsv ./data/raw/2007-2008/
	touch ./data/raw/2007-2008/SRA_files_are_Downloaded.txt


# Rename files and move them to the proper directories
allFastq = $(fastq20042005) $(fastq20052006) $(fastq20072008) $(fastqRun1293) $(fastqRun1304) $(fastqRun1412)

$(fastqRun1293) : ./data/raw/Run_1293/SRA_files_are_Downloaded.txt
	python scripts/change_names_sra.py -s ./data/raw/Run_1293/ -f ./data/raw/Run_1293/ -k ./data/raw/Run_1293/Run-1293.info.tsv  -run	


$(fastqRun1304) : ./data/raw/Run_1304/SRA_files_are_Downloaded.txt
	python scripts/change_names_sra.py -s ./data/raw/Run_1304/ -f ./data/raw/Run_1304/ -k ./data/raw/Run_1304/Run-1304.info.tsv  -run	

$(fastqRun1412) : ./data/raw/Run_1412/SRA_files_are_Downloaded.txt
	python scripts/change_names_sra.py -s ./data/raw/Run_1412/ -f ./data/raw/Run_1412/ -k ./data/raw/Run_1412/Run-1412.info.tsv  -run	


$(fastq20042005) : ./data/raw/2004-2005/SRA_files_are_Downloaded.txt
	python scripts/change_names_sra.py -s data/raw/2004-2005/ -f data/raw/2004-2005/ -k data/raw/2004-2005/2004-2005.info.tsv  -run	

$(fastq20052006) : ./data/raw/2005-2006/SRA_files_are_Downloaded.txt
	python scripts/change_names_sra.py -s data/raw/2005-2006/ -f data/raw/2005-2006/ -k data/raw/2005-2006/2005-2006.info.tsv  -run	

$(fastq20072008) : ./data/raw/2007-2008/SRA_files_are_Downloaded.txt
	python scripts/change_names_sra.py -s data/raw/2007-2008/ -f data/raw/2007-2008/ -k data/raw/2007-2008/2007-2008.info.tsv  -run	


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
#The if statements don't run if the files already exist.
#This is a hack incase you want to rerun the secondary analysis only and
# already have the results of the primary analysis

./data/processed/Run_1293/Variants/all.sum.csv ./data/processed/Run_1293/deepSNV/all.coverage.csv ./data/processed/Run_1293/deepSNV/*.fa: $(fastqRun1293)
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1293/ -o ./data/processed/Run_1293/ -r ./data/reference/Brisbane_seq_untranslated -p bris -d two.sided -m fisher -a 0.9

./data/processed/Run_1304/Variants/all.sum.csv ./data/processed/Run_1304/deepSNV/all.coverage.csv ./data/processed/Run_1304/deepSNV/*.fa: $(fastqRun1304)
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1304/ -o ./data/processed/Run_1304/ -r ./data/reference/Brisbane_seq_untranslated -p Bris -d two.sided -m fisher -a 0.9

./data/processed/2007-2008/Variants/all.sum.csv ./data/processed/2007-2008/deepSNV/all.coverage.csv ./data/processed/2007-2008/deepSNV/*.fa: $(fastq20072008)
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/2007-2008/ -o ./data/processed/2007-2008/ -r ./data/reference/Brisbane_seq_untranslated -p Brisbane -d two.sided -m fisher -a 0.9

./data/processed/Run_1412/Variants/all.sum.csv ./data/processed/Run_1412/deepSNV/all.coverage.csv ./data/processed/Run_1412/deepSNV/*.fa: $(fastqRun1412)
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/Run_1412/ -o ./data/processed/Run_1412/ -r ./data/reference/CalH3N2_untranslated -p CAcontrol -d two.sided -m fisher -a 0.9

./data/processed/2004-2005/Variants/all.sum.csv ./data/processed/2004-2005/deepSNV/all.coverage.csv ./data/processed/2004-2005/deepSNV/*.fa: $(fastq20042005) 
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/2004-2005/ -o ./data/processed/2004-2005/ -r ./data/reference/CalH3N2_untranslated -p Cal_H3N2 -d two.sided -m fisher -a 0.9

./data/processed/2005-2006/Variants/all.sum.csv ./data/processed/2005-2006/deepSNV/all.coverage.csv ./data/processed/2005-2006/deepSNV/*.fa: $(fastq20052006)
	python ~/variant_pipeline/bin/variantPipeline.py -i ./data/raw/2005-2006/ -o ./data/processed/2005-2006/ -r ./data/reference/CalH3N2_untranslated -p Cal_H3N2 -d two.sided -m fisher -a 0.9


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
./data/processed/Run_1293/HA.fa ./data/processed/Run_1293/NR.fa: ./data/processed/Run_1293/deepSNV/*.fa ./data/processed/concat_pos_bris.csv
	python ./scripts/consensus.pipe.py ./data/processed/Run_1293/deepSNV/ ./data/reference/Brisbane_H3N2_plasmids.fa ./data/processed/concat_pos_bris.csv $(musclePath)

./data/processed/Run_1304/HA.fa ./data/processed/Run_1304/NR.fa: ./data/processed/Run_1304/deepSNV/*.fa ./data/processed/concat_pos_bris.csv
	python ./scripts/consensus.pipe.py ./data/processed/Run_1304/deepSNV/ ./data/reference/Brisbane_H3N2_plasmids.fa ./data/processed/concat_pos_bris.csv $(musclePath)

./data/processed/Run_1412/HA.fa ./data/processed/Run_1412/NR.fa : ./data/processed/Run_1412/deepSNV/*.fa ./data/processed/concat_pos_CalH3N2.csv
	cp -v ./data/processed/2005-2006/deepSNV/536.removed.fa ./data/processed/Run_1412/deepSNV/536.removed.fa
	python ./scripts/consensus.pipe.py ./data/processed/Run_1412/deepSNV/ ./data/reference/CalH3N2SeqCont.fa ./data/processed/concat_pos_CalH3N2.csv $(musclePath)

./data/processed/2004-2005/HA.fa ./data/processed/2004-2005/NR.fa : ./data/processed/2004-2005/deepSNV/*.fa ./data/processed/concat_pos_CalH3N2.csv
	python ./scripts/consensus.pipe.py data/processed/2004-2005/deepSNV/ ./data/reference/CalH3N2SeqCont.fa ./data/processed/concat_pos_CalH3N2.csv $(musclePath)

./data/processed/2005-2006/HA.fa ./data/processed/2005-2006/NR.fa : ./data/processed/2005-2006/deepSNV/*.fa ./data/processed/concat_pos_CalH3N2.csv
	python ./scripts/consensus.pipe.py data/processed/2005-2006/deepSNV/ ./data/reference/CalH3N2SeqCont.fa ./data/processed/concat_pos_CalH3N2.csv $(musclePath)

# get the start and stop of each segment in the concatenated genome.
./data/processed/concat_pos_bris.csv: ./data/processed/Run_1293/deepSNV/all.coverage.csv ./data/processed/Run_1304/deepSNV/all.coverage.csv
	Rscript --vanilla scripts/get_concat_pos.R ./data/processed/Run_1293/deepSNV/all.coverage.csv ./data/processed/Run_1304/deepSNV/all.coverage.csv ./data/processed/concat_pos_bris.csv

./data/processed/concat_pos_CalH3N2.csv: ./data/processed/Run_1412/deepSNV/all.coverage.csv data/processed/2004-2005/deepSNV/all.coverage.csv  data/processed/2005-2006/deepSNV/all.coverage.csv
	Rscript --vanilla scripts/get_concat_pos.R ./data/processed/Run_1412/deepSNV/all.coverage.csv data/processed/2004-2005/deepSNV/all.coverage.csv  data/processed/2005-2006/deepSNV/all.coverage.csv ./data/processed/concat_pos_CalH3N2.csv
