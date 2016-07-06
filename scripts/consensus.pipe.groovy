
load './fasta.config.groovy'

parse= {
	doc "Take the concatenated consensus fasta file from deepSNV and deconcatenate it using the segmented positions in from the coverage file"
	output.dir = "${MAIN_DIR}/parsed_fa"
	filter("parsed"){
		exec "${SCRIPTS}/parse_consensus.py $input ${CONCAT_CSV} $output"
	}
}

coding = {
	doc "Trim the parsed sequences to just the coding regions defined by a separate fasta file which contains the coding region for each segment."
	output.dir = "${MAIN_DIR}/coding_fa"
	filter("coding"){
		exec "${SCRIPTS}/trim_to_coding.py ~/muscle3.8.31/ $input ${CODING_FA} -out_fa $output"
	}

}

concatenate = {
	doc "Concatenate desired segments from each sample into one file and use meta data to name the sequences in that file accordingly."
	produce("${MAIN_DIR}/*.fa"){
	exec "${SCRIPTS}/concat_seg.py ${MAIN_DIR}/coding_fa/ HA ./data/raw/2007_2008.meta.HAgm.csv ${MAIN_DIR}/HA.fa"
	exec "${SCRIPTS}/concat_seg.py ${MAIN_DIR}/coding_fa/ NR ./data/raw/2007_2008.meta.HAgm.csv ${MAIN_DIR}/NR.fa"
	}
}

clean_up = {
	doc "clean up the files that were put in the working directory"
	exec "mv ./fasta.config.groovy ${MAIN_DIR}/fasta.config.groovy"
	exec "mv ./commandlog.txt ${MAIN_DIR}/consensuscommandlog.txt"
}

run{
	"%.fa"*[parse+coding] + concatenate + clean_up

}
