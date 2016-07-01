#!/Library/Frameworks/R.framework/Resources/bin Rscript
## The goal here is to get the concatenated start and stop positions for each segement. I'm getting this information from the coverage csv

# The script takes in a list of arguments. the first n are coverage files from the deepSNV output. The n+1 is the file that will hold the concatenated start and stop positions for the files
args = commandArgs(trailingOnly=TRUE)

require(plyr)
cov_files=length(args)-1 # The number of coverage files we are using

all_cov<-data.frame(coverage=c(), concat.pos=c(), chr=c(), chr.pos=c(), Id=c())

print("Reading coverage files")
for (i in 1:cov_files) {
  read_cov=read.csv(args[i],stringsAsFactors = F)

  if(!(all.equal(names(read_cov),c("coverage","concat.pos","chr","chr.pos","Id")))){
    stop(paste0(args[i]," is not a proper coverage file. Expected column names to match c(coverage,concat.pos,chr,chr.pos,Id).\n"),call.=FALSE)
  }
  all_cov=rbind(all_cov,read_cov)
}

concat.pos<-ddply(all_cov,~chr, summarize, start=min(concat.pos), stop=max(concat.pos))
concat.pos$chr[concat.pos$chr=="N_A"]=as.character('NR')

print("writing output")
con <- file(args[cov_files+1], open="wt") # The last argument
writeLines(paste("# this csv file was created on ",date(), " It represents the ends of the segements aligned to the brisbane plasmid control (years 2007-2008)"), con)
write.csv( concat.pos, con,row.names=F)
close(con)
