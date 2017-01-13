require(plyr)
require(knitr)
require(reshape2)
#require(cowplot)
require(ggdendro)
require(grid)


######## Procesing #########
read_rbind<-function(list,list2=NULL){
  out<-data.frame()
  for(i in 1:length(list)){
    print(paste0("reading in ",list[i]))
    
    x<-read.csv(list[i],stringsAsFactors = F)
    if(is.null(list2)==F){
      x$run=list2[i]
      print(paste0("appending run column: ",list2[i]))
    }
    out<-rbind(x,out)
  }
  return(out)
}


regions.bed <-read.csv("../data/processed/bis_difference.csv",stringsAsFactors = F,comment.char = "#") # these are the untranslated regions from the brisbane sequence

coding.cut<-function(x){ # a helper function to remove the variants that lie in the primer regions
  chr<-unique(x$chr)
  start<-regions.bed$off.5[match(x$chr,regions.bed$chr)]
  stop<-regions.bed$off.3[match(x$chr,regions.bed$chr)]
  
  subset(x,pos>start & pos<stop)
}


infer<-function(x){ # helper function that identifies the variants that need to be infered.
  
  x<-mutate(x,var=ref,freq.var=1-total.freq,mutation=paste0(chr,"_",ref,pos,var))
  if(dim(x)[1]>1){ # if there are multiple called variants here, then only take the top one as the infered consensus. The coverage and will be the same
    x<-x[1,]
  }
  return(x)
  #x[,c(5,7:16,19:21)]<-NA
  #print(x)
}

infer_all<-function(data.df,cut.low,cut.high){ 
  data.df=mutate(data.df,id_pos=paste(chr,pos,Id,run,sep="_"))
  x<-ddply(data.df,~id_pos,summarize,total.freq=sum(freq.var)) # what is the total frequency of all the variants found at this position.
  print(x)
  data.df=mutate(data.df, total.freq=x$total.freq[match(id_pos,x$id_pos)])
  
  to_infer<-subset(data.df,total.freq>cut.low & total.freq<cut.high)
  
  infered<-ddply(to_infer,~id_pos,infer)
  #print(infered)
  data.df<-rbind(data.df,infered)
  data.df<-subset(data.df,select=-c(id_pos))
}

processing<-function(data.df,meta.df,pval,phred,mapq,read_cut,recip=T,gc){
  data.df.cut<-subset(data.df,MapQ>mapq & Phred>phred & Read_pos <read_cut[2] & Read_pos>read_cut[1] & p.val<pval) # subset dataframe
  
  data.df.cut<-ddply(data.df.cut,~chr,coding.cut) # interogating only the sites within coding regions
  # Removed for now as it looks like we just have the coding regions here
  
  ### work with meta data and then join data frames
  #meta<-rename(meta.df,c("Sample"="Id"))
  
  data.df.cut<-join(data.df.cut,meta.df,by="Id",type='left')
  data.df.cut<-subset(data.df.cut,Copy_num>gc)
  
  #data.df.cut<-mutate(data.df.cut,Id=Id,Id=paste(HOUSE_ID,season,Id,VAX,onset,sep="."))
  
  if(recip==T){
    data.df.cut<-infer_all(data.df.cut,0.01,0.99) # infer minor variants with major frequency at or below 99%
  }
  return(data.df.cut)
}




##### 



##### join duplicates/ high quality #####

sift_dups<-function(df){
  mean_wanted<-c("p.val","freq.var","sigma2.freq.var","n.tst.fw","cov.tst.fw","n.tst.bw","cov.tst.bw","n.ctrl.fw","cov.ctrl.fw","n.ctrl.bw","cov.ctrl.bw","raw.p.val","MapQ","Read_pos","Phred","total.freq")
  if(dim(df)[1]==2){ #it's found in both duplicates
    x<-df[1,] # grab the first row
    check<-df[,names(df)[!(names(df) %in% c(mean_wanted,"run"))]] # the run should be different
    if(any(check[1,]!=check[2,],na.rm = T)){
      stop(paste("Meta data doesn't match for sample", unique(x$Id)))
    }
    x[mean_wanted]<-colMeans(df[mean_wanted])
    x$run<-paste(df$run[1],df$run[2],sep=".")
    return(x)
  } 
  try(if(dim(df)[1]>2) stop(paste0("This mutation is found more than once in this sample : ", unique(df$Id), " - ", unique(df$mutation))))
  # add if statement to through error if the dim is greater than 2 
}

join_dups<-function(df){
  ddply(df,~Id+mutation,sift_dups)
}

quality<-function(df,freq_cut){
  starting_samples<-unique(df$Id) # samples we start with
  #good<-subset(df,gc_ul>=1e5)
  # check only one run/sample
  runs<-ddply(df,~Id,summarize,runs=length(unique(run))) # Id is unique to the sample. How many times was each sample sequenced? duplicates were sequenced in separate runs
  doubles<-subset(runs,runs==2) # Theses were sequenced twice
  dups_ran<-subset(df,Id %in% doubles$Id)
  
  dups_good<-join_dups(dups_ran) # We only want those sequenced twice with high enough titers.
  
  good<-subset(df,Copy_num>1e5 & !(Id %in% dups_ran$Id)) # Just incase some with high titers were sequenced twice
  out<-rbind(good,dups_good)
  
  ## filter to coding regions
  
  
  
  out<-subset(out,freq.var>freq_cut[1] & freq.var<freq_cut[2]) # filter on frequency
  # Were any with high titers sequenced twice?
  
  high_twice<-subset(df,Copy_num>1e5 & (Id %in% dups_ran$Id))
  if(dim(high_twice)[1]>0){
    print(paste0("Sample ",as.character(unique(high_twice$Id)), " was sequenced twice even though the titer was ",  as.character(unique(high_twice$Copy_num)),". It was treated as duplicates in the analysis"))
  }
  print(paste0("We removed ",dim(df)[1]-dim(out)[1]," variants of ",dim(df)[1]," and ", dim(out)[1], " remain."))
  
  final_samples<-unique(out$Id)
  
  if(length(which(!(starting_samples %in% final_samples)))>0){
    samples=starting_samples[which(!(starting_samples %in% final_samples))]
    runs<-c()
    for( i in 1:length(samples)){
      sub<-subset(df,Id==samples[i],select=c(run))
      runs[i]=paste(unique(sub$run),sep=".")
    }
    missing_table<-data.frame(id = samples,
      Copy_num=unique(df$Copy_num[df$Id%in%samples]),
      runs= runs)
    print(paste0("After processing no variants were found the following samples: "))
    print(missing_table)
  }
  
  return(out)
}

slide_window<-function(freq.df,window,step){
  
  
  spots <- seq(from=0, to=1, by=step)
  result <- data.frame(mean.freq=rep(F,times=length(spots)),number.donor=rep(F,times=length(spots)),number.transmitted=rep(F,times=length(spots)),prob.transmitted=rep(F,times=length(spots)))
  for(i in 1:(length(spots)-1)){
    
    var.df<-subset(freq.df,freq1>spots[i] & freq1<(spots[i]+window))
    #result$mean.freq[i] <- mean(c(spots[i],(spots[i]+window)))
    result$mean.freq[i] <- mean(var.df$freq1)
    
    result$number.donor[i]<-dim(var.df)[1]
    result$number.transmitted[i]<-length(which(var.df$freq2>0.002))
    result$prob.transmitted[i]<-result$number.transmitted[i]/result$number.donor[i]
  }
  return(result)
  
}

slide_window.reciep<-function(window,step,freq.df){
  
  
  spots <- seq(from=0, to=1, by=step)
  result <- data.frame(mean.freq=rep(F,times=length(spots)),number.recipt=rep(F,times=length(spots)),number.transmitted=rep(F,times=length(spots)),prob.transmitted=rep(F,times=length(spots)))
  for(i in 1:(length(spots)-1)){
    result$mean.freq[i] <- mean(c(spots[i],(spots[i]+window)))
    var.df<-subset(freq.df,freq2>=spots[i] & freq1<=(spots[i]+window))
    result$number.donor[i]<-dim(var.df)[1]
    result$number.transmitted[i]<-length(which(var.df$freq1>0.002))
    result$prob.transmitted[i]<-result$number.transmitted[i]/result$number.donor[i]
    #result$concat.pos[i]<-mean(c.pos[spots[i]:(spots[i]+window)])
  }
  return(result)
}



slide_window.pos<-function(freq.df,window,step){ # ddply by segment 
  end<-max(freq.df$pos)+step
  spots <- seq(from=0, to=end, by=step)
  result <- data.frame(meaan.pos=rep(F,times=length(spots)),number.donor=rep(F,times=length(spots)),number.transmitted=rep(F,times=length(spots)),prob.transmitted=rep(F,times=length(spots)))
  for(i in 1:(length(spots)-1)){
    
    var.df<-subset(freq.df,pos>spots[i] & pos<(spots[i]+window))
    result$mean.pos[i] <- mean(c(spots[i],(spots[i]+window))) # just the middle of the window ( could be change to reflect the actual positions)
    #    result$mean.pos[i] <- mean(var.df$pos)
    
    result$number.donor[i]<-dim(var.df)[1]
    result$number.transmitted[i]<-length(which(var.df$freq2>0.002))
    result$prob.transmitted[i]<-result$number.transmitted[i]/result$number.donor[i]
  }
  return(result)
  
}




########## Data handelling #########

bin_freq<-function(df,col,binsize){
  total<-dim(df)[1]
  bins=seq(0,(max(df[,col])+binsize),by=binsize)
  count=c()
  freq=c()
  for(i in 1:length(bins)){
    count[i]=length(which(df[,col]>bins[i] & df[,col]<bins[i+1]))
    freq[i]=length(which(df[,col]>bins[i] & df[,col]<bins[i+1]))/total
  }
  as.data.frame(cbind(bins,count,freq))
}

######## Coverage #########

slideFunct <- function(data, window, step){ #dapted from http://coleoguy.blogspot.com/2014/04/sliding-window-analysis.html
  coverage<-data$coverage
  concat.pos<-data$concat.pos
  # c.pos<-data$concat.pos
  total <- length(coverage)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- data.frame(mean=rep(F,times=length(spots)),concat.pos=rep(F,times=length(spots)))
  for(i in 1:length(spots)){
    result$mean[i] <- mean(coverage[spots[i]:(spots[i]+window)])
    result$concat.pos[i]<-mean(concat.pos[spots[i]:(spots[i]+window)])
    #result$concat.pos[i]<-mean(c.pos[spots[i]:(spots[i]+window)])
  }
  return(result)
}


cov_plot<-function(cov.df,title,window,step){
  
  cov.slid.df<-ddply(cov.df,~Id+chr,function(x) slideFunct(x,window = window,step=step))
  
  x.labels<-ddply(cov.slid.df,~chr,summarize,concat.pos=concat.pos[which(abs(concat.pos-mean(concat.pos))==(min(abs(concat.pos-mean(concat.pos)))))]) # No idea if this will get what I want but heres to hoping!
  
  #of course sometimes there are 2 good choices I'll take the first one
  x.labels<-ddply(x.labels,~chr,function(x) return(x[1,]))
  
  x.labels$chr[x.labels$chr %in% c("NR","N_A")]<-"NA"
  
  cov.plot<-ggplot(cov.slid.df, #subset(cov.slid.df,!(Sample%in%c("90","91","93"))),
                   mapping=aes(x=as.factor(concat.pos),
                               y=mean))+geom_boxplot(fill="white")
  
  cov.plot<-cov.plot+ggtitle(title)+ylab("Read depth")+scale_x_discrete(labels = x.labels$chr,breaks=x.labels$concat.pos)+xlab("Concatenated Genome Position")
  cov.plot<-cov.plot+theme(axis.title.y = element_text(vjust=1.2))
  cov.plot<-cov.plot+theme(legend.position="none")  + scale_y_continuous(limits=c(0,88000))
  return(cov.plot)
}

