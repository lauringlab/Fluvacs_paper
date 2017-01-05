require(plyr)
require(knitr)
require(reshape2)
#require(cowplot)
require(ggdendro)
require(grid)


######## Procesing #########


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
  data.df=mutate(data.df,id_pos=paste(chr,pos,Id,sep="_"))
  x<-ddply(data.df,~id_pos,summarize,total.freq=sum(freq.var)) # what is the total frequency of all the variants found at this position.
  data.df=mutate(data.df, total.freq=x$total.freq[match(id_pos,x$id_pos)])
  
  to_infer<-subset(data.df,total.freq>cut.low & total.freq<cut.high)
  
  infered<-ddply(to_infer,~id_pos,infer)
  #print(infered)
  data.df<-rbind(data.df,infered)
  data.df<-subset(data.df,select=-c(id_pos))
}

processing<-function(data.df,meta.df,pval,phred,mapq,read_cut,recip=T){
  data.df.cut<-subset(data.df,MapQ>mapq & Phred>phred & Read_pos <read_cut[2] & Read_pos>read_cut[1] & p.val<pval) # subset dataframe
  
  data.df.cut<-ddply(data.df.cut,~chr,coding.cut) # interogating only the sites within coding regions
  # Removed for now as it looks like we just have the coding regions here
  
  ### work with meta data and then join data frames
  #meta<-rename(meta.df,c("Sample"="Id"))
  
  meta<-meta.df
  data.df.cut<-join(data.df.cut,meta,by="Id",type='left')
  
  #data.df.cut<-mutate(data.df.cut,Lauring_Id=Id,Id=paste(HOUSE_ID,season,Id,VAX,onset,sep="."))
  
  if(recip==T){
    data.df.cut<-infer_all(data.df.cut,0.01,0.99) # infer minor variants with major frequency at or below 99%
  }
  return(data.df.cut)
}




##### 



##### join duplicates/ high quality #####


join_dups<-function(data1.df,data2.df){
  data1.df<-mutate(data1.df,samp_mut=paste0(Id,mutation))
  data2.df<-mutate(data2.df,samp_mut=paste0(Id,mutation))
  #data2.df<-rename(data2.df,c("freq.var"="freq.var2","p.val"="p.val2","MapQ"="MapQ2","Phred"="Phred2","Read_pos"="Read_pos2","sigma2.freq.var"="sigma2.freq.var2","n.tst.fw"="n.tst.fw2","cov.tst.fw"="cov.tst.fw2","n.tst.bw"="n.tst.bw2","cov.tst.bw"="cov.tst.bw2","n.ctrl.fw"="n.ctrl.fw2","cov.ctrl.fw"="cov.ctrl.fw2","n.ctrl.bw"="n.ctrl.bw2","cov.ctrl.bw"="cov.ctrl.bw2","raw.p.val"="raw.p.val2"))
  dups.df<-merge(data1.df,data2.df,type='inner',by=c("samp_mut","chr","pos","ref","var","Id","mutation", "Copy_num","Vax","season", "Intervention","HAI.Uruguay.preseason","HAI.Uruaguay.30.post.vax", "HAI.WI.vax.preseason","HAI.WI.30.post.vax","NAI.WI.vax.preseason", "NAI.WI.30.post.vax","Day.of.Infection.sample.collected","collection_date")) # inner so that only shared rows are kept
  #dups.df<-mutate(dups.df,freq.var=rowMeans(cbind(freq.var.x,freq.var.y)))
  #print(names(dups.df))
}


mean_dups<-function(data.df,ref.names){
  
  for(i in 1:length(ref.names)){
    matches<-which(grepl(paste0("^",ref.names[i]),names(data.df))==T) # the ^ makes sure the names starts with the matched name 
    #print(ref.names[matches])
    if(length(matches)==2){
      data.df<-mutate(data.df, new_col = rowMeans(cbind(data.df[,matches[1]],data.df[,matches[2]])))
      data.df<-rename(data.df,c("new_col"=ref.names[i]))
    }
  }
  return(data.df)
  
}

high_qual<-function(data1.df,dups.df,titer){
  htiter<-subset(data1.df,Copy_num>=1e5)
  dups.df<-mean_dups(dups.df,names(data1.df))
  dups.df<-subset(dups.df,select=c(names(data1.df)))
  #print(dups.df$mutation[dups.df$Id=="2"])
  all.df<-rbind(htiter,dups.df)
  subset(all.df,Copy_num>=titer)
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

