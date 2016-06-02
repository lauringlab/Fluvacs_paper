require(plyr)
require(knitr)
require(reshape2)
#require(cowplot)
require(ggdendro)
require(grid)

####### HEATMAP ######

make_heat_map<-function(data.df,s,title){ # s=season
  data.df<-subset(data.df,season==s| is.na(season))
  mut_count <- ddply(data.df, ~mutation, summarize, sample_count = length(mutation))
  
  pop_mut_count <- subset(mut_count, sample_count > 0)
  
  popular_muts <- subset(data.df, mutation %in% pop_mut_count$mutation)
  
  
  
  mut_table <- dcast(popular_muts, mutation ~ Id, value.var = "freq.var")
  mut_table[is.na(mut_table)] <- 0
  
  pat <- "([A-Z]+[.,1,2]?)_[A,T,C,G](\\d+)[A,T,C,G]"
  # pat<-'[A,T,C,G](\d+)[A,T,C,G]' if the mutation doesn't include a segment
  mut_table <- mutate(mut_table, pos = sub(pat, "\\2", mutation, perl = T), segment = sub(pat,
                                                                                          "\\1", mutation, perl = T))  # make \\1 if no segment
  
  mut_table$segment <- factor(mut_table$segment, levels = c("PB2", "PB1", "PA","HA", "NP", "NR", "M", "NS"))
  
  mut_table <- mut_table[order(c(mut_table$segment), as.numeric(mut_table$pos)),]
  
  
  
  x <- as.matrix(subset(mut_table, select = -c(pos, mutation, segment)))
  # dd.col <- as.dendrogram(hclust(dist(x))) # in our case this would cluster
  # the mutations
  
  dd.row <- as.dendrogram(hclust(as.dist(1 - cor(x))))
  #dd.row <- as.dendrogram(hclust(dist(x = x,method="manhattan"))) This way needs works
  row.ord <- order.dendrogram(dd.row)
  mut_table.no.meta <- subset(mut_table, select = -c(pos, mutation, segment))
  xx <- mut_table.no.meta[, row.ord]
  xx_names <- list(mut_table$mutation, names(xx))
  df <- as.data.frame(xx)
  colnames(df) <- xx_names[[2]]
  df$mutation <- xx_names[[1]]
  df$mutation <- with(df, factor(mutation, levels = mutation, ordered = TRUE))
  
  mdf <- melt(df, id.vars = "mutation")
  mdf<-mutate(mdf,segment = sub(pat, "\\1", mutation, perl = T))
  #print(head(mdf))
  mdf$value<-mdf$value*100
  
  mdf$value[(mdf$value==0)]<- 0.1
  mdf$value<-log10(mdf$value)
  #print(head(mdf))
  ddata_x <- dendro_data(dd.row)
  #ddata_y <- dendro_data(dd.col)
  
  ### Set up a blank theme
  theme_none <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour=NA),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank()
    #axis.ticks.length = element_blank()
  )
  
  
  my_colours<-colorRampPalette(c("blue4","yellow"))(n=200)
  breaks=c(seq(-1,0,length=1),seq(0,2,length=200))
  
  
  p1 <- ggplot(mdf, aes(x=variable, y=mutation)) +
    geom_tile(aes(fill=value))+
    scale_fill_gradientn(colours = my_colours, breaks=breaks)+#),trans="log10", limits=c(0.01,1))
    xlab("")+theme(axis.text.x=element_text(angle=90,vjust=0.6,size=8))+theme(legend.position="none")
  
  
  
  
  p1<-p1+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  p2 <- ggplot(segment(ddata_x)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +  theme_none + theme(axis.title.x=element_blank())+ggtitle(title)
  
  grid.newpage()
  print(p2, vp=viewport(0.94, 0.2, x=0.465, y=0.9))
  print(p1, vp=viewport(0.9, 0.85, x=0.45, y=0.4))
}
######## Procesing #########


regions.bed <-read.csv("../data/pr8_noprimingsites.bed.csv",stringsAsFactors = F,comment.char = "#") # these sites should be the same for Pr8 and these viruses but I'll check that later

primer.cut<-function(x){ # a helper function to remove the variants that lie in the primer regions
  chr<-unique(x$chr)
  start<-regions.bed$start[match(x$chr,regions.bed$chr)]
  stop<-regions.bed$stop[match(x$chr,regions.bed$chr)]
  
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

infer_all<-function(data.df,cut.low,cut.high){ # maybe coverage will be needed here as well.
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
  
  data.df.cut<-ddply(data.df.cut,~chr,primer.cut) # interogating only the sites within primer location 
  
  ### work with meta data and then join data frames
  meta<-rename(meta.df,c("Sample"="Id"))
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
  dups.df<-merge(data1.df,data2.df,type='inner',by=c("samp_mut","chr","pos","ref","var","Id","mutation","season","Copy_num","Vax"))
  #dups.df<-mutate(dups.df,freq.var=rowMeans(cbind(freq.var.x,freq.var.y)))
  #print(names(dups.df))
}


mean_dups<-function(data.df,ref.names){
  
  for(i in 1:length(ref.names)){
    matches<-which(grepl(paste0("^",ref.names[i]),names(data.df))==T) # the ^ makes sure the names starts with the matched name 
    #print(ref.names[i])
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
  all.df<-rbind(htiter,dups.df)
  subset(all.df,Copy_num>=titer)
}

###############################################################
########### Transmission functions ############################
###############################################################

# dist_tp<-function(data.df){ # implement with ddply to split by season and run
#   mut_table <- dcast(data.df, mutation ~ Id, value.var = "freq.var")
#   mut_table[is.na(mut_table)] <- 0
#   out<-data.frame(dist=c(),pair=c(),Ids=c())
#   for( i in 2:(dim(mut_table)[2]-1)){ # the first column is the mutation name
#     samp_home<-strsplit(names(mut_table)[i],split = ".",fixed = T)[[1]][1] 
#     samp_season<-strsplit(names(mut_table)[i],split = ".",fixed = T)[[1]][2]
#     samp_onset<-as.Date(strsplit(names(mut_table)[i],split = ".",fixed = T)[[1]][5])
#     samp_samp<-strsplit(names(mut_table)[i],split = ".",fixed = T)[[1]][3]
#     for( j in (i+1):dim(mut_table)[2]){
#       test_home<-strsplit(names(mut_table)[j],split = ".",fixed = T)[[1]][1]
#       test_season<-strsplit(names(mut_table)[j],split = ".",fixed = T)[[1]][2]
#       test_onset<-as.Date(strsplit(names(mut_table)[j],split = ".",fixed = T)[[1]][5])
#       test_samp<-strsplit(names(mut_table)[j],split = ".",fixed = T)[[1]][3]
#       if(samp_home==test_home & samp_season==test_season & abs(samp_onset-test_onset)<=7){
#         x<-mut_table[,c(i,j)]
#         l1<-dist(t(x),method = "manhattan")[1]
#         #print(l1)
#         pair="Household"
#         Id1=samp_samp
#         Id2=test_samp
#       }else {  # Think about excluding in home cases more than a week apart with else if(samp_home!=test_home)
#         x<-mut_table[,c(i,j)]
#         l1<-dist(t(x),method = "manhattan")[1]
#         #print(l1)
#         pair="Community"
#         Id1=samp_samp
#         Id2=test_samp
#       }
#       out.loop<-data.frame(dist=l1,pair=pair,Ids=paste(Id1,Id2,sep="_"))
#       out<-rbind(out,out.loop)
#     }
#   }
# return(out)
# } 






freq_tp<-function(data.df){
  mut_table <- dcast(data.df, mutation ~ Id, value.var = "freq.var")
  mut_table[is.na(mut_table)] <- 0
  #mut_table<-mut_table[,2:dim(mut_table)[2]] # remove the plasmid control
  
  freqs<-data.frame(mutation=c(),freq1=c(),freq2=c(),pair=c(),Id1=c(),Id2=c())
  for( i in 2:(dim(mut_table)[2]-1)){ # the first column is the mutation name
    samp_home<-strsplit(names(mut_table)[i],split = ".",fixed = T)[[1]][1] 
    samp_season<-strsplit(names(mut_table)[i],split = ".",fixed = T)[[1]][2]
    samp_onset<-as.Date(strsplit(names(mut_table)[i],split = ".",fixed = T)[[1]][5])
    samp_samp<-strsplit(names(mut_table)[i],split = ".",fixed = T)[[1]][3]
    for( j in (i+1):dim(mut_table)[2]){
      test_home<-strsplit(names(mut_table)[j],split = ".",fixed = T)[[1]][1]
      test_season<-strsplit(names(mut_table)[j],split = ".",fixed = T)[[1]][2]
      test_onset<-as.Date(strsplit(names(mut_table)[j],split = ".",fixed = T)[[1]][5])
      test_samp<-strsplit(names(mut_table)[j],split = ".",fixed = T)[[1]][3]
      if(samp_home==test_home & samp_season==test_season & abs(samp_onset-test_onset)<=7){
        x<-mut_table[,c(1,i,j)]
        x<-x[rowSums(x[,c(2,3)])!=0,] # elimate the cases where no vairant is found
        x$pair="Household"
        x$Id1=samp_samp
        x$Id2=test_samp
        if(samp_onset<=test_onset){
          names(x)<-c("mutation", "freq1","freq2","pair","Id1","Id2")
        }else{
          names(x)<-c("mutation","freq2","freq1","pair","Id2","Id1")
        }
        freqs<-rbind(freqs,x)
      }else if(samp_home!=test_home & samp_season==test_season){
        x<-mut_table[,c(1,i,j)]
        x<-x[rowSums(x[,c(2,3)])!=0,] # elimate the cases where no vairant is found
        x$pair="Community"
        x$Id1=samp_samp
        x$Id2=test_samp
        if(samp_onset<=test_onset){
          names(x)<-c("mutation", "freq1","freq2","pair","Id1","Id2")
        }else{
          names(x)<-c("mutation","freq2","freq1","pair","Id2","Id1")
        }
        freqs<-rbind(freqs,x)
        
      }
    }
  }
  pat="([A-Z]+[1-2]?)_([A|C|T|G]{1})([0-9]+)([A|C|T|G]{1})"
  freqs=mutate(freqs,chr=as.character(sub(pat,"\\1",mutation)),ref = sub(pat,"\\2",mutation),pos=as.numeric(sub(pat,"\\3",mutation)),var=sub(pat,"\\4",mutation))  
  freqs$chr[is.na(freqs$chr)]<-as.character("NR")
  return(freqs)
}

# freq_randompairs<-function(data.df){
#   mut_table <- dcast(data.df, mutation ~ Id, value.var = "freq.var")
#   mut_table[is.na(mut_table)] <- 0
#   #mut_table<-mut_table[,2:dim(mut_table)[2]] # remove the plasmid control
#   
#   freqs<-data.frame(freq1=c(),freq2=c())
#   pairs<-1
#   for( i in 2:(dim(mut_table)[2]-1)){ # the first column is the mutation name
#     samp_home<-strsplit(names(mut_table)[i],split = ".",fixed = T)[[1]][1]
#     samp_season<-strsplit(names(mut_table)[i],split = ".",fixed = T)[[1]][2]
#     samp_onset<-strsplit(names(mut_table)[i],split = ".",fixed = T)[[1]][5]
#     for( j in (i+1):dim(mut_table)[2]){
#       test_home<-strsplit(names(mut_table)[j],split = ".",fixed = T)[[1]][1]
#       test_season<-strsplit(names(mut_table)[j],split = ".",fixed = T)[[1]][2]
#       test_onset<-as.Date(strsplit(names(mut_table)[j],split = ".",fixed = T)[[1]][5])
#       if(samp_home!=test_home & samp_season==test_season){
#         x<-mut_table[,c(i,j)]
#         x<-x[rowSums(x)!=0,]
#         if(samp_onset<test_onset){
#           names(x)<-c("freq1","freq2")
#         }else{
#           names(x)<-c("freq2","freq1")
#         }
#         freqs<-rbind(freqs,x)
#         pairs<-pairs+1
#       }
#     }
#   }
#   return(freqs)
# }








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

