######## Procesing #########
processing<-function(data.df,meta.df,pval,phred,mapq,read_cut){
  data.df.cut<-subset(data.df,MapQ>mapq & Phred>phred & Read_pos <read_cut[2] & Read_pos>read_cut[1] & p.val<pval) # subset dataframe
  #meta<-subset(meta,select=c(LAURING_ID,HOUSE_ID,AGE,onset,collect,VAX,season,FINAL_RESULT,gc_ul))
  ### work with meta data and then join data frames
  meta<-rename(meta.df,c("Sample"="Id")) 
  
  data.df.cut<-join(data.df.cut,meta,by="Id",type="inner")
  data.df.cut<-mutate(data.df.cut,Lauring_Id=as.numeric(as.character(Id)),Id=paste(season,Id,Vax,sep="."))
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
    matches<-which(grepl(paste0("^",ref.names[i]),names(data.df))==T)
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

########## ############





####### HEATMAP ######

make_heat_map<-function(data.df){
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
  print(head(mdf))
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
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +  theme_none + theme(axis.title.x=element_blank())
  
  grid.newpage()
  print(p2, vp=viewport(0.94, 0.2, x=0.465, y=0.9))
  print(p1, vp=viewport(0.9, 0.85, x=0.45, y=0.4))
}


