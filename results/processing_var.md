---
title: "Fluvacs Processing"
author: "JT McCrone"
date: "July 7, 2016"
output: html_document
---

This document will run the first analysis on the variant calls. It will take in the csv files from the variant calling pipeline and filter out low quality calls as well as handle any samples done in duplicate. This script will then output 3 csv files. One that contains all minor quality variants. This will be used to make the figures. The second will just include quality calls from the HA segement. This will be used to identify ptutative antigenic variants. These go in the results directory as they are probably interesting for others as tables. The third is an improved metadata csv which includes information regarding the geometric mean of antibody titer and the day of collection in year.decimal format.





```r
meta.2007.08 <- read.csv("../data/raw/all.meta.csv", stringsAsFactors = F)
```


```r
var.2007.8 <- read.csv("../data/processed/Run_1293/Variants/all.sum.csv", stringsAsFactors = F)
x <- read.csv("../data/processed/Run_1304/Variants/all.sum.csv", stringsAsFactors = F)  # the rest of these samples
var.2007.8 <- rbind(var.2007.8, x)  # combine both runs
# other.seasons<-read.csv('../data/processed/Run_1412/Variants/all.sum.csv',stringsAsFactors
# = F)


# var.2004.5.df<-processing(data.df = other.seasons,meta.df =
# titer.2004.5,pval = 0.01,phred = 35,mapq = 30,read_cut = c(32,94))
# var.2004.5.df<-subset(var.2004.5.df,season=='04-05')
# var.2005.6.df<-processing(data.df = other.seasons,meta.df =
# titer.2005.6,pval = 0.01,phred = 35,mapq = 30,read_cut = c(32,94))
# var.2005.6.df<-subset(var.2005.6.df,season=='05-06')
var.2007.8.df <- processing(data.df = var.2007.8, meta.df = meta.2007.08, pval = 0.01, 
    phred = 35, mapq = 30, read_cut = c(32, 94))


# all.df<-rbind(var.2004.5.df,var.2005.6.df)
# all.df<-rbind(all.df,var.2007.8.df) all.df<-subset(all.df,freq.var>=0.01)
# ##### Now for the duplicate runs #######
# var.2004.5.df2<-read.csv('../data/processed/2004_2005/Variants/all.sum.csv',stringsAsFactors
# = F) var.2004.5.df2<-processing(data.df = var.2004.5.df2,meta.df =
# titer.2004.5,pval = 0.01,phred = 35,mapq = 30,read_cut = c(32,94))
# var.2005.6.df2<-read.csv('../data/processed/2005-2006/Variants/all.sum.csv',stringsAsFactors
# = F) var.2005.6.df2<-processing(data.df = var.2005.6.df2,meta.df =
# titer.2005.6,pval = 0.01,phred = 35,mapq = 30,read_cut = c(32,94))

var.2007.8.df2 <- read.csv("../data/processed/2007-2008/Variants/all.sum.csv", 
    stringsAsFactors = F)
var.2007.8.df2 <- processing(data.df = var.2007.8.df2, meta.df = meta.2007.08, 
    pval = 0.01, phred = 35, mapq = 30, read_cut = c(32, 94))

### Join duplicates #####

# dups.2004.5<-join_dups(data1.df = var.2004.5.df,data2.df = var.2004.5.df2)
# dups.2005.6<-join_dups(data1.df = var.2005.6.df,data2.df = var.2005.6.df2)
dups.2007.8 <- join_dups(data1.df = var.2007.8.df, data2.df = var.2007.8.df2)

##### Merge duplicates with intial data that was >1e5

# qual.2004.5<-high_qual(data1.df = var.2004.5.df,dups.df =
# dups.2004.5,titer = 1e3) # only duplicates above 1e3 kept
# qual.2005.6<-high_qual(data1.df = var.2005.6.df,dups.df =
# dups.2005.6,titer = 1e3)
qual.2007.8 <- high_qual(data1.df = var.2007.8.df, dups.df = dups.2007.8, titer = 1000)
# all.qual.df<-rbind(qual.2004.5,qual.2005.6)
# all.qual.df<-rbind(all.qual.df,qual.2007.8)
# all.qual.df<-subset(all.qual.df,freq.var>=0.01 & freq.var<=0.99)

var.2007.8.df <- subset(var.2007.8.df, freq.var >= 0.01 & freq.var <= 0.99)
qual.2007.8 <- subset(qual.2007.8, freq.var >= 0.01 & freq.var <= 0.99)

bris.bed <- read.csv("../data/processed/bis_difference.csv", stringsAsFactors = F, 
    comment.char = "#")

coding.adjust <- function(x) {
    chr <- unique(x$chr)
    start <- bris.bed$off.5[match(x$chr, bris.bed$chr)]
    
    mutate(x, coding.pos = pos - start)
}

qual.2007.8 <- ddply(qual.2007.8, ~chr, coding.adjust)
```


```r
gm_mean = function(x, na.rm = TRUE, zero.propagate = FALSE) {
    # from
    # http://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
    if (any(x < 0, na.rm = TRUE)) {
        return(NaN)
    }
    if (zero.propagate) {
        if (any(x == 0, na.rm = TRUE)) {
            return(0)
        }
        exp(mean(log(x), na.rm = na.rm))
    } else {
        exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
    }
}

# var.2007.8.df<-mutate(var.2007.8.df,responder.HA=HAI.WI.30.post.vax>40,responder.NA=NAI.WI.30.post.vax>40,HAI.geo=HAI.WI.30.post.vax>gm_mean(HAI.WI.30.post.vax),NAI.geo=NAI.WI.30.post.vax>gm_mean(NAI.WI.30.post.vax),responder.both=(responder.NA&responder.HA==T),
# geom.both=(HAI.geo&NAI.geo==T))

qual.2007.8 <- mutate(qual.2007.8, responder.HA = HAI.WI.30.post.vax > 40, responder.NA = NAI.WI.30.post.vax > 
    40, HAI.geo = HAI.WI.30.post.vax > gm_mean(HAI.WI.30.post.vax), NAI.geo = NAI.WI.30.post.vax > 
    gm_mean(NAI.WI.30.post.vax), responder.both = (responder.NA & responder.HA == 
    T), geom.both = (HAI.geo & NAI.geo == T))
minor = subset(qual.2007.8, freq.var <= 0.5)


### For consensus seqeunce meta data

HAI.gm <- ddply(meta.2007.08, ~season, summarise, gm_mean = gm_mean(HAI.post.vax))
NAI.gm <- ddply(meta.2007.08, ~season, summarise, gm_mean = gm_mean(NAI.post.vax))

HAI.cut <- function(x) {
    season <- unique(x$season)
    gm <- HAI.gm$gm_mean[which(HAI.gm$season == season)]
    # print(gm)
    mutate(x, HAI.geo = HAI.post.vax > gm)
}

NAI.cut <- function(x) {
    season <- unique(x$season)
    gm <- NAI.gm$gm_mean[which(NAI.gm$season == season)]
    # print(gm)
    mutate(x, NAI.geo = NAI.post.vax > gm)
}
meta.df <- ddply(meta.2007.08, ~season, HAI.cut)
meta.df <- ddply(meta.df, ~season, NAI.cut)

# meta.df<-mutate(meta.2007.08,HAI.geo=HAI.WI.30.post.vax>gm_mean(HAI.WI.30.post.vax),NAI.geo=NAI.WI.30.post.vax>gm_mean(NAI.WI.30.post.vax))
meta.df$collection_date <- as.Date(meta.df$collection_date, format = "%d-%b-%y")
require(lubridate)
```

```
## Loading required package: lubridate
```

```
## Loading required package: methods
```

```
## 
## Attaching package: 'lubridate'
```

```
## The following object is masked from 'package:plyr':
## 
##     here
```

```
## The following object is masked from 'package:base':
## 
##     date
```

```r
meta.df$collection_date <- decimal_date(meta.df$collection_date)
```

writing the outputs


```r
## Getting data for HA ha<-subset(qual.2007.8,chr=='HA')
ha <- subset(minor, chr == "HA")

dim(ha)
```

```

## [1] 38 45
```

```r
ha <- ha[order(ha$pos), ]
unique(ha$pos)
```

```
##  [1]  222  234  314  440  700  701  717  758  761  791  883  894  899  909
## [15] 1044 1077 1144 1150 1214 1226 1258 1265 1298 1325 1365 1433 1436 1487
## [29] 1505 1532 1569 1589 1595 1643
```

```r
write.csv(x = ha, file = "../results/2007-2008.HA.csv", row.names = F)
write.csv(x = minor, file = "../results/2007-2008.wg.csv", row.names = F)
write.csv(meta.df, "../data/raw/2007_2008.meta.HAgm.csv")
```


