---
title: "Fluvacs Processing"
author: "JT McCrone"
date: "July 7, 2016"
output: html_document
---

This document will run the first analysis on the variant calls. It will take in the csv files from the variant calling pipeline and filter out low quality calls as well as handle any samples done in duplicate. This script will then output 3 csv files. One that contains all minor quality variants. This will be used to make the figures. The second will just include quality calls from the HA segement. This will be used to identify putative antigenic variants. These go in the results directory as they are probably interesting for others as tables. The third is an improved metadata csv which includes information regarding the geometric mean of antibody titer and the day of collection in year.decimal format.


Read in the meta data file

```r
meta.all <- read.csv("../data/raw/all.meta.csv", stringsAsFactors = F)
```

Here I read in the variants calls, filter the calls for quality, apply the duplicate analysis on the appropriate sample, and trim for coding regions.


```r
var.2007.8 <- read.csv("../data/processed/Run_1293/Variants/all.sum.csv", stringsAsFactors = F)
x <- read.csv("../data/processed/Run_1304/Variants/all.sum.csv", stringsAsFactors = F)  # the rest of these samples
var.2007.8 <- rbind(var.2007.8, x)  # combine both runs

# Processing involves infering reciprocal variants

var.2007.8.df <- processing(data.df = var.2007.8, meta.df = meta.all, pval = 0.01, 
    phred = 35, mapq = 30, read_cut = c(32, 94))


other.seasons <- read.csv("../data/processed/Run_1412/Variants/all.sum.csv", 
    stringsAsFactors = F)


other.seasons.df <- processing(data.df = other.seasons, meta.df = meta.all, 
    pval = 0.01, phred = 35, mapq = 30, read_cut = c(32, 94))

var.2004.5.df <- subset(other.seasons.df, season == "04-05")
var.2005.6.df <- subset(other.seasons.df, season == "05-06")

# ##### Now for the duplicate runs #######
var.2004.5.df2 <- read.csv("../data/processed/2004-2005/Variants/all.sum.csv", 
    stringsAsFactors = F)
var.2004.5.df2 <- processing(data.df = var.2004.5.df2, meta.df = meta.all, pval = 0.01, 
    phred = 35, mapq = 30, read_cut = c(32, 94))
# 
var.2005.6.df2 <- read.csv("../data/processed/2005-2006/Variants/all.sum.csv", 
    stringsAsFactors = F)
var.2005.6.df2 <- mutate(var.2005.6.df2, Id = gsub(536, 530, Id))  # correct labling error in this run. There is no sample 536
var.2005.6.df2 <- processing(data.df = var.2005.6.df2, meta.df = meta.all, pval = 0.01, 
    phred = 35, mapq = 30, read_cut = c(32, 94))


var.2007.8.df2 <- read.csv("../data/processed/2007-2008/Variants/all.sum.csv", 
    stringsAsFactors = F)
var.2007.8.df2 <- processing(data.df = var.2007.8.df2, meta.df = meta.all, pval = 0.01, 
    phred = 35, mapq = 30, read_cut = c(32, 94))

### Join duplicates #####

dups.2004.5 <- join_dups(data1.df = var.2004.5.df, data2.df = var.2004.5.df2)
dups.2005.6 <- join_dups(data1.df = var.2005.6.df, data2.df = var.2005.6.df2)
dups.2007.8 <- join_dups(data1.df = var.2007.8.df, data2.df = var.2007.8.df2)

##### Merge duplicates with intial data that was >1e5

qual.2004.5 <- high_qual(data1.df = var.2004.5.df, dups.df = dups.2004.5, titer = 1000)  # only duplicates above 1e3 kept
qual.2005.6 <- high_qual(data1.df = var.2005.6.df, dups.df = dups.2005.6, titer = 1000)
qual.2007.8 <- high_qual(data1.df = var.2007.8.df, dups.df = dups.2007.8, titer = 1000)
# all.qual.df<-rbind(qual.2004.5,qual.2005.6)
# all.qual.df<-rbind(all.qual.df,qual.2007.8)
# all.qual.df<-subset(all.qual.df,freq.var>=0.01 & freq.var<=0.99)

# Filter for frequency.
var.2007.8.df <- subset(var.2007.8.df, freq.var >= 0.01 & freq.var <= 0.99)
qual.2007.8 <- subset(qual.2007.8, freq.var >= 0.01 & freq.var <= 0.99)
qual.2004.5 <- subset(qual.2004.5, freq.var >= 0.01 & freq.var <= 0.99)
qual.2005.6 <- subset(qual.2005.6, freq.var >= 0.01 & freq.var <= 0.99)




bris.bed <- read.csv("../data/processed/bis_difference.csv", stringsAsFactors = F, 
    comment.char = "#")

coding.adjust.bris <- function(x) {
    chr <- unique(x$chr)
    start <- bris.bed$off.5[match(x$chr, bris.bed$chr)]
    
    mutate(x, coding.pos = pos - start)
}


# Just the variants in the coding regions.
qual.2007.8 <- ddply(qual.2007.8, ~chr, coding.adjust.bris)


cal.bed <- read.csv("../data/processed/CalH3N2_difference.csv")
coding.adjust.cal <- function(x) {
    chr <- unique(x$chr)
    start <- cal.bed$off.5[match(x$chr, cal.bed$chr)]
    
    mutate(x, coding.pos = pos - start)
}

qual.2004.5 <- ddply(qual.2004.5, ~chr, coding.adjust.cal)
qual.2005.6 <- ddply(qual.2005.6, ~chr, coding.adjust.cal)
```

Initially we partitioned the data relative to the geometric mean of the antibody titer for a given season. In the final anlysis we used a titer of 40 as the cutoff. These cut offs were identical except for in 2005-2006 in which the mean was 30. This affected 2 samples.

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

# qual.2007.8<-mutate(qual.2007.8,responder.HA=HAI.WI.30.post.vax>40,responder.NA=NAI.WI.30.post.vax>40,HAI.geo=HAI.WI.30.post.vax>gm_mean(HAI.WI.30.post.vax),NAI.geo=NAI.WI.30.post.vax>gm_mean(NAI.WI.30.post.vax),responder.both=(responder.NA&responder.HA==T),
# geom.both=(HAI.geo&NAI.geo==T))
minor.07.08 = subset(qual.2007.8, freq.var <= 0.5)

minor.04.05 = subset(qual.2004.5, freq.var <= 0.5)

minor.05.06 = subset(qual.2005.6, freq.var <= 0.5)



### For consensus seqeunce meta data

HAI.gm <- ddply(meta.all, ~season, summarise, gm_mean = gm_mean(HAI.post.vax))
NAI.gm <- ddply(meta.all, ~season, summarise, gm_mean = gm_mean(NAI.post.vax))

HAI.cut <- function(x) {
    season <- unique(x$season)
    gm <- HAI.gm$gm_mean[which(HAI.gm$season == season)]
    print(gm)
    mutate(x, HAI.geo = HAI.post.vax > gm, HAI.cut = HAI.post.vax >= 40)
}

NAI.cut <- function(x) {
    season <- unique(x$season)
    gm <- NAI.gm$gm_mean[which(NAI.gm$season == season)]
    print(gm)
    mutate(x, NAI.geo = NAI.post.vax > gm, NAI.cut = NAI.post.vax >= 40)
}
meta.df <- ddply(meta.all, ~season, HAI.cut)
```

```
## [1] 31.21756
## [1] 30.6433
## [1] 63.16057
```

```r
meta.df <- ddply(meta.df, ~season, NAI.cut)
```

```
## [1] 1
## [1] 1
## [1] 31.33869
```

```r
minor.04.05 <- ddply(minor.04.05, ~season, HAI.cut)
```

```
## [1] 31.21756
```

```r
minor.04.05 <- ddply(minor.04.05, ~season, NAI.cut)
```

```
## [1] 1
```

```r
minor.05.06 <- ddply(minor.05.06, ~season, HAI.cut)
```

```
## [1] 30.6433
```

```r
minor.05.06 <- ddply(minor.05.06, ~season, NAI.cut)
```

```
## [1] 1
```

```r
minor.07.08 <- ddply(minor.07.08, ~season, HAI.cut)
```

```
## [1] 63.16057
```

```r
minor.07.08 <- ddply(minor.07.08, ~season, NAI.cut)
```

```
## [1] 31.33869
```

```r
# meta.df<-mutate(meta.all,HAI.geo=HAI.WI.30.post.vax>gm_mean(HAI.WI.30.post.vax),NAI.geo=NAI.WI.30.post.vax>gm_mean(NAI.WI.30.post.vax))
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
ha.07.08 <- subset(minor.07.08, chr == "HA")
ha.04.05 <- subset(minor.04.05, chr == "HA")
ha.05.06 <- subset(minor.05.06, chr == "HA")



# dim(ha) ha<-ha[order(ha$pos),] unique(ha$pos)
write.csv(x = ha.07.08, file = "../results/2007-2008.HA.csv", row.names = F)
write.csv(x = ha.04.05, file = "../results/2004-2005.HA.csv", row.names = F)
write.csv(x = ha.05.06, file = "../results/2005-2006.HA.csv", row.names = F)


write.csv(x = minor.07.08, file = "../results/2007-2008.wg.csv", row.names = F)
write.csv(x = minor.04.05, file = "../results/2004-2005.wg.csv", row.names = F)
write.csv(x = minor.05.06, file = "../results/2005-2006.wg.csv", row.names = F)
write.csv(meta.df, "../data/raw/meta.all.HAgm.csv")
```


