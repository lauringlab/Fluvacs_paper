FluVacs Figures - Draft
================
Kari and JT
May 2, 2016

-   [Initial data processing](#initial-data-processing)
-   [Frequency of Variants](#frequency-of-variants)
-   [Diversity measures](#diversity-measures)
    -   [SNV count per sample](#snv-count-per-sample)
    -   [Shannon's Entropy](#shannons-entropy)
-   [Heat map](#heat-map)
-   [Antigenic Sites](#antigenic-sites)
    -   [All data](#all-data-2)
    -   [High quality](#high-quality-2)

Initial data processing
=======================

``` r
# read in the csvs and add a season column to use later
titer.2004.5 <- read.csv("../Titers_status_2004-2005.csv", stringsAsFactors = F)
titer.2004.5$season <- "04-05"
titer.2005.6 <- read.csv("../Titers_status_2005-2006.csv", stringsAsFactors = F)
titer.2005.6$season <- "05-06"
titer.2007.8 <- read.csv("../Titers_status_2007-2008.csv", stringsAsFactors = F)
titer.2007.8$season <- "07-08"
```

``` r
var.2007.8 <- read.csv("../data/processed/Run_1293/Variants/all.sum.csv", stringsAsFactors = F)
x <- read.csv("../data/processed/Run_1304/Variants/all.sum.csv", stringsAsFactors = F)  # the rest of these samples
var.2007.8 <- rbind(var.2007.8, x)  # combine both runs
other.seasons <- read.csv("../data/processed/Run_1412/Variants/all.sum.csv", 
    stringsAsFactors = F)


var.2004.5.df <- processing(data.df = other.seasons, meta.df = titer.2004.5, 
    pval = 0.01, phred = 35, mapq = 30, read_cut = c(32, 94))

var.2004.5.df <- subset(var.2004.5.df, season == "04-05")
var.2005.6.df <- processing(data.df = other.seasons, meta.df = titer.2005.6, 
    pval = 0.01, phred = 35, mapq = 30, read_cut = c(32, 94))
var.2005.6.df <- subset(var.2004.5.df, season == "05-06")

var.2007.8.df <- processing(data.df = var.2007.8, meta.df = titer.2007.8, pval = 0.01, 
    phred = 35, mapq = 30, read_cut = c(32, 94))

all.df <- rbind(var.2004.5.df, var.2005.6.df)
all.df <- rbind(all.df, var.2007.8.df)

all.df <- subset(all.df, freq.var > 0.005)
##### Now for the duplicate runs #######

var.2004.5.df2 <- read.csv("../data/processed/2004_2005/Variants/all.sum.csv", 
    stringsAsFactors = F)
var.2004.5.df2 <- processing(data.df = var.2004.5.df2, meta.df = titer.2004.5, 
    pval = 0.01, phred = 35, mapq = 30, read_cut = c(32, 94))

var.2005.6.df2 <- read.csv("../data/processed/2005-2006/Variants/all.sum.csv", 
    stringsAsFactors = F)
var.2005.6.df2 <- processing(data.df = var.2005.6.df2, meta.df = titer.2005.6, 
    pval = 0.01, phred = 35, mapq = 30, read_cut = c(32, 94))

var.2007.8.df2 <- read.csv("../data/processed/2007-2008/Variants/all.sum.csv", 
    stringsAsFactors = F)
var.2007.8.df2 <- processing(data.df = var.2007.8.df2, meta.df = titer.2007.8, 
    pval = 0.01, phred = 35, mapq = 30, read_cut = c(32, 94))

### Join duplicates #####

dups.2004.5 <- join_dups(data1.df = var.2004.5.df, data2.df = var.2004.5.df2)
dups.2005.6 <- join_dups(data1.df = var.2005.6.df, data2.df = var.2005.6.df2)
dups.2007.8 <- join_dups(data1.df = var.2007.8.df, data2.df = var.2007.8.df2)

##### Merge duplicates with intial data that was >1e5

qual.2004.5 <- high_qual(data1.df = var.2004.5.df, dups.df = dups.2004.5, titer = 1000)  # only duplicates above 1e3 kept
qual.2005.6 <- high_qual(data1.df = var.2005.6.df, dups.df = dups.2005.6, titer = 1000)
qual.2007.8 <- high_qual(data1.df = var.2007.8.df, dups.df = dups.2007.8, titer = 1000)
all.qual.df <- rbind(qual.2004.5, qual.2005.6)
all.qual.df <- rbind(all.qual.df, qual.2007.8)

all.qual.df <- subset(all.qual.df, freq.var > 0.005)
```

So after all that work we are left with all the samples (all seasons) in all.df and only those variants from high quality data in all.qual.df. This includes samples with &gt;10<sup>5</sup> genomes per microliter or &gt; 10<sup>3</sup> genomes per microliter but were sequenced in duplicate.

> The notes about NAs being produced is a little concerning. I'll look into it later. My hunch is some of the samples from that were sequenced were not included in the meta data

Frequency of Variants
=====================

``` r
ggplot(subset(all.df, freq.var < 0.9), aes(x = freq.var)) + geom_histogram(color = "white") + 
    scale_x_log10() + ggtitle("All samples")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

``` r
ggplot(subset(all.qual.df, freq.var < 0.9), aes(x = freq.var)) + geom_histogram(color = "white") + 
    scale_x_log10() + ggtitle("Only high quality samples")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-2-2.png" style="display: block; margin: auto;" />

Diversity measures
==================

SNV count per sample
--------------------

``` r
count_muts <- function(data.df) {
    ddply(data.df, ~Lauring_Id + Vax + season, function(x) dim(x)[1])
}

all.snv <- count_muts(all.df)
qual.snv <- count_muts(all.qual.df)
```

### Whole genome

#### all data

``` r
ggplot(data = subset(all.snv, !is.na(Vax)), aes(y = V1, x = Vax)) + geom_dotplot(stackdir = "center", 
    binaxis = "y") + ylab("Number of SNV")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

#### High quality

``` r
ggplot(data = subset(qual.snv, !is.na(Vax)), aes(y = V1, x = Vax)) + geom_dotplot(stackdir = "center", 
    binaxis = "y") + ylab("Number of SNV")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

### HA NA

#### all data

``` r
all.snv <- count_muts(subset(all.df, chr %in% c("HA", "N_A", "NR")))
ggplot(data = subset(all.snv, !is.na(Vax)), aes(y = V1, x = Vax)) + geom_dotplot(stackdir = "center", 
    binaxis = "y") + ylab("Number of SNV")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

#### High quality

``` r
qual.snv <- count_muts(subset(all.qual.df, chr %in% c("HA", "N_A", "NR")))

ggplot(data = subset(qual.snv, !is.na(Vax)), aes(y = V1, x = Vax)) + geom_dotplot(stackdir = "center", 
    binaxis = "y") + ylab("Number of SNV")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

### Anitigenic sites

``` r
Ha_antigenic <- function(df) {
    non_coding <- 29
    
    df <- mutate(df, coding_pos = pos - non_coding, AA_pos = (coding_pos - 1)%/%3 + 
        1)
    
    
    df$antigenic <- F
    antigenic_sites <- c(122, 124, 126, 130, 131, 132, 133, 135, 137, 138, 140, 
        142, 143, 144, 145, 146, 150, 152, 168, 128, 129, 155, 156, 157, 158, 
        159, 160, 163, 164, 165, 186, 187, 188, 189, 190, 192, 193, 194, 196, 
        197, 198, 44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 
        280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312, 96, 
        102, 103, 117, 121, 167, 170, 171, 172, 173, 174, 175, 176, 177, 179, 
        182, 201, 203, 207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 
        226, 227, 228, 229, 230, 238, 240, 242, 246, 247, 248, 57, 59, 62, 63, 
        67, 75, 78, 80, 81, 82, 83, 86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 
        265)
    # antigenic sites: Site A:
    # 122,124,126,130,131,132,133,135,137,138,140,142,143,144,145,146,150,152,168
    # Site B:
    # 128,129,155,156,157,158,159,160,163,164,165,186,187,188,189,190,192,193,194,196,197,198
    # Site C:
    # 44,45,46,47,48,50,51,53,54,273,275,276,278,279,280,294,297,299,300,304,305,307,308,309,310,311,312
    # Site D:
    # 96,102,103,117,121,167,170,171,172,173,174,175,176,177,179,182,201,203,207,208,209,212,213,214,215,216,217,218,219,226,227,228,229,230,238,240,242,246,247,248
    # Site E:
    # 57,59,62,63,67,75,78,80,81,82,83,86,87,88,91,92,94,109,260,261,262,265 All
    # sites:
    # 122,124,126,130,131,132,133,135,137,138,140,142,143,144,145,146,150,152,168,128,129,155,156,157,158,159,160,163,164,165,186,187,188,189,190,192,193,194,196,197,198,44,45,46,47,48,50,51,53,54,273,275,276,278,279,280,294,297,299,300,304,305,307,308,309,310,311,312,96,102,103,117,121,167,170,171,172,173,174,175,176,177,179,182,201,203,207,208,209,212,213,214,215,216,217,218,219,226,227,228,229,230,238,240,242,246,247,248,57,59,62,63,67,75,78,80,81,82,83,86,87,88,91,92,94,109,260,261,262,265
    # Cluster Transitions: 145,155,156,158,159,189,193
    
    
    
    df$antigenic[df$chr == "HA" & df$AA_pos %in% antigenic_sites] <- T
    return(df)
}


all.df <- Ha_antigenic(all.df)

all.df.counts <- count_muts(subset(all.df, antigenic == T))

ggplot(data = subset(all.df.counts, !is.na(Vax)), aes(y = V1, x = Vax)) + geom_dotplot(stackdir = "center", 
    binaxis = "y", dotsize = 0.8) + ylab("Number of SNV") + ggtitle("all samples")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

``` r
all.qual.df <- Ha_antigenic(all.qual.df)

all.qual.df.counts <- count_muts(subset(all.qual.df, antigenic == T))
ggplot(data = subset(all.qual.df.counts, !is.na(Vax)), aes(y = V1, x = Vax)) + 
    geom_dotplot(stackdir = "center", binaxis = "y", dotsize = 0.8) + ylab("Number of SNV") + 
    ggtitle("quality samples")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-7-2.png" style="display: block; margin: auto;" />

``` r
all.df.counts <- count_muts(subset(all.df, antigenic == T & freq.var < 0.9))

ggplot(data = subset(all.df.counts), aes(y = V1, x = Vax)) + geom_dotplot(stackdir = "center", 
    binaxis = "y") + ylab("Number of SNV") + ggtitle("all samples frequency<90%")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-7-3.png" style="display: block; margin: auto;" />

``` r
all.qual.df.counts <- count_muts(subset(all.qual.df, antigenic == T & freq.var < 
    0.9))
ggplot(data = subset(all.qual.df.counts), aes(y = V1, x = Vax)) + geom_dotplot(stackdir = "center", 
    binaxis = "y") + ylab("Number of SNV") + ggtitle("quality samples (frequency<90%)")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-7-4.png" style="display: block; margin: auto;" />

Shannon's Entropy
-----------------

For the entropy calculations I am assuming 39000 potential variants. This is not right but probably close enough. It is just used to normalize the data. In the future we can use the fasta reference files to get the exact number.

### Whole genome

``` r
possible_vars <- 3900  # fix to be accurate

H <- function(x) {
    x <- subset(x, freq.var < 0.9)
    H_pos <- ddply(x, ~chr + pos, summarize, wt = 1 - sum(freq.var), H = -sum(c(freq.var * 
        log(freq.var), wt * log(wt))))
    sum(H_pos$H)/(possible_vars/3)
}

all.H <- ddply(all.df, ~Lauring_Id + Vax + season, H)

ggplot(data = subset(all.H), aes(y = V1, x = Vax)) + geom_dotplot(stackdir = "center", 
    binaxis = "y") + ylab("Shannon's Entropy") + scale_y_log10() + ggtitle("all samples")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

``` r
qual.H <- ddply(all.qual.df, ~Lauring_Id + Vax + season, H)
ggplot(data = subset(qual.H), aes(y = V1, x = Vax)) + geom_dotplot(stackdir = "center", 
    binaxis = "y") + ylab("Shannon's Entropy") + scale_y_log10() + ggtitle("quality samples")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-8-2.png" style="display: block; margin: auto;" />

### HA NA

``` r
all.HA.NA.H <- ddply(subset(all.df, chr %in% c("HA", "N_A", "NR")), ~Lauring_Id + 
    Vax + season, H)

ggplot(data = subset(all.HA.NA.H), aes(y = V1, x = Vax)) + geom_dotplot(stackdir = "center", 
    binaxis = "y") + ylab("Shannon's Entropy") + scale_y_log10() + ggtitle("all samples")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

``` r
qual.HA.NA.H <- ddply(subset(all.qual.df, chr %in% c("HA", "N_A", "NR")), ~Lauring_Id + 
    Vax + season, H)
ggplot(data = subset(qual.H), aes(y = V1, x = Vax)) + geom_dotplot(stackdir = "center", 
    binaxis = "y") + ylab("Shannon's Entropy") + scale_y_log10() + ggtitle("quality samples")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-9-2.png" style="display: block; margin: auto;" />

### Anitigenic sites

``` r
all.ant.H <- ddply(subset(all.df, antigenic == T), ~Lauring_Id + Vax + season, 
    H)

ggplot(data = subset(all.ant.H), aes(y = V1, x = Vax)) + geom_dotplot(stackdir = "center", 
    binaxis = "y") + ylab("Shannon's Entropy") + scale_y_log10() + ggtitle("all samples")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

``` r
qual.ant.H <- ddply(subset(all.qual.df, antigenic == T), ~Lauring_Id + Vax + 
    season, H)

ggplot(data = subset(qual.ant.H), aes(y = V1, x = Vax)) + geom_dotplot(stackdir = "center", 
    binaxis = "y") + ylab("Shannon's Entropy") + scale_y_log10() + ggtitle("quality samples")
```

    ## `stat_bindot()` using `bins = 30`. Pick better value with `binwidth`.

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-10-2.png" style="display: block; margin: auto;" />

Heat map
========

Theses are made with the quality samples

``` r
make_heat_map(subset(all.qual.df, season == "07-08"))
```

    ## Aggregation function missing: defaulting to length

    ##    mutation      variable value segment
    ## 1 PB2_A182G 07-08.79.TRUE    -1     PB2
    ## 2 PB2_G186A 07-08.79.TRUE    -1     PB2
    ## 3 PB2_T201C 07-08.79.TRUE    -1     PB2
    ## 4 PB2_G204A 07-08.79.TRUE    -1     PB2
    ## 5 PB2_G209A 07-08.79.TRUE    -1     PB2
    ## 6 PB2_G220A 07-08.79.TRUE    -1     PB2

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

``` r
# make_heat_map(subset(all.qual.df,season=='05-06'))
make_heat_map(subset(all.qual.df, season == "04-05"))
```

    ##    mutation       variable     value segment
    ## 1  PB2_G96A 04-05.402.TRUE -1.000000     PB2
    ## 2 PB2_T249C 04-05.402.TRUE -1.000000     PB2
    ## 3 PB2_A316G 04-05.402.TRUE  1.998882     PB2
    ## 4 PB2_C462T 04-05.402.TRUE -1.000000     PB2
    ## 5 PB2_C465T 04-05.402.TRUE  1.999526     PB2
    ## 6 PB2_A486G 04-05.402.TRUE -1.000000     PB2

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-11-2.png" style="display: block; margin: auto;" />

Antigenic Sites
===============

### All data

``` r
kable(subset(all.df, antigenic == T, select = c(Lauring_Id, mutation, freq.var)))
```

|       |  Lauring\_Id| mutation  |   freq.var|
|-------|------------:|:----------|----------:|
| 3     |          412| HA\_T460G |  0.9978714|
| 4     |          412| HA\_T516C |  0.9986327|
| 5     |          412| HA\_A610G |  0.9979880|
| 7     |          412| HA\_T755C |  0.0051081|
| 9     |          412| HA\_C863T |  0.0473373|
| 96    |          407| HA\_T213C |  0.0197986|
| 97    |          407| HA\_T460G |  0.9965571|
| 98    |          407| HA\_T516C |  0.9979601|
| 99    |          407| HA\_G544A |  0.9973478|
| 101   |          407| HA\_A610G |  0.9972039|
| 105   |          407| HA\_G715A |  0.9976775|
| 748   |          426| HA\_T460G |  0.9988866|
| 749   |          426| HA\_G497A |  0.9983497|
| 750   |          426| HA\_T516C |  0.9995262|
| 751   |          426| HA\_A610G |  0.9988186|
| 753   |          426| HA\_A717G |  0.0344117|
| 988   |          404| HA\_T460G |  0.9971313|
| 989   |          404| HA\_G497A |  0.9976937|
| 990   |          404| HA\_T516C |  0.9975737|
| 993   |          404| HA\_A610G |  0.9938350|
| 1092  |          411| HA\_T460G |  0.9970203|
| 1093  |          411| HA\_G497A |  0.9974661|
| 1094  |          411| HA\_T504C |  0.0198027|
| 1095  |          411| HA\_T516C |  0.9982140|
| 1096  |          411| HA\_A610G |  0.9972855|
| 1156  |          409| HA\_C483T |  0.9939058|
| 1157  |          409| HA\_T516C |  0.9986246|
| 1158  |          409| HA\_C597T |  0.9941361|
| 1159  |          409| HA\_A610G |  0.9983770|
| 1649  |          401| HA\_T460G |  0.9968130|
| 1650  |          401| HA\_T516C |  0.9976624|
| 1651  |          401| HA\_A610G |  0.9975520|
| 2057  |          402| HA\_T173C |  0.0082386|
| 2058  |          402| HA\_C178G |  0.0138956|
| 2062  |          402| HA\_T460G |  0.9983584|
| 2064  |          402| HA\_A610G |  0.9979201|
| 2069  |          402| HA\_G765A |  0.9972809|
| 2165  |          403| HA\_C483T |  0.9977536|
| 2166  |          403| HA\_T516C |  0.9986142|
| 2167  |          403| HA\_C597T |  0.9960192|
| 2168  |          403| HA\_A610G |  0.9977149|
| 2473  |          414| HA\_G449A |  0.0087982|
| 2474  |          414| HA\_C483T |  0.9987189|
| 2475  |          414| HA\_T508C |  0.0062232|
| 2476  |          414| HA\_T516C |  0.9987315|
| 2477  |          414| HA\_T593C |  0.0059936|
| 2478  |          414| HA\_C597T |  0.9972632|
| 2479  |          414| HA\_A610G |  0.9981937|
| 2482  |          414| HA\_T755C |  0.0064331|
| 3254  |          415| HA\_A288G |  0.0058375|
| 3255  |          415| HA\_T460G |  0.9986629|
| 3256  |          415| HA\_G497A |  0.9988006|
| 3258  |          415| HA\_T516C |  0.9989183|
| 3259  |          415| HA\_A610G |  0.9983850|
| 3262  |          415| HA\_G814A |  0.0057951|
| 3852  |          406| HA\_A199C |  0.0135042|
| 3854  |          406| HA\_T460G |  0.9980997|
| 3855  |          406| HA\_G497A |  0.9986187|
| 3856  |          406| HA\_T516C |  0.9987882|
| 3857  |          406| HA\_A610G |  0.9978076|
| 4240  |          400| HA\_G180A |  0.9952773|
| 4242  |          400| HA\_T460G |  0.9978954|
| 4243  |          400| HA\_T516C |  0.9982730|
| 4244  |          400| HA\_A610G |  0.9971264|
| 4324  |          427| HA\_T460G |  0.9971152|
| 4325  |          427| HA\_T516C |  0.9996469|
| 4327  |          427| HA\_A610G |  0.9971178|
| 4330  |          427| HA\_A630G |  0.0275211|
| 4333  |          427| HA\_G714A |  0.0353287|
| 4334  |          427| HA\_G715A |  0.9550557|
| 4338  |          427| HA\_C811A |  0.0323465|
| 1413  |           20| HA\_G335A |  0.9977312|
| 1414  |           20| HA\_C440T |  0.0071559|
| 1415  |           20| HA\_A565C |  0.9985224|
| 1508  |           58| HA\_G335A |  0.9963896|
| 1510  |           58| HA\_A565C |  0.9976782|
| 1628  |           15| HA\_G335A |  0.9971824|
| 1629  |           15| HA\_A565C |  0.9983310|
| 16951 |           28| HA\_G335A |  0.9970985|
| 16971 |           28| HA\_A565C |  0.9981969|
| 1767  |           71| HA\_G335A |  0.9954906|
| 1768  |           71| HA\_A565C |  0.9969662|
| 1833  |           10| HA\_G335A |  0.9961708|
| 1835  |           10| HA\_A565C |  0.9977655|
| 1837  |           10| HA\_A652G |  0.0088028|
| 1950  |           55| HA\_G335A |  0.9973494|
| 1951  |           55| HA\_C411T |  0.0680573|
| 1953  |           55| HA\_A565C |  0.9984146|
| 2045  |           40| HA\_G335A |  0.9973318|
| 2046  |           40| HA\_G434A |  0.0231497|
| 2047  |           40| HA\_A565C |  0.9985938|
| 2053  |           40| HA\_G854A |  0.9979359|
| 20581 |           40| HA\_G958A |  0.0455798|
| 21561 |           67| HA\_G335A |  0.9956266|
| 21571 |           67| HA\_A565C |  0.9972097|
| 21581 |           67| HA\_G610A |  0.0080955|
| 22461 |           64| HA\_G335A |  0.9961368|
| 22471 |           64| HA\_A565C |  0.9982134|
| 2317  |           74| HA\_G335A |  0.9939727|
| 2318  |           74| HA\_A565C |  0.9978975|
| 2394  |           27| HA\_G335A |  0.9972820|
| 2395  |           27| HA\_A565C |  0.9978890|
| 2398  |           27| HA\_G854A |  0.9986577|
| 2454  |           31| HA\_G335A |  0.9982555|
| 2455  |           31| HA\_A565C |  0.9988404|
| 25131 |           66| HA\_G335A |  0.9972223|
| 25141 |           66| HA\_A565C |  0.9979457|
| 25881 |           17| HA\_G335A |  0.9926969|
| 25891 |           17| HA\_A565C |  0.9920638|
| 2650  |           52| HA\_G335A |  0.9974599|
| 2651  |           52| HA\_A565C |  0.9983659|
| 2653  |           52| HA\_G685A |  0.0087493|
| 2730  |           49| HA\_G335A |  0.9973413|
| 2731  |           49| HA\_A565C |  0.9981817|
| 2796  |           59| HA\_G335A |  0.9969594|
| 2797  |           59| HA\_A565C |  0.9977509|
| 2911  |           34| HA\_G335A |  0.9964811|
| 2912  |           34| HA\_A565C |  0.9980709|
| 2981  |           70| HA\_G335A |  0.9968609|
| 2982  |           70| HA\_T518C |  0.9977758|
| 2983  |           70| HA\_A565G |  0.9981506|
| 2985  |           70| HA\_A717G |  0.0247364|
| 2987  |           70| HA\_C865T |  0.0162663|
| 3098  |           61| HA\_G335A |  0.9964441|
| 3099  |           61| HA\_A565C |  0.9978608|
| 3172  |           41| HA\_G335A |  0.9978889|
| 3173  |           41| HA\_A565C |  0.9984537|
| 3176  |           41| HA\_G854A |  0.9985950|
| 3241  |           30| HA\_G335A |  0.9973585|
| 3242  |           30| HA\_T462C |  0.0075009|
| 3243  |           30| HA\_A565C |  0.9983373|
| 33331 |            3| HA\_G335A |  0.9969682|
| 33341 |            3| HA\_A565C |  0.9982774|
| 33361 |            3| HA\_G671A |  0.9839711|
| 3408  |           73| HA\_G335A |  0.9964744|
| 3409  |           73| HA\_A565C |  0.9970282|
| 3479  |           57| HA\_G335A |  0.9974149|
| 3480  |           57| HA\_A565C |  0.9983615|
| 3547  |            6| HA\_C164A |  0.0399776|
| 3548  |            6| HA\_T228C |  0.0362521|
| 3549  |            6| HA\_C285A |  0.0431470|
| 3550  |            6| HA\_G335A |  0.9973246|
| 3551  |            6| HA\_A565C |  0.9981962|
| 3552  |            6| HA\_T594C |  0.0182470|
| 3553  |            6| HA\_A611G |  0.0058156|
| 3555  |            6| HA\_A651G |  0.0075603|
| 3651  |           56| HA\_G335A |  0.9967348|
| 3652  |           56| HA\_A565C |  0.9981518|
| 3748  |           72| HA\_A199G |  0.0116255|
| 3749  |           72| HA\_G335A |  0.9968578|
| 3750  |           72| HA\_C420T |  0.0135588|
| 3751  |           72| HA\_C548T |  0.0080640|
| 3752  |           72| HA\_A565C |  0.9981816|
| 38611 |           76| HA\_G335A |  0.9972253|
| 38621 |           76| HA\_G428A |  0.0184960|
| 38631 |           76| HA\_A565C |  0.9980131|
| 38641 |           76| HA\_T606G |  0.0073840|
| 4009  |            8| HA\_G335A |  0.9972604|
| 4010  |            8| HA\_A565C |  0.9979842|
| 4078  |           19| HA\_G335A |  0.9977268|
| 4079  |           19| HA\_A565C |  0.9986290|
| 4142  |           75| HA\_G335A |  0.9971999|
| 4145  |           75| HA\_A565C |  0.9978136|
| 4230  |           68| HA\_T311C |  0.0110569|
| 4231  |           68| HA\_G335A |  0.9969388|
| 4232  |           68| HA\_A565C |  0.9983324|
| 4233  |           68| HA\_G592A |  0.0192200|
| 43401 |           26| HA\_G335A |  0.9978690|
| 43411 |           26| HA\_A565C |  0.9983244|
| 43431 |           26| HA\_A672G |  0.1156970|
| 44131 |            4| HA\_G335A |  0.9958891|
| 44141 |            4| HA\_A565C |  0.9978261|
| 4481  |           16| HA\_G335A |  0.9946479|
| 4482  |           16| HA\_A565C |  0.9960869|
| 4545  |           18| HA\_G335A |  0.9884812|
| 4546  |           18| HA\_A565C |  0.9956490|
| 4608  |           78| HA\_G335A |  0.9967651|
| 4609  |           78| HA\_A565C |  0.9979265|
| 4671  |            7| HA\_G335A |  0.9968246|
| 4672  |            7| HA\_A565C |  0.9980902|
| 4734  |            5| HA\_G335A |  0.9973509|
| 4736  |            5| HA\_A565C |  0.9977732|
| 4830  |           51| HA\_G335A |  0.9963630|
| 4831  |           51| HA\_A565C |  0.9981302|
| 4834  |           51| HA\_C865T |  0.9970511|
| 4904  |           NA| HA\_G335A |  0.9977533|
| 4905  |           NA| HA\_A412G |  0.9980955|
| 4906  |           NA| HA\_C447T |  0.9990917|
| 4907  |           NA| HA\_G460T |  0.9987126|
| 4908  |           NA| HA\_T467A |  0.9983806|
| 4909  |           NA| HA\_G605T |  0.9985123|
| 4911  |           NA| HA\_G715A |  0.9986315|
| 4915  |           NA| HA\_C954T |  0.9985809|
| 5045  |           24| HA\_G335A |  0.9976960|
| 5046  |           24| HA\_A565C |  0.9987801|
| 5121  |           25| HA\_G335A |  0.9981932|
| 5122  |           25| HA\_G497A |  0.9993614|
| 5123  |           25| HA\_A565C |  0.9986573|
| 5184  |           46| HA\_G335A |  0.9975693|
| 5185  |           46| HA\_A565C |  0.9984604|
| 5245  |           33| HA\_G335A |  0.9979487|
| 5246  |           33| HA\_A565C |  0.9983568|
| 5317  |           50| HA\_A315G |  0.9984685|
| 5318  |           50| HA\_G335A |  0.9974576|
| 5319  |           50| HA\_A565C |  0.9984302|
| 5380  |            2| HA\_G335A |  0.9962115|
| 5381  |            2| HA\_A565C |  0.9980567|
| 5443  |           39| HA\_G335A |  0.9964492|
| 5444  |           39| HA\_A565C |  0.9974803|
| 5509  |           29| HA\_G335A |  0.9974372|
| 5510  |           29| HA\_A565C |  0.9984071|
| 5575  |           54| HA\_G335A |  0.9964840|
| 5576  |           54| HA\_A399G |  0.0082548|
| 5577  |           54| HA\_A565C |  0.9981359|
| 5579  |           54| HA\_G671A |  0.9815564|
| 5641  |           23| HA\_G335A |  0.9977300|
| 5642  |           23| HA\_G419A |  0.0082951|
| 5643  |           23| HA\_A565C |  0.9986095|
| 5648  |           23| HA\_C959T |  0.0392059|
| 5730  |           77| HA\_G335A |  0.9954350|
| 5731  |           77| HA\_A399G |  0.0101209|
| 5732  |           77| HA\_A565C |  0.9946714|
| 5846  |           22| HA\_G335A |  0.9973459|
| 5847  |           22| HA\_A565C |  0.9986670|
| 5912  |           43| HA\_G335A |  0.9973633|
| 5913  |           43| HA\_A565C |  0.9974667|
| 5915  |           43| HA\_G674A |  0.0054497|
| 5980  |           79| HA\_G335A |  0.9958628|
| 5981  |           79| HA\_T518C |  0.9978475|
| 5984  |           79| HA\_C709A |  0.0073681|
| 5987  |           79| HA\_C950A |  0.0098849|
| 6106  |           42| HA\_G335A |  0.9972996|
| 6107  |           42| HA\_T459C |  0.0100849|
| 6108  |           42| HA\_A565C |  0.9985562|
| 6172  |           47| HA\_G335A |  0.9979190|
| 6173  |           47| HA\_A565C |  0.9985586|
| 6239  |           80| HA\_G335A |  0.9978706|
| 6240  |           80| HA\_G464A |  0.0051679|
| 6241  |           80| HA\_G485A |  0.0082520|
| 6242  |           80| HA\_T504C |  0.0063246|
| 6243  |           80| HA\_A565C |  0.9988508|
| 6352  |           38| HA\_T276C |  0.0361201|
| 6353  |           38| HA\_C333A |  0.0286524|
| 6354  |           38| HA\_G335A |  0.9989949|
| 6355  |           38| HA\_G494A |  0.0290098|
| 6356  |           38| HA\_A565C |  0.9994678|
| 6429  |           14| HA\_G335A |  0.9979290|
| 6431  |           14| HA\_A565C |  0.9987825|
| 6434  |           14| HA\_C709T |  0.0073778|
| 6436  |           14| HA\_C809T |  0.0067067|
| 6437  |           14| HA\_G854A |  0.1048269|
| 6558  |           93| HA\_G335A |  0.9983432|
| 6560  |           93| HA\_A565C |  0.9991695|
| 6564  |           93| HA\_C954T |  0.9970880|
| 6633  |           91| HA\_G335A |  0.9991905|
| 6634  |           91| HA\_A565C |  0.9996051|
| 6692  |           12| HA\_G335A |  0.9979693|
| 6694  |           12| HA\_A565C |  0.9988072|
| 6695  |           12| HA\_G603A |  0.9925620|
| 6757  |          103| HA\_G335A |  0.9984687|
| 6758  |          103| HA\_A565C |  0.9986705|
| 6760  |          103| HA\_G715A |  0.0058576|
| 6832  |           36| HA\_G335A |  0.9983351|
| 6834  |           36| HA\_A565C |  0.9990739|
| 6835  |           36| HA\_T585A |  0.0222295|
| 6839  |           36| HA\_T918C |  0.0063328|
| 6934  |           85| HA\_G205A |  0.0149806|
| 6935  |           85| HA\_G335A |  0.9983505|
| 6936  |           85| HA\_A565C |  0.9989570|
| 7027  |           89| HA\_G335A |  0.9974933|
| 7029  |           89| HA\_A565C |  0.9986828|
| 7103  |          106| HA\_G335A |  0.9984824|
| 7104  |          106| HA\_A565C |  0.9988101|
| 7169  |           83| HA\_G335A |  0.9989679|
| 7170  |           83| HA\_A565C |  0.9992604|
| 7172  |           83| HA\_T663C |  0.9981590|
| 7175  |           83| HA\_C942T |  0.9986578|
| 7252  |          104| HA\_G335A |  0.9977110|
| 7254  |          104| HA\_A565C |  0.9987507|
| 7366  |           NA| HA\_G335A |  0.9984826|
| 7367  |           NA| HA\_A412G |  0.9981518|
| 7368  |           NA| HA\_C447T |  0.9992102|
| 7369  |           NA| HA\_G460T |  0.9990123|
| 7370  |           NA| HA\_T467A |  0.9990066|
| 7371  |           NA| HA\_G605T |  0.9989709|
| 7373  |           NA| HA\_G715A |  0.9988844|
| 7377  |           NA| HA\_C954T |  0.9987671|
| 7511  |           98| HA\_G335A |  0.9981452|
| 7514  |           98| HA\_A565C |  0.9989137|
| 7594  |           86| HA\_G335A |  0.9977122|
| 7595  |           86| HA\_A565C |  0.9982039|
| 7672  |           82| HA\_G335A |  0.9986986|
| 7673  |           82| HA\_A565C |  0.9987050|
| 7677  |           82| HA\_G854A |  0.0083207|
| 7679  |           82| HA\_A910G |  0.0114305|
| 7749  |           88| HA\_C214T |  0.9978096|
| 7750  |           88| HA\_G262A |  0.0062458|
| 7751  |           88| HA\_G335A |  0.9983515|
| 7752  |           88| HA\_A565C |  0.9989711|
| 7756  |           88| HA\_A754G |  0.0066067|
| 7757  |           88| HA\_C865T |  0.0137764|
| 7840  |          105| HA\_G335A |  0.9982814|
| 7841  |          105| HA\_T509C |  0.0099163|
| 7842  |          105| HA\_A565C |  0.9989095|
| 7845  |          105| HA\_T743C |  0.0071381|
| 7948  |          100| HA\_G335A |  0.9983947|
| 7949  |          100| HA\_T518C |  0.9957963|
| 8029  |           11| HA\_G335A |  0.9977941|
| 8030  |           11| HA\_A565C |  0.9987086|
| 8031  |           11| HA\_G603A |  0.9980866|
| 8093  |           84| HA\_G335A |  0.9987552|
| 8094  |           84| HA\_A565C |  0.9994419|
| 8176  |           37| HA\_G335A |  0.9986999|
| 8177  |           37| HA\_A565C |  0.9989790|
| 8268  |           44| HA\_G335A |  0.9985839|
| 8269  |           44| HA\_A565C |  0.9988419|
| 8272  |           44| HA\_G671A |  0.9626091|
| 8347  |           96| HA\_G335A |  0.9989645|
| 8348  |           96| HA\_A565C |  0.9989827|
| 8432  |           81| HA\_G335A |  0.9980835|
| 8433  |           81| HA\_A565C |  0.9985278|
| 8504  |           13| HA\_G335A |  0.9977982|
| 8506  |           13| HA\_A565C |  0.9987868|
| 8509  |           13| HA\_G854A |  0.4659828|
| 8599  |           92| HA\_G335A |  0.9983227|
| 8600  |           92| HA\_A565C |  0.9990205|
| 8681  |          101| HA\_G335A |  0.9976464|
| 8682  |          101| HA\_T518C |  0.9985878|
| 8767  |           90| HA\_G335A |  0.9963350|
| 8768  |           90| HA\_A565C |  0.9991500|

### High quality

``` r
kable(subset(all.qual.df, antigenic == T, select = c(Lauring_Id, mutation, freq.var)))
```

|       |  Lauring\_Id| mutation  |   freq.var|
|-------|------------:|:----------|----------:|
| 988   |          404| HA\_T460G |  0.9971313|
| 989   |          404| HA\_G497A |  0.9976937|
| 990   |          404| HA\_T516C |  0.9975737|
| 993   |          404| HA\_A610G |  0.9938350|
| 1649  |          401| HA\_T460G |  0.9968130|
| 1650  |          401| HA\_T516C |  0.9976624|
| 1651  |          401| HA\_A610G |  0.9975520|
| 192   |          400| HA\_A610G |  0.9977904|
| 195   |          400| HA\_G180A |  0.9961603|
| 197   |          400| HA\_T460G |  0.9985080|
| 198   |          400| HA\_T516C |  0.9986825|
| 257   |          402| HA\_A610G |  0.9982760|
| 259   |          402| HA\_G765A |  0.9977265|
| 262   |          402| HA\_T460G |  0.9987133|
| 283   |          403| HA\_A610G |  0.9982779|
| 288   |          403| HA\_C483T |  0.9983432|
| 289   |          403| HA\_C597T |  0.9967209|
| 293   |          403| HA\_T516C |  0.9988845|
| 363   |          406| HA\_A610G |  0.9982362|
| 366   |          406| HA\_G497A |  0.9984981|
| 368   |          406| HA\_T460G |  0.9985823|
| 369   |          406| HA\_T516C |  0.9991032|
| 404   |          409| HA\_A610G |  0.9986096|
| 409   |          409| HA\_C483T |  0.9967033|
| 410   |          409| HA\_C597T |  0.9962182|
| 414   |          409| HA\_T516C |  0.9989004|
| 483   |          414| HA\_A610G |  0.9984139|
| 487   |          414| HA\_C483T |  0.9988346|
| 488   |          414| HA\_C597T |  0.9965694|
| 492   |          414| HA\_T516C |  0.9989146|
| 563   |          426| HA\_A610G |  0.9942401|
| 566   |          426| HA\_G497A |  0.7869945|
| 567   |          426| HA\_T460G |  0.9940300|
| 568   |          426| HA\_T516C |  0.9796461|
| 605   |          427| HA\_A610G |  0.9982441|
| 607   |          427| HA\_A630G |  0.5127296|
| 610   |          427| HA\_C811A |  0.5151538|
| 613   |          427| HA\_G714A |  0.5166147|
| 615   |          427| HA\_T460G |  0.9985207|
| 616   |          427| HA\_T516C |  0.9979014|
| 1499  |           34| HA\_G335A |  0.9964811|
| 1500  |           34| HA\_A565C |  0.9980709|
| 72    |           10| HA\_A565C |  0.9983945|
| 83    |           10| HA\_G335A |  0.9967105|
| 148   |          100| HA\_G335A |  0.9978070|
| 150   |          100| HA\_T518C |  0.9957175|
| 22410 |          101| HA\_G335A |  0.9975342|
| 22610 |          101| HA\_T518C |  0.9981249|
| 29010 |          103| HA\_A565C |  0.9988935|
| 29810 |          103| HA\_G335A |  0.9979138|
| 29910 |          103| HA\_G715A |  0.0061154|
| 35910 |          104| HA\_A565C |  0.9987385|
| 36610 |          104| HA\_G335A |  0.9975562|
| 4301  |          105| HA\_A565C |  0.9989337|
| 4381  |          105| HA\_G335A |  0.9978429|
| 4931  |          106| HA\_A565C |  0.9986871|
| 5001  |          106| HA\_G335A |  0.9979818|
| 5581  |           11| HA\_A565C |  0.9988134|
| 5651  |           11| HA\_G335A |  0.9976365|
| 5661  |           11| HA\_G603A |  0.9979748|
| 6191  |           12| HA\_A565C |  0.9989583|
| 6281  |           12| HA\_G335A |  0.9978017|
| 6291  |           12| HA\_G603A |  0.9929640|
| 6821  |           14| HA\_A565C |  0.9987207|
| 6891  |           14| HA\_G335A |  0.9977344|
| 6901  |           14| HA\_G854A |  0.0788728|
| 745   |           15| HA\_A565C |  0.9986656|
| 752   |           15| HA\_G335A |  0.9973627|
| 806   |           16| HA\_A565C |  0.9974008|
| 815   |           16| HA\_G335A |  0.9959440|
| 872   |           17| HA\_A565C |  0.9954355|
| 879   |           17| HA\_G335A |  0.9947886|
| 930   |           18| HA\_A565C |  0.9971292|
| 937   |           18| HA\_G335A |  0.9928576|
| 9911  |           19| HA\_A565C |  0.9987964|
| 9981  |           19| HA\_G335A |  0.9974281|
| 10521 |            2| HA\_A565C |  0.9985182|
| 10601 |            2| HA\_G335A |  0.9967804|
| 1113  |           22| HA\_A565C |  0.9988689|
| 1122  |           22| HA\_G335A |  0.9973136|
| 1179  |           25| HA\_A565C |  0.9988729|
| 1186  |           25| HA\_G335A |  0.9977840|
| 1187  |           25| HA\_G497A |  0.9976218|
| 1242  |           26| HA\_A565C |  0.9986016|
| 1243  |           26| HA\_A672G |  0.1305439|
| 1250  |           26| HA\_G335A |  0.9966141|
| 1306  |           27| HA\_A565C |  0.9984427|
| 1314  |           27| HA\_G335A |  0.9969219|
| 1315  |           27| HA\_G854A |  0.9979336|
| 1367  |           28| HA\_A565C |  0.9986050|
| 1374  |           28| HA\_G335A |  0.9974500|
| 1428  |           29| HA\_A565C |  0.9986560|
| 1435  |           29| HA\_G335A |  0.9971020|
| 1494  |            3| HA\_A565C |  0.9986414|
| 15001 |            3| HA\_G335A |  0.9969509|
| 15011 |            3| HA\_G671A |  0.9851019|
| 15631 |           31| HA\_A565C |  0.9989570|
| 1570  |           31| HA\_G335A |  0.9979034|
| 1621  |           33| HA\_A565C |  0.9985810|
| 1629  |           33| HA\_G335A |  0.9976974|
| 16831 |           34| HA\_A565C |  0.9984352|
| 16901 |           34| HA\_G335A |  0.9968830|
| 1751  |           39| HA\_A565C |  0.9982252|
| 1758  |           39| HA\_G335A |  0.9969295|
| 1817  |            4| HA\_A565C |  0.9982301|
| 1824  |            4| HA\_G335A |  0.9962522|
| 1881  |           41| HA\_A565C |  0.9987356|
| 1890  |           41| HA\_G335A |  0.9976531|
| 1891  |           41| HA\_G854A |  0.9985191|
| 1949  |           42| HA\_A565C |  0.9984543|
| 1957  |           42| HA\_G335A |  0.9967091|
| 1960  |           42| HA\_T459C |  0.0097621|
| 2014  |           43| HA\_A565C |  0.9980756|
| 2021  |           43| HA\_G335A |  0.9972581|
| 2079  |           46| HA\_A565C |  0.9986484|
| 2087  |           46| HA\_G335A |  0.9975462|
| 2139  |           47| HA\_A565C |  0.9988333|
| 2148  |           47| HA\_G335A |  0.9977795|
| 2203  |           49| HA\_A565C |  0.9983454|
| 2211  |           49| HA\_G335A |  0.9973142|
| 2269  |            5| HA\_A565C |  0.9982648|
| 2277  |            5| HA\_G335A |  0.9970790|
| 2330  |           50| HA\_A315G |  0.9983631|
| 2331  |           50| HA\_A565C |  0.9984787|
| 2338  |           50| HA\_G335A |  0.9971437|
| 2392  |           51| HA\_A565C |  0.9982992|
| 2397  |           51| HA\_C865T |  0.9969113|
| 2400  |           51| HA\_G335A |  0.9966877|
| 2463  |           54| HA\_A399G |  0.0080412|
| 2464  |           54| HA\_A565C |  0.9985944|
| 2470  |           54| HA\_G335A |  0.9967786|
| 2471  |           54| HA\_G671A |  0.9839943|
| 2529  |           55| HA\_A565C |  0.9986810|
| 2533  |           55| HA\_C411T |  0.0658146|
| 2540  |           55| HA\_G335A |  0.9972418|
| 2595  |           56| HA\_A565C |  0.9985008|
| 2602  |           56| HA\_G335A |  0.9972228|
| 2656  |           57| HA\_A565C |  0.9985654|
| 2663  |           57| HA\_G335A |  0.9973826|
| 2717  |            6| HA\_A565C |  0.9985249|
| 2721  |            6| HA\_C285A |  0.0395353|
| 2726  |            6| HA\_G335A |  0.9973278|
| 2783  |           61| HA\_A565C |  0.9984241|
| 2792  |           61| HA\_G335A |  0.9967259|
| 2851  |           64| HA\_A565C |  0.9986240|
| 2861  |           64| HA\_G335A |  0.9967796|
| 2921  |           66| HA\_A565C |  0.9983496|
| 2929  |           66| HA\_G335A |  0.9973845|
| 2988  |           67| HA\_A565C |  0.9980889|
| 2997  |           67| HA\_G335A |  0.9964781|
| 3053  |            7| HA\_A565C |  0.9984670|
| 3061  |            7| HA\_G335A |  0.9970704|
| 3115  |           71| HA\_A565C |  0.9981030|
| 3122  |           71| HA\_G335A |  0.9965689|
| 3179  |           72| HA\_A565C |  0.9982216|
| 3187  |           72| HA\_G335A |  0.9971008|
| 3243  |           73| HA\_A565C |  0.9981769|
| 3250  |           73| HA\_G335A |  0.9967526|
| 3307  |           74| HA\_A565C |  0.9982042|
| 3314  |           74| HA\_G335A |  0.9941117|
| 3374  |           76| HA\_A565C |  0.9985004|
| 3383  |           76| HA\_G335A |  0.9975395|
| 3443  |           78| HA\_A565C |  0.9984367|
| 3450  |           78| HA\_G335A |  0.9972318|
| 3512  |           79| HA\_G335A |  0.9966175|
| 3514  |           79| HA\_T518C |  0.9978933|
| 3582  |            8| HA\_A565C |  0.9983838|
| 3590  |            8| HA\_G335A |  0.9967431|
| 3651  |           80| HA\_A565C |  0.9990122|
| 3658  |           80| HA\_G335A |  0.9976453|
| 3713  |           81| HA\_A565C |  0.9987757|
| 3721  |           81| HA\_G335A |  0.9975695|
| 3779  |           82| HA\_A565C |  0.9987927|
| 3787  |           82| HA\_G335A |  0.9981104|
| 3837  |           84| HA\_A565C |  0.9993071|
| 3844  |           84| HA\_G335A |  0.9981885|
| 3900  |           85| HA\_A565C |  0.9990354|
| 3909  |           85| HA\_G335A |  0.9980108|
| 3967  |           86| HA\_A565C |  0.9985693|
| 3973  |           86| HA\_G335A |  0.9975667|
| 4029  |           88| HA\_A565C |  0.9991538|
| 4033  |           88| HA\_C214T |  0.9980784|
| 4036  |           88| HA\_C865T |  0.0136091|
| 4039  |           88| HA\_G335A |  0.9980286|
| 4093  |           89| HA\_A565C |  0.9988542|
| 4101  |           89| HA\_G335A |  0.9975472|
| 4158  |           98| HA\_A565C |  0.9989341|
| 4165  |           98| HA\_G335A |  0.9976566|
