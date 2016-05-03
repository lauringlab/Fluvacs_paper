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
var.2005.6.df <- subset(var.2005.6.df, season == "05-06")

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
    binaxis = "y", dotsize = 0.7) + ylab("Number of SNV")
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
make_heat_map(subset(all.qual.df, season == "05-06"))
```

    ## Aggregation function missing: defaulting to length

    ##    mutation       variable value segment
    ## 1  PB2_G96A 05-06.509.TRUE    -1     PB2
    ## 2 PB2_C191T 05-06.509.TRUE    -1     PB2
    ## 3 PB2_A225G 05-06.509.TRUE    -1     PB2
    ## 4 PB2_G236A 05-06.509.TRUE    -1     PB2
    ## 5 PB2_A264G 05-06.509.TRUE    -1     PB2
    ## 6 PB2_G282A 05-06.509.TRUE    -1     PB2

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-11-2.png" style="display: block; margin: auto;" />

``` r
make_heat_map(subset(all.qual.df, season == "04-05"))
```

    ##    mutation       variable     value segment
    ## 1  PB2_G96A 04-05.402.TRUE -1.000000     PB2
    ## 2 PB2_T249C 04-05.402.TRUE -1.000000     PB2
    ## 3 PB2_A316G 04-05.402.TRUE  1.998882     PB2
    ## 4 PB2_C462T 04-05.402.TRUE -1.000000     PB2
    ## 5 PB2_C465T 04-05.402.TRUE  1.999526     PB2
    ## 6 PB2_A486G 04-05.402.TRUE -1.000000     PB2

<img src="FluVacs_Figures_files/figure-markdown_github/unnamed-chunk-11-3.png" style="display: block; margin: auto;" />

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
| 227   |          426| HA\_T460G |  0.9988866|
| 228   |          426| HA\_G497A |  0.9983497|
| 229   |          426| HA\_T516C |  0.9995262|
| 230   |          426| HA\_A610G |  0.9988186|
| 232   |          426| HA\_A717G |  0.0344117|
| 289   |          404| HA\_T460G |  0.9971313|
| 290   |          404| HA\_G497A |  0.9976937|
| 291   |          404| HA\_T516C |  0.9975737|
| 294   |          404| HA\_A610G |  0.9938350|
| 393   |          411| HA\_T460G |  0.9970203|
| 394   |          411| HA\_G497A |  0.9974661|
| 395   |          411| HA\_T504C |  0.0198027|
| 396   |          411| HA\_T516C |  0.9982140|
| 397   |          411| HA\_A610G |  0.9972855|
| 457   |          409| HA\_C483T |  0.9939058|
| 458   |          409| HA\_T516C |  0.9986246|
| 459   |          409| HA\_C597T |  0.9941361|
| 460   |          409| HA\_A610G |  0.9983770|
| 557   |          401| HA\_T460G |  0.9968130|
| 558   |          401| HA\_T516C |  0.9976624|
| 559   |          401| HA\_A610G |  0.9975520|
| 643   |          402| HA\_T173C |  0.0082386|
| 644   |          402| HA\_C178G |  0.0138956|
| 648   |          402| HA\_T460G |  0.9983584|
| 650   |          402| HA\_A610G |  0.9979201|
| 655   |          402| HA\_G765A |  0.9972809|
| 751   |          403| HA\_C483T |  0.9977536|
| 752   |          403| HA\_T516C |  0.9986142|
| 753   |          403| HA\_C597T |  0.9960192|
| 754   |          403| HA\_A610G |  0.9977149|
| 871   |          414| HA\_G449A |  0.0087982|
| 872   |          414| HA\_C483T |  0.9987189|
| 873   |          414| HA\_T508C |  0.0062232|
| 874   |          414| HA\_T516C |  0.9987315|
| 875   |          414| HA\_T593C |  0.0059936|
| 876   |          414| HA\_C597T |  0.9972632|
| 877   |          414| HA\_A610G |  0.9981937|
| 880   |          414| HA\_T755C |  0.0064331|
| 1003  |          415| HA\_A288G |  0.0058375|
| 1004  |          415| HA\_T460G |  0.9986629|
| 1005  |          415| HA\_G497A |  0.9988006|
| 1007  |          415| HA\_T516C |  0.9989183|
| 1008  |          415| HA\_A610G |  0.9983850|
| 1011  |          415| HA\_G814A |  0.0057951|
| 1089  |          406| HA\_A199C |  0.0135042|
| 1091  |          406| HA\_T460G |  0.9980997|
| 1092  |          406| HA\_G497A |  0.9986187|
| 1093  |          406| HA\_T516C |  0.9987882|
| 1094  |          406| HA\_A610G |  0.9978076|
| 1176  |          400| HA\_G180A |  0.9952773|
| 1178  |          400| HA\_T460G |  0.9978954|
| 1179  |          400| HA\_T516C |  0.9982730|
| 1180  |          400| HA\_A610G |  0.9971264|
| 1260  |          427| HA\_T460G |  0.9971152|
| 1261  |          427| HA\_T516C |  0.9996469|
| 1263  |          427| HA\_A610G |  0.9971178|
| 1266  |          427| HA\_A630G |  0.0275211|
| 1269  |          427| HA\_G714A |  0.0353287|
| 1270  |          427| HA\_G715A |  0.9550557|
| 1274  |          427| HA\_C811A |  0.0323465|
| 1413  |          525| HA\_G419A |  0.9609481|
| 1414  |          525| HA\_T460G |  0.9712407|
| 1415  |          525| HA\_C483T |  0.0278879|
| 1416  |          525| HA\_T516C |  0.9984902|
| 1417  |          525| HA\_C597T |  0.0250218|
| 1418  |          525| HA\_A610G |  0.9971511|
| 1421  |          525| HA\_A630G |  0.9714209|
| 1423  |          525| HA\_G714A |  0.9720066|
| 1426  |          525| HA\_C811A |  0.9769375|
| 1578  |          504| HA\_A417G |  0.9971369|
| 1579  |          504| HA\_T460G |  0.9977722|
| 1580  |          504| HA\_T516C |  0.9982391|
| 1581  |          504| HA\_T517C |  0.9975879|
| 1582  |          504| HA\_A610G |  0.9979108|
| 1583  |          504| HA\_A611G |  0.0130542|
| 1586  |          504| HA\_A630G |  0.9986690|
| 1588  |          504| HA\_G714A |  0.9974301|
| 1591  |          504| HA\_C811A |  0.9979700|
| 1744  |          509| HA\_G405A |  0.0161785|
| 1745  |          509| HA\_T406C |  0.0196468|
| 1746  |          509| HA\_T460G |  0.9988192|
| 1747  |          509| HA\_T516C |  0.9991434|
| 1750  |          509| HA\_A610G |  0.9916580|
| 1753  |          509| HA\_A630G |  0.9986798|
| 1755  |          509| HA\_C711A |  0.0096927|
| 1756  |          509| HA\_G714A |  0.9977869|
| 1758  |          509| HA\_C811A |  0.9985649|
| 1934  |          501| HA\_G419A |  0.9964644|
| 1935  |          501| HA\_T460G |  0.9964523|
| 1936  |          501| HA\_T516C |  0.9977596|
| 1937  |          501| HA\_G540T |  0.0127358|
| 1938  |          501| HA\_A610G |  0.9969778|
| 1941  |          501| HA\_A630G |  0.9976991|
| 1943  |          501| HA\_G714A |  0.9976355|
| 1945  |          501| HA\_C811A |  0.9964795|
| 1946  |          501| HA\_A853G |  0.0117715|
| 2112  |          517| HA\_G262A |  0.0078170|
| 2113  |          517| HA\_A417G |  0.9982020|
| 2114  |          517| HA\_T460G |  0.9980793|
| 2115  |          517| HA\_T516C |  0.9983091|
| 2116  |          517| HA\_T517C |  0.9984622|
| 2117  |          517| HA\_A610G |  0.9983691|
| 2120  |          517| HA\_A630G |  0.9985089|
| 2122  |          517| HA\_T636C |  0.0102764|
| 2124  |          517| HA\_G671T |  0.0078651|
| 2125  |          517| HA\_G714A |  0.9975177|
| 2127  |          517| HA\_C811A |  0.9979849|
| 2129  |          517| HA\_G958A |  0.0881559|
| 2323  |          520| HA\_C393A |  0.9939498|
| 2324  |          520| HA\_T460G |  0.9969847|
| 2325  |          520| HA\_T516C |  0.9979719|
| 2326  |          520| HA\_A610G |  0.9978417|
| 2329  |          520| HA\_A630G |  0.9974505|
| 2331  |          520| HA\_G714A |  0.9977224|
| 2333  |          520| HA\_C811A |  0.9963940|
| 2507  |          511| HA\_T460G |  0.9967962|
| 2508  |          511| HA\_T516C |  0.9976284|
| 2509  |          511| HA\_A610G |  0.9969909|
| 2512  |          511| HA\_A630G |  0.9976982|
| 2515  |          511| HA\_G714A |  0.9974383|
| 2517  |          511| HA\_C811A |  0.9969660|
| 2652  |          512| HA\_G205A |  0.9932573|
| 2653  |          512| HA\_T460G |  0.9972270|
| 2654  |          512| HA\_T516C |  0.9988915|
| 2655  |          512| HA\_G544A |  0.0190153|
| 2656  |          512| HA\_G603T |  0.9945140|
| 2657  |          512| HA\_A610G |  0.9987225|
| 2660  |          512| HA\_A630G |  0.9970140|
| 2662  |          512| HA\_G714A |  0.9971530|
| 2664  |          512| HA\_A768G |  0.9948011|
| 2665  |          512| HA\_C811A |  0.9971420|
| 2828  |          528| HA\_A417G |  0.0767076|
| 2829  |          528| HA\_G419A |  0.7289497|
| 2830  |          528| HA\_T441C |  0.7056674|
| 2831  |          528| HA\_T460G |  0.9996533|
| 2832  |          528| HA\_T516C |  0.9996469|
| 2833  |          528| HA\_A610G |  0.9993957|
| 2836  |          528| HA\_A630G |  0.9997433|
| 2838  |          528| HA\_G714A |  0.9804632|
| 2840  |          528| HA\_C811A |  0.9995596|
| 3015  |          530| HA\_G419A |  0.9982332|
| 3016  |          530| HA\_T441C |  0.9966070|
| 3017  |          530| HA\_T460G |  0.9983211|
| 3018  |          530| HA\_T516C |  0.9987062|
| 3019  |          530| HA\_A610G |  0.9982769|
| 3022  |          530| HA\_A630G |  0.9989105|
| 3024  |          530| HA\_G714A |  0.9965330|
| 3026  |          530| HA\_C811A |  0.9976806|
| 3178  |          507| HA\_T460G |  0.9984070|
| 3179  |          507| HA\_T492C |  0.0134297|
| 3180  |          507| HA\_T516C |  0.9986697|
| 3181  |          507| HA\_A607G |  0.0232290|
| 3182  |          507| HA\_A610G |  0.9983556|
| 3185  |          507| HA\_A630G |  0.9987599|
| 3187  |          507| HA\_G714A |  0.9979765|
| 3189  |          507| HA\_C811A |  0.9980754|
| 3192  |          507| HA\_A910C |  0.0061613|
| 3368  |          510| HA\_T460G |  0.9970391|
| 3369  |          510| HA\_T516C |  0.9981613|
| 3371  |          510| HA\_A610G |  0.9979895|
| 3374  |          510| HA\_A630G |  0.9980323|
| 3376  |          510| HA\_G714A |  0.9971246|
| 3378  |          510| HA\_C811A |  0.9971579|
| 3523  |          506| HA\_T305C |  0.0119297|
| 3525  |          506| HA\_C393A |  0.9947674|
| 3526  |          506| HA\_T421C |  0.0106132|
| 3527  |          506| HA\_T460G |  0.9970598|
| 3530  |          506| HA\_T516C |  0.9982807|
| 3531  |          506| HA\_T549C |  0.0052135|
| 3532  |          506| HA\_A610G |  0.9974200|
| 3535  |          506| HA\_A630G |  0.9981211|
| 3539  |          506| HA\_G714A |  0.9969325|
| 3542  |          506| HA\_G769A |  0.0103862|
| 3543  |          506| HA\_C811A |  0.9975394|
| 3757  |          515| HA\_G354A |  0.9962710|
| 3758  |          515| HA\_T460G |  0.9989409|
| 3759  |          515| HA\_T477C |  0.9967306|
| 3760  |          515| HA\_T499C |  0.9968776|
| 3761  |          515| HA\_T516C |  0.9991458|
| 3762  |          515| HA\_T518C |  0.9971857|
| 3764  |          515| HA\_A610G |  0.9989148|
| 3767  |          515| HA\_A630G |  0.9975270|
| 3769  |          515| HA\_C711T |  0.9971129|
| 3770  |          515| HA\_G714A |  0.9973777|
| 3772  |          515| HA\_C811A |  0.9972853|
| 3865  |          508| HA\_G419A |  0.9986021|
| 3866  |          508| HA\_T460G |  0.9985142|
| 3867  |          508| HA\_T516C |  0.9989525|
| 3868  |          508| HA\_A610G |  0.9985918|
| 3871  |          508| HA\_A630G |  0.9987475|
| 3873  |          508| HA\_G714A |  0.9982015|
| 3875  |          508| HA\_C811A |  0.9986847|
| 4037  |          524| HA\_A443G |  0.0083427|
| 4038  |          524| HA\_T460G |  0.9984386|
| 4039  |          524| HA\_T499C |  0.9964607|
| 4040  |          524| HA\_T516C |  0.9988556|
| 4041  |          524| HA\_C558T |  0.0085430|
| 4042  |          524| HA\_A610G |  0.9983090|
| 4045  |          524| HA\_A630G |  0.9986854|
| 4047  |          524| HA\_A638G |  0.0071251|
| 4048  |          524| HA\_C711T |  0.9965614|
| 4049  |          524| HA\_G714A |  0.9984147|
| 4053  |          524| HA\_T807C |  0.9953625|
| 4054  |          524| HA\_C811A |  0.9985608|
| 4175  |          503| HA\_T460G |  0.9966631|
| 4177  |          503| HA\_T516C |  0.9972142|
| 4178  |          503| HA\_C531T |  0.0089009|
| 4179  |          503| HA\_A610G |  0.9967145|
| 4182  |          503| HA\_A630G |  0.9975460|
| 4184  |          503| HA\_G714A |  0.9965141|
| 4186  |          503| HA\_C811A |  0.9971348|
| 4336  |          529| HA\_G419A |  0.9981430|
| 4337  |          529| HA\_T441C |  0.9976882|
| 4338  |          529| HA\_T460G |  0.9983396|
| 4339  |          529| HA\_T516C |  0.9987139|
| 4340  |          529| HA\_A610G |  0.9981029|
| 4343  |          529| HA\_A630G |  0.9986309|
| 4345  |          529| HA\_G714A |  0.9978133|
| 4347  |          529| HA\_C811A |  0.9980774|
| 4348  |          529| HA\_T954C |  0.0105081|
| 4499  |           20| HA\_G335A |  0.9977312|
| 4500  |           20| HA\_C440T |  0.0071559|
| 4501  |           20| HA\_A565C |  0.9985224|
| 4594  |           58| HA\_G335A |  0.9963896|
| 4596  |           58| HA\_A565C |  0.9976782|
| 4714  |           15| HA\_G335A |  0.9971824|
| 4715  |           15| HA\_A565C |  0.9983310|
| 4781  |           28| HA\_G335A |  0.9970985|
| 4783  |           28| HA\_A565C |  0.9981969|
| 4853  |           71| HA\_G335A |  0.9954906|
| 4854  |           71| HA\_A565C |  0.9969662|
| 4919  |           10| HA\_G335A |  0.9961708|
| 4921  |           10| HA\_A565C |  0.9977655|
| 4923  |           10| HA\_A652G |  0.0088028|
| 5036  |           55| HA\_G335A |  0.9973494|
| 5037  |           55| HA\_C411T |  0.0680573|
| 5039  |           55| HA\_A565C |  0.9984146|
| 5131  |           40| HA\_G335A |  0.9973318|
| 5132  |           40| HA\_G434A |  0.0231497|
| 5133  |           40| HA\_A565C |  0.9985938|
| 5139  |           40| HA\_G854A |  0.9979359|
| 5144  |           40| HA\_G958A |  0.0455798|
| 5242  |           67| HA\_G335A |  0.9956266|
| 5243  |           67| HA\_A565C |  0.9972097|
| 5244  |           67| HA\_G610A |  0.0080955|
| 5332  |           64| HA\_G335A |  0.9961368|
| 5333  |           64| HA\_A565C |  0.9982134|
| 5403  |           74| HA\_G335A |  0.9939727|
| 5404  |           74| HA\_A565C |  0.9978975|
| 5480  |           27| HA\_G335A |  0.9972820|
| 5481  |           27| HA\_A565C |  0.9978890|
| 5484  |           27| HA\_G854A |  0.9986577|
| 5540  |           31| HA\_G335A |  0.9982555|
| 5541  |           31| HA\_A565C |  0.9988404|
| 5599  |           66| HA\_G335A |  0.9972223|
| 5600  |           66| HA\_A565C |  0.9979457|
| 5674  |           17| HA\_G335A |  0.9926969|
| 5675  |           17| HA\_A565C |  0.9920638|
| 5736  |           52| HA\_G335A |  0.9974599|
| 5737  |           52| HA\_A565C |  0.9983659|
| 5739  |           52| HA\_G685A |  0.0087493|
| 5816  |           49| HA\_G335A |  0.9973413|
| 5817  |           49| HA\_A565C |  0.9981817|
| 5882  |           59| HA\_G335A |  0.9969594|
| 5883  |           59| HA\_A565C |  0.9977509|
| 5997  |           34| HA\_G335A |  0.9964811|
| 5998  |           34| HA\_A565C |  0.9980709|
| 6067  |           70| HA\_G335A |  0.9968609|
| 6068  |           70| HA\_T518C |  0.9977758|
| 6069  |           70| HA\_A565G |  0.9981506|
| 6071  |           70| HA\_A717G |  0.0247364|
| 6073  |           70| HA\_C865T |  0.0162663|
| 6184  |           61| HA\_G335A |  0.9964441|
| 6185  |           61| HA\_A565C |  0.9978608|
| 6258  |           41| HA\_G335A |  0.9978889|
| 6259  |           41| HA\_A565C |  0.9984537|
| 6262  |           41| HA\_G854A |  0.9985950|
| 6327  |           30| HA\_G335A |  0.9973585|
| 6328  |           30| HA\_T462C |  0.0075009|
| 6329  |           30| HA\_A565C |  0.9983373|
| 6419  |            3| HA\_G335A |  0.9969682|
| 6420  |            3| HA\_A565C |  0.9982774|
| 6422  |            3| HA\_G671A |  0.9839711|
| 6494  |           73| HA\_G335A |  0.9964744|
| 6495  |           73| HA\_A565C |  0.9970282|
| 6565  |           57| HA\_G335A |  0.9974149|
| 6566  |           57| HA\_A565C |  0.9983615|
| 6633  |            6| HA\_C164A |  0.0399776|
| 6634  |            6| HA\_T228C |  0.0362521|
| 6635  |            6| HA\_C285A |  0.0431470|
| 6636  |            6| HA\_G335A |  0.9973246|
| 6637  |            6| HA\_A565C |  0.9981962|
| 6638  |            6| HA\_T594C |  0.0182470|
| 6639  |            6| HA\_A611G |  0.0058156|
| 6641  |            6| HA\_A651G |  0.0075603|
| 6737  |           56| HA\_G335A |  0.9967348|
| 6738  |           56| HA\_A565C |  0.9981518|
| 6834  |           72| HA\_A199G |  0.0116255|
| 6835  |           72| HA\_G335A |  0.9968578|
| 6836  |           72| HA\_C420T |  0.0135588|
| 6837  |           72| HA\_C548T |  0.0080640|
| 6838  |           72| HA\_A565C |  0.9981816|
| 6947  |           76| HA\_G335A |  0.9972253|
| 6948  |           76| HA\_G428A |  0.0184960|
| 6949  |           76| HA\_A565C |  0.9980131|
| 6950  |           76| HA\_T606G |  0.0073840|
| 7095  |            8| HA\_G335A |  0.9972604|
| 7096  |            8| HA\_A565C |  0.9979842|
| 7164  |           19| HA\_G335A |  0.9977268|
| 7165  |           19| HA\_A565C |  0.9986290|
| 7228  |           75| HA\_G335A |  0.9971999|
| 7231  |           75| HA\_A565C |  0.9978136|
| 7316  |           68| HA\_T311C |  0.0110569|
| 7317  |           68| HA\_G335A |  0.9969388|
| 7318  |           68| HA\_A565C |  0.9983324|
| 7319  |           68| HA\_G592A |  0.0192200|
| 7426  |           26| HA\_G335A |  0.9978690|
| 7427  |           26| HA\_A565C |  0.9983244|
| 7429  |           26| HA\_A672G |  0.1156970|
| 7499  |            4| HA\_G335A |  0.9958891|
| 7500  |            4| HA\_A565C |  0.9978261|
| 7567  |           16| HA\_G335A |  0.9946479|
| 7568  |           16| HA\_A565C |  0.9960869|
| 7631  |           18| HA\_G335A |  0.9884812|
| 7632  |           18| HA\_A565C |  0.9956490|
| 7694  |           78| HA\_G335A |  0.9967651|
| 7695  |           78| HA\_A565C |  0.9979265|
| 7757  |            7| HA\_G335A |  0.9968246|
| 7758  |            7| HA\_A565C |  0.9980902|
| 7820  |            5| HA\_G335A |  0.9973509|
| 7822  |            5| HA\_A565C |  0.9977732|
| 7916  |           51| HA\_G335A |  0.9963630|
| 7917  |           51| HA\_A565C |  0.9981302|
| 7920  |           51| HA\_C865T |  0.9970511|
| 7989  |           24| HA\_G335A |  0.9976960|
| 7990  |           24| HA\_A565C |  0.9987801|
| 8065  |           25| HA\_G335A |  0.9981932|
| 8066  |           25| HA\_G497A |  0.9993614|
| 8067  |           25| HA\_A565C |  0.9986573|
| 8128  |           46| HA\_G335A |  0.9975693|
| 8129  |           46| HA\_A565C |  0.9984604|
| 8189  |           33| HA\_G335A |  0.9979487|
| 8190  |           33| HA\_A565C |  0.9983568|
| 8261  |           50| HA\_A315G |  0.9984685|
| 8262  |           50| HA\_G335A |  0.9974576|
| 8263  |           50| HA\_A565C |  0.9984302|
| 8324  |            2| HA\_G335A |  0.9962115|
| 8325  |            2| HA\_A565C |  0.9980567|
| 8387  |           39| HA\_G335A |  0.9964492|
| 8388  |           39| HA\_A565C |  0.9974803|
| 8453  |           29| HA\_G335A |  0.9974372|
| 8454  |           29| HA\_A565C |  0.9984071|
| 8519  |           54| HA\_G335A |  0.9964840|
| 8520  |           54| HA\_A399G |  0.0082548|
| 8521  |           54| HA\_A565C |  0.9981359|
| 8523  |           54| HA\_G671A |  0.9815564|
| 8585  |           23| HA\_G335A |  0.9977300|
| 8586  |           23| HA\_G419A |  0.0082951|
| 8587  |           23| HA\_A565C |  0.9986095|
| 8592  |           23| HA\_C959T |  0.0392059|
| 8674  |           77| HA\_G335A |  0.9954350|
| 8675  |           77| HA\_A399G |  0.0101209|
| 8676  |           77| HA\_A565C |  0.9946714|
| 8790  |           22| HA\_G335A |  0.9973459|
| 8791  |           22| HA\_A565C |  0.9986670|
| 8856  |           43| HA\_G335A |  0.9973633|
| 8857  |           43| HA\_A565C |  0.9974667|
| 8859  |           43| HA\_G674A |  0.0054497|
| 8924  |           79| HA\_G335A |  0.9958628|
| 8925  |           79| HA\_T518C |  0.9978475|
| 8928  |           79| HA\_C709A |  0.0073681|
| 8931  |           79| HA\_C950A |  0.0098849|
| 9050  |           42| HA\_G335A |  0.9972996|
| 9051  |           42| HA\_T459C |  0.0100849|
| 9052  |           42| HA\_A565C |  0.9985562|
| 9116  |           47| HA\_G335A |  0.9979190|
| 9117  |           47| HA\_A565C |  0.9985586|
| 9183  |           80| HA\_G335A |  0.9978706|
| 9184  |           80| HA\_G464A |  0.0051679|
| 9185  |           80| HA\_G485A |  0.0082520|
| 9186  |           80| HA\_T504C |  0.0063246|
| 9187  |           80| HA\_A565C |  0.9988508|
| 9296  |           38| HA\_T276C |  0.0361201|
| 9297  |           38| HA\_C333A |  0.0286524|
| 9298  |           38| HA\_G335A |  0.9989949|
| 9299  |           38| HA\_G494A |  0.0290098|
| 9300  |           38| HA\_A565C |  0.9994678|
| 9373  |           14| HA\_G335A |  0.9979290|
| 9375  |           14| HA\_A565C |  0.9987825|
| 9378  |           14| HA\_C709T |  0.0073778|
| 9380  |           14| HA\_C809T |  0.0067067|
| 9381  |           14| HA\_G854A |  0.1048269|
| 9502  |           93| HA\_G335A |  0.9983432|
| 9504  |           93| HA\_A565C |  0.9991695|
| 9508  |           93| HA\_C954T |  0.9970880|
| 9577  |           91| HA\_G335A |  0.9991905|
| 9578  |           91| HA\_A565C |  0.9996051|
| 9636  |           12| HA\_G335A |  0.9979693|
| 9638  |           12| HA\_A565C |  0.9988072|
| 9639  |           12| HA\_G603A |  0.9925620|
| 9701  |          103| HA\_G335A |  0.9984687|
| 9702  |          103| HA\_A565C |  0.9986705|
| 9704  |          103| HA\_G715A |  0.0058576|
| 9776  |           36| HA\_G335A |  0.9983351|
| 9778  |           36| HA\_A565C |  0.9990739|
| 9779  |           36| HA\_T585A |  0.0222295|
| 9783  |           36| HA\_T918C |  0.0063328|
| 9878  |           85| HA\_G205A |  0.0149806|
| 9879  |           85| HA\_G335A |  0.9983505|
| 9880  |           85| HA\_A565C |  0.9989570|
| 9971  |           89| HA\_G335A |  0.9974933|
| 9973  |           89| HA\_A565C |  0.9986828|
| 10047 |          106| HA\_G335A |  0.9984824|
| 10048 |          106| HA\_A565C |  0.9988101|
| 10113 |           83| HA\_G335A |  0.9989679|
| 10114 |           83| HA\_A565C |  0.9992604|
| 10116 |           83| HA\_T663C |  0.9981590|
| 10119 |           83| HA\_C942T |  0.9986578|
| 10196 |          104| HA\_G335A |  0.9977110|
| 10198 |          104| HA\_A565C |  0.9987507|
| 10309 |           98| HA\_G335A |  0.9981452|
| 10312 |           98| HA\_A565C |  0.9989137|
| 10392 |           86| HA\_G335A |  0.9977122|
| 10393 |           86| HA\_A565C |  0.9982039|
| 10470 |           82| HA\_G335A |  0.9986986|
| 10471 |           82| HA\_A565C |  0.9987050|
| 10475 |           82| HA\_G854A |  0.0083207|
| 10477 |           82| HA\_A910G |  0.0114305|
| 10547 |           88| HA\_C214T |  0.9978096|
| 10548 |           88| HA\_G262A |  0.0062458|
| 10549 |           88| HA\_G335A |  0.9983515|
| 10550 |           88| HA\_A565C |  0.9989711|
| 10554 |           88| HA\_A754G |  0.0066067|
| 10555 |           88| HA\_C865T |  0.0137764|
| 10638 |          105| HA\_G335A |  0.9982814|
| 10639 |          105| HA\_T509C |  0.0099163|
| 10640 |          105| HA\_A565C |  0.9989095|
| 10643 |          105| HA\_T743C |  0.0071381|
| 10746 |          100| HA\_G335A |  0.9983947|
| 10747 |          100| HA\_T518C |  0.9957963|
| 10827 |           11| HA\_G335A |  0.9977941|
| 10828 |           11| HA\_A565C |  0.9987086|
| 10829 |           11| HA\_G603A |  0.9980866|
| 10891 |           84| HA\_G335A |  0.9987552|
| 10892 |           84| HA\_A565C |  0.9994419|
| 10974 |           37| HA\_G335A |  0.9986999|
| 10975 |           37| HA\_A565C |  0.9989790|
| 11066 |           44| HA\_G335A |  0.9985839|
| 11067 |           44| HA\_A565C |  0.9988419|
| 11070 |           44| HA\_G671A |  0.9626091|
| 11145 |           96| HA\_G335A |  0.9989645|
| 11146 |           96| HA\_A565C |  0.9989827|
| 11230 |           81| HA\_G335A |  0.9980835|
| 11231 |           81| HA\_A565C |  0.9985278|
| 11302 |           13| HA\_G335A |  0.9977982|
| 11304 |           13| HA\_A565C |  0.9987868|
| 11307 |           13| HA\_G854A |  0.4659828|
| 11397 |           92| HA\_G335A |  0.9983227|
| 11398 |           92| HA\_A565C |  0.9990205|
| 11479 |          101| HA\_G335A |  0.9976464|
| 11480 |          101| HA\_T518C |  0.9985878|
| 11565 |           90| HA\_G335A |  0.9963350|
| 11566 |           90| HA\_A565C |  0.9991500|

### High quality

``` r
kable(subset(all.qual.df, antigenic == T, select = c(Lauring_Id, mutation, freq.var)))
```

|       |  Lauring\_Id| mutation  |   freq.var|
|-------|------------:|:----------|----------:|
| 289   |          404| HA\_T460G |  0.9971313|
| 290   |          404| HA\_G497A |  0.9976937|
| 291   |          404| HA\_T516C |  0.9975737|
| 294   |          404| HA\_A610G |  0.9938350|
| 557   |          401| HA\_T460G |  0.9968130|
| 558   |          401| HA\_T516C |  0.9976624|
| 559   |          401| HA\_A610G |  0.9975520|
| 192   |          400| HA\_A610G |  0.9977904|
| 195   |          400| HA\_G180A |  0.9961603|
| 197   |          400| HA\_T460G |  0.9985080|
| 198   |          400| HA\_T516C |  0.9986825|
| 257   |          402| HA\_A610G |  0.9982760|
| 259   |          402| HA\_G765A |  0.9977265|
| 262   |          402| HA\_T460G |  0.9987133|
| 283   |          403| HA\_A610G |  0.9982779|
| 288   |          403| HA\_C483T |  0.9983432|
| 2891  |          403| HA\_C597T |  0.9967209|
| 2931  |          403| HA\_T516C |  0.9988845|
| 3631  |          406| HA\_A610G |  0.9982362|
| 3661  |          406| HA\_G497A |  0.9984981|
| 3681  |          406| HA\_T460G |  0.9985823|
| 3691  |          406| HA\_T516C |  0.9991032|
| 404   |          409| HA\_A610G |  0.9986096|
| 409   |          409| HA\_C483T |  0.9967033|
| 410   |          409| HA\_C597T |  0.9962182|
| 414   |          409| HA\_T516C |  0.9989004|
| 483   |          414| HA\_A610G |  0.9984139|
| 487   |          414| HA\_C483T |  0.9988346|
| 488   |          414| HA\_C597T |  0.9965694|
| 492   |          414| HA\_T516C |  0.9989146|
| 5631  |          426| HA\_A610G |  0.9942401|
| 5661  |          426| HA\_G497A |  0.7869945|
| 5671  |          426| HA\_T460G |  0.9940300|
| 5681  |          426| HA\_T516C |  0.9796461|
| 6051  |          427| HA\_A610G |  0.9982441|
| 6071  |          427| HA\_A630G |  0.5127296|
| 6101  |          427| HA\_C811A |  0.5151538|
| 6131  |          427| HA\_G714A |  0.5166147|
| 6151  |          427| HA\_T460G |  0.9985207|
| 6161  |          427| HA\_T516C |  0.9979014|
| 3322  |          509| HA\_G405A |  0.0161785|
| 3332  |          509| HA\_T406C |  0.0196468|
| 3342  |          509| HA\_T460G |  0.9988192|
| 3352  |          509| HA\_T516C |  0.9991434|
| 3382  |          509| HA\_A610G |  0.9916580|
| 3412  |          509| HA\_A630G |  0.9986798|
| 3432  |          509| HA\_C711A |  0.0096927|
| 3442  |          509| HA\_G714A |  0.9977869|
| 3462  |          509| HA\_C811A |  0.9985649|
| 1095  |          511| HA\_T460G |  0.9967962|
| 1096  |          511| HA\_T516C |  0.9976284|
| 1097  |          511| HA\_A610G |  0.9969909|
| 1100  |          511| HA\_A630G |  0.9976982|
| 1103  |          511| HA\_G714A |  0.9974383|
| 1105  |          511| HA\_C811A |  0.9969660|
| 1603  |          530| HA\_G419A |  0.9982332|
| 1604  |          530| HA\_T441C |  0.9966070|
| 1605  |          530| HA\_T460G |  0.9983211|
| 1606  |          530| HA\_T516C |  0.9987062|
| 1607  |          530| HA\_A610G |  0.9982769|
| 1610  |          530| HA\_A630G |  0.9989105|
| 1612  |          530| HA\_G714A |  0.9965330|
| 1614  |          530| HA\_C811A |  0.9976806|
| 1766  |          507| HA\_T460G |  0.9984070|
| 1767  |          507| HA\_T492C |  0.0134297|
| 1768  |          507| HA\_T516C |  0.9986697|
| 1769  |          507| HA\_A607G |  0.0232290|
| 1770  |          507| HA\_A610G |  0.9983556|
| 1773  |          507| HA\_A630G |  0.9987599|
| 1775  |          507| HA\_G714A |  0.9979765|
| 1777  |          507| HA\_C811A |  0.9980754|
| 1780  |          507| HA\_A910C |  0.0061613|
| 1956  |          510| HA\_T460G |  0.9970391|
| 1957  |          510| HA\_T516C |  0.9981613|
| 1959  |          510| HA\_A610G |  0.9979895|
| 1962  |          510| HA\_A630G |  0.9980323|
| 1964  |          510| HA\_G714A |  0.9971246|
| 1966  |          510| HA\_C811A |  0.9971579|
| 2763  |          503| HA\_T460G |  0.9966631|
| 2765  |          503| HA\_T516C |  0.9972142|
| 2766  |          503| HA\_C531T |  0.0089009|
| 2767  |          503| HA\_A610G |  0.9967145|
| 2770  |          503| HA\_A630G |  0.9975460|
| 2772  |          503| HA\_G714A |  0.9965141|
| 2774  |          503| HA\_C811A |  0.9971348|
| 2924  |          529| HA\_G419A |  0.9981430|
| 2925  |          529| HA\_T441C |  0.9976882|
| 2926  |          529| HA\_T460G |  0.9983396|
| 2927  |          529| HA\_T516C |  0.9987139|
| 2928  |          529| HA\_A610G |  0.9981029|
| 29311 |          529| HA\_A630G |  0.9986309|
| 2933  |          529| HA\_G714A |  0.9978133|
| 2935  |          529| HA\_C811A |  0.9980774|
| 2936  |          529| HA\_T954C |  0.0105081|
| 11741 |          504| HA\_A417G |  0.9976149|
| 11751 |          504| HA\_A610G |  0.9984881|
| 11771 |          504| HA\_A630G |  0.9991493|
| 11801 |          504| HA\_C811A |  0.9986808|
| 11831 |          504| HA\_G714A |  0.9977156|
| 11861 |          504| HA\_T460G |  0.9987647|
| 11871 |          504| HA\_T516C |  0.9987600|
| 11881 |          504| HA\_T517C |  0.9976600|
| 1305  |          508| HA\_A610G |  0.9987240|
| 1307  |          508| HA\_A630G |  0.9990099|
| 1310  |          508| HA\_C811A |  0.9988566|
| 1312  |          508| HA\_G419A |  0.9980482|
| 1314  |          508| HA\_G714A |  0.9984538|
| 1317  |          508| HA\_T460G |  0.9989013|
| 1318  |          508| HA\_T516C |  0.9991749|
| 1474  |          529| HA\_A610G |  0.9983249|
| 1476  |          529| HA\_A630G |  0.9988436|
| 1479  |          529| HA\_C811A |  0.9983761|
| 1481  |          529| HA\_G419A |  0.9982299|
| 1483  |          529| HA\_G714A |  0.9969878|
| 1486  |          529| HA\_T441C |  0.9979795|
| 1487  |          529| HA\_T460G |  0.9987829|
| 1488  |          529| HA\_T516C |  0.9989257|
| 14992 |           34| HA\_G335A |  0.9964811|
| 15002 |           34| HA\_A565C |  0.9980709|
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
| 4302  |          105| HA\_A565C |  0.9989337|
| 4382  |          105| HA\_G335A |  0.9978429|
| 4932  |          106| HA\_A565C |  0.9986871|
| 5002  |          106| HA\_G335A |  0.9979818|
| 5582  |           11| HA\_A565C |  0.9988134|
| 5652  |           11| HA\_G335A |  0.9976365|
| 5662  |           11| HA\_G603A |  0.9979748|
| 6192  |           12| HA\_A565C |  0.9989583|
| 6282  |           12| HA\_G335A |  0.9978017|
| 6292  |           12| HA\_G603A |  0.9929640|
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
| 991   |           19| HA\_A565C |  0.9987964|
| 998   |           19| HA\_G335A |  0.9974281|
| 1052  |            2| HA\_A565C |  0.9985182|
| 1060  |            2| HA\_G335A |  0.9967804|
| 11131 |           22| HA\_A565C |  0.9988689|
| 11221 |           22| HA\_G335A |  0.9973136|
| 11792 |           25| HA\_A565C |  0.9988729|
| 11862 |           25| HA\_G335A |  0.9977840|
| 11872 |           25| HA\_G497A |  0.9976218|
| 12421 |           26| HA\_A565C |  0.9986016|
| 12431 |           26| HA\_A672G |  0.1305439|
| 12501 |           26| HA\_G335A |  0.9966141|
| 13061 |           27| HA\_A565C |  0.9984427|
| 13141 |           27| HA\_G335A |  0.9969219|
| 13151 |           27| HA\_G854A |  0.9979336|
| 13671 |           28| HA\_A565C |  0.9986050|
| 13741 |           28| HA\_G335A |  0.9974500|
| 14281 |           29| HA\_A565C |  0.9986560|
| 14351 |           29| HA\_G335A |  0.9971020|
| 14941 |            3| HA\_A565C |  0.9986414|
| 15001 |            3| HA\_G335A |  0.9969509|
| 15011 |            3| HA\_G671A |  0.9851019|
| 15631 |           31| HA\_A565C |  0.9989570|
| 15701 |           31| HA\_G335A |  0.9979034|
| 16212 |           33| HA\_A565C |  0.9985810|
| 16292 |           33| HA\_G335A |  0.9976974|
| 16831 |           34| HA\_A565C |  0.9984352|
| 16901 |           34| HA\_G335A |  0.9968830|
| 17511 |           39| HA\_A565C |  0.9982252|
| 17581 |           39| HA\_G335A |  0.9969295|
| 18171 |            4| HA\_A565C |  0.9982301|
| 18241 |            4| HA\_G335A |  0.9962522|
| 18811 |           41| HA\_A565C |  0.9987356|
| 18901 |           41| HA\_G335A |  0.9976531|
| 18911 |           41| HA\_G854A |  0.9985191|
| 19491 |           42| HA\_A565C |  0.9984543|
| 19571 |           42| HA\_G335A |  0.9967091|
| 19601 |           42| HA\_T459C |  0.0097621|
| 20141 |           43| HA\_A565C |  0.9980756|
| 20211 |           43| HA\_G335A |  0.9972581|
| 20791 |           46| HA\_A565C |  0.9986484|
| 20871 |           46| HA\_G335A |  0.9975462|
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
| 27831 |           61| HA\_A565C |  0.9984241|
| 27921 |           61| HA\_G335A |  0.9967259|
| 28511 |           64| HA\_A565C |  0.9986240|
| 28611 |           64| HA\_G335A |  0.9967796|
| 29212 |           66| HA\_A565C |  0.9983496|
| 29291 |           66| HA\_G335A |  0.9973845|
| 29881 |           67| HA\_A565C |  0.9980889|
| 29971 |           67| HA\_G335A |  0.9964781|
| 30531 |            7| HA\_A565C |  0.9984670|
| 30612 |            7| HA\_G335A |  0.9970704|
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
| 35121 |           79| HA\_G335A |  0.9966175|
| 3514  |           79| HA\_T518C |  0.9978933|
| 35821 |            8| HA\_A565C |  0.9983838|
| 3590  |            8| HA\_G335A |  0.9967431|
| 36511 |           80| HA\_A565C |  0.9990122|
| 3658  |           80| HA\_G335A |  0.9976453|
| 3713  |           81| HA\_A565C |  0.9987757|
| 37211 |           81| HA\_G335A |  0.9975695|
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
| 41011 |           89| HA\_G335A |  0.9975472|
| 4158  |           98| HA\_A565C |  0.9989341|
| 4165  |           98| HA\_G335A |  0.9976566|
