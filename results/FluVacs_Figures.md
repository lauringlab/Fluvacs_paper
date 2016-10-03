---
title: "FluVacs Figures - Draft"
author: "JT"
date: "June 3, 2016"
output:
  html_document:
    toc: true
    toc_depth: 4
    toc_float: true
    
---




# Figure 4

High quality is >10^5^ or >10^3^ sequenced in duplicate.
Currently I am making these plots with just the 2007-2008 samples. Once I have the meta datafor the other seasons it will be trivial to add them (if we want to do that) 

## A) Variant Frequencies

These are the variants between 1-50%. each bin is 1% wide. The y axis is log scaled which in ggplot leads to a problem when there are no observations (lines below axis). I can remove these before publication.
<img src="figure/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```
## png 
##   2
```

<img src="figure/unnamed-chunk-3-2.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" style="display: block; margin: auto;" />

```
## png 
##   2
```

The x axis is so large because there are 2 mutations that are found in 55 and 35 samples. Both of these are infered minor variants. They may be artifacts of our analysis.




<img src="figure/unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```
## png 
##   2
```

<img src="figure/unnamed-chunk-4-2.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

```
## png 
##   2
```


There are very few iSNV in NA. 


# Table 3
Average iSNV/ segment

All data 

|chr |IIV           |LAIV          |PLACEBO       |
|:---|:-------------|:-------------|:-------------|
|HA  |1.43 ± (1)    |1.27 ± (0.5)  |1.27 ± (0)    |
|M   |1 ± (0)       |1.14 ± (0)    |1 ± (0)       |
|NP  |1.75 ± (1.25) |1.33 ± (0.75) |1.4 ± (1)     |
|NR  |1 ± (0)       |1.17 ± (0)    |1.25 ± (0.25) |
|NS  |1 ± (0)       |1.33 ± (0.75) |1.29 ± (0.5)  |
|PA  |1.8 ± (0.75)  |1.7 ± (1)     |1.65 ± (1)    |
|PB1 |1.8 ± (1)     |1.46 ± (1)    |1.45 ± (1)    |
|PB2 |2.1 ± (2)     |1.78 ± (1)    |1.89 ± (1)    |



#Supplemental Figure 1)

These are just from the first runs. I'm not including the duplicates here - should I? These have a sliding window of 100 with a step of 100 no overlap.

<img src="figure/coverage-1.png" title="plot of chunk coverage" alt="plot of chunk coverage" style="display: block; margin: auto;" /><img src="figure/coverage-2.png" title="plot of chunk coverage" alt="plot of chunk coverage" style="display: block; margin: auto;" /><img src="figure/coverage-3.png" title="plot of chunk coverage" alt="plot of chunk coverage" style="display: block; margin: auto;" />

```
## png 
##   2
```

```
## png 
##   2
```

```
## png 
##   2
```

If we plot on a log scale the bars are well above 0. 


```
## Scale for 'y' is already present. Adding another scale for 'y', which
## will replace the existing scale.
```

<img src="figure/unnamed-chunk-6-1.png" title="plot of chunk unnamed-chunk-6" alt="plot of chunk unnamed-chunk-6" style="display: block; margin: auto;" />





# Supplemental Figure 2)
<img src="figure/extra_iSNV-1.png" title="plot of chunk extra_iSNV" alt="plot of chunk extra_iSNV" style="display: block; margin: auto;" /><img src="figure/extra_iSNV-2.png" title="plot of chunk extra_iSNV" alt="plot of chunk extra_iSNV" style="display: block; margin: auto;" />

```
## png 
##   2
```

```
## png 
##   2
```


# Appendix
## Linear plots

<img src="figure/unnamed-chunk-8-1.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" style="display: block; margin: auto;" /><img src="figure/unnamed-chunk-8-2.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" style="display: block; margin: auto;" /><img src="figure/unnamed-chunk-8-3.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" style="display: block; margin: auto;" /><img src="figure/unnamed-chunk-8-4.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" style="display: block; margin: auto;" />

## Linear model with all variables

NB : I'm not sure what I'm doing here and if it's valid.


```
## 
## Call:
## lm(formula = lines.07$muts ~ lines.07$HAI.WI + lines.07$NAI.WI + 
##     lines.07$Copy_num + lines.07$dpi)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -5.029 -1.423  0.095  1.492  6.366 
## 
## Coefficients:
##                     Estimate Std. Error t value Pr(>|t|)    
## (Intercept)        4.458e+00  1.115e+00   4.000 0.000185 ***
## lines.07$HAI.WI    1.417e-03  8.037e-04   1.763 0.083211 .  
## lines.07$NAI.WI   -2.191e-03  4.628e-03  -0.473 0.637695    
## lines.07$Copy_num  5.292e-05  1.392e-05   3.802 0.000352 ***
## lines.07$dpi       3.643e-02  3.200e-01   0.114 0.909750    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 2.552 on 57 degrees of freedom
##   (2 observations deleted due to missingness)
## Multiple R-squared:  0.2451,	Adjusted R-squared:  0.1922 
## F-statistic: 4.628 on 4 and 57 DF,  p-value: 0.00264
```

I think this means that maybe the copy number is the only variable that affects the iSNV count. What if we remove the one outlier at 16.

```
## 
## Call:
## lm(formula = no_out.07$muts ~ no_out.07$HAI.WI + no_out.07$NAI.WI + 
##     no_out.07$Copy_num + no_out.07$dpi)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -5.029 -1.423  0.095  1.492  6.366 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)         4.458e+00  1.115e+00   4.000 0.000185 ***
## no_out.07$HAI.WI    1.417e-03  8.037e-04   1.763 0.083211 .  
## no_out.07$NAI.WI   -2.191e-03  4.628e-03  -0.473 0.637695    
## no_out.07$Copy_num  5.292e-05  1.392e-05   3.802 0.000352 ***
## no_out.07$dpi       3.643e-02  3.200e-01   0.114 0.909750    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 2.552 on 57 degrees of freedom
##   (2 observations deleted due to missingness)
## Multiple R-squared:  0.2451,	Adjusted R-squared:  0.1922 
## F-statistic: 4.628 on 4 and 57 DF,  p-value: 0.00264
```

Yep. I looks like that was driving things. 

I'll have to think about whether or not this analysis is justified, valid, and needed.










It looks like the iSNV count peaks around day 3 or 4. But I don't know if that's statistically significant or robust.

