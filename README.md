CDK-mediated phosphoryations of IDRs
================
J.M. Valverde & G. Dubra

This repository provides the source code and data for all the analyses
performed in the publication “A Cyclin dependent kinase-mediated
phosphorylation switch of disordered protein condensation” (Valverde,
Dubra et al., 2023). This current GitHub readme file was automaticaly
generated, using R Studio, from a R markdown (README.Rmd) where all the
code used in the paper is recompiled.

- [Main figures](#main-figures)
  - [Fig. 2](#fig-2)
    - [b](#b)
    - [d](#d)
    - [e](#e)
  - [Fig. 3](#fig-3)
    - [b](#b-1)
    - [c](#c)
    - [d](#d-1)
    - [e](#e-1)
    - [f](#f)
  - [Fig. 4](#fig-4)
    - [a](#a)
  - [Fig. 6](#fig-6)
    - [c](#c-1)
    - [d](#d-2)
  - [Fig. 7](#fig-7)
    - [c](#c-2)
- [Supplemental information](#supplemental-information)
  - [Comparing predictors](#comparing-predictors)
  - [Supp. Fig. 5](#supp-fig-5)
    - [a](#a-1)
    - [d](#d-3)
    - [e](#e-2)
  - [Supp. Fig. 7](#supp-fig-7)
    - [Deleted LR domain](#deleted-lr-domain)


# Main figures

## Fig. 2

### b

![](README_files/figure-gfm/embryo%20heatmap%20-%20High-time%20resolutionr-1.png)<!-- -->
\### c

![](README_files/figure-gfm/single%20psites%20-%20biological%20processes%20-1.png)<!-- -->![](README_files/figure-gfm/single%20psites%20-%20biological%20processes%20-2.png)<!-- -->![](README_files/figure-gfm/single%20psites%20-%20biological%20processes%20-3.png)<!-- -->![](README_files/figure-gfm/single%20psites%20-%20biological%20processes%20-4.png)<!-- -->

### d

![](README_files/figure-gfm/single%20psites%20-%20CDK1%20Y15-1.png)<!-- -->

### e

## Fig. 3

### b

#### Xenopus

![](README_files/figure-gfm/Expected%20Vs%20observed%20Phospho%20S/T%20in%20IDRs%20-%20xenopus-1.png)<!-- -->

#### Human

![](README_files/figure-gfm/Expected%20Vs%20observed%20Phospho%20S/T%20in%20IDRs%20-%20human-1.png)<!-- -->

### c

    ## 
    ##  Wilcoxon signed rank test with continuity correction
    ## 
    ## data:  psites_diso by variable
    ## V = 206396, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -1.211715 -1.052488
    ## sample estimates:
    ## (pseudo)median 
    ##      -1.131956

![](README_files/figure-gfm/Expected%20vs%20Observed%20phospho%20S/T%20distribution%20in%20xenopus-1.png)<!-- -->

    ## 
    ##  Wilcoxon signed rank test with continuity correction
    ## 
    ## data:  psites_diso by variable
    ## V = 22188, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -1.890087 -1.491565
    ## sample estimates:
    ## (pseudo)median 
    ##      -1.679388

![](README_files/figure-gfm/Expected%20vs%20Observed%20phospho%20S/T%20distribution%20in%20xenopus-2.png)<!-- -->

    ## 
    ##  Wilcoxon signed rank test with continuity correction
    ## 
    ## data:  psites_diso by variable
    ## V = 1340, p-value = 1.019e-14
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -2.734868 -1.577823
    ## sample estimates:
    ## (pseudo)median 
    ##      -2.131005

![](README_files/figure-gfm/Expected%20vs%20Observed%20phospho%20S/T%20distribution%20in%20xenopus-3.png)<!-- -->

### d

![](README_files/figure-gfm/Enrichment%20of%20Phospho%20Ser%20and%20Ther%20in%20IDRs%20-1.png)<!-- -->

### e

#### Yeast

![](README_files/figure-gfm/Holt%20disorder%20percentage-1.png)<!-- -->

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Percentage of disorder by target
    ## W = 167690, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  21.63845 31.47638
    ## sample estimates:
    ## difference in location 
    ##               26.64479

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Percentage of disorder by target
    ## W = 166252, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  23.21311 34.28578
    ## sample estimates:
    ## difference in location 
    ##               28.64125

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Percentage of disorder by target
    ## W = 165710, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  20.31818 30.31861
    ## sample estimates:
    ## difference in location 
    ##               25.25882

#### Human

![](README_files/figure-gfm/human%20disorder%20percentage-1.png)<!-- -->

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Percentage of disorder by target
    ## W = 7707108, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  18.79420 22.82949
    ## sample estimates:
    ## difference in location 
    ##                20.8597

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Percentage of disorder by target
    ## W = 7269014, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  16.62156 21.14240
    ## sample estimates:
    ## difference in location 
    ##               18.87039

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Percentage of disorder by target
    ## W = 7348738, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  16.82185 21.06136
    ## sample estimates:
    ## difference in location 
    ##               18.95022

#### Xenopus

![](README_files/figure-gfm/xenopus%20disorder%20percentage-1.png)<!-- -->

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Percentage of disorder by ANOVA
    ## W = 307032, p-value = 2.358e-13
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -10.849034  -6.152332
    ## sample estimates:
    ## difference in location 
    ##               -8.45844

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Percentage of disorder by ANOVA
    ## W = 278128, p-value = 4.832e-13
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -12.27983  -7.01995
    ## sample estimates:
    ## difference in location 
    ##              -9.613807

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Percentage of disorder by ANOVA
    ## W = 319366, p-value = 5.781e-10
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -10.862103  -5.641645
    ## sample estimates:
    ## difference in location 
    ##              -8.265359

### f

![](README_files/figure-gfm/disorder%20percentages%20for%20other%20kinases-1.png)<!-- -->

    ##            subset     n   mean     sd   min     Q1 median     Q3     max
    ## 1 Phosphoproteome 15329 22.767 24.195 0.000  3.333 13.687 36.158 100.000
    ## 2             CDK   656 43.040 25.881 0.000 24.184 40.389 61.776 100.000
    ## 3            MAPK   408 36.539 25.173 0.000 15.237 31.296 53.273  98.487
    ## 4            AURK   412 44.438 24.141 0.291 27.811 40.845 60.857 100.000
    ## 5             PLK   459 39.847 25.215 0.000 18.771 38.403 57.566 100.000
    ## 6             NEK    43 28.806 22.094 2.448 10.339 25.741 43.054  85.624
    ## 7            DYRK    42 41.869 25.545 0.239 22.507 40.408 59.029  95.763
    ##   percZero
    ## 1    9.485
    ## 2    0.762
    ## 3    1.225
    ## 4    0.000
    ## 5    0.654
    ## 6    0.000
    ## 7    0.000

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Percentage of disorder by subset
    ## Kruskal-Wallis chi-squared = 1064.1, df = 6, p-value < 2.2e-16

![](README_files/figure-gfm/post-hoc%20test%20for%20all%20kinases%20disorder%20percentage-1.png)<!-- -->

    ##            subset     n nvalid   mean     sd   min     Q1 median     Q3     max
    ## 1 Phosphoproteome 15329  15300 30.599 28.045 0.000  8.231 20.914 45.864 100.000
    ## 2             CDK   656    656 48.864 28.506 0.000 26.305 46.592 70.344 100.000
    ## 3            MAPK   408    408 41.036 28.074 0.000 15.497 35.394 61.269 100.000
    ## 4            AURK   412    412 50.031 27.638 0.000 28.291 44.873 72.685 100.000
    ## 5             PLK   459    459 46.872 29.169 0.346 21.257 45.455 69.699 100.000
    ## 6             NEK    43     43 33.073 26.195 0.000 10.570 30.476 50.330 100.000
    ## 7            DYRK    42     42 46.485 26.216 6.310 25.558 45.340 61.612  99.548
    ##   percZero
    ## 1    0.791
    ## 2    0.457
    ## 3    0.245
    ## 4    0.243
    ## 5    0.000
    ## 6    2.326
    ## 7    0.000

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Percentage of disorder by subset
    ## Kruskal-Wallis chi-squared = 680.57, df = 6, p-value < 2.2e-16

![](README_files/figure-gfm/post-hoc%20test%20for%20all%20kinases%20disorder%20percentage-2.png)<!-- -->

    ##            subset     n nvalid   mean     sd    min     Q1 median     Q3
    ## 1 Phosphoproteome 15329  15309 43.241 25.440  1.980 21.942 37.717 62.359
    ## 2             CDK   656    656 61.190 23.712 10.583 43.233 61.230 80.563
    ## 3            MAPK   408    408 54.506 23.997  7.042 33.354 54.216 76.560
    ## 4            AURK   412    412 62.686 21.933 10.738 45.588 62.748 80.284
    ## 5             PLK   459    459 59.330 24.128 10.406 38.498 60.535 79.393
    ## 6             NEK    43     43 49.955 22.278 14.161 29.756 52.847 63.672
    ## 7            DYRK    42     42 59.229 24.553 19.616 40.189 59.548 74.960
    ##       max
    ## 1 100.000
    ## 2 100.000
    ## 3 100.000
    ## 4 100.000
    ## 5 100.000
    ## 6  95.455
    ## 7 100.000

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Percentage of disorder by subset
    ## Kruskal-Wallis chi-squared = 751.44, df = 6, p-value < 2.2e-16

![](README_files/figure-gfm/post-hoc%20test%20for%20all%20kinases%20disorder%20percentage-3.png)<!-- -->

## Fig. 4

### a

![](README_files/figure-gfm/CDK%20targets%20and%20Dynamic%20phosphoproteins%20in%20MLOs-1.png)<!-- -->
\### b

#### Human

![](README_files/figure-gfm/IUpred%20scores%20-%20human%20-1.png)<!-- -->![](README_files/figure-gfm/IUpred%20scores%20-%20human%20-2.png)<!-- -->![](README_files/figure-gfm/IUpred%20scores%20-%20human%20-3.png)<!-- -->![](README_files/figure-gfm/IUpred%20scores%20-%20human%20-4.png)<!-- -->![](README_files/figure-gfm/IUpred%20scores%20-%20human%20-5.png)<!-- -->![](README_files/figure-gfm/IUpred%20scores%20-%20human%20-6.png)<!-- -->![](README_files/figure-gfm/IUpred%20scores%20-%20human%20-7.png)<!-- -->

#### Xenopus

![](README_files/figure-gfm/IUpred%20scores%20-%20xenopus%20-1.png)<!-- -->![](README_files/figure-gfm/IUpred%20scores%20-%20xenopus%20-2.png)<!-- -->![](README_files/figure-gfm/IUpred%20scores%20-%20xenopus%20-3.png)<!-- -->![](README_files/figure-gfm/IUpred%20scores%20-%20xenopus%20-4.png)<!-- -->![](README_files/figure-gfm/IUpred%20scores%20-%20xenopus%20-5.png)<!-- -->![](README_files/figure-gfm/IUpred%20scores%20-%20xenopus%20-6.png)<!-- -->

## Fig. 6

### c

#### Interphase

![](README_files/figure-gfm/FRAP%20modelling%20-%20interphase%20-1.png)<!-- -->

#### Interphase (high expression)

![](README_files/figure-gfm/FRAP%20modelling%20-%20interphase%20(High%20expression)%20-1.png)<!-- -->

#### Mitosis

![](README_files/figure-gfm/FRAP%20modelling%20-%20mitosis-1.png)<!-- -->

### d

![](README_files/figure-gfm/FRAP%20quantification%20-1.png)<!-- -->

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  thalf by group
    ## W = 18.5, p-value = 0.1232
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -13.370047   1.150055
    ## sample estimates:
    ## difference in location 
    ##              -4.977137

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  thalf by group
    ## W = 0, p-value = 0.01421
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -23.110004  -6.819955
    ## sample estimates:
    ## difference in location 
    ##              -15.94311

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  thalf by group
    ## W = 15, p-value = 0.008276
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -15.209962  -5.260037
    ## sample estimates:
    ## difference in location 
    ##              -11.23145

![](README_files/figure-gfm/FRAP%20quantification%20-2.png)<!-- -->

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  percent by group
    ## W = 73, p-value = 0.005126
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  11.17000 45.76001
    ## sample estimates:
    ## difference in location 
    ##               27.58075

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  percent by group
    ## W = 24, p-value = 0.01421
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  17.58995 51.27007
    ## sample estimates:
    ## difference in location 
    ##                30.0655

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  percent by group
    ## W = 66, p-value = 0.5885
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -10.63997  17.64000
    ## sample estimates:
    ## difference in location 
    ##               2.790034

## Fig. 7

### c

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Foci_Nucleus by Treat_collapsed
    ## Kruskal-Wallis chi-squared = 205.15, df = 5, p-value < 2.2e-16

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
    ## 
    ## data:  Full_length_quant$Foci_Nucleus and Full_length_quant$Treat_collapsed 
    ## 
    ##         DMSO L- DMSO L+ OA L-   OA L+   PA L- 
    ## DMSO L+ 2.1e-14 -       -       -       -     
    ## OA L-   5.3e-11 0.0230  -       -       -     
    ## OA L+   8.6e-14 0.3149  0.1860  -       -     
    ## PA L-   7.0e-05 < 2e-16 < 2e-16 < 2e-16 -     
    ## PA L+   0.6424  2.0e-15 1.8e-11 2.1e-14 0.0011
    ## 
    ## P value adjustment method: BH

![](README_files/figure-gfm/optodroplet%20quantification%20-%20FL%20-1.png)<!-- -->

# Supplemental information

## Comparing predictors

This block compares the residues predicted as to be disordered across
the entire proteomes of yeast and human using different disorders. A
Venn diagram is generated to illustrate the overlap between different
methods

![](README_files/figure-gfm/compare%20index%20-1.png)<!-- -->![](README_files/figure-gfm/compare%20index%20-2.png)<!-- -->

## Supp. Fig. 5

### a

#### Yeast

![](README_files/figure-gfm/yeast%20composition%20-1.png)<!-- -->

#### Human

![](README_files/figure-gfm/human%20composition%20-1.png)<!-- -->

#### Xenopus

![](README_files/figure-gfm/xenopus%20composition%20-1.png)<!-- -->

### d

![](README_files/figure-gfm/Lollipop%20plot%20for%20all%20human%20kinases%20-1.png)<!-- -->

### e

Calculate the porcentage of disorder of yeast and human CDK targets vs
the rest of the phosphoproteome with all available predictors

![](README_files/figure-gfm/human%20disorder%20percentage%20for%20all%20predictors%20-1.png)<!-- -->![](README_files/figure-gfm/human%20disorder%20percentage%20for%20all%20predictors%20-2.png)<!-- -->

## Supp. Fig. 7

### Deleted LR domain

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Foci_Nucleus by Treat_collapsed
    ## Kruskal-Wallis chi-squared = 204.5, df = 5, p-value < 2.2e-16

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
    ## 
    ## data:  dLR_quant$Foci_Nucleus and dLR_quant$Treat_collapsed 
    ## 
    ##         DMSO L- DMSO L+ OA L-   OA L+   PA L- 
    ## DMSO L+ < 2e-16 -       -       -       -     
    ## OA L-   2.3e-05 1.5e-10 -       -       -     
    ## OA L+   < 2e-16 0.8143  2.9e-09 -       -     
    ## PA L-   0.7618  < 2e-16 2.9e-05 < 2e-16 -     
    ## PA L+   0.3916  < 2e-16 0.0014  8.9e-16 0.5839
    ## 
    ## P value adjustment method: BH

![](README_files/figure-gfm/optodroplet%20quantification%20-%20dLR%20-1.png)<!-- -->

### Deleted N-terminus (Repeats + LR)

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Foci_Nucleus by Treat_collapsed
    ## Kruskal-Wallis chi-squared = 204.5, df = 5, p-value < 2.2e-16

    ## 
    ##  Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
    ## 
    ## data:  Repeats_LR_quant$Foci_Nucleus and Repeats_LR_quant$Treat_collapsed 
    ## 
    ##         DMSO L- DMSO L+ OA L-   OA L+   PA L-  
    ## DMSO L+ 8.0e-09 -       -       -       -      
    ## OA L-   0.032   4.4e-05 -       -       -      
    ## OA L+   1.4e-09 0.341   3.3e-06 -       -      
    ## PA L-   3.3e-06 < 2e-16 3.0e-11 < 2e-16 -      
    ## PA L+   0.776   1.4e-09 0.014   2.8e-10 7.7e-06
    ## 
    ## P value adjustment method: BH

![](README_files/figure-gfm/optodroplet%20quantification%20-%20Deleted%20N-terminus%20(Repeats%20+%20LR)%20-1.png)<!-- -->
