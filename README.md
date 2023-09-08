CDK-mediated phosphoryations of IDRs
================
J.M. Valverde & G. Dubra

- [Supplemental information](#supplemental-information)
  - [Comparing predictors](#comparing-predictors)
  - [Supp. Fig. 5](#supp-fig-5)

This repository provides the source code and data for all the analyses
performed in the publication “A Cyclin dependent kinase-mediated
phosphorylation switch of disordered protein condensation” (Valverde,
Dubra et al., 2023). This current GitHub readme file was automaticaly
generated, using R Studio, from a R markdown (README.Rmd) where all the
code used in the paper is recompiled.

\#Main figures \##Fig. 3

\###b

#### Xenopus

![](README_files/figure-gfm/Expected%20Vs%20observed%20Phospho%20S/T%20in%20IDRs%20-%20xenopus-1.png)<!-- -->

#### Human

![](README_files/figure-gfm/Expected%20Vs%20observed%20Phospho%20S/T%20in%20IDRs%20-%20human-1.png)<!-- -->

\###c

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

\###d

![](README_files/figure-gfm/Enrichment%20of%20Phospho%20Ser%20and%20Ther%20in%20IDRs%20-1.png)<!-- -->

\###e

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

\####Human

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

\####Xenopus

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

\###f

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

# Supplemental information

## Comparing predictors

This block compares the residues predicted as to be disordered across
the entire proteomes of yeast and human using different disorders. A
Venn diagram is generated to illustrate the overlap between different
methods

![](README_files/figure-gfm/compare%20index%20-1.png)<!-- -->![](README_files/figure-gfm/compare%20index%20-2.png)<!-- -->

## Supp. Fig. 5

\###a \#### Yeast

![](README_files/figure-gfm/yeast%20composition%20-1.png)<!-- -->

#### Human

![](README_files/figure-gfm/human%20composition%20-1.png)<!-- -->

#### Xenopus

![](README_files/figure-gfm/xenopus%20composition%20-1.png)<!-- -->

\###d

![](README_files/figure-gfm/Lollipop%20plot%20for%20all%20human%20kinases%20-1.png)<!-- -->

\###e

Calculate the porcentage of disorder of human CDK targets vs the rest of
the phosphoproteome with all available predictors

![](README_files/figure-gfm/human%20disorder%20percentage%20for%20all%20predictors%20-1.png)<!-- -->![](README_files/figure-gfm/human%20disorder%20percentage%20for%20all%20predictors%20-2.png)<!-- -->

``` r
#######################################  HUMAN (Lila subset) ###############################################################

#########################################################################################################################################################


phosphoDiso_ST <- list()
phosphoDiso_obs <- numeric()
phosphoDiso_expct <- numeric()
phosphoDiso_expct_prob <- numeric()
for (i in 1:nrow(uniprot_all_predictions_phospho)) {
  # diso_fraction <- length(all_predictions_phospho[i,"disordered"][[1]])/nchar(all_predictions_phospho[i,"sequence"])
  # phosphoDiso_expct_uni[i] <- length(all_predictions_phospho[i,"psites"][[1]])*diso_fraction
  # Fraction of Ser And Thr that fall in disorder region
  TStotalIndexes <- as.numeric(gregexpr("S|T", uniprot_all_predictions_phospho[i,"sequence"])[[1]]) 
  TSinDiso_count <- sum(TStotalIndexes %in% uniprot_all_predictions_phospho[i,"disordered"][[1]])
  TSinDiso_fraction <- TSinDiso_count/length(TStotalIndexes)
  phosphoDiso_ST[[i]] <- TStotalIndexes
  phosphoDiso_expct_prob[i] <- TSinDiso_fraction
  phosphoDiso_expct[i] <- uniprot_all_predictions_phospho[i,"psites_count"][[1]]*TSinDiso_fraction
  phosphoDiso_obs[i] <- sum(uniprot_all_predictions_phospho[i,"psites"][[1]] %in% uniprot_all_predictions_phospho[i,"disordered"][[1]])
}

uniprot_all_predictions_phospho$ST_residues <- phosphoDiso_ST
uniprot_all_predictions_phospho$psites_obsv_diso <- phosphoDiso_obs
uniprot_all_predictions_phospho$psites_expct_diso <- phosphoDiso_expct
uniprot_all_predictions_phospho$psites_expct_diso_prob <- phosphoDiso_expct_prob



# only the embrios
# uniprot_all_predictions_phospho <- subset(uniprot_all_predictions_phospho, mpi.accesion %in% c("rif1","mcm4","cdc6","elys"))
# Non_Dynamic_embryoPsites_remapped_to_uniprot <- read_delim("utrech/db/uniprot/NonDynamic_embryoPsites_remapped_to_uniprot.tab", 
#                                                            delim = "\t", escape_double = FALSE, 
#                                                            trim_ws = TRUE)
# Non_Dynamic_embryoPsites_remapped_to_uniprot$non.dynamic_uniprot.psite <- lapply(Non_Dynamic_embryoPsites_remapped_to_uniprot$non.dynamic_uniprot.psite, function(x){ return(as.numeric(strsplit(x,";")[[1]]))})
# 
# uniprot_all_predictions_phospho <- merge(uniprot_all_predictions_phospho,Non_Dynamic_embryoPsites_remapped_to_uniprot,by.x = "mpi.accesion", by.y = "non.dynamic_mpi.accesion" )
# 
# 
# plotList <- IUpredScoresPlotGenerator(uniprot_all_predictions_phospho,id_col="ID",sites_col="non.dynamic_uniprot.psite",subset_sites_col="uniprot.psite",sequence_col="sequence")


# embrios + extracts
uniprot_all_predictions_phospho <- subset(uniprot_all_predictions_phospho, mpi.accesion %in% c("coil","npm","ki67","nup53","nup98","tp53b","nucl"))
Non_Dynamic_embryoPsites_remapped_to_uniprot <- read_delim("utrech/db/uniprot/NonDynamic_extractAndEmbryoPsites_remapped_to_uniprot", 
                                                           delim = "\t", escape_double = FALSE, 
                                                           trim_ws = TRUE)
Non_Dynamic_embryoPsites_remapped_to_uniprot$non.dynamic_uniprot.psite <- lapply(Non_Dynamic_embryoPsites_remapped_to_uniprot$non.dynamic_uniprot.psite, function(x){ return(as.numeric(strsplit(x,";")[[1]]))})

uniprot_all_predictions_phospho <- merge(uniprot_all_predictions_phospho,Non_Dynamic_embryoPsites_remapped_to_uniprot,by.x = "mpi.accesion", by.y = "non.dynamic_mpi.accesion",all.y = T)


plotList <- IUpredScoresPlotGenerator(uniprot_all_predictions_phospho,id_col="ID",sites_col="non.dynamic_uniprot.psite",subset_sites_col="uniprot.psite",sequence_col="sequence")



# 
# pdf("IUpredScores.pdf",width = 15,height = 3)
#  for (plot in plotList) {
#    print(plot)
#  }
# dev.off()

lilaSetXenopus <- IUpredScoresPlotGenerator(subset(all_predictions_phospho,ID %in% c("coil","npm","ki67","tp53b","nup53","nup98","nucl")),sites_col = "psites")

# JM_highly_phospho_mlos <- c("cndd3","dnli1","ube4b","tsc2","at2b1","rptor","caf1b","pcm1","abcf1","gemi5","cq028","tdrkh","tdrd6","ctr9","sf3b1","rbp2","armc9","chsp1","dkc1","eif3a","tacc3")
# 
# plotList_JM_highly_phospho_mlos <- IUpredScoresPlotGenerator(subset(all_predictions_phospho,ID %in% JM_highly_phospho_mlos))

pdf("../exportImages/pdfs/suppFig5/IUpredScores_plotList_JM_highly_phospho_mlos",width = 15,height = 3)
for (plot in plotList_JM_highly_phospho_mlos) {
  print(plot)
}
dev.off()


####################################################################################
#######################                    YEAST             #######################
####################################################################################
#CDK targets are the intersection of Holt and ubersax, but the phosphoproteome is only holt since it is the only set that contain residue level information.

phosphoDiso_ST <- list()
phosphoDiso_obs <- numeric()
phosphoDiso_expct <- numeric()
phosphoDiso_expct_prob <- numeric()
psites_count <- numeric()
#I select only those entries with at least one holt site
yeast_binomial_data <- subset(yeastIntersect_All_data, !isEmpty(site) & !is.na(site))

for (i in 1:nrow(yeast_binomial_data)) {
  TStotalIndexes <- as.numeric(gregexpr("S|T", yeast_binomial_data[i,"Sequence"])[[1]]) 
  disoIndexes_numeric <- as.numeric(strsplit(yeast_binomial_data[i,"iupl_disoIndexes"],",")[[1]])
  TSinDiso_count <- sum(TStotalIndexes %in% disoIndexes_numeric)
  TSinDiso_fraction <- TSinDiso_count/length(TStotalIndexes)
  phosphoDiso_ST[[i]] <- TStotalIndexes
  phosphoDiso_expct_prob[i] <- TSinDiso_fraction
  psites_count[i] <- length(yeast_binomial_data[i,"site"][[1]])
  phosphoDiso_expct[i] <- length(yeast_binomial_data[i,"site"][[1]])*TSinDiso_fraction
  phosphoDiso_obs[i] <- sum(yeast_binomial_data[i,"site"][[1]] %in% disoIndexes_numeric)
}

yeast_binomial_data$ST_residues <- phosphoDiso_ST
yeast_binomial_data$psites_obsv_diso <- phosphoDiso_obs
yeast_binomial_data$psites_expct_diso <- phosphoDiso_expct
yeast_binomial_data$psites_expct_diso_prob <- phosphoDiso_expct_prob
yeast_binomial_data$psites_count <- psites_count

yeast_binomial_data <- as.data.table(yeast_binomial_data)
yeast_binomial_data[,binom := purrr::pmap(.(x=psites_obsv_diso, n=psites_count, p=psites_expct_diso_prob), binom.test, alternative="greater")]
yeast_binomial_data[,binom_p := mapply("[[", binom, "p.value", SIMPLIFY = T)]
yeast_binomial_data[,binom_q := mapply(p.adjust, binom_p)]
yeast_binomial_data[,binom_sig := factor(ifelse(binom_q < 0.05, ifelse(binom_q < 0.01, "1% FDR", "5% FDR"), "n.s."), levels=c("n.s.","5% FDR", "1% FDR")) ]




ggplot(yeast_binomial_data) + 
  geom_point(aes(x=psites_obsv_diso,y=psites_expct_diso, colour = binom_sig,shape=target_intersect),size=2,alpha=0.80)+
  geom_abline(color="darkslategrey",slope = 1,size=0.5,linetype = "dashed")+
  ggpubr::theme_classic2() + 
  theme(text = element_text(size=15),legend.position = c(0.32,0.60),legend.box.just = "left",legend.box.margin = margin(2, 2, 2, 2),legend.box.background = element_rect(color="darkslategrey"),legend.title = element_text(size = 13)) +
  guides(color=guide_legend(title="Statistical significance")) +
  scale_x_continuous(limits = c(0, 60),breaks = c(seq(0, 60, by = 10)))+ xlab("Observed phospho S/T in IDR") +
  scale_y_continuous(limits = c(0, 60),breaks = c(seq(0, 60, by = 10)))+ ylab("Expected phospho S/T in IDR") + 
  scale_colour_manual(values = c(pal_jco()(10)[3],"#ffdd15ff",pal_jco()(10)[4]))+
  scale_shape_manual(values = c(16,4,4))
```
