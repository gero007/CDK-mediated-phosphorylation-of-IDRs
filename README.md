CDK-mediated phosphoryations of IDRs
================

This repository provides the source code and data for all the analyses performed in the publication "A Cyclin dependent kinase-mediated phosphorylation switch of disordered protein condensation" (Valverde, Dubra et al., 2023)


Authrs: J.M. Valverde & G. Dubra

- [Human](#human)
- [Fig. 3d](#fig-3d)
- [Fig. 3e](#fig-3e)
  - [Yeast datasets intersection](#yeast-datasets-intersection)
  - [Human](#human-1)
  - [Xenopus](#xenopus)
  - [Violin plots for the disorder percentage in all kinases’ targets
    for the three main predictors used in this
    paper](#violin-plots-for-the-disorder-percentage-in-all-kinases-targets-for-the-three-main-predictors-used-in-this-paper)
  - [Post-hoc test hypothesis for all possible
    comparisons](#post-hoc-test-hypothesis-for-all-possible-comparisons)
- [Supplemental information](#supplemental-information)
  - [Comparing predictors](#comparing-predictors)
  - [Supp. Fig. 5a](#supp-fig-5a)
    - [Yeast](#yeast)
    - [Human](#human-2)
    - [Xenopus](#xenopus-1)
- [Supp. Fig. 5d](#supp-fig-5d)
- [Supp. Fig. 5e](#supp-fig-5e)

\#Fig. 3b

\## Xenopus

``` r
all_protein_id_phospho <- read_lines("phosphorylation_data/Valverde_data/ids/xenopus_allphosphosites_mpi.txt")
cluster_A <- read_lines("phosphorylation_data/Valverde_data/ids/xenopusclusterA_mpi.txt")
cluster_B <- read_lines("phosphorylation_data/Valverde_data/ids/xenopusclusterB_mpi.txt")
cluster_C <- read_lines("phosphorylation_data/Valverde_data/ids/xenopusclusterC_mpi.txt")
cluster_D <- read_lines("phosphorylation_data/Valverde_data/ids/xenopusclusterD_mpi.txt")
significant_anova <- Reduce(union, list(cluster_A,cluster_B,cluster_C,cluster_D))
human_CDK1targets <- read_lines("phosphorylation_data/Valverde_data/ids/humanCDKtargets_mpi.txt")
xenopus_ANOVA_data <- read_delim("phosphorylation_data/Valverde_data/Xen_phospho_AnovaPhosphosites_curated.txt","\t", escape_double = FALSE, trim_ws = TRUE)
xenopus_extract_ANOVA_data <- read_delim("phosphorylation_data/Valverde_data/Xen_Extracts_phospho_AnovaPhosphosites_curated.txt","\t", escape_double = FALSE, trim_ws = TRUE)


disoPath <- "disorder_data/local_predictions/IUpred_run_xenopus/"
file.names <- dir(disoPath, pattern =".iupred")
ids <-unlist(strsplit(x = file.names,split = ".iupred",fixed = T))
all_predictions <- data.frame(ids)
colnames(all_predictions) <- "ID"
all_predictions$ID <- as.character(all_predictions$ID)
positions <- list()
scores <- list()
disordered <- list()

for (n in 1:length(file.names)) {
  aux_table <- read_delim(paste(disoPath,file.names[n],sep = ""),"\t", escape_double = FALSE, col_names = FALSE,col_types = cols(X1 = col_integer(),X2 = col_character(),X3 = col_double()),comment = "#", trim_ws = TRUE)
  aux_table <- as.data.frame(aux_table)
  colnames(aux_table) <- c("POS","RES","IUPRED_SCORE")
  aux_table$IUPRED_DISO <- aux_table$IUPRED_SCORE>=0.5
  all_predictions[n,"sequence"] <- paste(aux_table$RES,collapse = "")
  positions[[n]] <- as.numeric(aux_table$POS)
  scores[[n]] <- as.numeric(aux_table$IUPRED_SCORE)
  # all_predictions[n,"phospho"] <- NULL
  all_predictions[n,"threshold"] <- 0.5
  disordered[[n]] <- which(aux_table$IUPRED_DISO)
}
all_predictions$positions <-positions
all_predictions$IUPredScores <-scores
all_predictions$disordered <-disordered

# remove unannotated proteins
# all_protein_id_phospho <- all_protein_id_phospho[all_protein_id_phospho %in% all_predictions$ID]
all_predictions_phospho <- subset(all_predictions, ID %in% all_protein_id_phospho)
all_predictions_phospho$length <- nchar(all_predictions_phospho$sequence)

# Introduce phosphosites (from all_proteins variable) using the apply function
phosphosites <- read_delim("phosphorylation_data/Valverde_data/all_phosphosites_xenopus_nr.tab", "\t", escape_double = FALSE, col_types = cols(`Leading proteins` = col_skip(), Protein = col_skip()), trim_ws = TRUE)

# Add phosphosites from the extracts!!! uncomment if wanted. 
phosphosites_extracts <- read_delim("phosphorylation_data/Valverde_data/all_phosphosites_xenopus_extracts_nr.tab", "\t", escape_double = FALSE, trim_ws = TRUE)
phosphosites <- rbind(phosphosites,phosphosites_extracts)


phosphosites <- phosphosites %>% dplyr::rename(ID=Proteins,psites=`Positions within proteins`,seqWindow=`Sequence window`,UID=`Unique identifier`) %>% group_by(ID) %>% summarise_at(c("psites","seqWindow","UID"),function(x){paste(x, collapse=",")})
# The data has been opened by excel and some proteins have their names modified to be dates.


phosphosites$psites <- lapply(phosphosites$psites, function(x){ return(as.numeric(strsplit(x,",")[[1]]))})
phosphosites$psites_count <- sapply(phosphosites$psites, length)
phosphosites$UID <- lapply(phosphosites$UID, function(x){ return(as.character(strsplit(x,",")[[1]]))})

xenopus_dynamic <- xenopus_ANOVA_data$UID
# Add phosphosites from the extracts!!! uncomment if wanted.
xenopus_dynamic <- c(xenopus_dynamic,xenopus_extract_ANOVA_data$UID)
#

phosphosites$anova_psites <- apply(phosphosites, 1, function(x){x$psites[x$UID %in% xenopus_dynamic]})
all_predictions_phospho <- merge.data.frame(all_predictions_phospho,phosphosites,by = "ID")


all_predictions_phospho <- all_predictions_phospho %>% mutate(anova_sig=case_when(
  ID %in% significant_anova ~ "Dynamic",
  TRUE ~ "Non dynamic"
))


all_predictions_phospho <- all_predictions_phospho %>% mutate(hCDK1target=case_when(
  ID %in% human_CDK1targets ~ "Human CDK1 target",
  TRUE ~ "other"
))


phosphoDiso_ST <- list()
phosphoDiso_obs <- numeric()
phosphoDiso_expct <- numeric()
phosphoDiso_expct_prob <- numeric()
for (i in 1:nrow(all_predictions_phospho)) {
  TStotalIndexes <- as.numeric(gregexpr("S|T", all_predictions_phospho[i,"sequence"])[[1]]) 
  TSinDiso_count <- sum(TStotalIndexes %in% all_predictions_phospho[i,"disordered"][[1]])
  TSinDiso_fraction <- TSinDiso_count/length(TStotalIndexes)
  phosphoDiso_ST[[i]] <- TStotalIndexes
  phosphoDiso_expct_prob[i] <- TSinDiso_fraction
  phosphoDiso_expct[i] <- all_predictions_phospho[i,"psites_count"][[1]]*TSinDiso_fraction
  phosphoDiso_obs[i] <- sum(all_predictions_phospho[i,"psites"][[1]] %in% all_predictions_phospho[i,"disordered"][[1]])
}

all_predictions_phospho$ST_residues <- phosphoDiso_ST
all_predictions_phospho$psites_obsv_diso <- phosphoDiso_obs
all_predictions_phospho$psites_expct_diso <- phosphoDiso_expct
all_predictions_phospho$psites_expct_diso_prob <- phosphoDiso_expct_prob


# Binomial test:
# Calculate percentage of disordered region from all proteins
all_predictions_phospho <- as.data.table(all_predictions_phospho)
all_predictions_phospho[,binom := purrr::pmap(.(psites_obsv_diso, psites_count, psites_expct_diso_prob), binom.test, alternative="greater")]
all_predictions_phospho[,binom_p := mapply("[[", binom, "p.value", SIMPLIFY = T)]
all_predictions_phospho[,binom_q := mapply(p.adjust, binom_p)]
all_predictions_phospho[,binom_sig := factor(ifelse(binom_q < 0.05, ifelse(binom_q < 0.01, "1% FDR", "5% FDR"), "n.s."), levels=c("n.s.","5% FDR", "1% FDR")) ]



ggplot(all_predictions_phospho) + 
  geom_point(aes(x=psites_obsv_diso,y=psites_expct_diso, colour = binom_sig,shape=anova_sig),size=2,alpha=0.80)+
  geom_abline(color="darkslategrey",slope = 1,size=0.5,linetype = "dashed")+
  ggpubr::theme_classic2() + 
  theme(text = element_text(size=15),legend.position = c(0.32,0.60),legend.box.just = "left",legend.box.margin = margin(2, 2, 2, 2),legend.box.background = element_rect(color="darkslategrey"),legend.title = element_text(size = 13)) +
  guides(color=guide_legend(title="Statistical significance")) +
  scale_x_continuous(limits = c(0, 30),breaks = c(seq(0, 30, by = 5)))+ xlab("Observed phospho S/T in IDR") +
  scale_y_continuous(limits = c(0, 30),breaks = c(seq(0, 30, by = 5)))+ ylab("Expected phospho S/T in IDR") + 
  scale_colour_manual(values = c(pal_jco()(10)[3],"#ffdd15ff",pal_jco()(10)[4]))+
  scale_shape_manual(values = c(16,4,4))
```

## Human

``` r
disoPath <- "disorder_data/local_predictions/IUpred_run_human/"
file.names <- dir(disoPath, pattern =".iupred")
ids <-unlist(strsplit(x = file.names,split = ".iupred",fixed = T))
all_predictions <- data.frame(ids)
colnames(all_predictions) <- "ID"
all_predictions$ID <- as.character(all_predictions$ID)
positions <- list()
scores <- list()
disordered <- list()

for (n in 1:length(file.names)) {
  aux_table <- read_delim(paste(disoPath,file.names[n],sep = ""),"\t", escape_double = FALSE, col_names = FALSE,col_types = cols(X1 = col_integer(),X2 = col_character(),X3 = col_double()),comment = "#", trim_ws = TRUE)
  aux_table <- as.data.frame(aux_table)
  colnames(aux_table) <- c("POS","RES","IUPRED_SCORE")
  aux_table$IUPRED_DISO <- aux_table$IUPRED_SCORE>=0.5
  all_predictions[n,"sequence"] <- paste(aux_table$RES,collapse = "")
  positions[[n]] <- as.numeric(aux_table$POS)
  scores[[n]] <- as.numeric(aux_table$IUPRED_SCORE)
  # all_predictions[n,"phospho"] <- NULL
  # all_predictions[n,"threshold"] <- 0.5
  disordered[[n]] <- which(aux_table$IUPRED_DISO)
}
all_predictions$positions <-positions
all_predictions$IUPredScores <-scores
all_predictions$disordered <-disordered

# remove unannotated proteins
# all_protein_id_phospho <- all_protein_id_phospho[all_protein_id_phospho %in% all_predictions$ID]
all_predictions$length <- nchar(all_predictions$sequence)





human_data <- read_delim("phosphorylation_data/phosphositePlus_data/human_data_curated_V2_cleaned.tab", 
                         "\t", escape_double = FALSE, col_types = cols(MOD_RSD = col_character()), 
                         trim_ws = TRUE)

human_data$psites_CDK1 <- lapply(human_data$MOD_RSD, function(x){ return(as.numeric(strsplit(x,",")[[1]]))})
human_data$target <- rep("Cdk1 target",nrow(human_data))
human_data$target <- as.factor(human_data$target)

human_universe_data <- read_delim("phosphorylation_data/phosphositePlus_data/Phosphorylation_site_dataset", 
                                  "\t", escape_double = FALSE, col_types = cols(HU_CHR_LOC = col_skip(), 
                                                                                SITE_GRP_ID = col_skip(), MW_kD = col_skip(), 
                                                                                DOMAIN = col_skip(), `SITE_+/-7_AA` = col_skip(), 
                                                                                LT_LIT = col_skip(), MS_LIT = col_skip(), 
                                                                                MS_CST = col_skip(), `CST_CAT#` = col_skip()), 
                                  trim_ws = TRUE)

human_universe_data<-rename(human_universe_data,c(`ACC#`=ACC_ID))
# Select human proteins
human_universe_data <- subset(human_universe_data, ORGANISM == "human")
# Remove all information related to isoforms
human_universe_data <- subset(human_universe_data, !grepl("-",`ACC#`))
human_universe_data <- subset(human_universe_data, !grepl(" iso[0-9]",`PROTEIN`))
# Select targets with pS or pT
human_universe_data<-subset(human_universe_data,(substr(MOD_RSD,1,1)=="S"|substr(MOD_RSD,1,1)=="T"))
# format the MOD_RSD column and group by gene/protein/uniprot
human_universe_data <- human_universe_data %>% mutate(MOD_RSD=substr(MOD_RSD,1,nchar(MOD_RSD)-2)) 
human_universe_data <- human_universe_data %>% group_by(`ACC#`,GENE,PROTEIN) %>% summarise_at("MOD_RSD",function(x){paste(substr(x,2,2000), collapse=",")})
# Generate the psite column, with the vectors containing psites (for the contingency table analysis)
human_universe_data <- as.data.frame(human_universe_data)
human_universe_data$psites <- lapply(human_universe_data$MOD_RSD, function(x){ return(as.numeric(strsplit(x,",")[[1]]))})

# merge the tables and mark CDK1 targets and non CDK1 targets. Only by ACC, protein names in human data are still with  the isoform nomenclature
human_data <- merge.data.frame(x = human_data,y = human_universe_data,by = c("ACC#"),all = T,suffixes = c("_CDK1","_ALL"))
#remove PROTEIN and GENE comumns from x, and rename the y columns
human_data$GENE_CDK1<-NULL
human_data$PROTEIN_CDK1<-NULL
human_data <- rename(human_data,c(PROTEIN=PROTEIN_ALL,GENE=GENE_ALL))
# Adding the target category "Non Cdk1 target"
levels(human_data$target) <- c("Cdk1 target","Non Cdk1 target")
human_data$target[is.na(human_data$target)]<-"Non Cdk1 target"

human_data<-merge.data.frame(human_data,all_predictions,by.x = "ACC#",by.y = "ID")
human_data$psites_count <- sapply(human_data$psites, length)
human_data$psites_CDK1_count <- sapply(human_data$psites_CDK1, length)


phosphoDiso_ST <- list()
phosphoDiso_obs <- numeric()
phosphoDiso_expct <- numeric()
phosphoDiso_expct_prob <- numeric()
for (i in 1:nrow(human_data)) {
  TStotalIndexes <- as.numeric(gregexpr("S|T", human_data[i,"sequence"])[[1]]) 
  TSinDiso_count <- sum(TStotalIndexes %in% human_data[i,"disordered"][[1]])
  TSinDiso_fraction <- TSinDiso_count/length(TStotalIndexes)
  phosphoDiso_ST[[i]] <- TStotalIndexes
  phosphoDiso_expct_prob[i] <- TSinDiso_fraction
  phosphoDiso_expct[i] <- human_data[i,"psites_count"][[1]]*TSinDiso_fraction
  phosphoDiso_obs[i] <- sum(human_data[i,"psites"][[1]] %in% human_data[i,"disordered"][[1]])
}

human_data$ST_residues <- phosphoDiso_ST
human_data$psites_obsv_diso <- phosphoDiso_obs
human_data$psites_expct_diso <- phosphoDiso_expct
human_data$psites_expct_diso_prob <- phosphoDiso_expct_prob


human_data <- as.data.table(human_data)
human_data[,binom := purrr::pmap(.(psites_obsv_diso, psites_count, psites_expct_diso_prob), binom.test, alternative="greater")]
human_data[,binom_p := mapply("[[", binom, "p.value", SIMPLIFY = T)]
human_data[,binom_q := mapply(p.adjust, binom_p)]
human_data[,binom_sig := factor(ifelse(binom_q < 0.05, ifelse(binom_q < 0.01, "1% FDR", "5% FDR"), "n.s."), levels=c("n.s.","5% FDR", "1% FDR")) ]



ggplot(human_data) + 
  geom_point(aes(x=psites_obsv_diso,y=psites_expct_diso, colour = binom_sig,shape=target),size=2,alpha=0.80)+
  geom_abline(color="darkslategrey",slope = 1,size=0.5,linetype = "dashed")+
  ggpubr::theme_classic2() + 
  theme(text = element_text(size=15),legend.position = c(0.32,0.60),legend.box.just = "left",legend.box.margin = margin(2, 2, 2, 2),legend.box.background = element_rect(color="darkslategrey"),legend.title = element_text(size = 13)) +
  guides(color=guide_legend(title="Statistical significance")) +
  scale_x_continuous(limits = c(0, 60),breaks = c(seq(0, 60, by = 10)))+ xlab("Observed phospho S/T in IDR") +
  scale_y_continuous(limits = c(0, 60),breaks = c(seq(0, 60, by = 10)))+ ylab("Expected phospho S/T in IDR") + 
  scale_colour_manual(values = c(pal_jco()(10)[3],"#ffdd15ff",pal_jco()(10)[4]))+
  scale_shape_manual(values = c(16,4,4))
```

\#Fig. 3c

``` r
melted_all_predictions_phospho <- all_predictions_phospho %>% melt(id.vars=c("ID","anova_sig","hCDK1target","length"),value.name = "psites_diso",measure.vars=c("psites_expct_diso","psites_obsv_diso"))



# ALL
wilcox.test(psites_diso~variable,data=melted_all_predictions_phospho,paired=T,estimate =T,conf.int=T)
ggplot(melted_all_predictions_phospho) + 
  geom_boxplot(aes(y=psites_diso,x=variable, fill = variable,color = variable),outlier.shape = NA)+
  ggpubr::theme_classic2()  + 
  theme(text = element_text(size=20),legend.position = "none",axis.ticks.x = element_blank()) +
  geom_segment(aes(x = 1, y = 16.1, xend = 2, yend = 16.1)) + annotate(geom="text", x=1.5, y=16.3, label="***",size=10) +
  guides(color=guide_legend(title="Statistical significance")) +
  scale_x_discrete(labels = c("Expected","Observed")) + xlab(element_blank()) +
  scale_y_continuous(limits = c(0, 31),breaks = c(seq(0, 30, by = 2)),expand = c(0.05,0.05))+ ylab("Phospho S/T in IDR") + 
  scale_colour_manual(values = pal_jco()(10)[c(7,10)]) +
  scale_fill_manual(values = pal_jco()(10)[c(2,5)])

# ANOVA +
wilcox.test(psites_diso~variable,data=subset(melted_all_predictions_phospho,anova_sig == "Dynamic"),paired=T,estimate =T,conf.int=T)
ggplot(subset(melted_all_predictions_phospho,anova_sig == "Dynamic")) + 
  geom_boxplot(aes(y=psites_diso,x=variable, fill = variable,color = variable),outlier.shape = NA)+
  ggpubr::theme_classic2()  + 
  theme(text = element_text(size=20),legend.position = "none",axis.ticks.x = element_blank()) +
  geom_segment(aes(x = 1, y = 16.1, xend = 2, yend = 16.1)) + annotate(geom="text", x=1.5, y=16.3, label="***",size=10) +
  guides(color=guide_legend(title="Statistical significance")) +
  scale_x_discrete(labels = c("Expected","Observed")) + xlab(element_blank()) +
  scale_y_continuous(limits = c(0, 31),breaks = c(seq(0, 30, by = 2)),expand = c(0.05,0.05))+ ylab("Phospho S/T in IDR") + 
  scale_colour_manual(values = pal_jco()(10)[c(7,10)]) +
  scale_fill_manual(values = pal_jco()(10)[c(2,5)])

# HumanCDK targets in Dynamic
wilcox.test(psites_diso~variable,data=subset(melted_all_predictions_phospho,anova_sig == "Dynamic" & hCDK1target == "Human CDK1 target"),paired=T,estimate =T,conf.int=T)
ggplot(subset(melted_all_predictions_phospho,hCDK1target == "Human CDK1 target" & anova_sig == "Dynamic")) + 
  geom_boxplot(aes(y=psites_diso,x=variable, fill = variable,color = variable),outlier.shape = NA)+
  ggpubr::theme_classic2()  + 
  theme(text = element_text(size=20),legend.position = "none",axis.ticks.x = element_blank()) +
  geom_segment(aes(x = 1, y = 16.1, xend = 2, yend = 16.1)) + annotate(geom="text", x=1.5, y=16.3, label="***",size=10) +
  guides(color=guide_legend(title="Statistical significance")) +
  scale_x_discrete(labels = c("Expected","Observed")) + xlab(element_blank()) +
  scale_y_continuous(limits = c(0, 31),breaks = c(seq(0, 30, by = 2)),expand = c(0.05,0.05))+ ylab("Phospho S/T in IDR") + 
  scale_colour_manual(values = pal_jco()(10)[c(7,10)]) +
  scale_fill_manual(values = pal_jco()(10)[c(2,5)])
```

# Fig. 3d

![](index_V2_files/figure-gfm/Enrichment%20of%20Phospho%20Ser%20and%20Ther%20in%20IDRs%20-1.png)<!-- -->

# Fig. 3e

## Yeast datasets intersection

![](index_V2_files/figure-gfm/Holt%20disorder%20percentage-1.png)<!-- -->![](index_V2_files/figure-gfm/Holt%20disorder%20percentage-2.png)<!-- -->

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

## Human

![](index_V2_files/figure-gfm/human%20disorder%20percentage-1.png)<!-- -->

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Percentage of disorder by target
    ## W = 7707796, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  18.78995 22.82557
    ## sample estimates:
    ## difference in location 
    ##               20.85548

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Percentage of disorder by target
    ## W = 7269682, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  16.61683 21.13845
    ## sample estimates:
    ## difference in location 
    ##               18.86599

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Percentage of disorder by target
    ## W = 7349442, p-value < 2.2e-16
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  16.81769 21.05814
    ## sample estimates:
    ## difference in location 
    ##               18.94645

## Xenopus

![](index_V2_files/figure-gfm/xenopus%20disorder%20percentage-1.png)<!-- -->

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

\#Fig. 3f

## Violin plots for the disorder percentage in all kinases’ targets for the three main predictors used in this paper

![](index_V2_files/figure-gfm/disorder%20percentages%20for%20other%20kinases-1.png)<!-- -->

## Post-hoc test hypothesis for all possible comparisons

    ##            subset     n   mean     sd   min     Q1 median     Q3     max
    ## 1 Phosphoproteome 15332 22.775 24.207 0.000  3.333 13.688 36.166 100.000
    ## 2             CDK   656 43.040 25.881 0.000 24.184 40.389 61.776 100.000
    ## 3            MAPK   408 36.539 25.173 0.000 15.237 31.296 53.273  98.487
    ## 4            AURK   412 44.438 24.141 0.291 27.811 40.845 60.857 100.000
    ## 5             PLK   459 39.847 25.215 0.000 18.771 38.403 57.566 100.000
    ## 6             NEK    43 28.806 22.094 2.448 10.339 25.741 43.054  85.624
    ## 7            DYRK    42 41.869 25.545 0.239 22.507 40.408 59.029  95.763
    ##   percZero
    ## 1    9.483
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
    ## Kruskal-Wallis chi-squared = 1063.5, df = 6, p-value < 2.2e-16

![](index_V2_files/figure-gfm/post-hoc%20test%20for%20all%20kinases%20disorder%20percentage-1.png)<!-- -->

    ##            subset     n nvalid   mean     sd   min     Q1 median     Q3     max
    ## 1 Phosphoproteome 15332  15303 30.607 28.055 0.000  8.230 20.918 45.883 100.000
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
    ## Kruskal-Wallis chi-squared = 680.09, df = 6, p-value < 2.2e-16

![](index_V2_files/figure-gfm/post-hoc%20test%20for%20all%20kinases%20disorder%20percentage-2.png)<!-- -->

    ##            subset     n nvalid   mean     sd    min     Q1 median     Q3
    ## 1 Phosphoproteome 15332  15312 43.246 25.447  1.980 21.942 37.723 62.377
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
    ## Kruskal-Wallis chi-squared = 750.94, df = 6, p-value < 2.2e-16

![](index_V2_files/figure-gfm/post-hoc%20test%20for%20all%20kinases%20disorder%20percentage-3.png)<!-- -->

# Supplemental information

## Comparing predictors

This block compares the residues predicted as to be disordered across
the entire proteomes of yeast and human using different disorders. A
Venn diagram is generated to illustrate the overlap between different
methods

![](index_V2_files/figure-gfm/compare%20index%20-1.png)<!-- -->![](index_V2_files/figure-gfm/compare%20index%20-2.png)<!-- -->

## Supp. Fig. 5a

The compositional bias of IDRs is calculated with 3 different
predictores for each of the organisms considered in this paper.

### Yeast

![](index_V2_files/figure-gfm/yeast%20composition%20-1.png)<!-- -->

### Human

![](index_V2_files/figure-gfm/human%20composition%20-1.png)<!-- -->

### Xenopus

![](index_V2_files/figure-gfm/xenopus%20composition%20-1.png)<!-- -->

# Supp. Fig. 5d

![](index_V2_files/figure-gfm/Lollipop%20plot%20for%20all%20human%20kinases%20-1.png)<!-- -->

# Supp. Fig. 5e

Calculate the porcentage of disorder of human CDK targets vs the rest of
the phosphoproteome with all available predictors

![](index_V2_files/figure-gfm/human%20disorder%20percentage%20for%20all%20predictors%20-1.png)<!-- -->

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
