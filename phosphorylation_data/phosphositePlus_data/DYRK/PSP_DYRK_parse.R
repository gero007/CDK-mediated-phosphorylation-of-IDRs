library(readr)
library(dplyr)
library(tidyr)
# PSP Search Result Legend	 
# COLUMN HEADER	DESCRIPTION
# PROTEIN	Primary name in PhosphoSitePlusÂ® (PSP)
# GENE	Unique ID (UID) assigned to a gene by the HUGO Nomenclature Committee and adapted by HUPO for the cognate protein.
# ORGANISM	The species from which either the kinase and/or substrate proteins have been derived
# ACC_ID	Primary accession ID in PSP
# MW_(DA)	Molecular weight in Daltons.
# MOD_RSD	The location within the protein sequence of the amino acid residue that is posttranslationally modified (modsite).
# SITE_+/-7_AA	The modsite plus flanking sequence (+/- 7 rsds.). The modsite, as well as other residues within the flanking sequence that are known to be posttranslationally modified, are displayed in lower-case letters.
# SITE_GRP_ID	Unique identifier of the modification site and its homologous sites in all proteoforms and species
# DOMAIN	The Pfam residue in which the modsite is located.
# METHOD	The experimental method in which kinase-substrate relationships were reported. Vivo = in vivo: determined from reactions within living cells, cell cultures or organisms. Vitro = in vitro: determined from reactions outside of living cellular structures.
# LTP_LIT	The number of literature records using low-throughput (LTP) experimental techniques in which the modsite has been reported. LTP results may be more reliable than MS2 results.
# MS2_LIT	The number of literature records using tandem mass spectrometry (MS2) experimental techniques in which the modsite has been observed.
# MS2_CST	The number of curation sets (CS) derived from MS2 experiments performed at Cell Signaling Technology (CST) in which the modsite has been observed.
# CST_CAT#	The catalog number(s) of CST antibodies specific for the associated modsite.

#there's overlapping between the targets. I create 13 differents DF and I will merge them

phosphosite_DYRK1A_target_HS <- read_delim("PSP/DYRK/dyrk1a.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_DYRK1B_target_HS <- read_delim("PSP/DYRK/dyrk1b.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_DYRK2_target_HS <- read_delim("PSP/DYRK/dyrk2.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_DYRK3_target_HS <- read_delim("PSP/DYRK/dyrk3.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_DYRK4_target_HS <- read_delim("PSP/DYRK/dyrk4.tab","\t", escape_double = FALSE, trim_ws = TRUE)




DF_list <- list(phosphosite_DYRK1A_target_HS,
                phosphosite_DYRK1B_target_HS,
                phosphosite_DYRK2_target_HS,
                phosphosite_DYRK3_target_HS,
                phosphosite_DYRK4_target_HS
                )

phosphosite_DYRKFamilyTarget_HS <- Reduce(
  function(x, y, ...) merge(x, y, by = c("GENE","PROTEIN","ACC#","MOD_RSD"),all = T, ...),
  DF_list
)



dyrk_data <- phosphosite_DYRKFamilyTarget_HS %>% unite(col = "KINASE",
                                                       DYRK1A.KINASE,
                                                       DYRK1B.KINASE,
                                                       DYRK2.KINASE,
                                                       DYRK3.KINASE,
                                                       DYRK4.KINASE,
                                                       sep = ",",remove = T,na.rm=T )


dyrk_data <- dyrk_data %>% group_by(`ACC#`,GENE,PROTEIN) %>% summarise_at(c("MOD_RSD","KINASE"),function(x){paste(x, collapse=",")})

dyrk_data$DYRK_MOD_RSD <- sapply(dyrk_data$MOD_RSD, function(x){ return(paste(substr(strsplit(x,",")[[1]],start = 2,stop = 2000),collapse = ","))})
dyrk_data$DYRK_KINASE <- sapply(dyrk_data$KINASE, function(x){ return(paste(unique(strsplit(x,",")[[1]]),collapse = ","))})

dyrk_data$MOD_RSD <- NULL
dyrk_data$KINASE <- NULL
# It's simpler to eliminate the gene and protein names for the merge. slight changes in the naming make the marge inconsistent and duplicate the entries
dyrk_data$GENE <- NULL
dyrk_data$PROTEIN <- NULL


write.table(dyrk_data,file = "PSP/DYRK/PSP_DYRK_Target_HS.tab",quote = F,sep = "\t",row.names = F)
