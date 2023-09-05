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

phosphosite_PLK1_target_HS <- read_delim("PSP/PLK/plk1.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_PLK2_target_HS <- read_delim("PSP/PLK/plk2.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_PLK3_target_HS <- read_delim("PSP/PLK/plk3.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_PLK4_target_HS <- read_delim("PSP/PLK/plk4.tab","\t", escape_double = FALSE, trim_ws = TRUE)




DF_list <- list(phosphosite_PLK1_target_HS,
                phosphosite_PLK2_target_HS,
                phosphosite_PLK3_target_HS,
                phosphosite_PLK4_target_HS
                )

phosphosite_PLKFamilyTarget_HS <- Reduce(
  function(x, y, ...) merge(x, y, by = c("GENE","PROTEIN","ACC#","MOD_RSD"),all = T, ...),
  DF_list
)



plk_data <- phosphosite_PLKFamilyTarget_HS %>% unite(col = "KINASE",
                                                       PLK1.KINASE,
                                                       PLK2.KINASE,
                                                       PLK3.KINASE,
                                                       PLK4.KINASE,
                                                       sep = ",",remove = T,na.rm=T )


plk_data <- plk_data %>% group_by(`ACC#`,GENE,PROTEIN) %>% summarise_at(c("MOD_RSD","KINASE"),function(x){paste(x, collapse=",")})

plk_data$PLK_MOD_RSD <- sapply(plk_data$MOD_RSD, function(x){ return(paste(substr(strsplit(x,",")[[1]],start = 2,stop = 2000),collapse = ","))})
plk_data$PLK_KINASE <- sapply(plk_data$KINASE, function(x){ return(paste(unique(strsplit(x,",")[[1]]),collapse = ","))})

plk_data$MOD_RSD <- NULL
plk_data$KINASE <- NULL

#Add targets from Kettenbach et al. form Gerber's lab

# gerber_list <- read_lines("PSP/PLK/gerber_plk_list")
# 
# gerber_matrix <- str_split_fixed(gerber_list,"_\\(",2)
# 
# gerber_matrix[,2] <- trimws(gerber_matrix[,2],which = "right",whitespace = "\\)")
# gerber_data <- as.data.frame(gerber_matrix) %>% group_by(V1) %>% summarise_at(c("V2"),function(x){paste(x, collapse=":")})
# 
# colnames(gerber_data)<-c("ACC#","gerber_sites")
# gerber_data$gerber_sites<-str_replace_all(gerber_data$gerber_sites, ",", "-")
# gerber_data$gerber_sites<-str_replace_all(gerber_data$gerber_sites, ":", ";")
# 
# gerber_data<-merge.data.frame(gerber_data,plk_data,all = T,by = "ACC#")
# gerber_data$PLK_MOD_RSD<-str_replace_all(gerber_data$PLK_MOD_RSD, ",", ";")
# 
# write.table(gerber_data,file = "PSP/PLK/gerber_data.tab",quote = F,sep = "\t",row.names = F)

write.table(plk_data,file = "PSP/PLK//PSP_PLK_Target_HS.tab",quote = F,sep = "\t",row.names = F)
# algunos valores tienen espacios, por eso el gsub
# human_data$method <- sapply(human_data$method, function(x){ return(paste(unique(gsub(" ", "", strsplit(x,",")[[1]], fixed = TRUE)),collapse = ","))})



# The data is exported and mannually curated afterwards. Psites from bibliography were added. All the information is then deposited in the file human_data_curated.tab