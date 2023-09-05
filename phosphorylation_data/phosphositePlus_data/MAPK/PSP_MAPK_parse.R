
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

phosphosite_ERK1_target_HS <- read_delim("PSP/MAPK/erk1.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_ERK2_target_HS <- read_delim("PSP/MAPK/erk2.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_ERK3_target_HS <- read_delim("PSP/MAPK/erk3.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_ERK5_target_HS <- read_delim("PSP/MAPK/erk5.tab","\t", escape_double = FALSE, trim_ws = TRUE)

phosphosite_P38A_target_HS <- read_delim("PSP/MAPK/p38a.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_P38B_target_HS <- read_delim("PSP/MAPK/p38b.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_P38D_target_HS <- read_delim("PSP/MAPK/p38d.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_P38G_target_HS <- read_delim("PSP/MAPK/p38g.tab","\t", escape_double = FALSE, trim_ws = TRUE)

phosphosite_JNK1_target_HS <- read_delim("PSP/MAPK/jnk1.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_JNK1iso2_target_HS <- read_delim("PSP/MAPK/jnk1_iso2.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_JNK2_target_HS <- read_delim("PSP/MAPK/jnk2.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_JNK2iso2_target_HS <- read_delim("PSP/MAPK/jnk2_iso2.tab","\t", escape_double = FALSE, trim_ws = TRUE)
phosphosite_JNK3_target_HS <- read_delim("PSP/MAPK/jnk3.tab","\t", escape_double = FALSE, trim_ws = TRUE)

DF_list <- list(phosphosite_ERK1_target_HS,
                phosphosite_ERK2_target_HS,
                phosphosite_ERK3_target_HS,
                phosphosite_ERK5_target_HS,
                phosphosite_P38A_target_HS,
                phosphosite_P38B_target_HS,
                phosphosite_P38D_target_HS,
                phosphosite_P38G_target_HS,
                phosphosite_JNK1_target_HS,
                phosphosite_JNK1iso2_target_HS,
                phosphosite_JNK2_target_HS,
                phosphosite_JNK2iso2_target_HS,
                phosphosite_JNK3_target_HS
                )

phosphosite_MAPKFamilyTarget_HS <- Reduce(
  function(x, y, ...) merge(x, y, by = c("GENE","PROTEIN","ACC#","MOD_RSD"),all = T, ...),
  DF_list
)



mapk_data <- phosphosite_MAPKFamilyTarget_HS %>% unite(col = "KINASE",
                                                       ERK1.KINASE,
                                                       ERK2.KINASE,
                                                       ERK3.KINASE,
                                                       ERK5.KINASE,
                                                       P38A.KINASE,
                                                       P38B.KINASE,
                                                       P38D.KINASE,
                                                       P38G.KINASE,
                                                       JNK1.KINASE,
                                                       JNK1_iso2.KINASE,
                                                       JNK2.KINASE,
                                                       JNK2_iso2.KINASE,
                                                       JNK3.KINASE,
                                                       sep = ",",remove = T,na.rm=T )


mapk_data <- mapk_data %>% group_by(`ACC#`,GENE,PROTEIN) %>% summarise_at(c("MOD_RSD","KINASE"),function(x){paste(x, collapse=",")})

mapk_data$MAPK_MOD_RSD <- sapply(mapk_data$MOD_RSD, function(x){ return(paste(substr(strsplit(x,",")[[1]],start = 2,stop = 2000),collapse = ","))})
mapk_data$MAPK_KINASE <- sapply(mapk_data$KINASE, function(x){ return(paste(unique(strsplit(x,",")[[1]]),collapse = ","))})

mapk_data$MOD_RSD <- NULL
mapk_data$KINASE <- NULL

write.table(mapk_data,file = "PSP/MAPK/PSP_MAPK_Target_HS.tab",quote = F,sep = "\t",row.names = F)
# algunos valores tienen espacios, por eso el gsub
# human_data$method <- sapply(human_data$method, function(x){ return(paste(unique(gsub(" ", "", strsplit(x,",")[[1]], fixed = TRUE)),collapse = ","))})



# The data is exported and mannually curated afterwards. Psites from bibliography were added. All the information is then deposited in the file human_data_curated.tab