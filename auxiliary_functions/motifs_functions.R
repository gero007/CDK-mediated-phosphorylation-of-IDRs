
#Function to calculate the information content
informationContent_calculator <- function(sequences,backgroundProb){
  if (var(nchar(sequences))!=0){
    stop("Sequences with different lenghts cannot be aligned")
  } else {
    
    windowsLength <- unique(nchar(sequences))
    sequencesCount <- length(sequences)
    aux_matrix <- matrix(as.factor(unlist(strsplit(sequences,""))),byrow = T,nrow=sequencesCount,ncol = windowsLength)
    # not working now. Should return a count sumarrry of each aminoacide
    out_matrix <- matrix(nrow = 20, ncol = windowsLength)
    rownames(out_matrix) <- AA_ALPHABET[1:20]
    
    for (i in 1:ncol(aux_matrix)) {
      
      colCounts <- summary(as.factor(aux_matrix[,i]))
      # colCounts <- str_split_fixed(aux_matrix[,i],":( )*",2)
      # aux_names <- colCounts[,1]
      # colCounts <- as.numeric(trimws(colCounts[,2]))
      # names(colCounts) <- aux_names
      colCounts <- colCounts[(names(colCounts)!="_" & names(colCounts)!="")]
      ObservedProb <- colCounts/sum(colCounts)
      for (aminoacid in row.names(out_matrix)) {
        
        out_matrix[aminoacid,i] <- ObservedProb[aminoacid]*log2(ObservedProb[aminoacid]/backgroundProb[aminoacid])
        
      }
      
    }
  }
  out_matrix[is.na(out_matrix)] <- 0
  norm_out_matrix <- apply(out_matrix, 2, function(x){x/sum(abs(x))})
  return(norm_out_matrix)
}


#=====================| Count without overlapps |=====================

countMotifs <- function(motif_list, sequences){
  plk <- grepl(motif_list[["Plk"]],sequences)
  aurora <- grepl(motif_list[["Aurora A/B"]],sequences)
  ck <- grepl(motif_list[["Ck1"]],sequences) | grepl(motif_list[["Ck2"]],sequences)
  ddk <- grepl(motif_list[["Cdc7"]],sequences)
  cdk_min <- grepl(motif_list[["Cdk minimal"]],sequences)
  cdk_full <- grepl(motif_list[["Cdk full"]],sequences)
  total <- grepl(motif_list[["Total"]],sequences)
  
  output <- list(plk,aurora,ck,ddk,cdk_min,cdk_full,total)
  names(output) <- c("plk","aurora","ck","ddk","cdk_min","cdk_full","total")
  return(output)
}