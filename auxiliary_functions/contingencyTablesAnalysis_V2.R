
# Calculate the contingency table for S and T diso vs struct and Phospho vs Non Phospho

getStratContingencyArray <- function(df,sequence_col,diso_index_col,psites_col){
  
  df <- as.data.frame(df)
  indexST <- gregexpr("(S|T)",df[[sequence_col]])
  
  stratContingencyArray <- array(dim = c(2,2,nrow(df)))
  for (i in 1:nrow(df)) {
    
    
    indexDiso <- as.numeric(base::strsplit(df[[diso_index_col]][i],",")[[1]])
    indexPhosphoST <- df[[psites_col]][[i]]
    indexNonPhosphoST <- setdiff(as.numeric(indexST[[i]]),indexPhosphoST)
    
    ndp <- length(which(indexPhosphoST %in% indexDiso))
    nop <- length(which(!indexPhosphoST %in% indexDiso))
    ndnp <- length(which(indexNonPhosphoST %in% indexDiso))
    nonp<- length(which(!indexNonPhosphoST %in% indexDiso))
    # stratContingencyArray[1,1,i]<-ndp
    # stratContingencyArray[1,2,i]<-nop
    # stratContingencyArray[2,1,i]<-ndnp
    # stratContingencyArray[2,2,i]<-nonp
    stratContingencyArray[1,1,i]<-ndp
    stratContingencyArray[1,2,i]<-ndnp
    stratContingencyArray[2,1,i]<-nop
    stratContingencyArray[2,2,i]<-nonp
    # indexProt <- 1:nchar(df[i,sequence_col])
    # indexDiso <- disorderedRegions[[i]]
    # indexOrd <- setdiff(indexProt,indexDiso)
    # indexPhospho <- as.numeric(TPphosphoSites[[i]])
    # indexNonPhospho <- setdiff(indexProt,indexPhospho)
    # 
    # ndp <- length(which(indexPhospho %in% indexDiso))
    # nop <- length(which(indexPhospho %in% indexOrd))
    # ndnp <- length(which(indexNonPhospho %in% indexDiso))
    # nonp<- length(which(indexNonPhospho %in% indexOrd))
    # stratContingencyArray[1,1,i]<-ndp
    # stratContingencyArray[1,2,i]<-nop
    # stratContingencyArray[2,1,i]<-ndnp
    # stratContingencyArray[2,2,i]<-nonp
    
  }
return(stratContingencyArray)
}



# Calculate the contingency table for Psites diso vs struct and CDK vs other Kinases

getStratContingencyArrayPsites <- function(df,sequence_col,diso_index_col,cdk1_psites_col,psites_col,min_num_psites=2){
  
  # indexST <- gregexpr("(S|T)",df[,sequence_col])
  
  stratContingencyArray <- array(dim = c(2,2,nrow(df)))
  for (i in 1:nrow(df)) {
    
    
    indexDiso <- as.numeric(strsplit(df[i,diso_index_col],",")[[1]])
    indexCDKPhospho <- df[i,cdk1_psites_col][[1]]
    indexAllPhospho <- df[i,psites_col][[1]]
    indexOthersPhospho <- setdiff(indexAllPhospho,indexCDKPhospho)
    
    nCdkDiso <- sum(indexCDKPhospho %in% indexDiso)
    nCdkStruct <- sum(!indexCDKPhospho %in% indexDiso)
    nOtherDiso <- sum(indexOthersPhospho %in% indexDiso)
    nOtherStruct<- sum(!indexOthersPhospho %in% indexDiso)
    if(length(indexAllPhospho)>=min_num_psites){
      stratContingencyArray[1,1,i]<-nCdkDiso
      stratContingencyArray[1,2,i]<-nCdkStruct
      stratContingencyArray[2,1,i]<-nOtherDiso
      stratContingencyArray[2,2,i]<-nOtherStruct
    } else {stratContingencyArray[,,i] <- rep(0,4) }

    
  }
  filteredStratContingencyArray <- numeric()
  for(contArray in stratContingencyArray){
    if(!is.null(contArray)){
      filteredStratContingencyArray<-c(filteredStratContingencyArray,contArray)}
  }
  filteredStratContingencyArray <- array(filteredStratContingencyArray,dim = c(2,2,length(filteredStratContingencyArray)/4))
  return(filteredStratContingencyArray)
}





