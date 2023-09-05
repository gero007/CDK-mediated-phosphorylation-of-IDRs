########################################
############## Mobi-DB  ################
########################################

getMobiDB_PredDiso <- function(mobiDB_df,accession_column,sequence_col){
  
  output_data_frame = mobiDB_df
  # Generate a list with one  entry per proteins containing one named matrix with the disordered regions
  
  
  disorder <- mobiDB_df$disorder$predictors
  names(disorder) <- mobiDB_df[,accession_column]
  predictors_list <- lapply(disorder, function(x){
    m <- x$regions
    names(m) <- x$method
    return(m)
  })
  # explore all the available predictors
  available_predictors <- unique(unlist(sapply(predictors_list, names)))
  
  # For each predictor I create three columns for storing the results
  for (predictor in available_predictors) {
    output_data_frame[,paste(predictor,"perc",sep =  "_")] <- rep(NA,nrow(output_data_frame))
    output_data_frame[,paste(predictor,"disoIndexes",sep = "_")] <- rep(NA,nrow(output_data_frame))
    output_data_frame[,paste(predictor,"disoSeq",sep = "_")] <- rep(NA,nrow(output_data_frame))
  }
  # Make the operations and store the results in the data frame
  # predictors list and output_data_frame are sorted in the same way, i'll iterate using an incrementator for simplicity
  for (i in 1:length(predictors_list)) {
    # print(i)
    for (n in 1:length(predictors_list[[i]])) {
      pred_name <- names(predictors_list[[i]][n])
      mat <- predictors_list[[i]][[n]]
      # print(mat)
      # testing if there are disorder regions
      if (class(mat)[1]=="matrix") {
        
        #Percentage of disorder
        nDiso<-sum(as.numeric(mat[,2])-as.numeric(mat[,1])+1)
        proteinLength <- output_data_frame$Length[[i]]
        # print(paste(pred_name,"perc",sep = "_"))
        # print(output_data_frame[i,paste(pred_name,"perc",sep = "_")])
        # print((nDiso/proteinLength)*100)
        output_data_frame[i,paste(pred_name,"perc",sep = "_")] <- (nDiso/proteinLength)*100
        
        #Disorder indexes
        disoIdx <-numeric()
        disoIdx <- apply(mat, 1, function(x){disoIdx<-c(disoIdx,x[1]:x[2])
        return(disoIdx)})
        disoIdx <- unlist(disoIdx)
        # print(paste(disoIdx,collapse = ","))
        output_data_frame[i,paste(pred_name,"disoIndexes",sep = "_")] <- paste(disoIdx,collapse = ",")
        
        #Disorder sequence
        disoSeq <- character()
        for (row in 1:nrow(mat)) {
          startIdx <- as.numeric(mat[row,1])
          endIdx <- as.numeric(mat[row,2])
          
          disoSeq <- paste(disoSeq,substr(output_data_frame[i,sequence_col],startIdx,endIdx),sep="")
        }
        output_data_frame[i,paste(pred_name,"disoSeq",sep = "_")] <- disoSeq
      } else {
        output_data_frame[i,paste(pred_name,"perc",sep = "_")] <- 0
        output_data_frame[i,paste(pred_name,"disoIndexes",sep = "_")] <- "0"
        output_data_frame[i,paste(pred_name,"disoSeq",sep = "_")] <- ""
      }
    }
  }
  
  return(output_data_frame)
}

########################################
############# SPOT Disorder ############
########################################


getSpotDisoList <- function(disoPath){ 
  output_list <- list()
  for (disoFile in dir(path = disoPath,pattern = "*.spotd")) {
    id=strsplit(disoFile,split = "\\.")[[1]][1]
    # print(id)
    # print(paste(disoPath,disoFile,sep = ""))
    auxdf <- read.delim(paste(disoPath,disoFile,sep = ""), header=FALSE, comment.char="#")
    disoRunningLenghts <- rle(as.vector(auxdf$V4))
    start=1
    disoMatrix <- matrix(ncol = 4,nrow = length(disoRunningLenghts$lengths))
    for (i in 1:length(disoRunningLenghts$lengths)) {
      end=start+disoRunningLenghts$lengths[i]-1
      disoMatrix[i,] <- as.character(c(start,end,disoRunningLenghts$values[i],disoRunningLenghts$lengths[i]))
      # print(c(start,end,test$values[i],test$lengths[i]))
      start=end+1
    }
    output_list[[id]] <- disoMatrix
  }
  # print(disoMatrix)
  # print(which(targetDf$acc==id))
  return(output_list)
}


getSpotPredDiso <- function(df,disoPath,accession_col,length_col=NULL,sequence_col){
  spotDisorderList <- getSpotDisoList(disoPath)
  output_data_frame = df
  
  perc <- numeric()
  disoIdxList <- c()
  disoSeqList <- c()
  for (i in 1:nrow(df)) {
    mat<-spotDisorderList[[df[i,accession_col]]]
    if (is.null(mat)){
      nDiso<-NA
      perc[i] <- nDiso
      disoIdx <- NA
      disoIdxList[i]<-disoIdx
      disoSeq <- NA
      disoSeqList[i]<-disoSeq
    } else {
      # print(i)
      # print(accession_col)
      #select only disordered regions
      mat <- matrix(mat[mat[,3]=="D",],ncol=4)
      # print(df$acc[[i]])
      # print(mat)
      if (nrow(mat)!=0) {
        #calculate the number of disoredered AAs
        nDiso<-sum(as.numeric(mat[,4]))
        # print(nDiso)
        if (!is.null(length_col)){
          proteinLength <- as.numeric(df[i,length_col])
        } else {proteinLength <- nchar(df[i,sequence_col])}
        # print(proteinLength)
        perc[i] <- (nDiso/proteinLength)*100
        # print((nDiso/proteinLength)*100)
        disoIdx <-c()
        disoIdx <-apply(mat, 1, function(x){disoIdx<-c(disoIdx,x[1]:x[2])
        return(disoIdx)})
        disoIdx <- paste(unlist(disoIdx),collapse = ",")
        disoIdxList[i]<-disoIdx
        disoSeq <- character()
        for (row in 1:nrow(mat)) {
          startIdx <- as.numeric(mat[row,1])
          endIdx <- as.numeric(mat[row,2])
          disoSeq <- paste(disoSeq,substr(df[i,sequence_col],startIdx,endIdx),sep = "")
        }
        # print(disoSeq)
        disoSeqList[i]<-disoSeq
        
      } else {
        nDiso<-0
        perc[i] <- nDiso
        disoIdx <- "0"
        disoIdxList[i]<-disoIdx
        disoSeq <- ""
        disoSeqList[i]<-disoSeq
      }}
  }
  output_data_frame$spot_perc <- as.numeric(perc)
  output_data_frame$spot_disoIndexes <- disoIdxList
  output_data_frame$spot_disoSeq <- disoSeqList
  return(output_data_frame)
}

########################################
################ VSL2B ##################
########################################


getVSLDisoMatrix <- function(disoPath){ 
  output_list <- list()
  for (disoFile in dir(path = disoPath,pattern = "*.vsl2b")) {
    id=strsplit(disoFile,split = "\\.")[[1]][1]
    # print(id)
    # print(paste(disoPath,disoFile,sep = ""))
    lines <- readLines(paste(disoPath,disoFile,sep = ""))
    lines <- lines[grepl(pattern = "[0-9]-[0-9]",lines)]
    indexes <- as.numeric(unlist(strsplit(lines,"-")))
    disoMatrix <- matrix(indexes, ncol = 2, byrow = T)
    # disoMatrix <- cbind(disoMatrix,rep("D",nrow(disoMatrix))) 
    output_list[[id]] <- disoMatrix
  }
  return(output_list)
}


getVSLPredDiso <- function(df,disoPath,accession_col,length_col=NULL,sequence_col){
  vslDisorderList <- getVSLDisoMatrix(disoPath)
  output_data_frame = df
  
  perc <- numeric()
  disoIdxList <- c()
  disoSeqList <- c()
  for (i in 1:nrow(df)) {
    mat<-vslDisorderList[[df[i,accession_col]]]
    if (is.null(mat)){
      nDiso<-NA
      perc[i] <- nDiso
      disoIdx <- NA
      disoIdxList[i]<-disoIdx
      disoSeq <- NA
      disoSeqList[i]<-disoSeq
    } else {
      # print(i)
      # print(accession_col)
      #select only disordered regions
      # print(df$acc[[i]])
      # print(mat)
      if (nrow(mat)!=0) {
        #calculate the number of disoredered AAs
        nDiso<-sum(mat[,2]-mat[,1]+1)
        # print(nDiso)
        if (!is.null(length_col)){
          proteinLength <- as.numeric(df[i,length_col])
        } else {proteinLength <- nchar(df[i,sequence_col])}
        # print(proteinLength)
        perc[i] <- (nDiso/proteinLength)*100
        # print((nDiso/proteinLength)*100)
        disoIdx <-c()
        disoIdx <-apply(mat, 1, function(x){disoIdx<-c(disoIdx,x[1]:x[2])
        return(disoIdx)})
        disoIdx <- paste(unlist(disoIdx),collapse = ",")
        disoIdxList[i]<-disoIdx
        disoSeq <- character()
        for (row in 1:nrow(mat)) {
          startIdx <- as.numeric(mat[row,1])
          endIdx <- as.numeric(mat[row,2])
          disoSeq <- paste(disoSeq,substr(df[i,sequence_col],startIdx,endIdx),sep = "")
        }
        # print(disoSeq)
        disoSeqList[i]<-disoSeq
        
      } else {
        nDiso<-0
        perc[i] <- nDiso
        disoIdx <- "0"
        disoIdxList[i]<-disoIdx
        disoSeq <- ""
        disoSeqList[i]<-disoSeq
      }}
  }
  output_data_frame$vsl_perc <- as.numeric(perc)
  output_data_frame$vsl_disoIndexes <- disoIdxList
  output_data_frame$vsl_disoSeq <- disoSeqList
  return(output_data_frame)
}

########################################
################ IUpred ##################
########################################


getIUpredDisoList <- function(disoPath){ 
  output_list <- list()
  for (disoFile in dir(path = disoPath,pattern = "*.iupred")) {
    id=strsplit(disoFile,split = "\\.")[[1]][1]
    # print(id)
    # print(paste(disoPath,disoFile,sep = ""))
    auxdf <- read.delim(paste(disoPath,disoFile,sep = ""), header=FALSE, comment.char="#")
    auxdf$V4 <- cut(auxdf$V3,breaks = c(-0.1,0.5,1),labels = c("O","D"))
    disoRunningLenghts <- rle(as.vector(auxdf$V4))
    start=1
    disoMatrix <- matrix(ncol = 4,nrow = length(disoRunningLenghts$lengths))
    for (i in 1:length(disoRunningLenghts$lengths)) {
      end=start+disoRunningLenghts$lengths[i]-1
      disoMatrix[i,] <- as.character(c(start,end,disoRunningLenghts$values[i],disoRunningLenghts$lengths[i]))
      # print(c(start,end,test$values[i],test$lengths[i]))
      start=end+1
    }
    output_list[[id]] <- disoMatrix
  }
  # print(disoMatrix)
  # print(which(targetDf$acc==id))
  return(output_list)
}



getIUpredPredDiso <- function(df,disoPath,accession_col,length_col=NULL,sequence_col){
  IUpredDisorderList <- getIUpredDisoList(disoPath)
  output_data_frame = df
  
  perc <- numeric()
  disoIdxList <- c()
  disoSeqList <- c()
  for (i in 1:nrow(df)) {
    mat<-IUpredDisorderList[[df[i,accession_col]]]
    if (is.null(mat)){
      nDiso<-NA
      perc[i] <- nDiso
      disoIdx <- NA
      disoIdxList[i]<-disoIdx
      disoSeq <- NA
      disoSeqList[i]<-disoSeq
    } else {
      # print(i)
      # print(accession_col)
      #select only disordered regions
      mat <- matrix(mat[mat[,3]=="D",],ncol=4)
      # print(df$acc[[i]])
      # print(mat)
      if (nrow(mat)!=0) {
        #calculate the number of disoredered AAs
        nDiso<-sum(as.numeric(mat[,4]))
        # print(nDiso)
        if (!is.null(length_col)){
          proteinLength <- as.numeric(df[i,length_col])
        } else {proteinLength <- nchar(df[i,sequence_col])}
        # print(proteinLength)
        perc[i] <- (nDiso/proteinLength)*100
        # print((nDiso/proteinLength)*100)
        disoIdx <-c()
        disoIdx <-apply(mat, 1, function(x){disoIdx<-c(disoIdx,x[1]:x[2])
        return(disoIdx)})
        disoIdx <- paste(unlist(disoIdx),collapse = ",")
        disoIdxList[i]<-disoIdx
        disoSeq <- character()
        for (row in 1:nrow(mat)) {
          startIdx <- as.numeric(mat[row,1])
          endIdx <- as.numeric(mat[row,2])
          disoSeq <- paste(disoSeq,substr(df[i,sequence_col],startIdx,endIdx),sep = "")
        }
        # print(disoSeq)
        disoSeqList[i]<-disoSeq
        
      } else {
        nDiso<-0
        perc[i] <- nDiso
        disoIdx <- "0"
        disoIdxList[i]<-disoIdx
        disoSeq <- ""
        disoSeqList[i]<-disoSeq
      }}
  }
  output_data_frame$iupl_perc <- as.numeric(perc)
  output_data_frame$iupl_disoIndexes <- disoIdxList
  output_data_frame$iupl_disoSeq <- disoSeqList
  return(output_data_frame)
}