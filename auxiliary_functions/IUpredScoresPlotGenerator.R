IUpredScoresPlotGenerator <- function(dataframe,id_col="ID",sites_col="psites",subset_sites_col="anova_psites",sequence_col="sequence",highlightSP=FALSE){
  plotList <- list()
  for(i in 1:nrow(dataframe)){
    dataframe <-  as.data.frame(dataframe)
    plotTittle <- dataframe[i,id_col]
    aux_df <-  as.data.frame(cbind(dataframe[i,"positions"][[1]],dataframe[i,"IUPredScores"][[1]]))
    names(aux_df) <- c("positions","IUPredScores")
    breaksXaxis <- c(seq(0, nrow(aux_df), by = as.integer(nrow(aux_df)/10)))
    if (nrow(aux_df) < 100) {
      breaksXaxis <- c(seq(0, nrow(aux_df), by = 20))
    } else if (nrow(aux_df) > 100  & nrow(aux_df) <= 500) {
      breaksXaxis <- c(seq(0, nrow(aux_df), by = 50))
    } else if (nrow(aux_df) > 500  & nrow(aux_df) <= 1000) {
      breaksXaxis <- c(seq(0, nrow(aux_df), by = 100))
    } else if (nrow(aux_df) > 1000  & nrow(aux_df) <= 2000) {
      breaksXaxis <- c(seq(0, nrow(aux_df), by = 200))
    } else {
      breaksXaxis <- c(seq(0, nrow(aux_df), by = 500))
    }
    
    if(highlightSP) {
      
      STP_sites <- as.numeric(gregexpr("[S|T]P", as.data.frame(dataframe)[i,sequence_col])[[1]])
      if(STP_sites[1]==-1){STP_sites <- numeric(0)}
      
      plotList[[i]] <- ggplot(aux_df,aes(x=positions,y=IUPredScores)) +
        geom_col(aes(fill=IUPredScores >= 0.5 ,color=IUPredScores >= 0.5)) +
        ggpubr::theme_classic2() + 
        theme(text = element_text(size=19),legend.position = "none") +
        scale_x_continuous(limits = c(0, nrow(aux_df)+1),breaks = breaksXaxis,expand = c(0.005,0.005))+ xlab("Positions") +
        scale_y_continuous(limits = c(-0.1, 1),breaks = c(seq(0,1,by=0.25)),expand = c(0.05,0.05)) + ylab("Score") +
        scale_fill_manual(values = c(pal_jco()(10)[3],"#f6921eff")) +
        scale_color_manual(values = c(pal_jco()(10)[3],"#f6921eff")) +
        annotate("point",x = dataframe$ST_residues[[i]],y=rep(0,length(dataframe$ST_residues[[i]])),shape="|",color=pal_jco()(10)[1],fill=pal_jco()(10)[1],size=8) +
        annotate("point",x = STP_sites,y=rep(0,length(STP_sites)),shape="|",color="deeppink4",fill="deeppink4",size=8) +
        annotate("point",x = dataframe[i,sites_col][[1]],y=rep(-0.08,length(dataframe[i,sites_col][[1]])),shape=21,color=pal_jco()(10)[8],fill="#ffdd15ff",size=6) +
        annotate("point",x = dataframe[i,subset_sites_col][[1]],y=rep(-0.08,length(dataframe[i,subset_sites_col][[1]])),shape=21,color=pal_jco()(10)[8],fill="olivedrab3",size=6) +
        geom_hline(yintercept = 0.5,size=0.5,linetype = "dashed",color = "darkslategrey") +
        ggtitle(plotTittle)
    } else {
    
      plotList[[i]] <- ggplot(aux_df,aes(x=positions,y=IUPredScores)) +
        geom_col(aes(fill=IUPredScores >= 0.5 ,color=IUPredScores >= 0.5)) +
        ggpubr::theme_classic2() + 
        theme(text = element_text(size=19),legend.position = "none") +
        scale_x_continuous(limits = c(0, nrow(aux_df)+1),breaks = breaksXaxis,expand = c(0.005,0.005))+ xlab("Positions") +
        scale_y_continuous(limits = c(-0.1, 1),breaks = c(seq(0,1,by=0.25)),expand = c(0.05,0.05)) + ylab("Score") +
        scale_fill_manual(values = c(pal_jco()(10)[3],"#f6921eff")) +
        scale_color_manual(values = c(pal_jco()(10)[3],"#f6921eff")) +
        annotate("point",x = dataframe$ST_residues[[i]],y=rep(0,length(dataframe$ST_residues[[i]])),shape="|",color=pal_jco()(10)[1],fill=pal_jco()(10)[1],size=8) +
        annotate("point",x = dataframe[i,sites_col][[1]],y=rep(-0.08,length(dataframe[i,sites_col][[1]])),shape=21,color=pal_jco()(10)[8],fill="#ffdd15ff",size=6) +
        annotate("point",x = dataframe[i,subset_sites_col][[1]],y=rep(-0.08,length(dataframe[i,subset_sites_col][[1]])),shape=21,color=pal_jco()(10)[8],fill="olivedrab3",size=6) +
        geom_hline(yintercept = 0.5,size=0.5,linetype = "dashed",color = "darkslategrey") +
        ggtitle(plotTittle)
    }
    

  }
  names(plotList)<-dataframe[,id_col]
  return(plotList)
}



#_________________________________________________TEST FOR DIFFERENT KINASES__________________________________________

IUpredScoresPlotGenerator_AllKinases <- function(dataframe,id_col="ID",sites_CDK="psite_CDK1",sites_MAPK="psite_MAPK",sites_AURK="psite_aurk",sequence_col="sequence"){
  plotList <- list()
  for(i in 1:nrow(dataframe)){
    dataframe <-  as.data.frame(dataframe)
    plotTittle <- dataframe[i,id_col]
    aux_df <-  as.data.frame(cbind(dataframe[i,"positions"][[1]],dataframe[i,"IUPredScores"][[1]]))
    names(aux_df) <- c("positions","IUPredScores")
    breaksXaxis <- c(seq(0, nrow(aux_df), by = as.integer(nrow(aux_df)/10)))
    if (nrow(aux_df) < 100) {
      breaksXaxis <- c(seq(0, nrow(aux_df), by = 20))
    } else if (nrow(aux_df) > 100  & nrow(aux_df) <= 500) {
      breaksXaxis <- c(seq(0, nrow(aux_df), by = 50))
    } else if (nrow(aux_df) > 500  & nrow(aux_df) <= 1000) {
      breaksXaxis <- c(seq(0, nrow(aux_df), by = 100))
    } else if (nrow(aux_df) > 1000  & nrow(aux_df) <= 2000) {
      breaksXaxis <- c(seq(0, nrow(aux_df), by = 200))
    } else {
      breaksXaxis <- c(seq(0, nrow(aux_df), by = 500))
    }
    
    plotList[[i]] <- ggplot(aux_df,aes(x=positions,y=IUPredScores)) +
      geom_col(aes(fill=IUPredScores >= 0.5 ,color=IUPredScores >= 0.5)) +
      ggpubr::theme_classic2() + 
      theme(text = element_text(size=19),legend.position = "none") +
      scale_x_continuous(limits = c(0, nrow(aux_df)+1),breaks = breaksXaxis,expand = c(0.005,0.005))+ xlab("Positions") +
      scale_y_continuous(limits = c(-0.25, 1),breaks = c(seq(0,1,by=0.25)),expand = c(0.05,0.05)) + ylab("Score") +
      scale_fill_manual(values = c(pal_jco()(10)[3],"#f6921eff")) +
      scale_color_manual(values = c(pal_jco()(10)[3],"#f6921eff")) +
      annotate("point",x = dataframe$ST_residues[[i]],y=rep(0,length(dataframe$ST_residues[[i]])),shape="|",color=pal_jco()(10)[1],fill=pal_jco()(10)[1],size=10) +
      annotate("point",x = dataframe[i,sites_MAPK][[1]],y=rep(-0.08,length(dataframe[i,sites_MAPK])),shape="|",color=pal_jco()(10)[1],fill=pal_jco()(10)[1],size=10) +
      annotate("point",x = dataframe[i,sites_AURK][[1]],y=rep(-0.16,length(dataframe[i,sites_AURK][[1]])),shape="|",color=pal_jco()(10)[1],fill=pal_jco()(10)[1],size=10) +
      annotate("point",x = dataframe[i,sites_CDK][[1]],y=rep(-0.08,length(dataframe[i,sites_CDK][[1]])),shape=21,color=pal_jco()(10)[8],fill="#ffdd15ff",size=6) +
      annotate("point",x = dataframe[i,sites_MAPK][[1]],y=rep(-0.16,length(dataframe[i,sites_MAPK][[1]])),shape=21,color=pal_jco()(10)[8],fill="olivedrab3",size=6) +
      annotate("point",x = dataframe[i,sites_AURK][[1]],y=rep(-0.24,length(dataframe[i,sites_AURK][[1]])),shape=21,color=pal_jco()(10)[8],fill="firebrick3",size=6) +
      geom_hline(yintercept = 0.5,size=0.5,linetype = "dashed",color = "darkslategrey") +
      ggtitle(plotTittle)

    
    
  }
  names(plotList)<-dataframe[,id_col]
  return(plotList)
}
