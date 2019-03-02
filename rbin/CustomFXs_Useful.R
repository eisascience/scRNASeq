

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

MDSmyDF <- function(dfx, labelsDF, factorV, title = "MDS Plot", col_vector){
  
  # mds <- (as.data.frame(dfx)  %>%
  #        cbind(cmdscale(as.matrix(dist(dfx) ))) )  
  
  # dfx = as.data.frame(t(assay(dds)+1))
  # labelsDF = dge$samples[,c("Group", "GroupDate", "Day")]
  #factorV = labelsDF[, colorbyN]
  
  if(length(factorV) == nrow(dfx)){
    mds <- cmdscale(as.matrix(dist(dfx)))
    colnames(mds) <- c("MDS1", "MDS2")
    
    mds <- cbind(mds, labelsDF) #samples as rows add a cols of labels
    
    
    p1 <- ggplot(mds, aes(x = MDS1, y = MDS2)) +
      theme_bw() +
      geom_hline(yintercept = 0, color = "gray70") +
      geom_vline(xintercept = 0, color = "gray70") +
      geom_point(aes(color = factorV), alpha = 0.8, size = 3) +
      coord_cartesian(xlim = c(min(mds[,1])-5, max(mds[,1])+5)) + 
      scale_color_manual(values = col_vector, name="Samples")
    
    # the graphic with ggrepel
    p1 + geom_text_repel(aes(y = MDS2 + 0.25), label = factorV) +  
      ggtitle(paste("MDS of:",title ,sep=" "))+ 
      theme(plot.title = element_text(hjust = 0.5))
    
  } else {
    print("rotate t() your dataX")
  }
  
  
}

PCAmyDF <- function (dfx, labels, factorV, title = "PCA Plot", scale, center, col_vector, namePointBL = F) {
  # dfx = log.cpm.one; scale = F; center = F; namePointBL = F
  # factor = dge$samples$GroupDate
  # title = "PCA Plot"
  if(class(labels) == "function") {
    print("no labels, using factor as names")
    labels = as.character(factorV)
  }
  if(length(as.character(factorV)) != length(labels)) {
    print("labels length != factor length, using factor as names")
    labels = as.character(factorV)
  }
  
  dfx.pca <- prcomp(t(dfx), scale.=scale, center = center)
  
  MinXY <- min(c(round(min(dfx.pca$rotation[,1:2]) - abs(min(dfx.pca$rotation[,1:2])*0.5),1), -1) )
  MaxXY <- max(c(round(max(dfx.pca$rotation[,1:2]) + abs(max(dfx.pca$rotation[,1:2])*0.5),1),  1) )
  
  
  if(namePointBL){
    autoplot(dfx.pca) +
      theme_bw() +
      geom_point(aes(color = factorV), alpha = 0.8, size = 3) + 
      scale_color_manual(values = col_vector, name="Samples") + 
      geom_text_repel(aes(y = PC2, label = labels))  +  
      ggtitle(paste("PCA of:",title ,sep=" ")) + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      xlim(c(MinXY,MaxXY)) + ylim(c(MinXY,MaxXY))
  } else {
    autoplot(dfx.pca) +
      theme_bw() +
      geom_point(aes(color = factorV), alpha = 0.8, size = 3) + 
      scale_color_manual(values = col_vector, name="Samples")  +  
      ggtitle(paste("PCA of:",title ,sep=" "))+ 
      theme(plot.title = element_text(hjust = 0.5)) + 
      xlim(c(MinXY,MaxXY)) + ylim(c(MinXY,MaxXY))
  }
  
  
  
  
}

transposedt <- function(dt, varlabel="myVar") {
  require(data.table)
  dtrows = names(dt)
  dtcols = as.list(c(dt[,1]))
  dtt = transpose(dt)
  dtt[, eval(varlabel) := dtrows]
  setnames(dtt, old = names(dtt), new = c(dtcols[[1]], eval(varlabel)))
  dtt = dtt[-1,]
  setcolorder(dtt, c(eval(varlabel), names(dtt)[1:(ncol(dtt) - 1)]))
  return(dtt)
}