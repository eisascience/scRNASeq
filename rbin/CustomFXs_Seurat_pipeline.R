






CTL_Immune_GeneList <- function(QuickGO.path="./data/QuickGO"){
  
  SGS.LS <- list()
  SGS.LS$HighlyActivated <- c("TNFRSF9", "NFKBID", "CCL4", "CCL4L2", "IRF8", "IxNG", "CD83", "CD82", "PLEK", "RGCC") #"ENSMMUG00000013779"
  SGS.LS$LessActivated <- c("LTB", "IL7R", "PTPRCAP", "GZMK")
  SGS.LS$Pheno1 <- c("MT1M", "MGEA5")
  SGS.LS$CD8Canonical <- c("CD8A", "CD8B", "CD3G", "CD3D", "CD3E") 
  SGS.LS$BCellCanonical <- c("CD19", "MS4A1", "CD79A") 
  SGS.LS$CD4Canonical <- c("CD4", "IL7R")
  SGS.LS$NKCanonical <- c("GZMK", "LYAR", "NKG7", "GZMA", "GNLY", "FGFBP2", "FCGR3A", "CCL4", "GZMH")
  
  AnnotationFiles <- list.files(QuickGO.path, pattern = "QuickGO", full.names = T)
  GeneLists <- list()
  
  for(AnnFile in AnnotationFiles){
    # AnnFile = AnnotationFiles[1]
    tempName = gsub(".txt", "", gsub("_","",gsub("annotations-", "", gsub(".tsv","",gsub("-","",basename(AnnFile))))))
    GeneLists$Extra[[tempName]] <-  read.table(
      AnnFile,
      sep="\t", header=TRUE, row.names = NULL, fill = TRUE )
  }
  
  SGS.LS$QuickGOgenes <- data.table::rbindlist(lapply(GeneLists$Extra, function(setX){
    subset(setX, TAXON.ID == 9606)[,c("GO.NAME", "SYMBOL")] #10090 = mouse ; 9606 = human
  }))
  
  return(SGS.LS)
}



MakeSerObjs_10XFolders <- function(counts.path = NULL, 
                                   min.cells = 0, 
                                   min.genes = 0,
                                   ProjName="10X",
                                   save.path = NULL,
                                   returnList=F){
  
  require(Seurat)
  if(returnList) TempLS <- list()

  if(!is.null(counts.path)){
    if(is.null(save.path)) save.path <- counts.path
    
    exp.dirs <- list.dirs(counts.path, recursive = F, full.names = T)
    
    if(returnList) TempLS$exp.dirs <- exp.dirs
    
    for(xN in 1:length(exp.dirs)){
      # xN = 1
      # use the list.files below to make sure the expected files are there....
      #list.files(exp.dirs[[xN]], full.names = T, recursive = T)
      
      if(returnList) TempLS$SeuratObjs <- list()
      
      print(exp.dirs[[xN]])
      print("Reading in 10X folder...")
      Seurat10X  <- Read10X(data.dir = exp.dirs[[xN]])
      
      print("Converting to Seurat Obj....")
      SeuratObjs <- CreateSeuratObject(raw.data = Seurat10X,
                                       min.cells = min.cells,   #genes expressed in >= 5 cells
                                       min.genes = min.genes, #Keep all cells with at least 200 detected genes
                                       project = paste(ProjName, xN, sep="_expN"))
      
      
      print("Making matrix sparse...")
      SeuratObjs <- MakeSparse(object = SeuratObjs)
      
      if(returnList) TempLS$SeuratObjs[[basename(exp.dirs[[xN]])]] <- SeuratObjs
      
      print("saving...")
      saveRDS(SeuratObjs, 
              paste(save.path, "/SeuratObj_", ProjName, "_", xN , ".rds", sep=""))
      
      
    }; remove(xN)
    
    if(returnList) return(TempLS)
  }
  
}



PreProcess_SerObjFolders <- function(SerObj.path = NULL,
                                   ProjName="10X",
                                   save.path = NULL, save.fig.path = NULL, 
                                   returnList=F, save.fig = T,
                                   ENSMB.tag="ENSMM", RhesusConvDavid = F,
                                   nUMI.high = 20000, nGene.high = 3000, pMito.high = 0.15, 
                                   nUMI.low = 0.99, nGene.low = 200, pMito.low = -Inf,
                                   RhesusConvDavid.path = "./data/Rhesus/David6.8_ConvertedRhesus_ENSMMUG.txt",
                                   fvg.x.low.cutoff = 0.01, fvg.x.high.cutoff = 4.5, fvg.y.cutoff = 1.5,
                                   KeepGene.LS =NULL, 
                                   nDimPCA=15,
                                   returnList=F){
  
  require(Seurat)
  if(returnList) TempLS <- list()

  if(!is.null(SerObj.path)){
    if(is.null(save.path)) {
      save.path <- paste(SerObj.path, "/SerProc", sep="")
      if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
    }
    if(is.null(save.fig.path)) {
      save.fig.path <- paste(SerObj.path, "/SerProcFigs", sep="")
      if(!dir.exists(save.fig.path)) dir.create(save.fig.path, recursive = T)
    }

    
    
    all_RDS  <- list.files(SerObj.path, full.names = T, pattern = ".rds")
    SeurObj_RDS <-  all_RDS[grep("SeuratObj.rds", all_RDS)]
    
    
    for(xN in 1:length(SerObj.files)){
      # xN=1
      if(returnList) TempLS$SeuratObjs <- list()
      
      SeuratObjs <- readRDS(SerObj.files[xN])
      
      mito.genes <- grep(pattern = "^MT-", x = rownames(x = SeuratObjs@raw.data), value = TRUE)
      length(mito.genes)
      
      percent.mito <- Matrix::colSums(SeuratObjs@raw.data[mito.genes, ]) / Matrix::colSums(SeuratObjs@raw.data)
      
      SeuratObjs <- AddMetaData(object = SeuratObjs,
                                             metadata = percent.mito,
                                             col.name = "percent.mito")
      
      TotalPerGeneExpressed      <- rowSums(SeuratObjs@raw.data)
      TotalPerGeneExpressed.perc <- round(TotalPerGeneExpressed/sum(TotalPerGeneExpressed)*100, 5)
      
      TotalPerCellExpressed <- colSums(SeuratObjs@raw.data)
      TotalPerCellExpressed.perc <- round(TotalPerCellExpressed/sum(TotalPerCellExpressed)*100, 5)
      
      SeuratObjs <- AddMetaData(object = SeuratObjs,
                              metadata = TotalPerCellExpressed,
                              col.name = "SumTotGeneExprPerCell")
      SeuratObjs <- AddMetaData(object = SeuratObjs,
                              metadata = TotalPerCellExpressed.perc,
                              col.name = "PercentSumTotGeneExprPerCell")
      

      
      if(save.fig) png(filename =  paste(save.fig.path, "/ViolinPlotTrio_preFilt.png", sep=""), width = 15, height = 10, units = "in", res=200)
      VlnPlot(object = SeuratObjs,
              features.plot = c("nGene", "nUMI", "percent.mito", "PercentSumTotGeneExprPerCell"),
              nCol = 3, cols.use = col_vector, x.lab.rot=T, size.x.use = 11)
      if(save.fig) dev.off()
      
  
      #Genes that dont map to a specific name
      noGeneSYM <- rownames(SeuratObjs@raw.data)[grepl(ENSMB.tag, rownames(SeuratObjs@raw.data))]
      
      length(noGeneSYM)
      
      # write.table(noGeneSYM, 
      #             "./10X/Rhesus_ENSMMUG.csv",
      #             sep=", ", , row.names = F, quote = F,
      #             col.names = F)
      
      if(save.fig)  
      if(RhesusConvDavid){
        
        David6.8ConvTable <- data.frame(read.csv(RhesusConvDavid.path, sep = "\t", header = T))
        rownames(David6.8ConvTable) <- David6.8ConvTable$From
        David6.8ConvTable <- David6.8ConvTable[noGeneSYM, ]
        length(unique(noGeneSYM)); length((noGeneSYM))
        rownames(David6.8ConvTable) <- noGeneSYM
        
        David6.8ConvTable$Final <- as.character(David6.8ConvTable$To)
        
        David6.8ConvTable$Final[which(is.na(David6.8ConvTable$To))] <- rownames(David6.8ConvTable)[which(is.na(David6.8ConvTable$To))]
        
        rownames(SeuratObjs@raw.data)[grepl("ENSMM", rownames(SeuratObjs@raw.data))] <- David6.8ConvTable$Final
        
        length(rownames(SeuratObjs@raw.data)); length(unique(rownames(SeuratObjs@raw.data)))
        
        
        duplicatedGeneNames <- names(which(table(rownames(SeuratObjs@raw.data))>1))
        
        #change the second duplace name to add a .2
        #perhaps can avg the expr?
        for(geneN in duplicatedGeneNames){
          rownames(SeuratObjs@raw.data)[which(rownames(SeuratObjs@raw.data)==geneN)[2]] <- paste(geneN, ".2", sep="")
          
        }
        
        rownames(SeuratObjs@data) <- rownames(SeuratObjs@raw.data)
      }
      
      
      
      SeuratObjs <- FilterCells(object = SeuratObjs,
                                             subset.names = c("nUMI", "nGene", "percent.mito"),
                                             low.thresholds = c(nUMI.low,   nGene.low,     pMito.low),
                                             high.thresholds = c(nUMI.high, nGene.high,    pMito.high))
      
      if(save.fig) png(filename =  paste(save.fig.path, "/ViolinPlotTrio_postFilt.png", sep=""), width = 15, height = 10, units = "in", res=200)
      VlnPlot(object = SeuratObjs,
              features.plot = c("nGene", "nUMI", "percent.mito", "PercentSumTotGeneExprPerCell"),
              nCol = 3, cols.use = col_vector, x.lab.rot=T, size.x.use = 11)
      #CleaningLS$Figs$Combo$ViolinPlotTrio_postFilt <- recordPlot()
      if(save.fig) dev.off()
      
      
      SeuratObjs <- NormalizeData(object = SeuratObjs, 
                                               normalization.method = "LogNormalize", 
                                               scale.factor = 10000)
      
      
      
      SeuratObjs <- FindVariableGenes(object = SeuratObjs,
                                                   mean.function = ExpMean,
                                                   dispersion.function = LogVMR,
                                                   x.low.cutoff = fvg.x.low.cutoff, #X-axis function is the mean expression level
                                                   x.high.cutoff = fvg.x.high.cutoff, #based on plot viz
                                                   y.cutoff = fvg.y.cutoff, #Y-axis it is the log(Variance/mean)
                                                   num.bin = 40) #y.cutoff = 1 sd away from averge within a bin
      
      SeuratObjs <- ScaleData(object = SeuratObjs, vars.to.regress = c("nUMI", "percent.mito"))

      
      if(!is.null(KeepGene.LS)){
        length(SeuratObjs@var.genes)      
        SeuratObjs@var.genes <- unique(c(SeuratObjs@var.genes, as.character(unlist(KeepGene.LS))))
        length(SeuratObjs@var.genes)
      }
      
      

      
      SeuratObjs <- RunPCA(object = SeuratObjs,
                                        pc.genes = SeuratObjs@var.genes,
                                        do.print = F)
      
      SeuratObjs <- ProjectPCA(object = SeuratObjs, do.print = FALSE)
      
      
      if(save.fig) png(filename =  paste(save.fig.path, "/PCElbowPlot.png", sep=""), width = 10, height = 10, units = "in", res=200)
      PCElbowPlot(SeuratObjs)
      if(save.fig) dev.off()
      
      if(save.fig) png(filename =  paste(save.fig.path, "/PCAHeatmap.png", sep=""), width = 10, height = 10, units = "in", res=200)
      PCHeatmap(object = SeuratObjs,
                pc.use = 1:5,
                cells.use = 200,
                do.balanced = TRUE,
                label.columns = FALSE,
                use.full = FALSE)
      if(save.fig) dev.off()
      
      
      SeuratObjs <- FindClusters(object = SeuratObjs, 
                               reduction.type = "pca", 
                               dims.use = 1:nDimPCA, 
                               resolution = 0.6, 
                               print.output = 0, 
                               save.SNN = TRUE)
      
      SeuratObjs <- StashIdent(object = SeuratObjs, 
                             save.name = "Cluster_PCA_0.6")
      
      if(save.fig) png(filename =  paste(save.fig.path, "/PCAplot_clust.png", sep=""), width = 10, height = 10, units = "in", res=200)
      PCAPlot(object = SeuratObjs, dim.1 = 1, dim.2 = 2)
      if(save.fig) dev.off()
      
      
      
      SeuratObjs <- RunTSNE(object = SeuratObjs, 
                                         reduction.type = "pca", 
                                         dims.use = 1:nDimPCA)
      
      if(save.fig) png(filename =  paste(save.fig.path, "/TSNEplot_clust.png", sep=""), width = 10, height = 10, units = "in", res=200)
      TSNEPlot(object = SeuratObjs, do.label = TRUE)
      if(save.fig) dev.off()
      
      SeuratObjs <- RunUMAP(SeuratObjs, 
                                         cells.use = NULL, 
                                         dims.use = 1:nDimPCA, 
                                         reduction.use = "pca", 
                                         max.dim = 2L,
                                         reduction.name = "umap", 
                                         reduction.key = "UMAP", 
                                         n_neighbors = 30L,
                                         min_dist = 0.3, 
                                         metric = "correlation", 
                                         seed.use = 42)


      if(save.fig) png(filename =  paste(save.fig.path, "/UMAPplot_clust.png", sep=""), width = 10, height = 10, units = "in", res=200)
      DimPlot(object = SeuratObjs, reduction.use = 'umap')
      if(save.fig) dev.off()
      
      
      SeuratObjs <- SetAllIdent(object = SeuratObjs, id = "expN")
      
      if(save.fig) png(filename =  paste(save.fig.path, "/PCAplot_expN.png", sep=""), width = 10, height = 10, units = "in", res=200)
      PCAPlot(object = SeuratObjs, dim.1 = 1, dim.2 = 2)
      if(save.fig) dev.off()
      
      if(save.fig) png(filename =  paste(save.fig.path, "/TSNEplot_expN.png", sep=""), width = 10, height = 10, units = "in", res=200)
      TSNEPlot(object = SeuratObjs, do.label = TRUE)
      if(save.fig) dev.off()
      
      if(save.fig) png(filename =  paste(save.fig.path, "/UMAPplot_expN.png", sep=""), width = 10, height = 10, units = "in", res=200)
      DimPlot(object = SeuratObjs, reduction.use = 'umap')
      if(save.fig) dev.off()
      
      
      print("saving...")
      saveRDS(SeuratObjs, 
              paste(save.path, "/SeuratObj_", ProjName, "_", xN , ".rds", sep=""))
      
      if(returnList) TempLS$SeuratObjs[[basename(exp.dirs[[xN]])]] <- SeuratObjs
      

    }
    
    if(returnList) return(TempLS)
    

  }
}