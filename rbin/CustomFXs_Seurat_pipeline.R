

MakeSerObjs_10Xfolders <- function(counts.path = NULL, 
                                   min.cells = 0, 
                                   min.genes = 0,
                                   ProjName="10X",
                                   save.path = NULL){
  
  require(Seurat)

  if(!is.null(counts.path)){
    if(is.null(save.path)) save.path <- counts.path
    
    
    #counts.files   <- list.files(counts.path, full.names = T, recursive = T); counts.files
    exp.dirs <- list.dirs(counts.path, recursive = F, full.names = T)
    
    for(xN in 1:length(exp.dirs)){
      # xN = 1
      # use the list.files below to make sure the expected files are there....
      #list.files(exp.dirs[[xN]], full.names = T, recursive = T)
      
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
      
      print("saving...")
      saveRDS(SeuratObjs, 
              paste(save.path, "/SeuratObj_", ProjName, "_", xN , ".rds", sep=""))
      
      
    }; remove(xN)
    
  }
  
}

