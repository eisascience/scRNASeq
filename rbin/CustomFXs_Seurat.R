# Customized functions from or for Seurat analysis...


AddModuleScoreSindv <- function(Ser_object=NULL, geneSet=NULL){
  # compute individual scores of genes in a gene list 
  #      as oppose to Seurat default ... it too does it individually, but controls the naming
  
  # Also we scale the results between 0 and 1 for probability equivalence 
  #     the range01() is in the Useful.R function set
  
  
  for(geneN in geneSet){
    print(geneN)
    #geneN = geneSet[1]
    Ser_object <- AddModuleScore(Ser_object, 
                                 genes.list = geneN, 
                                 genes.pool = rownames(Ser_object@data), 
                                 n.bin = 25,
                                 seed.use = 1, 
                                 ctrl.size = 100, 
                                 use.k = FALSE, 
                                 enrich.name = geneN,
                                 random.seed = 1)
    
    Ser_object@meta.data[,paste(geneN, "1", sep="")] <- range01(Ser_object@meta.data[,paste(geneN, "1", sep="")])
    
  }
  #colnames(Ser_object@meta.data[,paste(geneSet, "1", sep="")]) <- geneSet
  return(Ser_object)
} 