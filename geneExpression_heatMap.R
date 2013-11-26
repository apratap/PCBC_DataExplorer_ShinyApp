#load the memoised version of pheatmap
source("memoised_pheatmap.R")


#gene expression heatmap logic
get_geneExpression_heatMap <- function(m,annotation,...){
  
  #change the data type to integer and 
  # NOT FOR NOW : transform the count to log 2 scale
  m <- apply(m,2,function(x) as.numeric(x))
  
  #removing those genes which dont vary much across the samples
  # so any gene with SD < .2 across the samples will be dropped 
  drop_genes <- which(apply(m,1,sd) < .2)
  m <-  m[-drop_genes,]
  
  #check if remaining genes can be clustered
  if(nrow(m) < 3){
    error_msg <- sprintf('After removing genes with SD < .20 %d genes remain \n 
                          More than 2 needed to cluster \n\n', nrow(m))
    stop(error_msg)
  }
  else{ #do the clustering and heatmap
    #scaling genes across experiments
    mat.scaled <- t(scale(t(m)))  
    memoised_pheatmap(mat.scaled,
                      scale="none",
                      annotation = annotation,
                      clustering_distance_rows = "correlation",
                      clustering_distance_cols = "correlation",
                      clustering_method = "average",
                      border_color = NA,
                      drawRowD = FALSE
    )
  }
  
}

