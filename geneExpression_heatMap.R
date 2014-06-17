#load the memoised version of pheatmap
source("memoised_pheatmap.R")


#gene expression heatmap logic
get_geneExpression_heatMap <- function(m, annotation = NA ,
                                       clustering_distance_rows = "correlation",
                                       clustering_distance_cols = "correlation",
                                       explicit_rownames = 'none', ...){
  if(nrow(m) <= 2){
    return(memoised_pheatmap(m, cluster_rows=FALSE,
                             scale="none",
                             annotation = annotation,
                             drawRowD = FALSE,
                             border_color = NA,
                             explicit_rownames = explicit_rownames,...))
  }
  if(nrow(m) >= 70){
    #removing those genes which dont vary much across the samples
    #so any gene with SD < .2 across the samples will be dropped 
    drop_genes <- which(apply(m,1,sd) < .2)
    #following step to remove the bug seen 
    #when m <-  m[-drop_genes,] is done directly and length(drop_genes) = 0
    if(length(drop_genes) != 0){
      m <-  m[-drop_genes,]  #filtering a mat , IMP
      #also remove the same from the explicit rownames as those genes are taken out in anycase
      explicit_rownames <- explicit_rownames[-drop_genes] #filtering a vector no , needed
    } 
  }
  
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
                      clustering_distance_rows = clustering_distance_rows,
                      clustering_distance_cols = clustering_distance_cols,
                      clustering_method = "average",
                      border_color = NA,
                      drawRowD = FALSE,
                      explicit_rownames = explicit_rownames,...)
  }
}

