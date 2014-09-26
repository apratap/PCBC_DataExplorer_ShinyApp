library(memoise)

#for faster rendering caching the computationally expensive functions
memoised_corAndPvalue <- memoise(function(...) corAndPvalue(...))


get_expMatrix_withcorrelated_genes <- function(geneIds, expMatrix, corThreshold){
  cat('Calculating correlated genes ....')  
  #expression matrix with selected genes
  m1 <- expMatrix[rownames(expMatrix) %in%  geneIds,]
  #expression matrix with which the selected genes will be correlated
  m2 <- expMatrix
  #calculate correlation
  res <- memoised_corAndPvalue(t(m1),t(m2),nThreads=4)
  cor <- round(res$cor,digits=3)
  cor <- abs(cor) >= corThreshold
  #columns of the cor matrix which have correlation with some gene > corThreshold
  cols_to_select <- apply(cor,2,any)
  correlated_genes <- colnames(cor)[cols_to_select]
  cat('Done','\n')
  expMatrix[rownames(expMatrix) %in% correlated_genes,]
}


#filter metadata
get_filtered_metadata <- function(input, metadata){
  filtered_metadata <- metadata
  
  if( length(input$linetype) != 0 ){
    filtered_metadata <- filter(filtered_metadata, Line_Type %in% input$linetype)
  }
  if( length(input$gene_combination) != 0 ){
    filtered_metadata <- filter(filtered_metadata, Reprogramming_Gene_Combination %in% input$gene_combination)  
  }
  if(length(input$vector_type) != 0){
    filtered_metadata <- filter(filtered_metadata, Reprogramming_Vector_Type %in% input$vector_type)
  }
  if(length(input$tissue_origin) != 0){
    filtered_metadata <- subset(filtered_metadata, Tissue_of_Origin %in% input$tissue_origin)
  }
  if(length(input$diff_state) != 0){
    filtered_metadata <- subset(filtered_metadata, Differentiation_State %in% input$diff_state)
  }
  if(length(input$cell_origin) != 0){
    filtered_metadata <- subset(filtered_metadata, Cell_Type_of_Origin %in% input$cell_origin)
  }
  filtered_metadata
}


#create the annotation data frame for the heatmap
get_filteredAnnotation <- function(input,metadata){
  if(length(input$heatmap_annotation_labels) == 0){
    stop('please select atleast one heatmap annotation variable \n\n')      
  }
  else{
    annotation <- metadata[,c(input$heatmap_annotation_labels),drop=F]
    rownames(annotation) <- metadata$Sample
    annotation
  }
}
