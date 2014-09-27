
#Define the server the logic
shinyServer(function(input,output,session){
  #get the list of user submitted genes
  user_submitted_geneList <- reactive({
    input$custom_search
    geneList <- isolate(input$custom_gene_list)
    geneList <- unlist(strsplit(geneList, split=c('[\\s+,\\n+\\r+)]'),perl=T))
    #conevert everything to upper case
    geneList <- toupper(geneList)
    geneList <- geneList[ !geneList == "" ] #remove the blank entries
    geneList
  })
  
  #get the list of user submitted genes
  user_submitted_miRNAlist <- reactive({
    input$custom_search
    miRNAlist <- isolate(input$custom_miRNA_list)
    miRNAlist  <- unlist(strsplit(miRNAlist,split=c('[\\s+,\\n+\\r+)]'),perl=T))
    #conevert everything to upper case
    miRNAlist  <- tolower(miRNAlist)
    miRNAlist  <- miRNAlist[ !miRNAlist == "" ] #remove the blank entries
    miRNAlist
  })
  
  #get list of miRNAs
  selected_miRNAs <- reactive({
    #get the list of geneIds that were selected by the user
    # + ones correlated with other genes (if corr option selected) 
    #this is the reason why not getting geneIds from selected_genes() as it wont have the correlated genes
    geneIds <- rownames(get_filtered_mRNA_matrix())
    #get miRNA targetting the selected genes
    selected_miRNAs <- filter(miRNA_to_genes, GeneID %in% geneIds)
    selected_miRNAs <- unique(paste(selected_miRNAs$miRNA1,selected_miRNAs$miRNA2,sep=','))
    selected_miRNAs
  })
  
  #get list of genes in current pathway or user entered list
  selected_genes <- reactive({
    if( input$genelist_type == 'custom_gene_list'  ){
      genes <- unique(user_submitted_geneList())
      miRNAs <- user_submitted_miRNAlist()
      #get miRNA targetting the selected genes
      selected_miRNAs_targetGenes <- filter(miRNA_to_genes, miRNAPrecursor %in% miRNAs | miRNA1 %in% miRNAs | 
                                            miRNA2 %in% miRNAs)
      selected_miRNAs_targetGenes <- unique(selected_miRNAs_targetGenes$GeneID)
      genes <- unique(c(genes, selected_miRNAs_targetGenes))
    } else if( input$genelist_type == 'precomputed_significant_geneList'){
      if(input$enrichedPathways == 'ALL'){
        genes_in_selected_GeneList <- sigGenes_lists[[input$selected_Significant_GeneList]]
        genes <- unique(genes_in_selected_GeneList)
      } else {
        #1. get a list of all genes in the selected enriched pathway
        #trimming the suffix : #pdj-
        pathway = gsub('#p.adj_.*','',input$enrichedPathways)
        genes_in_pathway <- MSigDB$C2.CP.KEGG[[pathway]]
        genes_in_selected_GeneList <- sigGenes_lists[[input$selected_Significant_GeneList]]
        genes <- intersect(genes_in_pathway, genes_in_selected_GeneList)
      }
   } else if( input$genelist_type == 'pathway'){
      genes <- as.character(unlist(pathways_list[input$selected_pathways]))
   } else  genes 
  })
  
  
  
  #get list of pathways enriched in the geneList selected by the user
  get_enrichedPathways <- reactive({
      #return the enriched pathway for a gene list
      #labels contain the pvalue of the FET test       
      precomputed_enrichedPathways_in_geneLists[[input$selected_Significant_GeneList]]
  })

  #update the enriched pathways for the user selected genelist
  observe({
    enriched_Pathways = sort(get_enrichedPathways())
    updateSelectInput(session = session,
                      inputId = "enrichedPathways",
                      label = sprintf('Enriched pathway/s: %d (?)', sum(! enriched_Pathways %in% c('NA','ALL'))),
                      choices = enriched_Pathways,
                      selected = enriched_Pathways[[1]]                  
    ) 
  })
  
  output$mRNA_compute_time <- renderPrint({
    print(mRNA_heatmap_compute_results$results$time)
  })
  
  get_filtered_mRNA_matrix <- reactive({
    #a.) subset on sample names based on user selected filters
    filtered_metadata <- get_filtered_metadata(input,combined_metadata)
    filtered_mRNA_NormCounts <- mRNA_NormCounts[,colnames(mRNA_NormCounts) %in% filtered_metadata$Sample ]
    #b.) subset based on selected genes
    selected_genesId <- convert_to_ensemblIds(selected_genes())
    if(input$incl_corr_genes == 'TRUE' & input$genelist_type == 'custom_gene_list'){ 
      filtered_mRNA_NormCounts <- get_expMatrix_withcorrelated_genes(selected_genesId, filtered_mRNA_NormCounts, input$corr_threshold)
    } else {  
      filtered_mRNA_NormCounts <- filtered_mRNA_NormCounts[rownames(filtered_mRNA_NormCounts) %in% selected_genesId,]
    }
    filtered_mRNA_NormCounts
  })
  
   #return the mRNA heatMap plot
   output$mRNA_heatMap <- renderPlot({  
     m <- get_filtered_mRNA_matrix()
     # zero variance filter
     rows_to_keep <- apply(m,1,var) > 0
     m <- m[rows_to_keep, ]
     m <- data.matrix(m)
     
     validate( need( ncol(m) != 0, "Filtered mRNA expression matrix contains 0 Samples") )
     validate( need( nrow(m) != 0, "Filtered mRNA expression matrix contains 0 genes") )
     validate( need(nrow(m) < 10000, "Filtered mRNA expression matrix contains > 10000 genes. MAX LIMIT 10,000 ") )
     fontsize_row=8
     fontsize_col=8
     if(nrow(m) > 100){ fontsize_row = 0 }
     if(ncol(m) > 50){ fontsize_col=0 }
     #convert ensembl ID's to gene name
     explicit_rownames = hg19_annot %>%
                               filter(ENSEMBL %in% rownames(m)) %>%
                               group_by(ENSEMBL) %>%
                               summarise(SYMBOL = unique(SYMBOL)[1])
     explicit_rownames <- explicit_rownames$SYMBOL
     #annotation
     filtered_metadata <- get_filtered_metadata(input, combined_metadata)
     annotation <- get_filteredAnnotation(input, filtered_metadata)
     
     withProgress(session, {
       setProgress(message = "clustering & rendering heatmap, please wait", 
                   detail = "This may take a few moments...")
       expHeatMap(m,annotation,
                  clustering_distance_rows = input$clustering_distance,
                  clustering_distance_cols = input$clustering_distance,
                  fontsize_col=fontsize_col, 
                  fontsize_row=fontsize_row,
                  scale=T,
                  clustering_method = input$clustering_method,
                  explicit_rownames = explicit_rownames)
     }) #END withProgress
   })
  
  output$microRNA_heatMap <- renderPlot({
    #get the microRNA expression matrix
    filtered_microRNA_NormCounts <- miRNA_normCounts[row.names(miRNA_normCounts) %in% selected_miRNAs(),]
    
    #subset on sample names based on user selected filters 
    filtered_metadata <- get_filtered_metadata(input,combined_metadata)
    filtered_microRNA_NormCounts <- filtered_microRNA_NormCounts[ ,colnames(filtered_microRNA_NormCounts) %in% filtered_metadata$Sample]
    
    #annotation <- get_filteredAnnotation(input,filtered_miRNA_metadata)
    m <- filtered_microRNA_NormCounts
    # zero variance filter
    rows_to_keep <- apply(m,1,var) > 0
    m <- m[rows_to_keep, ]
    m <- data.matrix(m)
    validate( need( nrow(m) != 0, "Filtered miRNA expression matrix contains 0 genes") )
    validate( need(nrow(m) < 10000, "Filtered miRNA expression matrix contains > 10000 genes. MAX LIMIT 10,000 ") )
    fontsize_row=4
    fontsize_col=8
    if(nrow(m) > 200){ fontsize_row = 0 }
    if(ncol(m) > 50){ fontsize_col=0 }
    annotation <- get_filteredAnnotation(input, filtered_metadata)
    
    withProgress(session, {
      setProgress(message = "clustering & rendering heatmap, please wait", 
                  detail = "This may take a few moments...")
      expHeatMap(m,annotation,
                 clustering_distance_rows = input$clustering_distance,
                 clustering_distance_cols = input$clustering_distance,
                 fontsize_col=fontsize_col, 
                 fontsize_row=fontsize_row,
                 scale=T,
                 clustering_method = input$clustering_method,
                 color=colorRampPalette(rev(brewer.pal(n = 7, name = "BrBG")))(100))
    }) #END withProgress
  
  })

  #get list of miRNAs
  selected_methProbes <- reactive({
    #get the list of geneIds that were selected by the user
    # + ones correlated with other genes (if corr option selected) 
    #this is the reason why not getting geneIds from selected_genes() as it wont have the correlated genes
     geneIds <- rownames(get_filtered_mRNA_matrix())
     #convert to entrezID
     entrez_geneIds <- convert_to_EntrezIds(geneIds)
     
     flt_res <- filter(meth_to_gene, entrezID %in% entrez_geneIds)
     selected_methProbes <- unique(flt_res$methProbe)
     selected_methProbes
  })
  
  
  output$methylation_heatMap <- renderPlot({
    #get the filtered methylation data
    flt_meth_data <- meth_data[row.names(meth_data) %in% selected_methProbes(),]
   
    #subset on sample names based on user selected filters 
    filtered_metadata <- get_filtered_metadata(input,combined_metadata)
    flt_meth_data <- flt_meth_data[ ,colnames(flt_meth_data) %in% filtered_metadata$Sample]
    
    m <- flt_meth_data
    validate( need( nrow(m) != 0, "Filtered methylation data matrix contains 0 genes") )
    
    # zero variance filter
    var_methProbe <- apply(m,1,var)
    rows_to_keep <- var_methProbe > .01
    m <- m[rows_to_keep, ]
    m <- data.matrix(m)
    
    annotation <- get_filteredAnnotation(input,filtered_metadata)
    validate( need( nrow(m) != 0, "Filtered methylation data matrix contains 0 genes") )
    validate( need(nrow(m) < 5000, "Filtered methylation data matrix > 5000 genes. MAX LIMIT 5,000 ") )
    fontsize_row=8
    fontsize_col=8
    if(nrow(m) > 100){ fontsize_row = 0 }
    if(ncol(m) > 50){ fontsize_col=0 }
    
    withProgress(session, {
      setProgress(message = "clustering & rendering heatmap, please wait", 
                  detail = "This may take a few moments...")
      expHeatMap(m,annotation,
                 clustering_distance_rows = input$clustering_distance,
                 clustering_distance_cols = input$clustering_distance,
                 fontsize_col=fontsize_col, 
                 fontsize_row=fontsize_row,
                 clustering_method = input$clustering_method)
    }) #END withProgress
  })

  #create a table with selected gene list and merge with some annotation
  output$geneExpTable <- renderDataTable({
    filtered_mRNA_NormCounts <- subset(mRNA_NormCounts, symbol %in% selected_genes())
    df <- merge(filtered_mRNA_NormCounts[,1:3], hg19_gene_annot, by.x='symbol',by.y='SYMBOL')
    df
  })

  output$mRNA_summary <- renderTable({
    summary <- data.frame('Category' =  c('#Uniq genes in current list/pathway', '#genes found with exp values', 
                                          '#samples'),
                          'Value'    =  c( length(selected_genes()), 
                                           nrow(mRNA_heatmap_compute_results$filtered_mRNANormCounts),
                                           as.integer(ncol(mRNA_heatmap_compute_results$filtered_mRNANormCounts)-3))
    )
  })

  #prepare data for download
  output$download_mRNAData <- downloadHandler(
    filename = function() { paste('PCBC_geneExpr_data.csv')},
    content  = function(file){
      #ordering the rows based on the clustering order as determined by heatmap clustering
#       row_order =  mRNA_heatmap_compute_results$results[[1]]$order
#       col_order =  mRNA_heatmap_compute_results$results[[2]]$order
#       df = subset(mRNA_NormCounts, symbol %in% selected_genes())
      #reorder rows based on clustering
#       df = df[row_order,]
      #reorder columns based on clustering
      #just reodering the cols after first three cols of annotation
#       ordered_cols_df <- df[,c(4:ncol(df))][,col_order]
#       df <- cbind(df[,c(1:3)], ordered_cols_df)
      write.csv(get_filtered_mRNA_matrix(),file,row.names=T, col.names=T)
    })

  #prepare data for download
  output$download_miRNAData <- downloadHandler(
    filename = function() { paste('PCBC_microRNAExpr_data.csv')},
    content  = function(file){
      #get the microRNA expression matrix
      filtered_microRNA_NormCounts <- miRNA_normCounts[row.names(miRNA_normCounts) %in% selected_miRNAs(),]
      #subset on sample names based on user selected filters 
      filtered_metadata <- get_filtered_metadata(input,combined_metadata)
      filtered_microRNA_NormCounts <- filtered_microRNA_NormCounts[ ,colnames(filtered_microRNA_NormCounts) %in% filtered_metadata$Sample]
      write.csv(filtered_microRNA_NormCounts,file,row.names=T, col.names=T)
    })

  #prepare data for download
  output$download_methylationData <- downloadHandler(
    filename = function() { paste('PCBC_methylation_data.csv')},
    content  = function(file){
      #get the methylation  matrix
      filtered_meth_data <- meth_data[row.names(meth_data) %in% selected_methProbes(),]
      #subset on sample names based on user selected filters 
      filtered_metadata <- get_filtered_metadata(input,combined_metadata)
      filtered_meth_data <- filtered_meth_data[ ,colnames(filtered_meth_data) %in% filtered_metadata$Sample]
      write.csv(filtered_meth_data,file,row.names=T, col.names=T)
    })


  output$microRNA_summary <- renderTable({
    summary <- data.frame('Category' =  c('#Uniq genes in current list/pathway', 
                                          '#Uniq miRNAs targetting these genes',
                                          '#Uniq miRNAs(with expression values) targetting in these genes',
                                          '#samples',
                                          'overall #uniq miRNAs with matching ensembl geneId'),
                          'Value'    =  c( length(selected_genes()), 
                                           microRNA_heatmap_compute_results$num_miRNA,
                                           nrow(microRNA_heatmap_compute_results$filtered_microRNANormCounts),
                                           as.integer(ncol(microRNA_heatmap_compute_results$filtered_microRNANormCounts)),
                                           length(unique(miRNA_to_genes$Pathway)))
    )
  })
  
  output$meth_data_notes <- reactive({global_meth_data_notes})
  output$mRNA_data_notes <- reactive({global_mRNA_data_notes})
  output$miRNA_data_notes <- reactive({global_miRNA_data_notes})

  output$topgene_linkOut <- reactive({
    prefix =  '<form action="https://toppgene.cchmc.org/CheckInput.action" method="post" target="_blank" display="inline">
    <input type="hidden" name="query" value="TOPPFUN">
    <input type="hidden" id="type" name="type" value="HGNC">
    <input type="hidden" name="training_set" id="training_set" value="'
    suffix = '">
    <input type="Submit" class="btn shiny-download-link shiny-bound-output", value="Enrichment Analysis in ToppGene">
    </form>'
    geneIds = rownames(get_filtered_mRNA_matrix())
    geneIds = convert_to_HUGOIds(geneIds)
    geneIds <- paste(geneIds,collapse=" ")
    #generate the HTML content
    htmlContent <- paste(c(prefix,geneIds,suffix), collapse="")
  })
  
  #reactive value to store precomputed shiny results of mRNA data
  mRNA_heatmap_compute_results <- reactiveValues() 
  

  mRNA_cache_time <- reactiveValues()
  output$mRNA_cache_time = renderPrint({
    print(mRNA_cache_time$time)
  })
  
  output$microRNA_compute_time = renderPrint({
    print(microRNA_heatmap_compute_results$time)
  })
  
  #reactive value to store precomputed shiny results of mRNA data
  microRNA_heatmap_compute_results <- reactiveValues()
  
  
  
})

#######
# TEST CODE
#######

#   #create summary table

# #gene list to display
# output$selected_genes <- renderPrint({
#   selected_genes <- selected_geneNormCounts()
#   selected_genes <- as.character(selected_genes$symbol)
#   print(selected_genes,quote=FALSE)
# })




#   get_matrix <- reactive({
#     # get the filtered geneExp counts 
#     m <- selected_geneNormCounts()
#     #add the row names
#     #PURE HACK : since many ensembly IDs have same gene names 
#     # and rownames(matrix) cant have duplicates
#     # forcing the heatmap to render explicity passed rownames
#     #rownames(m) <- m$gene_id
#     explicit_rownames <- as.vector(m$symbol)
#     #convert to matrix
#     m <- as.matrix(m)
#     # eliminate the first 3 cols to get rid of the annotation and convert to matrix
#     m <- m[,4:ncol(m)]
#     
#     m <- apply(m,2,as.numeric)
#     
#     #removing those genes which dont vary much across the samples
#     # so any gene with SD < .2 across the samples will be dropped 
#     drop_genes <- which(apply(m,1,sd) < .2)
#     #following step to remove the bug seen 
#     #when m <-  m[-drop_genes,] is done directly and length(drop_genes) = 0
#     if(length(drop_genes) != 0){
#       m <-  m[-drop_genes,]  #filtering a mat , IMP
#       #also remove the same from the explicit rownames as those genes are taken out in anycase
#       explicit_rownames <- explicit_rownames[-drop_genes] #filtering a vector no , needed
#     }
#     mat.scaled <- t(scale(t(m))) 
#   })

#testing interactive shiny heatmap
#   output$test_heatmap <- renderHeatmap(
#     get_matrix()
#   )

#    output$test <- renderText({
#      selected_genes()
# #      #print(paste( "samples:", length(selected_samples()) , sep=": "))
# #      paste( "genecounts dim:" , dim(selected_geneNormCounts()))
#    })

#function to render a dynamic dropdown on the UI
#    output$enrichedPathways <- renderUI({
#      enriched_Pathways = get_enrichedPathways()
#      selectInput("enrichedPathways",
#                  sprintf("Enriched Pathways: %d", sum(! enriched_Pathways %in% c('NA','ALL'))), 
#                  choices = sort(enriched_Pathways)
#      )
#    })


# output$mRNA_cached_heatMap <- renderImage({
#   #a.) subset based on genes found in a pathway or user defined list
#   filtered_mRNANormCounts <- subset(mRNA_NormCounts, symbol %in% selected_genes())
#   #b.) subset on sample names based on user selected filters + rebind the gene names (first 3 cols)
#   filtered_mRNA_metata <- get_filtered_metadata(input,mRNA_metadata)
#   filtered_mRNA_samples <- filtered_mRNA_metata$bamName
#   filtered_mRNANormCounts <- cbind( filtered_mRNANormCounts[,1:3],
#                                     filtered_mRNANormCounts[, names(filtered_mRNANormCounts) %in% filtered_mRNA_samples  ])
#   m <- filtered_mRNANormCounts
#   #add the row names
#   #PURE HACK : since many ensembly IDs have same gene names 
#   # and rownames(matrix) cant have duplicates
#   # forcing the heatmap to render explicity passed rownames
#   #rownames(m) <- m$gene_id
#   explicit_rownames <- as.vector(m$symbol)
#   #convert to matrix
#   m <- as.matrix(m, drop=FALSE)
#   # eliminate the first 3 cols to get rid of the annotation and convert to matrix
#   m <- m[,4:ncol(m)]
#   annotation <- get_filtered_genesAnnotation(input,filtered_mRNA_metata)
#   #create a md5 of matrix and annotation
#   md5=digest(c(m,annotation), algo='md5')
#   plot_file = paste0(cache_dir,'/',md5,'.png')
#   start_time = proc.time()
#   if ( ! file.exists(plot_file) ){
#     png(plot_file)
#     #png(plot_file,width=24, height=16, units="in",res=300)
#     mRNA_heatmap_compute_results$results <- get_geneExpression_heatMap(m,annotation,explicit_rownames = explicit_rownames)
#     dev.off()
#   }
#   mRNA_cache_time$time = proc.time() - start_time
#   list(src= plot_file)
# },deleteFile=FALSE)
# })
#   

  
 


