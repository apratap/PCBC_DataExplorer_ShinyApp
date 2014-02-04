#Define the server the logic
shinyServer(function(input,output,session){
  
  error <- reactiveValues( state='FALSE',
                           message='None')
  
  #Default Error state
  output$error <- renderText({error$state})
  
  #select the list of samples to show
  #Based on USER Selected filter options
  #directly obtained from filtered meta data
  filtered_sample_names <- reactive({
    get_filtered_metadata()$bamName
  })
  
  #filter the metadata based on the user selected options
  get_filtered_metadata <- reactive({
    filtered_metadata <- metadata
    #1. Filter based on user selected line type
    if( length(input$mod_linetype) != 0 ){
      filtered_metadata <- subset(filtered_metadata, mod_linetype %in% input$mod_linetype)
    }
    #2. filter based on user selected short differentiation name
    if( length(input$mod_diffnameshort) != 0 ){
      filtered_metadata <- subset(filtered_metadata, mod_diffnameshort %in% input$mod_diffnameshort)  
    }
    #3. filter based on cell origin
    if(length(input$cell_origin) != 0){
      filtered_metadata <- subset(filtered_metadata, mod_origcell %in% input$cell_origin)
    }
    #4. filter based on induction genes
    if(length(input$induction_genes) != 0){
      filtered_metadata <- subset(filtered_metadata, inductiongenes %in% input$induction_genes)
    }
    # converting to match the  sample names in the geneNorm counts matrix
    filtered_metadata$bamName <- gsub('-','.',as.vector(filtered_metadata$decoratedName))
    filtered_metadata
  })
  
  #####
  #get list of genes in current pathway or user entered list
  #####
  selected_genes <- reactive({
    if( input$genelist_type == 'custom_gene_list'  ){
      genes <- unique(user_submitted_geneList())
    }
    else if( input$genelist_type == 'precomputed_significant_geneList'){
      if(input$enrichedPathways == 'ALL'){
        genes_in_selected_GeneList <- sigGenes_lists[[input$selected_Significant_GeneList]]
        genes <- unique(genes_in_selected_GeneList)
      }
      else{
        #1. get a list of all genes in the selected enriched pathway
        #trimming the suffix : #pdj-
        pathway = gsub('#p.adj_.*','',input$enrichedPathways)
        genes_in_pathway <- MSigDB$C2.CP.KEGG[[pathway]]
        genes_in_selected_GeneList <- sigGenes_lists[[input$selected_Significant_GeneList]]
        genes <- intersect(genes_in_pathway, genes_in_selected_GeneList)
      }
    }
    else if( input$genelist_type == 'pathway'){
      pathway = input$selected_pathway
      #get the list of genes in the KEGG pathway
      genes <- MSigDB$C2.CP.KEGG[[pathway]]
    }
    if(length(genes) <= 2){
      error_msg = sprintf('too few genes %d to cluster', length(genes))
      stop(error_msg)
    }
    else{
      genes
    }
  })
  
  
  
    
  #########
  #get list of pathways enriched in the geneList selected by the user
  #########
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
  

  ###
  #subset the gene norm counts to select only the genes found 
  ###
  selected_geneNormCounts <- reactive({
    #a.) subset based on genes found in a pathway or user defined list
    filtered_geneNormCounts <- subset(geneNormCounts, symbol %in% selected_genes())
    #print(dim(filtered_geneNormCounts))
    #b.) subset on sample names based on user selected filters + rebind the gene names (first 3 cols)
    filtered_geneNormCounts <- cbind( filtered_geneNormCounts[,1:3],
                                      filtered_geneNormCounts[, names(filtered_geneNormCounts) %in% filtered_sample_names()  ])
  })
  
  
  #get the list of user submitted genes
  user_submitted_geneList <- reactive({
    geneList <- unlist(strsplit(input$custom_gene_list,split=c('[\\s+,\\n+\\r+)]'),perl=T))
    #conevert everything to upper case
    geneList <- toupper(geneList)
    geneList <- geneList[ !geneList == "" ] #remove the blank entries
    geneList
  })
    
  #create the annotation data frame for the heatmap
  get_filtered_genesAnnotation <- reactive({
    if(length(input$heatmap_annotation_labels) == 0){
      stop('please select atleast one heatmap annotation variable \n\n')      
    }
    else{
      filtered_metadata <- subset(get_filtered_metadata(), bamName %in% filtered_sample_names())
      annotation_cols <- heatmap_annotation_cols[input$heatmap_annotation_labels]
      annotation <- as.data.frame(filtered_metadata[,annotation_cols])
      names(annotation) <- input$heatmap_annotation_labels
      #assign the sample names to row names so that the heatmap function could use them for labelling
      rownames(annotation) <- filtered_metadata$bamName
      annotation
    }
  })

  #reactive value to store precomputed shiny results
  heatmap_compute <- reactiveValues() 
  
  
   #return the heatMap plot
   output$heatMap <- renderPlot({  
     # get the filtered geneExp counts 
     m <- selected_geneNormCounts()
     #add the row names
     #PURE HACK : since many ensembly IDs have same gene names 
     # and rownames(matrix) cant have duplicates
     # forcing the heatmap to render explicity passed rownames
     #rownames(m) <- m$gene_id
     explicit_rownames <- as.vector(m$symbol)
     #convert to matrix
     m <- as.matrix(m)
     # eliminate the first 3 cols to get rid of the annotation and convert to matrix
     m <- m[,4:ncol(m)]
     #m <- apply(m,2,as.numeric)
     #m <- log2(m+1)
     annotation <- get_filtered_genesAnnotation()
     #plot the heatmap
     heatmap_compute$results <- get_geneExpression_heatMap(m,annotation,explicit_rownames = explicit_rownames)  #available in geneExpression_heatMap.R
     
   })
  
  
  #create summary table
  output$summary <- renderTable({
     
     summary <- data.frame('Category' =  c('#Uniq genes in current list/pathway', '#genes found with exp values', 
                                           '#samples'),
                           'Value'    =  c( length(selected_genes()), nrow(selected_geneNormCounts()),
                                            as.integer(ncol(selected_geneNormCounts())-3))
                          )
   })
  
  #####
  #create a table with selected gene list and merge with some annotation
  #####
  output$geneExpTable <- renderDataTable({     
    df <- merge(selected_geneNormCounts()[,1:3], hg19_gene_annot, by.x='symbol',by.y='SYMBOL')
    df
    })
  
   #prepare data for download
   output$downloadData <- downloadHandler(
     filename = function() { paste('PCBC_geneExpr_data.csv')},
     content  = function(file){
       #ordering the rows based on the clustering order as determined by heatmap clustering
       row_order =  heatmap_compute$results[[1]]$order
       col_order =  heatmap_compute$results[[2]]$order
       df = selected_geneNormCounts()
       #reorder rows based on clustering
       df = df[row_order,]
       #reorder columns based on clustering
       #just reodering the cols after first three cols of annotation
       ordered_cols_df <- df[,c(4:ncol(df))][,col_order]
       df <- cbind(df[,c(1:3)], ordered_cols_df)
       write.csv(df,file,row.names=FALSE)
     })
    
    #gene list to display
    output$selected_genes <- renderPrint({
      selected_genes <- selected_geneNormCounts()
      selected_genes <- as.character(selected_genes$symbol)
      print(selected_genes,quote=FALSE)
    })
  
  
  output$topgene_linkOut <- reactive({
    
    prefix =  '<form action="http://toppgene.cchmc.org/CheckInput.action" method="post" target="_blank" display="inline">
               <input type="hidden" name="query" value="TOPPFUN">
               <input type="hidden" id="type" name="type" value="HGNC">
               <input type="hidden" name="training_set" id="training_set" value="'
    suffix = '">
              <input type="Submit" class="btn shiny-download-link shiny-bound-output", value="Enrichment Analysis in ToppGene">
              </form>'
    genes <- paste(selected_genes(),collapse=" ")

    #generate the HTML content
    htmlContent <- paste(c(prefix,genes,suffix), collapse="")
  })

  
  #######
  # TEST
  #######
  
  get_matrix <- reactive({
    # get the filtered geneExp counts 
    m <- selected_geneNormCounts()
    #add the row names
    #PURE HACK : since many ensembly IDs have same gene names 
    # and rownames(matrix) cant have duplicates
    # forcing the heatmap to render explicity passed rownames
    #rownames(m) <- m$gene_id
    explicit_rownames <- as.vector(m$symbol)
    #convert to matrix
    m <- as.matrix(m)
    # eliminate the first 3 cols to get rid of the annotation and convert to matrix
    m <- m[,4:ncol(m)]
    
    m <- apply(m,2,as.numeric)
    
    
    #removing those genes which dont vary much across the samples
    # so any gene with SD < .2 across the samples will be dropped 
    drop_genes <- which(apply(m,1,sd) < .2)
    #following step to remove the bug seen 
    #when m <-  m[-drop_genes,] is done directly and length(drop_genes) = 0
    if(length(drop_genes) != 0){
      m <-  m[-drop_genes,]  #filtering a mat , IMP
      #also remove the same from the explicit rownames as those genes are taken out in anycase
      explicit_rownames <- explicit_rownames[-drop_genes] #filtering a vector no , needed
    }
    mat.scaled <- t(scale(t(m))) 
  })
  
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
  
})