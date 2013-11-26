#load the modules
library("shiny")


#load custom heatmap source code
source("./geneExpression_heatMap.R")


#Define the server the logic
shinyServer(function(input,output,session){
  
  #to check when the plot is rendered
  plotRendered <- reactiveValues(flag=FALSE)
  
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
    if ( length(input$enrichedPathways) == 0 ){ 
    }
    if(input$enrichedPathways == 'NA' ){
      user_submitted_gene_set()
    }
    else if(input$enrichedPathways == 'ALL'){
      #get union of all the genes in the enriched pathways for this geneList
      enriched_pathways <- gsub('#p.adj_.*','', get_enrichedPathways())
      genes <- lapply(enriched_pathways,function(x) {MSigDB$C2.CP.KEGG[[x]]})
      Reduce(union,genes)
    }
    else{
      #1. get a list of all genes in the selected enriched pathway
      #trimming the suffix : #pdj-
      pathway = gsub('#p.adj_.*','',input$enrichedPathways)
      MSigDB$C2.CP.KEGG[[pathway]]
    }
  })
  
  
  #########
  #get list of pathways enriched in the geneList selected by the user
  #########
  get_enrichedPathways <- reactive({
    if(input$user_selected_geneList == "Custom gene list"){
      'NA'
    }
    else{
      #return the enriched pathway for a gene list
      #labels contain the pvalue of the FET test       
      precomputed_enrichedPathways_in_geneLists[[input$user_selected_geneList]]
    }
  })
  
  
  #update the enriched pathways for the user selected genelist
  observe({
    enriched_Pathways = sort(get_enrichedPathways())
    updateSelectInput(session = session,
                      inputId = "enrichedPathways",
                      label = sprintf('Enriched pathway/s: %d', sum(! enriched_Pathways %in% c('NA','ALL'))),
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
    #b.) subset on sample names based on user selected filters + rebind the gene names (first 3 cols)
    filtered_geneNormCounts <- cbind( filtered_geneNormCounts[,1:3],
                                      filtered_geneNormCounts[, names(filtered_geneNormCounts) %in% filtered_sample_names()  ])
  })
  
  
  #get the list of user submitted genes
  user_submitted_gene_set <- reactive({
    temp_geneL <- unlist(strsplit(input$custom_gene_list,split=c('[\\s+,\\n+\\r+)]'),perl=T))
    #conevert everything to upper case
    temp_geneL <- toupper(temp_geneL)
  })
  
  #create the annotation data frame for the heatmap
  get_filtered_genesAnnotation <- reactive({
    filtered_metadata <- subset(get_filtered_metadata(), bamName %in% filtered_sample_names())
    annotation_cols <- heatmap_annotation_cols[input$heatmap_annotation_labels]
    annotation <- as.data.frame(filtered_metadata[,annotation_cols])
    names(annotation) <- input$heatmap_annotation_labels
    #assign the sample names to row names so that the heatmap function could use them for labelling
    rownames(annotation) <- filtered_metadata$bamName
    annotation
  })
  
  
  #    #return the heatMap plot
  #    output$heatMap <- renderPlot({  
  #      plotRendered$set <- FALSE      
  #      # get the filtered geneExp counts 
  #      m <- as.matrix(selected_geneNormCounts())
  #      # eliminate the first 3 cols to get rid of the annotation and convert to matrix
  #      m <- m[,4:ncol(m)]
  #      annotation <- get_filtered_genesAnnotation()
  #      #plot the heatmap
  #      get_geneExpression_heatMap(m,annotation)  #available in geneExpression_heatMap.R
  #      plotRendered$set <- TRUE
  #    })
  
  
  #   #function to display if a plot is being rendered
  #   output$msg <- renderText({
  #     if(plotRendered$set == TRUE){
  #       ''
  #     }
  #     else if(plotRendered$set == FALSE){
  #       'rendering...'
  #     }
  #   })
  
  
  #   #create summary table
  #   output$summary <- renderTable({     
  #      summary <- data.frame('Category' =  c('#genes_currentPathway', '#genes_found_in_dataset'),
  #                            'Value'    =  c( length(selected_genes()), nrow(selected_geneNormCounts()))
  #                           )
  #    })
  
  
  
  #   #####
  #   #create final geneexpression data table
  #   #####
  #   output$geneExpTable <- renderDataTable({     
  #     selected_geneNormCounts()
  #     })
  
  
  
  #    #prepare data for download
  #    output$downloadData <- downloadHandler(
  #      filename = function() { paste('PCBC_geneExpr_data_for_',gsub(' ','_',input$user_selected_pathway),'_pathway.csv')},
  #      content  = function(file){
  #        write.csv(selected_geneNormCounts(),file)
  #      })
  
  
  #######
  # TEST
  #######
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