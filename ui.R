library(shinyIncubator)

#main UI code
shinyUI( fluidPage(
  
  sidebarLayout(
    ##################################
    #SIDE BAR PANEL FOR USER OPTIONS
    ##################################
    sidebarPanel(
      progressInit(),
      h4('1. Select a gene list'),
      tabsetPanel(
        id = 'genelist_type',
        #TAB PANEL 1 : custom gene list
        tabPanel('My Search',h5('1.a. Search on a custom gene list:'),
                 helpText("Accepts HUGO/Ensembl/Entrez gene ids"),
                 tags$textarea(id="custom_gene_list",rows=8,cols=200,paste0(sample_gene_list, collapse=', ')),
                 checkboxInput('incl_corr_genes', 'also include correlated genes', value = FALSE),
                 sliderInput('corr_threshold', label='Correlation Threshold', min=0.5, max=1.0, value=0.9, step=0.05),
                 br(),
                 h5('1.b. Add miRNA Targets:'),
                 helpText("Accepts mirbase ids"),
                 tags$textarea(id="custom_miRNA_list",rows=4,cols=200),
                 br(),
                 actionButton("custom_search", h4("Search")),
                 br(),
                 value='custom_gene_list'
        ), #END TAB PANEL 1
        # TAB PANEL 2 : select a pathway
        tabPanel(
          'Pathways',
          selectInput("selected_pathways",
                      "1.a. Select Pathway/s",
                      choices = names(pathways_list), selectize=T, multiple=T, width='400px',
                      selected = names(pathways_list)[c(1:2)]),
          br(), br(), br(), br(),
          value='pathway'
        ),  #END TAB PANEL 2
        #TAB PANEL 3 : precomputed sig gene list
        tabPanel(
          'Significant Genes',
          #tags$div(title="Significantly enriched gene lists as a result of pairwise comparison of all PCBC samples",
                   selectInput( "selected_Significant_GeneList",
                                "Precomputed Significant gene lists (?)", selectize=FALSE, 
                                choices = sort(names(precomputed_enrichedPathways_in_geneLists)) #loaded from getDATA.R
                              ),
          #dynamically updated with TOOL TIP
          tags$div(title='Enriched KEGG pathways in the selected gene list based on FET test (padj <.05)',
                   selectInput('enrichedPathways',
                               'Enriched Pathways (?)',
                               choices='ALL',
                               selectize=FALSE
                               )
                   ),
          value = 'precomputed_significant_geneList'
        ) #END TAB PANEL 3
      ),#END TABSET 
  
      br(),br(),
      #heatmap annotation labels
      selectInput('heatmap_annotation_labels', h4('2. Annotate Samples by:'),
                  choices  = colnames(combined_metadata)[-1],  #-1 to remove the first value "Sample"
                  selected='Differentiation_State'),
      br(),
      
      #FILTER OPTIONS
      h4('3. Filter samples by:'),
      #1. filter based on mod_linetype
      selectInput('linetype', h5('Line type'), choices=unique(combined_metadata$Line_Type),
                  selectize=T, multiple=T, selected=c('ESC','iPSC')),
      selectInput('gene_combination', h5('Reprogramming Gene Combination'), choices=unique(combined_metadata$Reprogramming_Gene_Combination),
                  selectize=T, multiple=T),
      selectInput('vector_type', h5('Reprogramming Vector Type'), choices=unique(combined_metadata$Reprogramming_Vector_Type),
                  selectize=T, multiple=T),
      selectInput('tissue_origin', h5('Tissue of Origin'), choices=unique(combined_metadata$Tissue_of_Origin),
                  selectize=T, multiple=T),
      selectInput('diff_state', h5('Differentiation State'), choices=unique(combined_metadata$Differentiation_State),
                  selectize=T, multiple=T),
      selectInput('cell_origin', h5('Cell Type of Origin'), choices=unique(combined_metadata$Cell_Type_of_Origin),
                  selectize=T, multiple=T),
      selectInput('gender', h5('Gender'), choices=unique(combined_metadata$Gender),
                  selectize=T, multiple=T),
      
      br(),
      h4('4. Heatmap Settings:'),
      #distance metric
      selectInput("clustering_distance","Distance Calculation",
                  choices=c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                  selectize=T, multiple=F, selected="correlation"),
      #linkage 
      selectInput("clustering_method","Clustering Method",
                  choices=c("ward", "single", "complete", "average", "mcquitty", "median", "centroid"),
                  selectize=T, multiple=F, selected="complete")
    ), # END sidebarpanel

  
    #####################
    #Main shiny panel
    #####################
    mainPanel(
        tabsetPanel(
          tabPanel("mRNA", 
                   plotOutput("mRNA_heatMap",height="700px",width="auto",hoverId=NULL),
                   br(),br(),br(),br(),
                   htmlOutput("topgene_linkOut"),
                   downloadButton('download_mRNAData','Download mRNA expression data'),
                   br(),br(),
                   htmlOutput('mRNA_data_notes')
          ),
          tabPanel("microRNA",
                   plotOutput("microRNA_heatMap",height="700px",width="auto",hoverId=NULL),
                   br(),br(),br(),br(),
                   downloadButton('download_miRNAData','Download microRNA expression data'),
                   br(),br(),
                   htmlOutput('miRNA_data_notes')
          ),
          tabPanel("methylation",
                   plotOutput("methylation_heatMap",height="700px",width="auto",hoverId=NULL),
                   br(),br(),br(),br(), 
                   downloadButton('download_methylationData','Download methylation data'),
                   br(),br(),
                   htmlOutput('meth_data_notes')
          )
      ) #END tabset panel
    )# END mainPanel 
) #END sidebarLayout
) #END fluidpage
) #END shinyUI
