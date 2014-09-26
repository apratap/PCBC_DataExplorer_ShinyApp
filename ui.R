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
                 helpText("Accepts HUGO/Ensembl/Entrez gene id's separated either by comma, space, line"),
                 tags$textarea(id="custom_gene_list",rows=8,cols=200,paste0(sample_gene_list, collapse=', ')),
                 checkboxInput('incl_corr_genes', 'also include correlated genes', value = FALSE),
                 sliderInput('corr_threshold', label='Correlation Threshold', min=0.5, max=1.0, value=0.9, step=0.05),
                 br(),
                 h5('1.b. Search on a custom miRNA list:'),
                 helpText("Accepts mirbase ids separated by comma, space,line, comma "),
                 tags$textarea(id="custom_miRNA_list",rows=4,cols=200),
                 br(),
                 value='custom_gene_list'
        ), #END TAB PANEL 1
        # TAB PANEL 2 : select a pathway
        tabPanel(
                'Pathways',
                selectInput("selected_pathways",
                            "Pathways",
                            choices = names(pathways_list), selectize=T, multiple=T, width='400px',
                            selected = names(pathways_list)[c(1:2)]),
                br(), br(), br(), br(), br(), br(),
                value='pathway'
        ),  #END TAB PANEL 2
        #TAB PANEL 3 : precomputed sig gene list
        tabPanel(
          'Sig gene lists',
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
  
      br(),
      #heatmap annotation labels
      selectInput('heatmap_annotation_labels', h4('2. Color heatmap by:'),
                   choices  = colnames(combined_metadata),selected='Differentiation_State'),
      br(),
      
      
      
      #FILTER OPTIONS
      h4('3. Filter samples by:'),
      #1. filter based on mod_linetype
      selectInput('linetype', h5('Line type'), choices=unique(combined_metadata$Line_Type),
                  selectize=T, multiple=T, selected=c('ESC')),
      br(),
      selectInput('gene_combination', h5('Reprogramming Gene Combination'), choices=unique(combined_metadata$Reprogramming_Gene_Combination),
                  selectize=T, multiple=T),
      br(),
      selectInput('vector_type', h5('Reprogramming Vector Type'), choices=unique(combined_metadata$Reprogramming_Vector_Type),
                  selectize=T, multiple=T),
      br(),
      selectInput('tissue_origin', h5('Tissue of Origin'), choices=unique(combined_metadata$Tissue_of_Origin),
                  selectize=T, multiple=T),
      br(),
      selectInput('diff_state', h5('Differentiation State'), choices=unique(combined_metadata$Differentiation_State),
                  selectize=T, multiple=T),
      br(),
      selectInput('cell_origin', h5('Cell Type of Origin'), choices=unique(combined_metadata$Cell_Type_of_Origin),
                  selectize=T, multiple=T),
      br(),
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
                   br(),
                   br(),
                   h5('summary'),
                   tableOutput("mRNA_summary"),
                   br(),
                   h5('compute time'),
                   verbatimTextOutput("mRNA_compute_time")
                   ),
#           tabPanel("mRNA (cached)",
#                    plotOutput("mRNA_cached_heatMap", height="auto", width="700px", hoverId=NULL ),
#                    verbatimTextOutput("mRNA_cache_time")
#           ),
          tabPanel("microRNA",
                   plotOutput("microRNA_heatMap",height="700px",width="auto",hoverId=NULL),
                   br(),
                   br(),
                   h5('summary'),
                   #tableOutput("microRNA_summary"),
                   br(),
                   h5('compute time'),
                   verbatimTextOutput('microRNA_compute_time')
          ),
          tabPanel("Explore Data", 
                   htmlOutput("topgene_linkOut"),
                   downloadButton('downloadData','Download Expression Data'),
                   br(), br(), br(),
                   dataTableOutput("geneExpTable"))
          ) #END tabset panel
    )# END mainPanel 
) #END sidebarLayout
) #END fluidpage
) #END shinyUI
