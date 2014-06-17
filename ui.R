


#main UI code
shinyUI( fluidPage(
  
  #Application title
  #titlePanel(""),
  
  sidebarLayout(
    ##################################
    #SIDE BAR PANEL FOR USER OPTIONS
    ##################################
    sidebarPanel(
      h4('1. Select a gene list'),
      tabsetPanel(
        id = 'genelist_type',
        #TAB PANEL 1 : custom gene list
        tabPanel('My gene list',h5('Search on a custom gene list:'),
                 tags$textarea(id="custom_gene_list",rows=8,cols=200,paste0(sample_gene_list, collapse=', ')),
                 helpText("Accepts HUGO gene names. Gene names may be separated by comma, space,line, comma "),
                 br(),
                 value='custom_gene_list'
        ), #END TAB PANEL 1
        # TAB PANEL 2 : select a pathway
        tabPanel(
                'Pathways',
                selectInput("selected_pathway",
                            "Pathways",
                            choices = sort(names(MSigDB$C2.CP.KEGG)), # global.R
                            selectize=FALSE
                            ),     
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
      br(),
      #heatmap annotation labels
      checkboxGroupInput('heatmap_annotation_labels', h4('2. Color heatmap by:'),
                         choices  = names(heatmap_annotation_cols),
                         selected = c('Diff Name')),
      br(),
  
      #FILTER OPTIONS
      h4('3. Filter samples by:'),
      #1. filter based on mod_linetype
      checkboxGroupInput('mod_linetype', h5('Line type'),choices=sort(unique(metadata$mod_linetype) ) ),
      #2. filter based on diff_short_name
      checkboxGroupInput('mod_diffnameshort',h5('Differentiation name'),choices= sort(unique(metadata$mod_diffnameshort) ) ),
      #3. filter based on cell origin
      checkboxGroupInput('cell_origin',h5('Cell origin'),choices=sort(unique(metadata$mod_origcell) ) ),
      #4. filter based on induction genes
      checkboxGroupInput('induction_genes',h5('Induction genes used'),choices=sort(unique(metadata$inductiongenes) ) )
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
                   h5('compute time'),
                   verbatimTextOutput("test")
                   ),
          tabPanel("mRNA (cached)",
                   plotOutput("mRNA_cached_heatMap", height="auto", width="700px", hoverId=NULL ),
                   verbatimTextOutput("mRNA_cache_time")
          ),
          tabPanel("miRNA (test)",
                   plotOutput("microRNA_heatMap",height="700px",width="auto",hoverId=NULL),
                   tableOutput("microRNA_summary")
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
