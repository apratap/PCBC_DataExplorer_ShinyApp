library("synapseClient")
library("gdata")
library("plyr")
library("org.Hs.eg.db")
library("shiny")
library("digest")

#load the memoised version of pheatmap
source("geneExpression_heatMap.R")

# #load the shiny based d3 app
# if (!require("devtools"))
#   install.packages("devtools")
# if(!require("heatmap"))
#   devtools::install_github("d3-heatmap", "jcheng5")
#library("heatmap")


#load the external files
# available through shiny
includeScript('css/tooltip.css')

#login to synapse
synapseLogin()


#create a dir to store plots if it doesnt exist
cache_dir <- ".plotcache"
dir.create(cache_dir, showWarnings = FALSE)
cache_dir <- normalizePath('.plotcache')
plot_cache_lookup <- list()

## WORKAROUND : mainly to avoid creating on the fly as org.Hs.eg.db on shiny server is old
#precomputed in precompute.R
hg19_gene_annot <- readRDS("precomputed_data/precomputed_hg19_gene_annot.RDS")


###
#get the MsigDB object
###
cat('Reading the MSIGDB object from synapse...')
MSIGDB_syn<-synGet("syn2227979")
load(MSIGDB_syn@filePath) #available as MSigDB R object
cat('..Done\n\n')


#####
#get the mRNA expression data
#####
source("mRNA_data_prep.R")


#####
#get the miRNA expression data
#####
source("miRNA_data_prep.R")




#filter metadata
get_filtered_metadata <- function(input, metadata){
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
}


#create the annotation data frame for the heatmap
get_filtered_genesAnnotation <- function(input,metadata){
  if(length(input$heatmap_annotation_labels) == 0){
    stop('please select atleast one heatmap annotation variable \n\n')      
  }
  else{
    annotation_cols <- heatmap_annotation_cols[input$heatmap_annotation_labels]
    annotation <- as.data.frame(metadata[,annotation_cols])
    names(annotation) <- input$heatmap_annotation_labels
    #assign the sample names to row names so that the heatmap function could use them for labelling
    rownames(annotation) <- metadata$bamName
    annotation
  }
}




#create heatmap annotation labels
heatmap_annotation_cols <- c('Gender'='donorsex','line type'='mod_linetype',
                             'Diff Name'='mod_diffnameshort','Cell Origin'='mod_origcell',
                             'Induction Genes' = 'inductiongenes')


#1. sex
sex <- unique(mRNA_metadata$donorsex)
sex <- sex[sex != "None"]



#sample gene list of the user input area
df <- read.table("precomputed_data/pre_selected_genelist.txt",sep="\t")
sample_gene_list <- as.character(unique(df$V5))

#########
#read the precomputed enriched pathway list
########
df_precomputed_enrichedPathways_in_geneLists = readRDS("precomputed_data/precomputed_enrichedPathways_in_geneLists.rds")
df_precomputed_enrichedPathways_in_geneLists$pathways_with_pvalue =  paste(df_precomputed_enrichedPathways_in_geneLists$pathways,
                                                                           '#p.adj_',
                                                                           format.pval(df_precomputed_enrichedPathways_in_geneLists$p.adj,digits=2),
                                                                           sep='')

#creating a list of list 
precomputed_enrichedPathways_in_geneLists = split(df_precomputed_enrichedPathways_in_geneLists$pathways_with_pvalue,
                                                  df_precomputed_enrichedPathways_in_geneLists$significant_gene_list_name)


#HACK
#For each geneList add another PATHWAY TYPE "ALL" which indicates use all the pathways for the shiny SERVER/UI
# in this case genes in all the enriched pathways will be shown on the heatmap
precomputed_enrichedPathways_in_geneLists <- lapply(precomputed_enrichedPathways_in_geneLists,function(x) { x[length(x)+1] = 'ALL'; x})



#add place holder for custom gene list
#precomputed_enrichedPathways_in_geneLists[['Custom gene list']] = 'NA'




