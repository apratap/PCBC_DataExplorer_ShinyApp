#load the memoised version of pheatmap
options(stringsAsFactors = F)
options(warn=-1)
library(devtools)
library("synapseClient")
library("gdata")
library("shiny")
library("digest")
library("dplyr")
library("org.Hs.eg.db")

#source the heatmap code
source_url("https://raw.githubusercontent.com/apratap/apRs/master/expression_heatmap.R")

#source generic heatmap functions
source_url("https://raw.githubusercontent.com/apratap/apRs/master/generic_annotation_functions.R")

#login to synapse
synapseLogin()

## WORKAROUND : mainly to avoid creating on the fly as org.Hs.eg.db on shiny server is old
#precomputed in precompute.R
hg19_gene_annot <- readRDS("precomputed_data/precomputed_hg19_gene_annot.RDS")

#get the MsigDB object
cat('Reading the MSIGDB object from synapse...')
MSIGDB_syn<-synGet("syn2227979")
load(MSIGDB_syn@filePath) #available as MSigDB R object
pathways_list <- c(MSigDB$C2.CP.BIOCARTA, MSigDB$C2.CP.KEGG, MSigDB$C2.CP.REACTOME)
cat('..Done\n\n')

#get the mRNA expression data
source("mRNA_data_prep.R")

#get the miRNA expression data
source("miRNA_data_prep.R")

#get the methylation data
source("methylation_data_prep.R")

#get the global functions
source("global_functions.R")

#prepare single global metadata
column_names <- c('Sample', colnames(mRNA_metadata)[-1])
colnames(mRNA_metadata) <- c(1:8)
colnames(miRNA_metadata) <- c(1:8)
colnames(meth_metadata) <- c(1:8)
combined_metadata <- rbind(mRNA_metadata, miRNA_metadata, meth_metadata, deparse.level = 0)
colnames(combined_metadata) <- gsub('\\s+','_',column_names,perl=T)


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


#####################
#OLD code
#####################

# #create a dir to store plots if it doesnt exist
# cache_dir <- ".plotcache"
# dir.create(cache_dir, showWarnings = FALSE)
# cache_dir <- normalizePath('.plotcache')
# plot_cache_lookup <- list()

# #load the shiny based d3 app
# if (!require("devtools"))
#   install.packages("devtools")
# if(!require("heatmap"))
#   devtools::install_github("d3-heatmap", "jcheng5")
#library("heatmap")

#load the external files
# available through shiny
#includeScript('css/tooltip.css')



