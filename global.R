library("synapseClient")
library("gdata")
library("plyr")
library("org.Hs.eg.db")
library("shiny")



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



# #login to synapse
synapseLogin()


## WORKAROUND : mainly to avoid creating on the fly as org.Hs.eg.db on shiny server is old
#precomputed in precompute.R
hg19_gene_annot <- readRDS("precomputed_hg19_gene_annot.RDS")


###
#get the MsigDB object
###
cat('Reading the MSIGDB object from synapse...')
MSIGDB_syn<-synGet("syn2227979")
load(MSIGDB_syn@filePath) #available as MSigDB R object
cat('..Done\n\n')


###
#get the PCBC samples gene normalized counts
###
cat('Reading the PCBC normalized geneExp data from Synapse')
syn_geneNormCounts <- synGet('syn2247799')
#read in the file
geneNormCounts <- read.table(syn_geneNormCounts@filePath,header=T,sep='\t')
#keep only uniq gene names
#might need to be improved
geneNormCounts <- geneNormCounts[!duplicated(geneNormCounts$symbol),]
cat('..Done\n\n')

######
#get the list siginificant genes from comparative analysis in synapse
#####
cat('Reading the precomputed significant genelist')
sigGenes_lists <- readRDS("precomputed_sigGenes_lists.rds")
cat('..Done\n\n')

###
#get the metadata from the synapse for PCBC samples
###
METADATA_ID <- 'syn2248030'
query <- sprintf('select * from entity where parentId=="%s"', METADATA_ID)
metadata <- synQuery(query)
cols_to_be_deleted = c('entity.benefactorId', 'entity.concreteType', 'entity.createdByPrincipalId', 
                       'entity.createdOn', 'entity.createdByPrincipalId', 'entity.id', 
                       'entity.modifiedOn', 'entity.modifiedByPrincipalId', 'entity.noteType', 
                       'entity.versionLabel', 'entity.versionComment', 'entity.versionNumber', 
                       'entity.parentId', 'entity.description', 'entity.eTag')
#delete the unwanted cols
metadata <- metadata[,!names(metadata) %in% cols_to_be_deleted]

#remove the prefix 'entity.' from the df col names
names(metadata) <- gsub('entity.','',names(metadata))


########
#create new higher order labels for heatmap labels
#groups are merged into higher order
########
#1. linetype
# merge everthing other than hESC , iPSC, None into tissue
metadata$mod_linetype <-  metadata$linetype
index_to_be_replaced <-  ! metadata$linetype %in% c('hESC','iPSC','None')
metadata$mod_linetype[index_to_be_replaced] = 'tissue'


#2. diffnameshort
metadata$mod_diffnameshort <- metadata$diffnameshort
#2A. EB-LF should be renamed EB and similarly SC-LF should be SC
metadata$mod_diffnameshort <- sub('SC-.*','SC',metadata$mod_diffnameshort)
metadata$mod_diffnameshort <- sub('EB-.*','EB',metadata$mod_diffnameshort)
#2B. anything except the following should be grouped into Other
# SC, EB, DE, ECTO, MESO, MESO-30, MESO-15
index_to_be_replaced <- ! metadata$mod_diffnameshort %in% c('SC', 'EB', 'DE', 'ECTO', 'MESO', 'MESO-30', 'MESO-15')
metadata$mod_diffnameshort[index_to_be_replaced] = 'Others'

#3. origcell
metadata$mod_origcell <- metadata$origcell
#3A. categorize UCB CD34+ into CD34+
metadata$mod_origcell <- sub('^.*CD34\\+$','CD34+',metadata$mod_origcell,perl=T)
#3B. categ bone marrow to BM MSC
metadata$mod_origcell <- sub('bone marrow','BM MSC',metadata$mod_origcell,perl=T)
#3C. adultHeart to tissue
metadata$mod_origcell <- sub('adultHeart','tissue',metadata$mod_origcell,perl=T)
#3C. categ the following into Others
Others_group <- c('embryo', 'primEndothel', 'skin', 'None')
index_to_be_replaced <- metadata$mod_origcell %in% Others_group
metadata$mod_origcell[index_to_be_replaced] <- 'Others'



##########
#create possible heatmap annotation labels
##########
heatmap_annotation_cols <- c('Gender'='donorsex','line type'='mod_linetype',
                             'Diff Name'='mod_diffnameshort','Cell Origin'='mod_origcell',
                             'Induction Genes' = 'inductiongenes')


#1. sex
sex <- unique(metadata$donorsex)
sex <- sex[sex != "None"]




###
#sample gene list of the user input area
###
sample_gene_list <- c("GP1BA","GP1BB", "EPO", "CD33", "TNF", "GP9", "ITGAM", "CD34",      
                      "CD36", "GP5", "ITGA4", "ITGA3", "KITLG", "ITGA2B",
                      "CCL18", "CCL3", "CACNG3", "AP2B1" )
sample_gene_list <- paste(sample_gene_list,",")


sigGenes_lists[["Mesoderm Day 5_vs_Endoderm Differentiated Cells"]]

#########
#read the precomputed enriched pathway list
########
df_precomputed_enrichedPathways_in_geneLists = readRDS("precomputed_enrichedPathways_in_geneLists.rds")
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




