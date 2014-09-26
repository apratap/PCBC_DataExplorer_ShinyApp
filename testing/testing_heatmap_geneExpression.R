source("http://bioconductor.org/biocLite.R")
biocLite("ggbio")


#load the modules
library(synapseClient)
library(gplots)
library(RColorBrewer)
library(ggbio)
library(preprocessCore)

#login to synapse
synapseLogin()


data("CRC", package = "biovizBase")

autoplot(hg19sub, layout="circle", fill="gray70")

data("CRC", package = "biovizBase")
p <- ggbio() + circle(hg19sub, geom = "ideo", fill = "gray70")
p  <- p + circle(hg19sub, geom = "scale", size = 2) 
p  <- p + circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
p

p <- ggbio(buffer=1) + circle(hg19sub, geom = "ideo", fill = "gray90")
p  <- p + circle(hg19sub, geom = "scale", size = 3) 
p  <- p + circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 4)
p
p
p

?layout_circle

?circle

# 
# 
# #get the metadata from the synapse about the above sampled
# #get the meta about PCBC samples
# METADATA_ID <- 'syn2024470'
# query <- sprintf('select * from entity where parentId=="%s"', METADATA_ID)
# metadata <- synQuery(query)
# dim(metadata)
# 
# cols_to_be_deleted = c('entity.benefactorId', 'entity.concreteType', 'entity.createdByPrincipalId', 
#                        'entity.createdOn', 'entity.createdByPrincipalId', 'entity.id', 
#                        'entity.modifiedOn', 'entity.modifiedByPrincipalId', 'entity.noteType', 
#                        'entity.versionLabel', 'entity.versionComment', 'entity.versionNumber', 
#                        'entity.parentId', 'entity.description', 'entity.eTag')
# metadata <- metadata[,!names(metadata) %in% cols_to_be_deleted]
# 
# #remove the prefix 'entity.' from the df col names
# names(metadata) <- gsub('entity.','',names(metadata))
# 
# 
# apply(metadata,2,unique)
# 
# #get the pathways from MSigDB
# #get the MsigDB object
# MSIGDB<-synGet("syn2227979")
# load(MSIGDB@filePath) #available as MSigDB R object
# 
# 
# 
# #select genes in a given pathway
# select_pathwayGenes <- MSigDB$C2.CP.KEGG$KEGG_HEMATOPOIETIC_CELL_LINEAGE
# selected_geneNormCounts <- subset(geneNormCounts, symbol %in% select_pathwayGenes)
# dim(selected_geneNormCounts)
# 
# dim(geneNormCounts)
# head(geneNormCounts)
# #filter based on samples selected
# dim(metadata)
# metadata$bamName <- gsub('-','.',metadata$bamName) #fix R reading error '-' is converted to .
# filtered_sample_names <- intersect(names(selected_geneNormCounts),metadata$bamName)
# selected_geneNormCounts <- selected_geneNormCounts[ , names(selected_geneNormCounts) %in% filtered_sample_names] 
# 
# selected_metadata <- subset(metadata, bamName %in% filtered_sample_names)
# selected_metadata <- metadata
# annotation <- data.frame(sex = selected_metadata$donorsex.cell.lines,
#                          level_1_diff_state = selected_metadata$grouplevel1differentiationstate,
#                          level_3_diff_state = selected_metadata$grouplevel3differentiationstate)
# 
# #assign the sample names to row names so that the heatmap function could use them for labelling
# rownames(annotation) <- selected_metadata$bamName
# 
# #temp
# rownames(annotation) <- metadata$bamName
# 
# 
# # eliminate the first 3 cols to get rid of the annotation
# m <- as.matrix(selected_geneNormCounts[4:ncol(selected_geneNormCounts)])
# 
# 
# #temp
# selected_rows <- sample(seq(1:nrow(geneNormCounts)),5000,replace=F)
# temp_geneCounts <-  geneNormCounts[selected_rows,]
# m <- as.matrix(temp_geneCounts[4:ncol(temp_geneCounts)])
# 
# #change the data type to integer
# m <- apply(m,2,as.numeric)
# 
# #removing those genes which dont vary much across the samples
# # so any gene with SD < .2 across the samples will be dropped 
# drop_genes <- which(apply(m,1,sd) < .2)
# m <-  m[-drop_genes,]
# 
# 
# #scaling genes across experiments
# mat.scaled <- t(scale(t(m)))
# 
# dim(m)
# source("~/dev//apRs/shiny_apps/heatMaps/memoised_pheatmap.R")
# 
# system.time(x <- cluster_mat(mat.scaled,"correlation",'average'))
# system.time(x <- memoised_cluster_mat(mat.scaled,"correlation",'average'))
# 
# forget(memoised_cluster_mat)
# png('test.png')
# system.time(heatmap(mat.scaled,
#                      annotation = annotation,
#                      scale="none",
#                      clustering_distance_rows = "manhattan",
#                      clustering_distance_cols = "manhattan",
#                      clustering_method = "average",
#                      fontsize_col = 5
#             )
# )
# system.time(heatmap(mat.scaled,keep.dendro=F))
# dev.off()
# X11.options()
# 
