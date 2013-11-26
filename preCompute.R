library("synapseClient")
library("gdata")
library("plyr")

# #login to synapse
synapseLogin()




#########
#Read the hg19 genes annotaiton and save a precomputed df
##########
k <- keys(org.Hs.eg.db,keytype="SYMBOL")
hg19_gene_annot <- select(org.Hs.eg.db, keys=k, columns=c("GENENAME","ALIAS"), keytype="SYMBOL")
hg19_gene_annot <- ddply(hg19_gene_annot,.variables=c('SYMBOL','GENENAME'), .fun = function(x) paste(x$ALIAS,collapse=', ') )
hg19_gene_annot['ALIAS'] <- hg19_gene_annot$V1
hg19_gene_annot$V1 <- NULL
saveRDS(hg19_gene_annot,"hg19_gene_annot.RDS")





####
#1. get the names of all the genes
###
#get the PCBC samples gene normalized counts
syn_geneNormCounts <- synGet('syn1968267')
#read in the file to get the names of allGenes
geneNormCounts <- read.table(syn_geneNormCounts@filePath,header=T,sep='\t')
allGenes <- toupper(unique(geneNormCounts$symbol))


######
#2. get the siginificant gene lists
#####
#get the list siginificant genes from comparative analysis in synapse
sigGenes_synId <- "syn2244163"
syn_sigGenes <- synGet(sigGenes_synId)
#read the file
sigGenes <- read.xls(syn_sigGenes@filePath)
#make sure all the gene names are upper case and remove the suffix
sigGenes$symbol <- toupper(sigGenes$symbol)
sigGenes_lists <- split(sigGenes$symbol,sigGenes$name)
names(sigGenes_lists) <- gsub('-fold2.0_adjp0.05','', names(sigGenes_lists) ) 
saveRDS(sigGenes_lists , file="precomputed_sigGenes_lists.rds")


####
#3. load the pathways info
####
MSIGDB<-synGet("syn2227979")
load(MSIGDB@filePath) #available as MSigDB R object




######
#4. FET test 
######
#which pathways are enriched in a particular gene list v/s others gene lists


#run the actual FET test
FET_test <- function(selected_pathway_genes,selected_geneList,allGenes){
  #create a empty 2,2 matrix
  mat2x2 = mat.or.vec(2,2)
  mat2x2[1,1] <- length(intersect(selected_pathway_genes,selected_geneList))
  mat2x2[1,2] <- length(setdiff(selected_geneList,selected_pathway_genes)) 
  mat2x2[2,1] <- length(setdiff(selected_pathway_genes,selected_geneList))
  mat2x2[2,2] <- length(allGenes) - (mat2x2[1,1] + mat2x2[1,2] + mat2x2[2,1] )
  
  x<-fisher.test(mat2x2)
  x$p.value
}



#returs a df of enriched pathways from a genelist 
#comparing against a selected pathway database
find_enrichedPathways <- function(geneList,selected_pathway_db){
  df = ldply(selected_pathway_db,.fun=FET_test,geneList,allGenes)
  names(df) = c('pathways','p.value')
  df$p.adj = p.adjust(df$p.value,method="hochberg")
  df
}


#run the FET test for each significant gene list testing for enriched pathways
df = ldply(sigGenes_lists,.fun=find_enrichedPathways,selected_pathway_db,.progress="text")
df$significant_gene_list_name = df$.id  
df$.id <- NULL

#keep only those pathways which have p.adj < .05
enrichedPathways_in_sigGenes_list <- subset(df, p.adj < .05)


#save the precomputed results
saveRDS(enrichedPathways_in_sigGenes_list, file="precomputed_enrichedPathways_in_geneLists.rds")




#######################
#######################
#######################

#testing
selected_pathway_db <- MSigDB$C2.CP.KEGG
selected_pathway_genes <- unique(selected_pathway_db[[2]])
selected_geneList <- unique(sigGenes_lists[[2]])

#test 1
FET_test(selected_pathway_genes,selected_geneList,allGenes)

#test 2
find_enrichedPathways(selected_geneList,selected_pathway_db)

#test 3
df = ldply(sigGenes_lists,.fun=find_enrichedPathways,selected_pathway_db,.progress="text")
