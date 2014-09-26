###
#get the PCBC samples geneExp normalized counts
###
cat('Reading the PCBC normalized mRNA Exp data from Synapse')
mRNA_NormCounts <- synGet('syn2701943')
#read in the file
mRNA_NormCounts <- read.delim(mRNA_NormCounts@filePath, header=T, sep='\t', as.is=T, stringsAsFactors = F, check.names=F)
rownames(mRNA_NormCounts) <- gsub('\\..*', '',mRNA_NormCounts$gene_id)
mRNA_NormCounts$symbol <- NULL
mRNA_NormCounts$gene_id <- NULL
mRNA_NormCounts$locus <- NULL
#apply(mRNA_NormCounts,2,class)
#mRNA_NormCounts <- as.data.frame(apply(mRNA_NormCounts,2,as.numeric))
#rownames(mRNA_NormCounts)
cat('..Done\n\n')



###
#get the metadata from the synapse for PCBC geneExp samples
###
mRNA_metadata <- synGet('syn2731147')
mRNA_metadata <- read.delim(mRNA_metadata@filePath, header=T, sep='\t',as.is=T, stringsAsFactors = F, check.names=F)
rownames(mRNA_metadata) <- mRNA_metadata[,'Decorated Name']
#keep only that metadata for samples which we have expression data
mRNA_metadata <- mRNA_metadata[rownames(mRNA_metadata) %in% colnames(mRNA_NormCounts),]


#get the list siginificant genes from comparative analysis in synapse
cat('Reading the precomputed significant genelist')
sigGenes_lists <- readRDS("precomputed_data/precomputed_sigGenes_lists.rds")
cat('..Done\n\n')
