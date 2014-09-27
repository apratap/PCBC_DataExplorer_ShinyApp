cat('Reading the PCBC methylation data....')
meth_data <- synGet('syn2731494')
meth_data <- read.delim(meth_data@filePath, header=T, sep='\t', as.is=T, stringsAsFactors = F, 
                        check.names=F)
rownames(meth_data) <- meth_data[,1]
meth_data[,1] <- NULL
cat('Done \n\n')


#meth to gene annotation
cat('Reading the PCBC methylation to genes mapping file....')
meth_to_gene <- synGet('syn2731501')
meth_to_gene <- read.delim(meth_to_gene@filePath, header=T, sep='\t', as.is=T, stringsAsFactors = F, 
                           check.names=F)
meth_to_gene$entrezID <-  as.character(meth_to_gene$entrezID)
cat('Done \n\n')

#miRNA metadata
cat('Reading the PCBC  methlyation metadata')
meth_metadata <- synGet('syn2731151') 
meth_metadata <- read.delim(meth_metadata@filePath, header=T, sep='\t',as.is=T, stringsAsFactors = F, check.names=F)
rownames(meth_metadata) <- meth_metadata$Sample
#keep only those samples that are present in the expression matrix
rows_to_keep <- rownames(meth_metadata) %in% colnames(meth_data)
meth_metadata <- meth_metadata[rows_to_keep, ]
cat('..Done\n\n')





