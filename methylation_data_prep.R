
#miRNA metadata
cat('Reading the PCBC  methlyation metadata')
meth_metadata <- synGet('syn2731151') 
meth_metadata <- read.delim(meth_metadata@filePath, header=T, sep='\t',as.is=T, stringsAsFactors = F, check.names=F)
rownames(meth_metadata) <- meth_metadata$Sample



#keep only those samples that are present in the expression matrix
#rows_to_keep <- rownames(miRNA_metadata) %in% colnames(miRNA_normCounts)
#miRNA_metadata <- miRNA_metadata[rows_to_keep, ]
cat('..Done\n\n')