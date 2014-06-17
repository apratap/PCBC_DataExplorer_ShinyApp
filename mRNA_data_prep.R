###
#get the PCBC samples gene normalized counts
###
cat('Reading the PCBC normalized mRNA Exp data from Synapse')
mRNA_NormCounts <- synGet('syn2247799')
#read in the file
mRNA_NormCounts <- read.table(mRNA_NormCounts@filePath,header=T,sep='\t', stringsAsFactors=FALSE)
#keep only uniq gene names #might need to be improved
mRNA_NormCounts <- mRNA_NormCounts[!duplicated(mRNA_NormCounts$symbol),]
cat('..Done\n\n')

######
#get the list siginificant genes from comparative analysis in synapse
#####
cat('Reading the precomputed significant genelist')
sigGenes_lists <- readRDS("precomputed_data/precomputed_sigGenes_lists.rds")
cat('..Done\n\n')

###
#get the metadata from the synapse for PCBC samples
###
mRNA_METADATA_ID <- 'syn2248030'
query <- sprintf('select * from entity where parentId=="%s"', mRNA_METADATA_ID)
mRNA_metadata <- synQuery(query)
cols_to_be_deleted = c('entity.benefactorId', 'entity.concreteType', 'entity.createdByPrincipalId', 
                       'entity.createdOn', 'entity.createdByPrincipalId', 'entity.id', 
                       'entity.modifiedOn', 'entity.modifiedByPrincipalId', 'entity.noteType', 
                       'entity.versionLabel', 'entity.versionComment', 'entity.versionNumber', 
                       'entity.parentId', 'entity.description', 'entity.eTag')
#delete the unwanted cols
mRNA_metadata <- mRNA_metadata[,!names(mRNA_metadata) %in% cols_to_be_deleted]

#remove the prefix 'entity.' from the df col names
names(mRNA_metadata) <- gsub('entity.','',names(mRNA_metadata))


########
#create new higher order labels for heatmap labels
#groups are merged into higher order
########
#1. linetype
# merge everthing other than hESC , iPSC, None into tissue
mRNA_metadata$mod_linetype <-  mRNA_metadata$linetype
index_to_be_replaced <-  ! mRNA_metadata$linetype %in% c('hESC','iPSC','None')
mRNA_metadata$mod_linetype[index_to_be_replaced] = 'tissue'


#2. diffnameshort
mRNA_metadata$mod_diffnameshort <- mRNA_metadata$diffnameshort
#2A. EB-LF should be renamed EB and similarly SC-LF should be SC
mRNA_metadata$mod_diffnameshort <- sub('SC-.*','SC',mRNA_metadata$mod_diffnameshort)
mRNA_metadata$mod_diffnameshort <- sub('EB-.*','EB',mRNA_metadata$mod_diffnameshort)
#2B. anything except the following should be grouped into Other
# SC, EB, DE, ECTO, MESO, MESO-30, MESO-15
index_to_be_replaced <- ! mRNA_metadata$mod_diffnameshort %in% c('SC', 'EB', 'DE', 'ECTO', 'MESO', 'MESO-30', 'MESO-15')
mRNA_metadata$mod_diffnameshort[index_to_be_replaced] = 'Others'

#3. origcell
mRNA_metadata$mod_origcell <- mRNA_metadata$origcell
#3A. categorize UCB CD34+ into CD34+
mRNA_metadata$mod_origcell <- sub('^.*CD34\\+$','CD34+',mRNA_metadata$mod_origcell,perl=T)
#3B. categ bone marrow to BM MSC
mRNA_metadata$mod_origcell <- sub('bone marrow','BM MSC',mRNA_metadata$mod_origcell,perl=T)
#3C. adultHeart to tissue
mRNA_metadata$mod_origcell <- sub('adultHeart','tissue',mRNA_metadata$mod_origcell,perl=T)
#3C. categ the following into Others
Others_group <- c('embryo', 'primEndothel', 'skin', 'None')
index_to_be_replaced <- mRNA_metadata$mod_origcell %in% Others_group
mRNA_metadata$mod_origcell[index_to_be_replaced] <- 'Others'