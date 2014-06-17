library("preprocessCore")


###
#get the PCBC samples raw miRNA counts
###
cat('Reading the PCBC raw miRNA Exp data from Synapse')
miRNA_rawCounts <- synGet('syn2247832')
#read in the file
miRNA_rawCounts <- read.table(miRNA_rawCounts@filePath,header=T,sep='\t')
#keep only uniq miRNA #might need to be improved
miRNA_rawCounts <- miRNA_rawCounts[!duplicated(miRNA_rawCounts$mir),]
rownames(miRNA_rawCounts) <- miRNA_rawCounts$mir
miRNA_rawCounts$mir <- NULL
cat('..Done\n\n')

#Filter miRNA that have 0 expression in > 50 percent samples
miRNA_to_keep = apply(miRNA_rawCounts,1,function(x) sum(is.na(x))/length(x) < .50)
miRNA_rawCounts <- miRNA_rawCounts[miRNA_to_keep,]
#subtitute NA values with 0
miRNA_rawCounts[is.na(miRNA_rawCounts)] = 0

#Quantile Normalize
miRNA_normCounts <- normalize.quantiles(as.matrix(miRNA_rawCounts))
colnames(miRNA_normCounts) <- colnames(miRNA_rawCounts)
#split the paired miRNA name and use the second name
temp_miRNAs_names <- as.data.frame(do.call('rbind',strsplit(rownames(miRNA_rawCounts),',')))
rownames(miRNA_normCounts) <- temp_miRNAs_names[,2]


#get the miRNA to genes mapping table from synapse
miRNA_to_genes <- synGet('syn2246991')
miRNA_to_genes <- read.table(miRNA_to_genes@filePath, header=T, stringsAsFactors=FALSE)

#length(unique(miRNA_to_genes$Pathway))
#sum(unique(row.names(miRNA_normCounts)) %in% unique(miRNA_to_genes$Pathway))


#sum(unique(miRNA_to_genes$GeneID) %in% unique(mRNA_NormCounts$mod_gene_id))




# gsub("(.*)","",mRNA_NormCounts$gene_id)
# mRNA_NormCounts$gene_id



                     
###
#get the metadata from the synapse for PCBC samples
###
# miRNA_METADATA_ID <- 'syn2248741'
# query <- sprintf('select * from entity where parentId=="%s"', miRNA_METADATA_ID)
# miRNA_metadata <- synQuery(query)
# cols_to_be_deleted = c('entity.benefactorId', 'entity.concreteType', 'entity.createdByPrincipalId', 
#                        'entity.createdOn', 'entity.createdByPrincipalId', 'entity.id', 
#                        'entity.modifiedOn', 'entity.modifiedByPrincipalId', 'entity.noteType', 
#                        'entity.versionLabel', 'entity.versionComment', 'entity.versionNumber', 
#                        'entity.parentId', 'entity.description', 'entity.eTag')
# #delete the unwanted cols
# miRNA_metadata <- miRNA_metadata[,!names(miRNA_metadata) %in% cols_to_be_deleted]
# 
# #remove the prefix 'entity.' from the df col names
# names(miRNA_metadata) <- gsub('entity.','',names(miRNA_metadata))



####
#TEMP : getting the metadata from an external file
####
miRNA_metadata_ID <- 'syn2501833'
miRNA_metadata <- synGet(miRNA_metadata_ID)
miRNA_metadata <- read.table(miRNA_metadata@filePath,sep=",",header=T, stringsAsFactors=F)

########
#create new higher order labels for heatmap labels
#groups are merged into higher order
########

#1. linetype
# merge everthing other than hESC , iPSC, None into tissue
miRNA_metadata$Cell.Line.Type[miRNA_metadata$Cell.Line.Type == 'ESC'] = 'hESC'
miRNA_metadata$mod_linetype <-  miRNA_metadata$Cell.Line.Type
#index_to_be_replaced <-  ! miRNA_metadata$linetype %in% c('hESC','iPSC','None')
index_to_be_replaced <-  ! miRNA_metadata$Cell.Line.Type %in% c('hESC','iPSC','None')
miRNA_metadata$mod_linetype[index_to_be_replaced] = 'tissue'


#2. diffnameshort
miRNA_metadata$mod_diffnameshort <- miRNA_metadata$Diffname.Short
#2A. EB-LF should be renamed EB and similarly SC-LF should be SC
miRNA_metadata$mod_diffnameshort <- sub('SC-.*','SC',miRNA_metadata$mod_diffnameshort)
miRNA_metadata$mod_diffnameshort <- sub('EB-.*','EB',miRNA_metadata$mod_diffnameshort)
#2B. anything except the following should be grouped into Other
# SC, EB, DE, ECTO, MESO, MESO-30, MESO-15
index_to_be_replaced <- ! miRNA_metadata$mod_diffnameshort %in% c('SC', 'EB', 'DE', 'ECTO', 'MESO', 'MESO-30', 'MESO-15')
miRNA_metadata$mod_diffnameshort[index_to_be_replaced] = 'Others'


#3. origcell
miRNA_metadata$mod_origcell <- miRNA_metadata$Cell.Type.of.Origin
table(miRNA_metadata$Cell.Type.of.Origin)
#miRNA_metadata$mod_origcell <- miRNA_metadata$origcell
#3A. categorize UCB CD34+ into CD34+
miRNA_metadata$mod_origcell <- sub('^.*CD34\\+$','CD34+',miRNA_metadata$mod_origcell,perl=T)
miRNA_metadata$mod_origcell <- sub('^.*CD34\\+.*$','CD34+',miRNA_metadata$mod_origcell,perl=T)
#3B. categ bone marrow to BM MSC
miRNA_metadata$mod_origcell <- sub('bone marrow','BM MSC',miRNA_metadata$mod_origcell,perl=T)
#3C. adultHeart to tissue
miRNA_metadata$mod_origcell <- sub('adultHeart','tissue',miRNA_metadata$mod_origcell,perl=T)
#3C. categ the following into Others
Others_group <- c('embryo', 'primEndothel', 'skin', 'None', 'cardiac endothelial cell',
                  'mesenchymal stem cell', 'cardiac myocyte', 'inner cell mass', 'mononuclear',
                  'pulmonary artery endothelium', 'N/A')
index_to_be_replaced <- miRNA_metadata$mod_origcell %in% Others_group
miRNA_metadata$mod_origcell[index_to_be_replaced] <- 'Others'

#4. Induction genes
miRNA_metadata$inductiongenes <- miRNA_metadata$Reprogramming.Gene.Combination

#5. decoratedName
miRNA_metadata$decoratedName <- gsub('-','.',as.vector(miRNA_metadata$uid))

#6. Gender
miRNA_metadata$donorsex <- miRNA_metadata$Gender


