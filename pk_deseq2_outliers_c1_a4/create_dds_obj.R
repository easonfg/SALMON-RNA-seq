rm(list=ls()) 

library("tximport")
library("readr")
library("tximportData")

dir <- '../quants'
filenames = c(list.files(dir))
filenames = filenames[-grep('A4', filenames)]
filenames = filenames[-grep('C1', filenames)]
filenames
# snames 

# label one group with condition '0' another with condition '1'
# c = startsWith(snames, 'S')
# samples = data.frame(snames, condition = as.character(as.integer(c)))

# make other conditions
conditions = c(paste0('Treatment', 1:3), paste0('Control', 2:4))
conditions

# snames = paste(1:20, '_quant', sep = '')
# snames = c(paste0('Treatment', 1:4), paste0('Control', 1:4))
snames = sapply(filenames, function(x) strsplit(x, '_')[[1]][1])
snames

samples = data.frame(snames, condition = factor(c(rep('Treatment', 3),
                                                rep('Control',3))))
print(samples)


rownames(samples) <- samples$snames
files <- file.path(dir, filenames, "quant.sf")
files

### make own gencode from gff ###
library('GenomicFeatures')
# custom_txDb = makeTxDbFromGFF(file = '/Users/henry/Desktop/salmon_tutorial/homo/gencode.v32.chr_patch_hapl_scaff.annotation.gtf')
# saveDb(custom_txDb, '/Users/henry/Desktop/salmon_tutorial/homo/txdb')
custom_txDb = loadDb('/Users/hhuang/Desktop/ketone_exp/ketone_DESeq2/mouse_txdb')
k <- keys(custom_txDb, keytype = "TXNAME")
tx2gene <- select(custom_txDb, k, "GENEID", "TXNAME")
#  k <- keys(custom_txDb, keytype = "EXONNAME")
# tx2gene <- select(custom_txDb, k, "GENEID",  "EXONNAME")
head(tx2gene)
### make own gencode from gff ###


txi <- tximport(files, type="salmon", tx2gene=tx2gene)
dim(txi$counts)


#### using R package annotationhub####################
## does work, nothing wrong with it, returned just little bit fewer genes than other method
# library(AnnotationHub)
# library(ensembldb)
# 
# # Connect to AnnotationHub
# ah <- AnnotationHub()
#
# # Create a transcript dataframe
# mouse_ens <- query(ah, c("Mus musculus", "EnsDb"))
# mouse_ens
# 
# # Extract annotations of interest
# mouse_ens <- mouse_ens[["AH78811"]]
# 
# txdb <- transcripts(mouse_ens, return.type = "data.frame") %>%
#   dplyr::select(tx_id_version, gene_id)
# head(txdb)
# txdb <- txdb[grep("ENSMUST", txdb$tx_id),]
# 
# # Create a gene-level dataframe
# genedb <- genes(mouse_ens, return.type = "data.frame")  %>%
#   dplyr::select(gene_id, symbol)
# 
# # Merge the two dataframes together
# library(dplyr)
# annotations <- inner_join(txdb, genedb)
# txi2 <- tximport(files, type="salmon", tx2gene=annotations)
# dim(txi2$counts)
#### using R package annotationhub####################


library("DESeq2")
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~ condition)

dds$snames
dds$condition
head(samples)


# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]

dds <- DESeq(dds)

# colData(dds)$cell_type = factor(c(rep('iPSC', 8), rep('Microglia', 12)))
# colData(dds)$treatment = factor(c(rep('No Treatment', 2),
#                                   rep('BHB', 2),
#                                   rep('LPS', 2),
#                                   rep('LPS_BHB',2),
#                                   rep('No Treatment', 3),
#                                   rep('BHB', 3),
#                                   rep('LPS', 3),
#                                   rep('LPS_BHB',3)))

saveRDS(dds, file = 'dds_object')

assay(dds)
dds$condition

