# Load libraries
library(AnnotationHub)
library(ensembldb)

# Connect to AnnotationHub
ah <- AnnotationHub()

# Explore all species information available
unique(ah$species) %>% View()

# Explore the types of Data Objects available
unique(ah$rdataclass) %>% View()

# Explore the Data Providers
unique(ah$dataprovider) %>% View()

human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
human_ens

# Extract annotations of interest
human_ens <- human_ens[["AH75011"]]

# Extract gene-level information
genes(human_ens, return.type = "data.frame") %>% View()

# Create a gene-level dataframe 
annotations_ahb <- genes(human_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, entrezid, gene_biotype) 

head(annotations_ahb)
dim(annotations_ahb)
library(purrr)
which(map(annotations_ahb$entrezid, length) > 1)

# Create a transcript dataframe
mouse_ens <- query(ah, c("Mus musculus", "EnsDb"))
mouse_ens

# Extract annotations of interest
mouse_ens <- mouse_ens[["AH78811"]]

txdb <- transcripts(mouse_ens, return.type = "data.frame") %>%
  dplyr::select(tx_id, gene_id)
head(txdb)
txdb <- txdb[grep("ENSMUST", txdb$tx_id),]

# Create a gene-level dataframe
genedb <- genes(mouse_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, symbol)

# Merge the two dataframes together
library(dplyr)
annotations <- inner_join(txdb, genedb)
head(annotations)
dim(annotations)

custom_txDb = loadDb('/Users/hhuang/Desktop/salmon_tutorial/homo/txdb')
custom_txDb
k <- keys(custom_txDb, keytype = "TXNAME")
tx2gene <- dplyr::select(custom_txDb, k, "GENEID", "TXNAME")

short.ensembl = sapply(strsplit(tx2gene$TXNAME, '\\.'), function(x) x[1])
length(intersect(short.ensembl, annotations$tx_id))
length(setdiff(short.ensembl, annotations$tx_id))
(setdiff(short.ensembl, annotations$tx_id))
length(setdiff(annotations$tx_id, short.ensembl))
(setdiff(annotations$tx_id, short.ensembl))
tx2gene$TXNAME
dim(tx2gene)
k
