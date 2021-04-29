rm(list=ls()) 

library("DESeq2")

dds = readRDS('dds_object')


### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
### Extract the rlog matrix from the object
rld_mat <- assay(rld) 
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)  
rld_cor
### Compute pairwise dist values
rld_dist <- dist(t(rld_mat)  )
rld_dist
### Load pheatmap package
library(pheatmap)
### Plot heatmap
# meta 
colData(dds)
sample.names = colData(dds)$snames
sample.names
sampletype = colData(dds)$condition
# Cells = sapply(as.vector(sampletype), function(x) strsplit(x, '_')[[1]][1])
# treatment = sapply(as.vector(sampletype), function(x) paste(strsplit(x, '_')[[1]][-1], collapse = '_'))

rownames(rld_cor) = sample.names
colnames(rld_cor) = sample.names

# meta = data.frame(row.names = sample.names, Cells = Cells, Treatment = treatment)
meta = data.frame(row.names = sample.names, Treatment = colData(dds)$condition)
pheatmap(rld_cor, annotation = meta)
# pheatmap(rld_dist, annotation = meta)

### Principal component plot of the samples###
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
library(ggplot2)
# ggplot(pcaData, aes(PC1, PC2, color=sub_condition, shape=treatment)) +
ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_text(aes(vjust = -1))+
  coord_fixed()


##PCA at lower PC
## The plotPCA() function will only return the values for PC1 and PC2.
## If you would like to explore the additional PCs in your data or if you would like to identify genes that contribute most to the PCs,
## you can use the prcomp() function. For example, to plot any of the PCs we could run the following code:
pca <- prcomp(t(rld_mat))
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x, sampletype)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = Cells, shape = treatment))


#### dispersion plot
# You expect your data to generally scatter around the curve,
# with the dispersion decreasing with increasing mean expression levels. If you see a cloud or different shapes,
# then you might want to explore your data more to see if you have contamination (mitochondrial, etc.) or outlier samples.
# Note how much shrinkage you get across the whole range of means in the plotDispEsts() plot for any experiment with low degrees of freedom
plotDispEsts(dds)


###Heatmap of the sample-to-sample distances##
# vsd <- vst(dds, blind=FALSE)
# sampleDists <- dist(t(assay(vsd)))
# library("RColorBrewer")
# sampleDistMatrix <- as.matrix(sampleDists)
# # sampleDistMatrix <- (sampleDists)
# colnames(sampleDistMatrix) = sample.names
# rownames(sampleDistMatrix) = sample.names
# # colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix, annotation = meta)
# pheatmap(sampleDistMatrix)
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)
# pheatmap(sampleDistMatrix,
#          cluster_rows=FALSE,
#          cluster_cols=FALSE)
###Heatmap of the sample-to-sample distances##

###Heatmap of the genes distances##
# geneDists <- dist((assay(vsd)))
###Heatmap of the genes distances##

### Principal component plot of the samples###
# rld <- rlog(dds, blind=TRUE)
# rld_PCA = plotPCA(rld, intgroup=c("sub_condition", "treatment"), returnData=TRUE)
# percentVar <- round(100 * attr(rld_PCA, "percentVar"))
# ggplot(rld_PCA, aes(PC1, PC2, color=sub_condition, shape=treatment)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   coord_fixed()





## correct normalization method
# rld <- rlog(dds, blind=TRUE)
# rld_mat <- assay(rld)
### Principal component plot of the samples###

## pearson correlation heatmap##
# cor.matrix <- data.frame(cor(counts(dds,normalized=TRUE)))
# cor.matrix <- data.frame(cor(rld_mat))
# pheatmap(cor.matrix, cluster_rows=FALSE, cluster_cols=FALSE)

### use plotly to not cluster
# library(plotly)
# ax_y <- list(
#   title = "",
#   zeroline = FALSE,
#   showline = FALSE,
#   showgrid = FALSE
# )
# 
# p = plot_ly(x = rownames(cor.matrix), y = colnames(cor.matrix), z = as.matrix(cor.matrix), 
#             colors = colorRamp(c('blue', 'white', "red")), type = "heatmap") %>%
#   layout(yaxis = ax_y)
# print(p)
# 
# ## correlation assesesment
# rowSums(cor.matrix[grep('Br12CD', rownames(cor.matrix)), grep('Br12CD', colnames(cor.matrix))])
# rowSums(cor.matrix[grep('Br12KD', rownames(cor.matrix)), grep('Br12KD', colnames(cor.matrix))])
# rowSums(cor.matrix[grep('Br26CD', rownames(cor.matrix)), grep('Br26CD', colnames(cor.matrix))])
# rowSums(cor.matrix[grep('Br26KD', rownames(cor.matrix)), grep('Br26KD', colnames(cor.matrix))])
# rowSums(cor.matrix[grep('LvrCD', rownames(cor.matrix)), grep('LvrCD', colnames(cor.matrix))])
# rowSums(cor.matrix[grep('LvrKD', rownames(cor.matrix)), grep('LvrKD', colnames(cor.matrix))])
# ## pearson correlation heatmap##
# 
# 
# 
# # keep <- rowSums(counts(dds)) >= 10
# # dds <- dds[keep,]
# 
# 
# ## look for outliers##
# threshold_n = 100
# norm_tpm = data.frame(sapply(seq(1,threshold_n), function(x) colSums(counts(dds,normalized=TRUE)>x)))
# norm_tpm$subject <- rownames(norm_tpm)
# colnames(norm_tpm) = as.vector(c(seq(1:threshold_n), 'subject'))
# norm_tpm$subject <- factor(norm_tpm$subject)
# 
# # library("tidyverse")
# # norm_tpm_tibble = as_tibble(norm_tpm)
# # df <- norm_tpm_tibble %>%
# #   select(1:threshold_n, subject) %>%
# #   gather(key = "variable", value = "value", -date)
# # 
# 
# ## normalized TPM threshold curve
# library(tidyr)
# sample_name = 'LvrCD'
# norm_tpm_long <- gather(norm_tpm, threshold, count, 1:threshold_n, factor_key=TRUE)
# norm_tpm_long$threshold = as.double(norm_tpm_long$threshold)
# sub_norm_tpm_long = norm_tpm_long[grep(sample_name, norm_tpm_long$subject),]
# outlier_assess = ggplot(sub_norm_tpm_long, aes(x = threshold, y = count)) + 
#   geom_line(aes(color = subject, linetype = subject))
# print(outlier_assess)
# 
# ## correlation scatterplot
# plot_data = data.frame(cbind(Br12CD1 = log2(counts(dds,normalized=TRUE))[,'LvrCD3_quant'],
#                              Br12CD2 = log2(counts(dds,normalized=TRUE))[,'LvrCD1_quant']))
# ggplot(plot_data, aes(x=Br12CD1, y=Br12CD2)) + geom_point()
# 
# 
# 
