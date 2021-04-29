rm(list=ls())
graphics.off()

library(gdata)
library(fgsea)
library(tidyverse)
library(stringr)
library(EnhancedVolcano)

run_fgsea = function(res, genemap, pathways.hallmark, graph_title){
  pathway.name = str_extract(names(pathways.hallmark)[1], '[^_]+')
  # browser()
  res_copy = res
  res_copy$rank_metric = res_copy$log2FoldChange * -log10(res_copy$pvalue)
  # res_copy2 = res_copy[,c('SYMBOL', 'log2FoldChange')]
  # res_copy2 = res_copy[,c('SYMBOL', 'rank_metric')]
  res_copy2 = res_copy[,c('SYMBOL', 'stat')]
  # res_copy2 = res_copy[,c('SYMBOL', 'pvalue')]
  res_copy2 = na.omit(res_copy2)
  
  #get ranks of genes
  ranks <- deframe(res_copy2)
  head(ranks, 20)
  ranks = ranks[order(ranks)]
  
  #get results
  fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
  
  #tidy up
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  # Show in a nice table:
  fgseaResTidy %>%
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
    arrange(padj) %>%
    DT::datatable()
  
  fgseaResTidy = fgseaResTidy[fgseaResTidy$padj < 0.05,]
  
  # browser()
  #Plot the normalized enrichment scores
  # quartz(graph_title)
  jpeg(paste('gsea_results/', graph_title, '_', pathway.name, '.jpeg', sep = ''),
       units="in", width=10, height=10, res=500)
  gg_obj = ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=graph_title) +
    theme_minimal()
  print(gg_obj)
  dev.off()
  # ggsave(paste('/Users/hhuang/Desktop/Microglia/microglia_deseq2/gsea_results/', graph_title, sep = ''))
  
  #What genes are in each of these pathways?
  # t = pathways.hallmark %>%
  #   enframe("pathway", "SYMBOL") %>%
  #   unnest() %>%
  #   inner_join(as_tibble(res_copy), cols=c("SYMBOL" ))
  
  # return(list(res = fgseaResTidy, genes = t))
}

add_symbol = function(res_copy, genemap){
  res_copy$ensembl = rownames(res_copy)
  idx <- match( res_copy$ensembl, genemap$ensembl_gene_id_version )
  res_copy$SYMBOL <- genemap$hgnc_symbol[ idx ]
  res_copy[is.na(res_copy$SYMBOL), 'SYMBOL'] = res_copy[is.na(res_copy$SYMBOL), 'ensembl']
  res_copy[res_copy$SYMBOL == '', 'SYMBOL'] = res_copy[res_copy$SYMBOL == '', 'ensembl']
  return(res_copy)
}


change_sym = function(name_ls, genemap){
  # name_ls = sapply( strsplit( name_ls, split="\\+" ), "[", 1 )
  
  idx <- match( name_ls, genemap$ensembl_gene_id_version )
  name_sym <- genemap$hgnc_symbol[ idx ]
  # res_copy2 = res_copy[,c('SYMBOL', 'stat')]
  # res_copy2 = na.omit(res_copy2)
  return(name_sym)
}

get_gene_name = function(epi){
  library('biomaRt')
  ensembl <- useMart('ENSEMBL_MART_ENSEMBL', dataset="hsapiens_gene_ensembl")
  annotation <- getBM(attributes=c("uniprotswissprot", "hgnc_symbol", "uniprot_gn_symbol"),
                      filters="uniprotswissprot",
                      values=levels(epi$UniProtIds), mart=ensembl)
  protID = levels(epi$UniProtIds)
  ## split the items with ; and add to list
  protID2 = c(protID[-grep(';', protID)],
              unlist(strsplit(protID[grep(';', protID)], ';')))
  idx <- match(protID2, annotation$uniprotswissprot )
  prot_gene_ID <- annotation$hgnc_symbol[ idx ]
  return(prot_gene_ID)
}



run_pipeline = function(condition1, condition2){
  
  colnames(dds) = colData(dds)$condition
  
  graph_title = paste(condition1, '_vs_', condition2, sep = '')
  print(graph_title)
  
  res_M_i_no_treatment <- results(dds, contrast= c('condition', condition1, condition2))
  unshrunken.res = res_M_i_no_treatment
  mcols(res_M_i_no_treatment)
  
  ### shrunken lfc. good for visualizations
  # shrunken.res <- lfcShrink(dds, contrast=c('condition', condition1, condition2), res=res_M_i_no_treatment)
  res_M_i_no_treatment <- lfcShrink(dds, contrast=c('condition', condition1, condition2), res=res_M_i_no_treatment)
  shrunken.res = res_M_i_no_treatment
  
  ## MA plots
  
  jpeg(paste('MAplots/', graph_title, '.jpeg', sep = ''),
       units="in", width=10, height=10, res=500)
  # plotMA(shrunken.res, ylim=c(-2,2))
  print(plotMA(shrunken.res))
  dev.off()
  # plotMA(unshrunken.res, ylim=c(-2,2))
  condition1
  condition2
  
  summary(unshrunken.res, alpha = 0.05)
  
  
  ## add symbol and filter out NA
  res_M_i_no_treatment = res_M_i_no_treatment[!is.na(res_M_i_no_treatment$padj),]
  #res_M_i_no_treatment = add_symbol(res_M_i_no_treatment, genemap)
  #dim(res_M_i_no_treatment)
  #(res_M_i_no_treatment)
  short_ensembl = str_extract(rownames(res_M_i_no_treatment), '[^.]+')
  short_ensembl
  res_M_i_no_treatment$ensembl = short_ensembl
  res_M_i_no_treatment
  
  # https://hbctraining.github.io/DGE_workshop_salmon/lessons/AnnotationDbi_lesson.html
  # https://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb
  library(org.Mm.eg.db)
  ens2symbol <- AnnotationDbi::select(org.Mm.eg.db,
                                      key= short_ensembl, 
                                      columns=c('ENTREZID', "SYMBOL"),
                                      keytype="ENSEMBL")
  
  dim(ens2symbol)
  length(short_ensembl)
  dim(res_M_i_no_treatment)
  
  res_M_i_no_treatment = inner_join(as_tibble(res_M_i_no_treatment), as_tibble(ens2symbol), by=c("ensembl"="ENSEMBL"))
  ##fill NA symbols with ensembl ids
  res_M_i_no_treatment[is.na(res_M_i_no_treatment$SYMBOL), 'SYMBOL'] = res_M_i_no_treatment[is.na(res_M_i_no_treatment$SYMBOL), 'ensembl']
  
  
  ### normalized count
  normalized_counts <- counts(dds, normalized=T) %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    mutate(gene = str_extract(gene, '[^.]+')) %>%
    as_tibble() %>%
    left_join(ens2symbol, by=c("gene" = "ENSEMBL"))

  ### heatmap of significant genes
  ### Set thresholds
  library(dplyr)
  padj.cutoff <- 0.05
  ### sig genes
  res_tableOE_tb <- res_M_i_no_treatment %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  sigOE <- res_tableOE_tb %>%
    filter(padj < padj.cutoff)
  sigOE
  
  write.csv(sigOE, 'sig_genes.csv')
  write.csv(res_M_i_no_treatment, 'all_genes.csv')

  rel.indx = c(grep(paste0(condition1, '(\\.|$)'), colnames(normalized_counts)), grep(paste0(condition2,'(\\.|$)'), colnames(normalized_counts)))
  norm_OEsig <- normalized_counts[,c(1,rel.indx)] %>%
    filter(gene %in% sigOE$ensembl)
  # sigOE
  # normalized_counts$gene
  ### Set a color palette
  library(RColorBrewer)
  heat_colors <- rev(brewer.pal(7, "RdYlBu"))

  ### Run pheatmap using the metadata data frame for the annotation
  library(pheatmap)

  jpeg(paste('sig_gene_heatmap/', graph_title, '.jpeg', sep = ''),
       units="in", width=10, height=10, res=500)
  pheatmap(norm_OEsig[,-1],
           color = heat_colors,
           cluster_rows = T,
           show_rownames = F,
           border_color = NA,
           fontsize = 10,
           scale = "row",
           fontsize_row = 10,
           height = 20)
  dev.off()

  ########### Over-representation analysis##############################
  ### clusterprofiler
  ## Create background dataset for hypergeometric testing using all genes tested for significance in the results
  allOE_genes <- as.character(res_M_i_no_treatment$ensembl)

  ## Extract significant results
  sigOE <- dplyr::filter(res_M_i_no_treatment, padj < 0.05)

  sigOE_genes <- as.character(sigOE$ensembl)

  ## Run GO enrichment analysis
  library(DOSE)
  library(pathview)
  library(clusterProfiler)
  ego <- enrichGO(gene = sigOE_genes,
                  universe = allOE_genes,
                  keyType = "ENSEMBL",
                  OrgDb = org.Mm.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)

  ## Dotplot
  jpeg(paste('over_representation_analysis/GO/dotplots/', graph_title, '.jpeg', sep = ''),
       units="in", width=10, height=10, res=500)
  print(dotplot(ego, showCategory=50))
  dev.off()
  ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
  jpeg(paste('over_representation_analysis/GO/enrichment_clusters/', graph_title, '.jpeg', sep = ''),
       units="in", width=10, height=10, res=500)
  print(emapplot(ego, showCategory = 50))
  dev.off()

  ## category netplot
  ## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
  OE_foldchanges <- sigOE$log2FoldChange

  names(OE_foldchanges) <- sigOE$ensembl

  ## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
  jpeg(paste('over_representation_analysis/GO/netplots/', graph_title, '.jpeg', sep = ''),
       units="in", width=10, height=10, res=500)
  print(cnetplot(ego,
           categorySize="pvalue",
           showCategory = 10,
           foldChange=OE_foldchanges,
           vertex.label.font=6))
  dev.off()


  ego <- enrichGO(gene = sigOE_genes,
                  universe = allOE_genes,
                  keyType = "ENSEMBL",
                  OrgDb = org.Mm.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)

  ## Dotplot
  jpeg(paste('over_representation_analysis/GO_mf/dotplots/', graph_title, '.jpeg', sep = ''),
       units="in", width=10, height=10, res=500)
  print(dotplot(ego, showCategory=50))
  dev.off()
  ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
  if (nrow(ego)>1){
    jpeg(paste('over_representation_analysis/GO_mf/enrichment_clusters/', graph_title, '.jpeg', sep = ''),
         units="in", width=10, height=10, res=500)
    print(emapplot(ego, showCategory = 50))
    dev.off()
  }

  ## category netplot
  ## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
  OE_foldchanges <- sigOE$log2FoldChange

  names(OE_foldchanges) <- sigOE$ensembl

  ## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
  if (nrow(ego)>0){
    jpeg(paste('over_representation_analysis/GO_mf/netplots/', graph_title, '.jpeg', sep = ''),
         units="in", width=10, height=10, res=500)
    print(cnetplot(ego,
                   categorySize="pvalue",
                   showCategory = 10,
                   foldChange=OE_foldchanges,
                   vertex.label.font=6))
    dev.off()
  }

  ## If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
  # OE_foldchanges <- ifelse(OE_foldchanges > 2, 2, OE_foldchanges)
  # OE_foldchanges <- ifelse(OE_foldchanges < -2, -2, OE_foldchanges)
  #
  # cnetplot(ego,
  #          categorySize="pvalue",
  #          showCategory = 5,
  #          foldChange=OE_foldchanges,
  #          vertex.label.font=6)



  ## Run KEGG enrichment analysis
  allOE_genes <- as.character(res_M_i_no_treatment$ENTREZID)
  ## Extract significant results
  sigOE <- dplyr::filter(res_M_i_no_treatment, padj < 0.05)
  sigOE_genes <- as.character(sigOE$ENTREZID)

  ekegg <- enrichKEGG(gene = sigOE_genes,
                  universe = allOE_genes,
                  organism = 'mmu', 
                  keyType = "ncbi-geneid",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)

  if (nrow(ekegg) > 0) {
    ## Dotplot
    jpeg(paste('over_representation_analysis/KEGG/dotplots/', graph_title, '.jpeg', sep = ''),
         units="in", width=10, height=10, res=500)
    print(dotplot(ekegg, showCategory=50))
    dev.off()
    ## Enrichmap clusters the 50 most significant (by padj) KEGG terms to visualize relationships between terms
    jpeg(paste('over_representation_analysis/KEGG/enrichment_clusters/', graph_title, '.jpeg', sep = ''),
         units="in", width=10, height=10, res=500)
    print(emapplot(ekegg, showCategory = 50))
    dev.off()

    ## category netplot
    ## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
    OE_foldchanges <- sigOE$log2FoldChange

    names(OE_foldchanges) <- sigOE$ensembl

    ## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
    jpeg(paste('over_representation_analysis/KEGG/netplots/', graph_title, '.jpeg', sep = ''),
         units="in", width=10, height=10, res=500)
    print(cnetplot(ekegg,
             categorySize="pvalue",
             showCategory = 10,
             foldChange=OE_foldchanges,
             vertex.label.font=6))
    dev.off()
  }
  ########### Over-representation analysis##############################

  ########### Functional class scoring##############################
  ## foldchange
  ## Remove any NA values (reduces the data by quite a bit)
  res_entrez <- dplyr::filter(res_M_i_no_treatment, ENTREZID != "NA")
  ## Remove any Entrez duplicates
  res_entrez <- res_entrez[which(duplicated(res_entrez$ENTREZID) == F), ]
  ## Extract the foldchanges
  foldchanges <- res_entrez$log2FoldChange
  ## Name each fold change with the corresponding Entrez ID
  names(foldchanges) <- res_entrez$ENTREZID
  ## Sort fold changes in decreasing order
  foldchanges <- sort(foldchanges, decreasing = TRUE)
  head(foldchanges)

  ## statistics
  stats.res_entrez <- dplyr::filter(res_M_i_no_treatment, ENTREZID != "NA")
  ## Remove any Entrez duplicates
  stats.res_entrez <- stats.res_entrez[which(duplicated(stats.res_entrez$ENTREZID) == F), ]
  ## Extract the foldchanges
  test.stats <- stats.res_entrez$stat
  ## Name each fold change with the corresponding Entrez ID
  names(test.stats) <- stats.res_entrez$ENTREZID
  ## Sort fold changes in decreasing order
  test.stats <- sort(test.stats, decreasing = TRUE)
  head(test.stats)


  ## GSEA using gene sets from KEGG pathways
  # gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
  gseaKEGG <- gseKEGG(geneList = test.stats, # ordered named vector of fold changes (Entrez IDs are the associated names)
                      organism = "mmu", # supported organisms listed below
                      nPerm = 1000, # default number permutations
                      minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                      pvalueCutoff = 0.05, # padj cutoff value
                      verbose = FALSE)

  ## Extract the GSEA results
  gseaKEGG_results <- gseaKEGG@result

  jpeg(paste('clusterProfiler_gsea/KEGG/', graph_title, '.jpeg', sep = ''),
       units="in", width=10, height=10, res=500)
  gg_obj = ggplot(gseaKEGG_results, aes(reorder(Description, NES), NES)) +
    geom_col(aes(fill=p.adjust<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=graph_title) +
    theme_minimal()
  print(gg_obj)
  dev.off()

  # ## Plot the GSEA plot for a single enriched pathway, `hsa03040`
  # gseaplot(gseaKEGG, geneSetID = gseaKEGG_results$ID[1])
  #
  # ## Use the Pathview R package to integrate the KEGG pathway data from clusterProfiler into pathway images:
  ### detach("package:dplyr", unload=TRUE) # first unload dplyr to avoid conflicts
  #
  # ## Output images for a single significant KEGG pathway
  # pathview(gene.data = foldchanges,
  #          pathway.id = gseaKEGG_results$ID[1],
  #          species = "mmu",
  #          limit = list(gene = 2, # value gives the max/min limit for foldchanges
  #                       cpd = 1))


  # GSEA using gene sets associated with BP Gene Ontology terms
  # gseaGO <- gseGO(geneList = foldchanges,
  gseaGO <- gseGO(geneList = test.stats,
                  OrgDb = org.Mm.eg.db,
                  ont = 'BP',
                  nPerm = 1000,
                  minGSSize = 20,
                  pvalueCutoff = 0.05,
                  verbose = FALSE)

  gseaGO_results <- gseaGO@result
  # View(gseaGO_results)

  jpeg(paste('clusterProfiler_gsea/GO/', graph_title, '.jpeg', sep = ''),
       units="in", width=10, height=10, res=500)
  gg_obj = ggplot(gseaGO_results, aes(reorder(Description, NES), NES)) +
    geom_col(aes(fill=p.adjust<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=graph_title) +
    theme_minimal()
  print(gg_obj)
  dev.off()

  ### other custom GMT files from MSigDB
  library(GSEABase)

  # Load in GMT file of gene sets (we downloaded from the Broad Institute for MSigDB)

  ## hallmark
  # c2 <- read.gmt('/Users/hhuang/Desktop/ketone_exp/ketone_DESeq2/msigDB/h.all.v7.2.entrez.gmt')
  # 
  # msig <- GSEA(foldchanges, TERM2GENE=c2, verbose=FALSE)
  # 
  # msig_df <- data.frame(msig)
  # msig_df = msig_df[msig_df$p.adjust<0.05,]
  # 
  # jpeg(paste('clusterProfiler_gsea/hallmark/', graph_title, '.jpeg', sep = ''),
  #      units="in", width=10, height=10, res=500)
  # gg_obj = ggplot(msig_df, aes(reorder(Description, NES), NES)) +
  #   geom_col(aes(fill=p.adjust<0.05)) +
  #   coord_flip() +
  #   labs(x="Pathway", y="Normalized Enrichment Score",
  #        title=graph_title) +
  #   theme_minimal()
  # print(gg_obj)
  # dev.off()
  ###########Functional class scoring##############################


  # ## SPIA ###
  # library(SPIA)
  # ## Significant genes is a vector of fold changes where the names are ENTREZ gene IDs.
  # ## The background set is a vector of all the genes represented on the platform.
  # background_entrez <- res_entrez$ENTREZID
  # sig_res_entrez <- res_entrez[which(res_entrez$padj < 0.05), ]
  # sig_entrez <- sig_res_entrez$log2FoldChange
  # names(sig_entrez) <- sig_res_entrez$ENTREZID
  # head(sig_entrez)
  #
  # spia_result <- spia(de=sig_entrez, all=background_entrez, organism="mmu")
  # # View(head(spia_result, n=20))
  # # plotP(spia_result, threshold=0.26)
  # # subset(spia_result, ID == "04060")
  #
  # spia_result = spia_result[spia_result$pGFdr < 0.05,]
  # jpeg(paste('SPIA/', graph_title, '.jpeg', sep = ''),
  #      units="in", width=10, height=10, res=500)
  # gg_obj = ggplot(spia_result, aes(reorder(Name, tA), tA)) +
  #   geom_col(aes(fill=pGFdr<0.05)) +
  #   coord_flip() +
  #   labs(x="Pathway", y="perturbation accumulation",
  #        title=graph_title) +
  #   theme_minimal()
  # print(gg_obj)
  # dev.off()
  ## SPIA ###

  # browser()
  #volcano plot
  jpeg(paste('volcano_plot/', graph_title, '.jpeg', sep = ''),
       units="in", width=10, height=10, res=500)

  # browser()
  res_M_i_no_treatment[order(res_M_i_no_treatment$log2FoldChange, decreasing = T),]
  volcano.plot = EnhancedVolcano(res_M_i_no_treatment,
                  # lab = res_M_i_no_treatment$SYMBOL,
                  lab = res_M_i_no_treatment$ensembl,
                  pCutoff = 0.05,
                  FCcutoff = 0,
                  x = 'log2FoldChange',
                  y = 'padj',
                  # legendLabels=c('Not sig.','Log (base 2) FC','padj',
                  #                'padj & Log (base 2) FC'),
                  legendPosition = 'bottom',
                  legendLabSize = 16,
                  legendIconSize = 5.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.2,
                  colConnectors = 'grey30')
  print(volcano.plot)
  dev.off()

  #volcano plot
  jpeg(paste('volcano_plot_symbol/', graph_title, '.jpeg', sep = ''),
       units="in", width=10, height=10, res=500)

  # browser()
  res_M_i_no_treatment[order(res_M_i_no_treatment$log2FoldChange, decreasing = T),]
  volcano.plot = EnhancedVolcano(res_M_i_no_treatment,
                                 lab = res_M_i_no_treatment$SYMBOL,
                                 # lab = res_M_i_no_treatment$ensembl,
                                 pCutoff = 0.05,
                                  FCcutoff = 0,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 # legendLabels=c('Not sig.','Log (base 2) FC','padj',
                                 #                'padj & Log (base 2) FC'),
                  legendPosition = 'bottom',
                  legendLabSize = 16,
                  legendIconSize = 5.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.2,
                  colConnectors = 'grey30')
  print(volcano.plot)
  dev.off()
  # 
  # # run_fgsea(res_M_i_no_treatment, genemap, geo.pathways, graph_title)
  # # run_fgsea(res_M_i_no_treatment, genemap, hallmark.pathways, graph_title)
  # # run_fgsea(res_M_i_no_treatment, genemap, kegg.pathways,graph_title)
}

# read dds data
library("DESeq2")
dds = readRDS('dds_object')

## load biomart for mapping between ensembl id to gene symbols
# library( "biomaRt" )
# ensembl = useMart( "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl" )
# genemap <- getBM( attributes = c("ensembl_gene_id_version", "hgnc_symbol"),
#                   mart = ensembl )

## load msigDB pathways
hallmark.pathways <- gmtPathways('/Users/hhuang/Desktop/ketone_exp/ketone_DESeq2/msigDB/h.all.v7.1.symbols.gmt')
kegg.pathways<- gmtPathways('/Users/hhuang/Desktop/ketone_exp/ketone_DESeq2/msigDB/c2.cp.kegg.v7.1.symbols.gmt')
geo.pathways <- gmtPathways('/Users/hhuang/Desktop/ketone_exp/ketone_DESeq2/msigDB/c5.all.v7.1.symbols.gmt')

# delete and make new directories
unlink( c('gsea_results', 'MAplots', 'sig_gene_heatmap', 'over_representation_analysis/GO/dotplots/', 
  'over_representation_analysis/GO/enrichment_clusters/', 'over_representation_analysis/GO/netplots/', 
  'over_representation_analysis/GO_mf/dotplots/', 'over_representation_analysis/GO_mf/enrichment_clusters/', 
  'over_representation_analysis/GO_mf/netplots/', 'over_representation_analysis/KEGG/dotplots/', 
  'over_representation_analysis/KEGG/enrichment_clusters/', 'over_representation_analysis/KEGG/netplots/', 
  'clusterProfiler_gsea/KEGG/', 'clusterProfiler_gsea/GO/', 'clusterProfiler_gsea/hallmark/', 'SPIA/', 
  'volcano_plot', 'volcano_plot_symbol'),
  recursive = T
)
sapply(c('gsea_results', 'MAplots', 'sig_gene_heatmap', 'over_representation_analysis/GO/dotplots/', 
  'over_representation_analysis/GO/enrichment_clusters/', 'over_representation_analysis/GO/netplots/', 
  'over_representation_analysis/GO_mf/dotplots/', 'over_representation_analysis/GO_mf/enrichment_clusters/', 
  'over_representation_analysis/GO_mf/netplots/', 'over_representation_analysis/KEGG/dotplots/', 
  'over_representation_analysis/KEGG/enrichment_clusters/', 'over_representation_analysis/KEGG/netplots/', 
  'clusterProfiler_gsea/KEGG/', 'clusterProfiler_gsea/GO/', 'clusterProfiler_gsea/hallmark/', 'SPIA/', 
  'volcano_plot', 'volcano_plot_symbol'),
  function(x) dir.create(x)
)

colData(dds)$condition
colData(dds)
run_pipeline('Treatment', 'Control')


# run_pipeline('iPSC_no', 'Microglia_no')

# assay(dds)
# colData(dds)
# assay(dds)[grep('ENSG00000070047', rownames(assay(dds))),]
# assay(dds)[grep('ENSG00000234287', rownames(assay(dds))),]
# 
# 
# dds.mtx = data.frame(assay(dds))
# dds.mtx
# short_ensembl = str_extract(rownames(dds.mtx), '[^.]+')
# dds.mtx$ensembl = short_ensembl
# 
# ens2symbol <- AnnotationDbi::select(org.Mm.eg.db,
#                                     key= short_ensembl, 
#                                     columns="SYMBOL",
#                                     keytype="ENSEMBL")
# 
# dds.mtx = inner_join(as_tibble(dds.mtx), as_tibble(ens2symbol), by=c("ensembl"="ENSEMBL"))
# dds.mtx[is.na(dds.mtx$SYMBOL), 'SYMBOL'] = dds.mtx[is.na(dds.mtx$SYMBOL), 'ensembl']
# dds.mtx$SYMBOL
