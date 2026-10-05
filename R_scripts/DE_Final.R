##################################
### DIFFERENTIAL ANALYSIS ########
##################################

# Loading Libraries 
library("tidyverse")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggrepel")
library("genefilter")
library("gplots")
library("dplyr")
library("reshape2")
library("ggplot2")
library("VennDiagram")
library("grid")
library("clusterProfiler")
library("org.Hs.eg.db")
library("ggalluvial")
library("patchwork")
library("here")
# TO RUN THIS SCRIPT YOU NEED TO USE FUNCTIONS INSIDE FILE functions_dds 
# RUN THAT FIRST 
#################################
# EDA-VARIABILITY  ##############
#################################

dds <- readRDS(here::here("data","FINAL_DDS.rds")
### Filtering  
keep <- rowSums(counts(dds) >= 10) >= min(table(dds$design))
dds_filtered <- dds[keep,]
design(dds_filtered) <- ~ donor + design   

# VSD TRANSFORMATION
vsd <- vst(dds_filtered, blind = FALSE)
head(assay(vsd),3)

# Clustering 
sampleDists <-dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$CC, vsd$design, vsd$donor, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$CC, vsd$design, vsd$donor, sep="-")

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

#Plot 
clustering_plot <- pheatmap(
  sampleDistMatrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  col = colors,
  fontsize = 14,          
  fontsize_col = 12,
  main = "Euclidean distance between samples",
  annotation_col = data.frame(Design = vsd$design,
                              row.names = colnames(sampleDistMatrix),
  annotation_legend = TRUE,
  legend = TRUE,
  silent = TRUE            

ggsave(
  filename = "clustering_plot.png",
  plot = clustering_plot$gtable,
  path = "results/DE_RESULTS",
  width = 16,
  height = 14,          
  dpi = 600,             # 600 for high resolution
  device = "png",
  bg = "white"           # better a white background for PNG
)


# PCA 
# ggplot
pcadata <- plotPCA(vsd, intgroup = c( "withEC", "duration"), returnData=TRUE)
percentVar <- round(100 * attr(pcadata, "percentVar"))
pcadata$sample_id <- sapply(strsplit(rownames(pcadata), "_"), function(x) x[2])
PCA_VST <- ggplot(pcadata, aes(PC1, PC2, color = duration, shape = withEC, label = sample_id)) + 
  geom_point(size = 4) +
  geom_text_repel(size = 4) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  ggtitle("PCA with VST data") +
  theme(
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 14),
    legend.key.size = unit(1.2, "cm"),      
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12)
  )

print(PCA_VST)

# PNG alta risoluzione (per slide/poster/PowerPoint)
ggsave(
  filename = "PCA_plot.png",
  plot = PCA_VST,
  path = "results/DE_RESULTS/",
  width = 16,
  height = 14,
  dpi = 600,
  device = "png",
  bg = "white"
)


####################################
##DIFFERENTIAL EXPRESSION ANALYSIS##
####################################
# WE RUN THE ANALYSIS USING ONLY THE PAIRED SAMPLE
# WE REMOVE THE THREE ADDITIONAL UNPAIRED SAMPLES AT 72 

# samples non paired
exclude <- c("CC09", "CC10", "CC12")
# mantaining only columns without those prefixes
dds_filtered <- dds_filtered[, !grepl(paste(exclude, collapse="|"),
                                      colnames(dds_filtered))]

# eliminate levels 
dds_filtered$donor  <- droplevels(dds_filtered$donor)
dds_filtered$design <- droplevels(dds_filtered$design)
table(dds_filtered$donor, dds_filtered$design)

# check
colnames(dds_filtered)
dim(dds_filtered)


# DEA AT DIFFERENT TIMEPOINTS 
design(dds_filtered) <- ~ donor + design
dds_final <-DESeq(dds_filtered)


## SAVE objects ##
saveRDS(dds_filtered,"results/DE_RESULTS/FINAL_DDS_LY.rds")
saveRDS(dds_final,"results/DE_RESULTS/DDS_RESULTS_LY.rds")

# EC vs CONTROL 
res_06_EC_CTR <- results(dds_final, contrast = c("design","06h_EC", "REST_noEC"))
res_24_EC_CTR <- results(dds_final, contrast = c("design","24h_EC", "REST_noEC"))
res_72_EC_CTR <- results(dds_final, contrast = c("design","72h_EC", "REST_noEC"))

# UP and DOWN regulated genes ( pval < 0.05, logfoldchange >< 2)

up_06_EC <- Up_genes(res_06_EC_CTR,2)
dw_06_EC <- Dw_genes(res_06_EC_CTR, 2)

up_24_EC <- Up_genes(res_24_EC_CTR,2)
dw_24_EC <- Dw_genes(res_24_EC_CTR,2)

up_72_EC <- Up_genes(res_72_EC_CTR,2)
dw_72_EC <- Dw_genes(res_72_EC_CTR,2)

### VENN DIAGRAMS ###
### useful for activation dynamics during time due to EC 
create_venn_3(up_06_EC,up_24_EC,up_72_EC,c("UP 06","UP 24", "UP 72"), "Upreguation wth EC")
create_venn_3(dw_06_EC,dw_24_EC,dw_72_EC,c("Dw 06","Dw 24", "Dw 72"), "downregulation wth EC")


### alluvia 
# alluvia plot function needs a df indicating UP, DW, NE genes
get_status <- function(res) {
  # Vector of "NE" (Not Expressed / Not Significant)
  status <- rep("NE", nrow(res))
  
  # Define UP: p-adj significant E log2FC positive
  idx_up <- which(res$padj < 0.05 & res$log2FoldChange > 2 & res$baseMean > 10 )
  status[idx_up] <- "UP"
  
  # Define DW: p-adj significative E log2FC negative
  idx_dw <- which(res$padj < 0.05 & res$log2FoldChange < -2 &  res$baseMean > 10 )
  status[idx_dw] <- "DW"
  
  return(status)
}

# DF creation 
alluvial_df <- data.frame(
  Gene = rownames(res_06_EC_CTR),
  R_06 = get_status(res_06_EC_CTR),
  R_24 = get_status(res_24_EC_CTR),
  R_72 = get_status(res_72_EC_CTR)
)

# Filtering 
df_summary_input <- alluvial_df %>%
  filter(R_06 != "NE" | R_24 != "NE" | R_72 != "NE") %>%
  group_by(R_06, R_24, R_72) %>%
  summarise(n = n(), .groups = "drop")


create_alluvia(df_summary_input,"LY activation LFC > 2, pval < 0.05 ")

#################################
#### DEA TIMEPOINTS vs CO-COLTURE
#################################

res_72 <- results(dds_final, contrast = c("design", "72h_EC", "72h_noEC"))
res_24 <- results(dds_final, contrast =c("design", "24h_EC", "24h_noEC"))
res_06 <- results(dds_final, contrast =c("design","06h_EC","06h_noEC"))

### ALLUVIA PLOT 
# get_status function 
get_status <- function(res) {
  # Creiamo un vettore di "NE" (Not Expressed / Not Significant)
  status <- rep("NE", nrow(res))
  
  # UP: p-adj significant E log2FC positive
  idx_up <- which(res$padj < 0.05 & res$log2FoldChange > 0.5)
  status[idx_up] <- "UP"
  
  # Definiamo i DW: p-adj significant E log2FC negative
  idx_dw <- which(res$padj < 0.05 & res$log2FoldChange < -0.5)
  status[idx_dw] <- "DW"
  
  return(status)
}

alluvial_df <- data.frame(
  Gene = rownames(res_06),
  R_06 = get_status(res_06),
  R_24 = get_status(res_24),
  R_72 = get_status(res_72)
)

df_summary_input <- alluvial_df %>%
  filter(R_06 != "NE" | R_24 != "NE" | R_72 != "NE") %>%
  group_by(R_06, R_24, R_72) %>%
  summarise(n = n(), .groups = "drop") # 'n' invece di 'Freq'


alluvia <- create_alluvia(df_summary_input,"LY activation")

ggsave(
  "results/DE_RESULTS/alluvia72.png",
  alluvia,
  width = 8,
  height = 6,
  dpi = 300
)

##################
## ANALYSIS AT 72h 
#################

res_72_shrunk <- lfcShrink(
  dds_final,
  contrast = c("design", "72h_EC", "72h_noEC"),
  type = "ashr" )


#### VOLCANO with LABELS
# Call
Volcano_dds(res_72_shrunk, "72h EC vs noEC", 1)

# poster 
p72 <- Volcano_dds(res_72_shrunk, "Coculture vs noCoculture", 1)
ggsave(
  "results/DE_RESULTS/volcano72.png",
  p72,
  width = 8,
  height = 6,
  dpi = 300
)



#####
Volcano_dds <- function(res, n_labels_up = 12, n_labels_down = 12) {
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ensembl") %>%
    mutate(
      symbol = mapIds(org.Hs.eg.db,
                      keys = ensembl,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first"),
      label = ifelse(!is.na(symbol), symbol, ensembl)
    )
  
  res_df$significance <- "NS"
  res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange >  1] <- "Up"
  res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"
  
  top_up   <- res_df %>% filter(significance == "Up") %>%
    arrange(padj, desc(log2FoldChange)) %>% head(n_labels_up)
  top_down <- res_df %>% filter(significance == "Down") %>%
    arrange(padj, log2FoldChange) %>% head(n_labels_down)
  to_label <- bind_rows(top_up, top_down)
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "#8B1E3F",
                                  "Down" = "#56B4E9",
                                  "NS" = "grey80")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05),  linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data          = to_label,
      aes(label     = label),
      size          = 3,
      fontface      = "bold",
      max.overlaps  = 20,
      box.padding   = 0.4,
      point.padding = 0.3,
      segment.color = "grey50",
      segment.size  = 0.3,
      nudge_x = ifelse(to_label$log2FoldChange > 0, 1, -1)
    ) +
    coord_cartesian(xlim = c(-7, 7), ylim = c(0, 55)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    labs(
      x = "log2 Fold Change Coculture vs no Coculture",
      y = "-log10 (adjusted p-value)"
    )
  return(p)
}

p2 <- Volcano_dds(res_72_shrunk)
ggsave(
  "results/DE_RESULTS/volcano72_try.png",
  p2,
  width = 8,
  height = 6,
  dpi = 300
)

#####

## HEATMAP BUT ONLY AT 72 
### up regulated genes 
up_genes <- subset(res_72, padj < 0.05 & log2FoldChange > 1)
up_genes <- up_genes[order(up_genes$stat, decreasing = TRUE), ]
up_genes$symbol <- mapIds(org.Hs.eg.db,
                          keys = rownames(up_genes),
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")
up_genes <- up_genes[!is.na(up_genes$symbol), ]

vsd <- vst(dds_filtered, blind = FALSE)
vsd_72h <- vsd[, colData(vsd)$design %in% c("72h_EC", "72h_noEC")]

# 2. top genes UP based on log2FoldChange
top_up_genes <- up_genes %>%
  as.data.frame() %>%
  rownames_to_column("ensembl") %>%
  arrange(desc(log2FoldChange)) %>%
  head(40) %>%          # top 40 
  pull(ensembl)

# Verify
cat("samples at 72h:", ncol(vsd_72h), "\n")
cat("Genes selected:", length(top_up_genes), "\n")


# function call 
hetmap_72_up <- Youwant_heatmap(
  vsd        = vsd_72h,
  gene_list  = top_up_genes,
  annotation = c("withEC", "donor"),
  order      = "design",
  title      = "Top 40 UP genes (by LFC) - 72h EC vs noEC"
)

ggsave(
  "results/DE_RESULTS/hetmap_72_up.png",
  hetmap_72_up,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)
### down regulated genes 
down_genes <- subset(res_72, padj < 0.05 & log2FoldChange < -1)
down_genes <- down_genes[order(down_genes$stat, decreasing = FALSE), ]
down_genes$symbol <- mapIds(org.Hs.eg.db,
                  keys=rownames(down_genes),
                  column="SYMBOL",
                  keytype="ENSEMBL",
                  multiVals="first")
down_genes <- down_genes[!is.na(down_genes$symbol),]


vsd_72h <- vsd[, colData(vsd)$design %in% c("72h_EC", "72h_noEC")]

# Selection top genes downregulated for log2FoldChange 
top_down_genes <- down_genes %>%
  as.data.frame() %>%
  rownames_to_column("ensembl") %>%
  arrange((log2FoldChange)) %>%
  head(40) %>%          #  top 40 per LFC
  pull(ensembl)

# Verify
cat("Campioni a 72h:", ncol(vsd_72h), "\n")
cat("Geni selezionati:", length(top_up_genes), "\n")

# Call for the function 
heatmap_72_dw <- Youwant_heatmap(
  vsd        = vsd_72h,
  gene_list  = top_down_genes,
  annotation = c("withEC", "donor"),
  order      = "design",
  title      = "Top 40 DOWN genes (by LFC) - 72h EC vs noEC"
)
ggsave(
  "results/DE_RESULTS/hetmap_72_dw.png",
  heatmap_72_dw,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)

#############
### GO ######
#############

# log2fold choiche
up_72 <- Up_genes(res_72, 1)
dw_72 <- Dw_genes(res_72, 1)

GO_72_Up <- run_GO_enrichment(up_72, " Up genes")
GO_72_Dw <- run_GO_enrichment(dw_72, "Down Genes")

#### GO WITH COMPARECLUSTER ####
#list with up and down regulated genes

gene_list_72 <- list(
  Up_regulated = up_72,
  Down_regulated = dw_72
)
orgdb = org.Hs.eg.db
# compareCluster per GO enrichment
cc_go_72 <- compareCluster(
  geneClusters = gene_list_72,
  fun = "enrichGO",
  OrgDb = orgdb,          # variabile orgdb
  keyType = "ENSEMBL",
  ont = "BP",             
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
# save for poster and results 
saveRDS(cc_go_72, "DE_RESULTS/cc_go_72.rds")

# Dotplot comparativo
go_72 <- dotplot(
  cc_go_72,
  showCategory = 10,
  label_format = 50,
  title = "GO Enrichment - Up vs Down regulated (72h)"
) +
  theme(
    axis.text.y = element_text(size = 10)
  )

