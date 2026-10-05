### ============
### FUNCTIONS
### ============

### ============
### HEATMAP
### ============
Youwant_heatmap <- function(vsd, gene_list, annotation, order = NULL, title) {
  
  # Remove NA rownames
  vsd <- vsd[!is.na(rownames(vsd)), ]
  
  # Filtering for present genes
  genes_present <- intersect(gene_list, rownames(assay(vsd)))
  cat("Requested genes:", length(gene_list), 
      "| genes found:", length(genes_present), "\n")
  
  mat <- assay(vsd)[genes_present, ]
  mat <- mat[complete.cases(mat), ]
  mat <- mat - rowMeans(mat)
  mat <- mat[matrixStats::rowVars(mat) > 0, ]
  
  # Annotation
  anno <- as.data.frame(colData(vsd)[, c(annotation, order), drop = FALSE])
  rownames(anno) <- colnames(mat)
  
  # ordering per variabile
  anno[[order]] <- factor(anno[[order]])
  ord <- order(anno[[order]])
  mat  <- mat[, ord]
  anno <- anno[ord, ]
  
  breaks <- seq(-3, 3, length.out = 101)
  
  pheatmap(
    mat,
    breaks         = breaks,
    annotation_col = anno,
    cluster_rows   = TRUE,
    cluster_cols   = FALSE,
    show_rownames  = FALSE,
    border_color   = NA,
    fontsize       = 10,
    main           = title
  )
}

  

### ============
### VENN DIAGRAMS
### ============

create_venn_with_percentages <- function(set1, set2, set3, 
                                         labels, 
                                         title,
                                         filename = NULL) {
  
  # Intersections
  A_only <- setdiff(set1, union(set2, set3))
  B_only <- setdiff(set2, union(set1, set3))
  C_only <- setdiff(set3, union(set1, set2))
  
  AB_only <- setdiff(intersect(set1, set2), set3)
  AC_only <- setdiff(intersect(set1, set3), set2)
  BC_only <- setdiff(intersect(set2, set3), set1)
  
  ABC <- Reduce(intersect, list(set1, set2, set3))
  
  union_all <- length(unique(c(set1, set2, set3)))
  
  # Percentages
  
  pA   <- length(A_only)  / union_all * 100
  pB   <- length(B_only)  / union_all * 100
  pC   <- length(C_only)  / union_all * 100
  pAB  <- length(AB_only) / union_all * 100
  pAC  <- length(AC_only) / union_all * 100
  pBC  <- length(BC_only) / union_all * 100
  pABC <- length(ABC)     / union_all * 100
  
  myCol <- brewer.pal(3, "Pastel2")
  
  # diagram
  venn_plot <- venn.diagram(
    x = list(set1, set2, set3),
    category.names = labels,
    fill = myCol,
    lwd = 2,
    lty = "blank",
    cex = 0.7,
    cat.cex = 0.4,
    fontface = "bold",
    filename = NULL 
  )
  
  # open to save the plot
  if(!is.null(filename)) {
    png(filename, width = 990, height = 990, res = 450)
  }
  
  # Draw
  grid.newpage()
  
  # Venn diagram
  grid.draw(venn_plot)
  
  # precentages in the plot
  grid.text(paste0("\n(", round(pA, 1), "%)"), x = 0.20, y = 0.70,gp = gpar(fontsize = 5))
  grid.text(paste0("\n(", round(pB, 1), "%)"), x = 0.80, y = 0.70,gp = gpar(fontsize = 5))
  grid.text(paste0("\n(", round(pC, 1), "%)"), x = 0.50, y = 0.15,gp = gpar(fontsize = 5))
  
  grid.text(paste0("\n(", round(pAB, 1), "%)"), x = 0.50, y = 0.80,gp = gpar(fontsize = 5))
  grid.text(paste0("\n(", round(pAC, 1), "%)"), x = 0.35, y = 0.40,gp = gpar(fontsize = 5))
  grid.text(paste0("\n(", round(pBC, 1), "%)"), x = 0.65, y = 0.40,gp = gpar(fontsize = 5))
  
  grid.text(paste0("\n(", round(pABC, 1), "%)"), x = 0.50, y = 0.55,gp = gpar(fontsize = 5))
  
  # Close visualisation
  if(!is.null(filename)) {
    dev.off()
  }
  
  # Print statistics
  cat("\n", title, "\n")
  cat("Total genes:", union_all, "\n")
  cat("Only 72h:", length(A_only), sprintf("(%.1f%%)", pA), "\n")
  cat("Only 24h:", length(B_only), sprintf("(%.1f%%)", pB), "\n")
  cat("Only 06h:", length(C_only), sprintf("(%.1f%%)", pC), "\n")
  cat("72h ∩ 24h:", length(AB_only), sprintf("(%.1f%%)", pAB), "\n")
  cat("72h ∩ 06h:", length(AC_only), sprintf("(%.1f%%)", pAC), "\n")
  cat("24h ∩ 06h:", length(BC_only), sprintf("(%.1f%%)", pBC), "\n")
  cat("All three:", length(ABC), sprintf("(%.1f%%)", pABC), "\n\n")
}


### ============
###  VENN DIAGRAM 
### ============

# Function for venn with TWO circles 

create_venn_2 <- function(set1, set2, 
                          labels, 
                          title,
                          filename = NULL) {
  
  # Intersections
  AB_only <- intersect(set1, set2)
  A_only <- setdiff(set1, intersect(set1, set2))
  B_only <- setdiff(set2, intersect(set1, set2))
  
  union_all <- length(unique(c(set1, set2)))
  
  # Percentages
  
  pA   <- length(A_only)  / union_all * 100
  pB   <- length(B_only)  / union_all * 100
  pAB  <- length(AB_only) / union_all * 100
  
  
  # diagram
  venn_plot <- venn.diagram(
    x = list(set1, set2),
    category.names = labels,
    fill = c("#4DAF4A", "#377EB8"),
    lwd = 2,
    lty = "blank",
    cex = 2,
    cat.cex = 1,
    fontface = "bold",
    filename = NULL )
  
  # open to save the plot
  if(!is.null(filename)) {
    png(filename, width = 990, height = 990, res = 450)
  }
  
  # Draw
  grid.newpage()
  
  # Venn diagram
  grid.draw(venn_plot)
  
  # precentages in the plot
  grid.text(paste0("\n(", round(pB, 1), "%)"), x = 0.08, y = 0.45,gp = gpar(fontsize = 9))
  grid.text(paste0("\n(", round(pA, 1), "%)"), x = 0.92, y = 0.45,gp = gpar(fontsize = 9))
  grid.text(paste0("\n(", round(pAB, 1), "%)"), x = 0.50, y = 0.40,gp = gpar(fontsize = 9))
  
  
  # Close visualisation
  if(!is.null(filename)) {
    dev.off()
  }
  
  # print Statstics
  cat("\n", title, "\n")
  cat("Total genes:", union_all, "\n")
  cat("Only EC:", length(A_only), sprintf("(%.1f%%)", pA), "\n")
  cat("Only noEC:", length(B_only), sprintf("(%.1f%%)", pB), "\n")
  cat("EC and noEC", length(AB_only), sprintf("(%.1f%%)", pAB), "\n")
}

### ==============================================================
### VENN DIAGRAM FOR UP AND DOWN REGULATED GENES ALL TOGETHER 
### ==============================================================

create_venn_3 <- function(set1, set2, set3, 
                          labels, 
                          title,
                          filename = NULL) {
  
  # Intersections
  A_only <- setdiff(set1, union(set2, set3))
  B_only <- setdiff(set2, union(set1, set3))
  C_only <- setdiff(set3, union(set1, set2))
  
  AB_only <- setdiff(intersect(set1, set2), set3)
  AC_only <- setdiff(intersect(set1, set3), set2)
  BC_only <- setdiff(intersect(set2, set3), set1)
  
  ABC <- Reduce(intersect, list(set1, set2, set3))
  
  union_all <- length(unique(c(set1, set2, set3)))
  
  # Percentages
  
  pA   <- length(A_only)  / union_all * 100
  pB   <- length(B_only)  / union_all * 100
  pC   <- length(C_only)  / union_all * 100
  pAB  <- length(AB_only) / union_all * 100
  pAC  <- length(AC_only) / union_all * 100
  pBC  <- length(BC_only) / union_all * 100
  pABC <- length(ABC)     / union_all * 100
  
  myCol <- brewer.pal(3, "Pastel2")
  
  # diagram
  venn_plot <- venn.diagram(
    x = list(set1, set2, set3),
    category.names = labels,
    fill = myCol,
    lwd = 2,
    lty = "blank",
    cex = 1.5,
    cat.cex = 1.5,
    fontface = "bold",
    filename = NULL 
  )
  
  # open to save the plot
  if(!is.null(filename)) {
    png(filename, width = 990, height = 990, res = 450)
  }
  
  # Draw
  grid.newpage()
  
  # Venn diagram
  grid.draw(venn_plot)
  
  # precentages in the plot
  grid.text(paste0("\n(", round(pA, 1), "%)"), x = 0.20, y = 0.70,gp = gpar(fontsize = 12))
  grid.text(paste0("\n(", round(pB, 1), "%)"), x = 0.80, y = 0.70,gp = gpar(fontsize = 12))
  grid.text(paste0("\n(", round(pC, 1), "%)"), x = 0.50, y = 0.15,gp = gpar(fontsize = 12))
  
  grid.text(paste0("\n(", round(pAB, 1), "%)"), x = 0.50, y = 0.80,gp = gpar(fontsize = 12))
  grid.text(paste0("\n(", round(pAC, 1), "%)"), x = 0.35, y = 0.40,gp = gpar(fontsize = 12))
  grid.text(paste0("\n(", round(pBC, 1), "%)"), x = 0.65, y = 0.40,gp = gpar(fontsize = 12))
  
  grid.text(paste0("\n(", round(pABC, 1), "%)"), x = 0.50, y = 0.55,gp = gpar(fontsize = 12))
  
  # Close visualisation
  if(!is.null(filename)) {
    dev.off()
  }
  
  # print statistics
  cat("\n", title, "\n")
  cat("Total genes:", union_all, "\n")
  cat("Only 72h:", length(A_only), sprintf("(%.1f%%)", pA), "\n")
  cat("Only 24h:", length(B_only), sprintf("(%.1f%%)", pB), "\n")
  cat("Only 06h:", length(C_only), sprintf("(%.1f%%)", pC), "\n")
  cat("72h ∩ 24h:", length(AB_only), sprintf("(%.1f%%)", pAB), "\n")
  cat("72h ∩ 06h:", length(AC_only), sprintf("(%.1f%%)", pAC), "\n")
  cat("24h ∩ 06h:", length(BC_only), sprintf("(%.1f%%)", pBC), "\n")
  cat("All three:", length(ABC), sprintf("(%.1f%%)", pABC), "\n\n")
}



### ==============
### VOLCANO PLOT 
### ==============

Volcano_dds <- function(res, title,thr_lgf,  n_labels_up = 12, n_labels_down = 12) {
  

  # df converstion and addition of symbols
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ensembl") %>%
    mutate(
      symbol = mapIds(org.Hs.eg.db,
                      keys = ensembl,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first"),
      # symbol, ifelse ensembl
      label = ifelse(!is.na(symbol), symbol, ensembl)
    )

  # Significance and classification 
  res_df$significance <- "NS"
  res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > thr_lgf] <- "Up"
  res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < - thr_lgf] <- "Down"
  
  # Seleziona i top geni da etichettare
  # top genes to lable 
  # top UP: most significant with higher LFC
  top_up <- res_df %>%
    filter(significance == "Up") %>%
    arrange(padj, desc(log2FoldChange)) %>%
    head(n_labels_up)
  
  # top DOWN: most significant with lower LFC
  top_down <- res_df %>%
    filter(significance == "Down") %>%
    arrange(padj, log2FoldChange) %>%
    head(n_labels_down)
  
  # Union of labelling genes
  to_label <- bind_rows(top_up, top_down)
  
  # Plot
  Volcano_plot <- ggplot(res_df, 
                         aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significance), 
               alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up"   = "#8B1E3F", 
                                  "Down" = "#56B4E9", 
                                  "NS"   = "grey80")) +
    geom_vline(xintercept = c(-0.5, 0.5), 
               linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), 
               linetype = "dashed", color = "grey40") +
    # labels with ggrepel
    geom_text_repel(
      data = to_label,
      aes(label = label),
      size          = 3,
      fontface      = "bold",
      max.overlaps  = 20,
      box.padding   = 0.4,
      point.padding = 0.3,
      segment.color = "grey50",
      segment.size  = 0.3,
      # Up and DOWN in different parts
      nudge_x = ifelse(to_label$log2FoldChange > 0, 0.5, -0.5)
    ) +
    coord_cartesian(
      xlim = c( -7, 7),
      ylim = c(0, 50)
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title  = paste("Volcano plot", title),
      x      = "log2 Fold Change",
      y      = "-log10 adjusted p-value",
      color  = "Regulation"
    )
  
  return(Volcano_plot)
}


### ==============
### UP AND DOWN REGULATED GENES
### ==============
# function for up and down regulated genes (you decide the threshold)
  Up_genes <- function(res_df, thr) {
    rownames(res_df)[
      !is.na(res_df$padj) & !is.na(res_df$log2FoldChange)&
        res_df$padj < 0.05 & res_df$log2FoldChange > thr]
  }
  Dw_genes <- function(res_df, thr) {
    rownames(res_df)[
      !is.na(res_df$padj) & !is.na(res_df$log2FoldChange) &
        res_df$padj < 0.05 & res_df$log2FoldChange < -thr
    ]
  }
  
  
  
### ==============
### GSEA 
### ==============
  
  
GSEA_PLOT <- function(dds_res, plot_title) {
    
  # Prepare ranked gene list
    res_gse <- dds_res %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      filter(!is.na(stat))
    
    
    geneList <- res_gse$stat
    names(geneList) <- res_gse$gene
    geneList <- sort(geneList, decreasing = TRUE)
    
    # GSEA son GO Biological Process
    gse_result <- gseGO(
      geneList      = geneList,
      ont           = "BP",
      keyType       = "ENSEMBL",
      OrgDb         = org.Hs.eg.db,
      minGSSize     = 10,
      maxGSSize     = 500,
      pvalueCutoff  = 0.1,
      eps           = 0,           
      scoreType     = "std",       # corrected for stat with values +/-
      pAdjustMethod = "BH",
      verbose       = FALSE
    )
    
    print(
      dotplot(gse_result,
              showCategory = 8,
              split = ".sign",
              font.size = 9,
              title = paste(plot_title, "(padj < 0.1)")) +
        facet_grid(. ~ .sign)
    )
    
    return(gse_result)
}
  
  
### ==============
### GO ANALYSIS 
### ==============
  
run_GO_enrichment <- function(gene_list,
                                title_prefix = "",
                                ont = "ALL",
                                show_n = 10,
                                orgdb = org.Hs.eg.db) {
    
    # Enrichment
    go_res <- enrichGO(
      gene = gene_list,
      OrgDb = orgdb,
      keyType = "ENSEMBL",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    
    # Barplot
    bar_plot <- barplot(
      go_res,
      color = "p.adjust",
      title = paste0("Top ", show_n, " GO Enrichment - ", title_prefix),
      showCategory = show_n,
      label_format = 50
    ) #+ theme_minimal(base_size = 9)
    
    print(bar_plot)
    
    # Dotplot
    dot_plot <- dotplot(
      go_res,
      color = "p.adjust",
      title = paste0("Top ", show_n, " GO Enrichment - ", title_prefix),
      showCategory = show_n,
      label_format = 50
    )
    
    print(dot_plot)
    
    return(go_res)
  }
  
### ==============
### ALLUVIA PLOT 
### ==============
# before running this function, data must be prepared in long format
  
  create_alluvia <- function(df_summary,title,filename = NULL)
  {
    
    count_06 <- as.data.frame(table(df_summary$R_06))
    count_24 <- as.data.frame(table(df_summary$R_24))
    count_72 <- as.data.frame(table(df_summary$R_72))
    
    colnames(count_06) <- c("stratum","n")
    colnames(count_24) <- c("stratum","n")
    colnames(count_72) <- c("stratum","n")
    
    count_06$time <- "06h"
    count_24$time <- "24h"
    count_72$time <- "72h"
    counts_all <- rbind(count_06, count_24, count_72)
    
    my_colors <- c(
      UP = "#8B1E3F",
      NE = "#fdf0d5",
      DW = "#56B4E9"
    )
    
   
    ggplot(as.data.frame(df_summary),
           aes(y = n,
               axis1 = R_06, axis2 = R_24, axis3 = R_72)) +
      geom_alluvium(aes(fill = R_06), width = 1/8) +
      
      # (stratum) colored by gene type
      geom_stratum(width = 1/8, aes(fill = after_stat(stratum)), alpha = 0.9) +
      
      # labels
      geom_text(stat = "stratum",
                aes(label = paste0(after_stat(stratum),
                                   "\n", after_stat(count))),
                size = 3,
                color = "black")  +
      scale_fill_manual(values = my_colors) +
      scale_x_discrete(limits = c("06h","24h","72h")) +
      guides(fill = "none") +
      theme_minimal(base_size = 10) +
      ggtitle(title)
    
  }
  
