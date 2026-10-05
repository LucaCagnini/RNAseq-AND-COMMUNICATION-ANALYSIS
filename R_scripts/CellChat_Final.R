########################
###CELL CHAT ANALYSIS###
########################

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
library("clusterProfiler")
library("org.Hs.eg.db")
library("readxl")
library("here")

# Installation CellChat from official repo
remotes::install_github("sqjin/CellChat")
library("CellChat")
library("patchwork")
library("readxl")
library("edgeR")

#####################
### LOAD DATASETS ###
#####################
# For cell to cell communciaiton anyalsis we need as many RDS objects as many are the RNAseq cell types we have
#in this care we have two cell types; lymphocytes (LY) and beta cells (EC)
# RNAseq lymphocytes
dds_LY  <- readRDS(here::here("data", "FINAL_DDS.rds"))
degs_LY <- readRDS(here::here("data", "DDS_RESULTS.rds"))

# RNAseq EndoC
dds_EC  <- readRDS(here::here("data", "DDS_EC"))
degs_EC <- readRDS(here::here("data", "DEGS_EC"))  

# we integrated the RNAseq data with data coming from a cytookine luminex panel, measuring the concentration of 71 cytokine 
#in the coculture media 
data    <- read.table(here::here("cytokine_array", "Cytokines_Concentration_Clean.txt"),
                      header = TRUE, sep = "\t")
Cy_list <- read_excel(here::here("cytokine_array", "Table_genes_repectors (1).xlsx"))
Cy_list <- as.data.frame(Cy_list)

df_luminex <- data %>%
  group_by(Cytokine, Condition) %>%
  summarise(
    Mean_Value = mean(Value, na.rm = TRUE),
    SD_Value = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )


###################################
### SCATTER PLOT CYTOKINE CROSSTALK 
###################################

##### NORMALISATION 
#separated normalisation for RNAseq from EC or from LY
dds_EC <- estimateSizeFactors(dds_EC)
dds_LY <- estimateSizeFactors(dds_LY)

counts_EC_norm <- counts(dds_EC, normalized = TRUE)
counts_LY_norm <- counts(dds_LY, normalized = TRUE)

common_genes <- intersect(rownames(counts_EC_norm), rownames(counts_LY_norm))
counts_EC_norm <- counts_EC_norm[common_genes, ]
counts_LY_norm <- counts_LY_norm[common_genes, ]
counts_norm_joint <- cbind(
  counts_EC_norm,
  counts_LY_norm
)

# meta data 
meta_EC <- as.data.frame(colData(dds_EC)) %>%
  rownames_to_column("sample_id") %>%
  mutate(cell_type = "Beta")

meta_LY <- as.data.frame(colData(dds_LY)) %>%
  rownames_to_column("sample_id") %>%
  mutate(cell_type = "Lymphocyte")

meta_joint <- bind_rows(meta_EC, meta_LY)

##### ABBUNDANCE FUNCTION
# function to calculate abbundance for timepoint

get_abundance_by_time <- function(counts_EC_matrix, counts_LY_matrix, time_label) {
  
  cols_EC <- grep(paste0("CD3_", time_label), colnames(counts_EC_matrix), value = TRUE)
  cols_LY <- grep(paste0(time_label, "_EC"),  colnames(counts_LY_matrix), value = TRUE)
  
  cat("Timepoint:", time_label,
      "| EC samples:", length(cols_EC),
      "| LY samples:", length(cols_LY), "\n")
  
  common_genes <- intersect(rownames(counts_EC_matrix), rownames(counts_LY_matrix))
  
  data.frame(
    ensembl_id   = common_genes,
    EC_Abundance = rowMeans(counts_EC_matrix[common_genes, cols_EC, drop = FALSE]),
    LY_Abundance = rowMeans(counts_LY_matrix[common_genes, cols_LY, drop = FALSE])
  )
}

df_ab_06h <- get_abundance_by_time(counts_EC_norm, counts_LY_norm, "6h")
df_ab_24h <- get_abundance_by_time(counts_EC_norm, counts_LY_norm, "24h")
df_ab_72h <- get_abundance_by_time(counts_EC_norm, counts_LY_norm, "72h")

##### DESEQ2 RESULTS 
res_EC     <- DESeq(dds_EC)
res_EC_06h <- results(res_EC, contrast = c("design", "CD3_6h",  "Endo"))
res_EC_24h <- results(res_EC, contrast = c("design", "CD3_24h", "Endo"))
res_EC_72h <- results(res_EC, contrast = c("design", "CD3_72h", "Endo"))

res_LY     <- DESeq(dds_LY)
res_LY_06h <- results(res_LY, contrast = c("design", "06h_EC", "06h_noEC"))
res_LY_24h <- results(res_LY, contrast = c("design", "24h_EC", "24h_noEC"))
res_LY_72h <- results(res_LY, contrast = c("design", "72h_EC", "72h_noEC"))

#### DATA PREPARATION 
# function for the creation of the dataset integrating logfoldchange ratio and abundance between EC and LY
prepare_ratio_data <- function(ly_res, be_res, abundance_df, lum_data, time_label, min_prot = 5) {
  
  # Mapping symbol if not present
  if (!"symbol" %in% colnames(abundance_df)) {
    abundance_df$symbol <- mapIds(org.Hs.eg.db,
                                  keys      = abundance_df$ensembl_id,
                                  column    = "SYMBOL",
                                  keytype   = "ENSEMBL",
                                  multiVals = "first")
  }
  #lymphocites 
  ly_df <- as.data.frame(ly_res) %>%
    rownames_to_column("ensembl_id") %>%
    dplyr::select(ensembl_id, Ly_LFC = log2FoldChange)
  # beta cells
  be_df <- as.data.frame(be_res) %>%
    rownames_to_column("ensembl_id") %>%
    dplyr::select(ensembl_id, Be_LFC = log2FoldChange)
  # cytokine 
  df <- Cy_list %>%
    dplyr::select(Gene = `CYT_GENE_SYMBOL_(Ligand)`, Luminex_Name = CYT_ALIAS) %>%
    distinct() %>%
    inner_join(abundance_df, by = c("Gene" = "symbol")) %>%
    inner_join(be_df, by = "ensembl_id") %>%
    inner_join(ly_df, by = "ensembl_id") %>%
    left_join(lum_data %>% filter(Condition == time_label),
              by = c("Luminex_Name" = "Cytokine"))
  
  # for time lable = 06, we don't have luminex data
  if (time_label != "6h") {
    df <- df %>% filter(!is.na(Mean_Value) & Mean_Value > min_prot)
  }
  
  df %>% mutate(
    Timepoint        = time_label,
    Expression_Ratio = log2((LY_Abundance + 1) / (EC_Abundance + 1)),
    Ly_LFC           = ifelse(is.na(Ly_LFC), 0, Ly_LFC),
    Be_LFC           = ifelse(is.na(Be_LFC), 0, Be_LFC),
    LFC_Ratio        = abs(Ly_LFC) - abs(Be_LFC)
  )
}

### ASSEMBLY
df_final_06h <- prepare_ratio_data(res_LY_06h, res_EC_06h, df_ab_06h, df_luminex, "6h")
df_final_24h <- prepare_ratio_data(res_LY_24h, res_EC_24h, df_ab_24h, df_luminex, "24h")
df_final_72h <- prepare_ratio_data(res_LY_72h, res_EC_72h, df_ab_72h, df_luminex, "72h")
df_master       <- bind_rows(df_final_06h, df_final_24h, df_final_72h)

##################################
## CYTOKINE CROSSTALK - 24h vs 72h 
##################################
# plot showing differences in therms of Expression-Ration and color timepoint in two timepoints

df_plot_complete <- bind_rows(df_final_24h, df_final_72h)
ggplot(df_plot_complete, aes(x = LFC_Ratio, y = Expression_Ratio, color = Timepoint)) +
  
  # Axes (0,0)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  
  # Points
  geom_point(data = function(x) filter(x, Mean_Value > 5 | abs(Expression_Ratio) > 2),
             aes(size = Mean_Value), alpha = 0.6) +
  
  # lables
  geom_text_repel(data = function(x) filter(x, Mean_Value > 5 | abs(Expression_Ratio) > 2),
                  aes(label = Gene),
                  size = 2.5, fontface = "bold", max.overlaps = 10) +
  
  
  # scales
  scale_size_continuous(range = c(2, 10), name = "Protein (pg/ml)") +
  scale_color_manual(values = c("24h" = "steelblue", "72h" = "darkred")) +
  
  labs(
    title = "Cellular Dominance Landscape: Expression vs. Reactivity",
    subtitle = "Y-axis: log2(Ly/Beta) Abundance | X-axis: Ly-Beta ΔLFC",
    x = "Differential Reactivity (Ly_LFC - Be_LFC)",
    y = "Abundance Dominance log2((LY_Ab + 1) / (EC_Ab + 1))"
  ) +
  theme_minimal() 

###########################################
# CYTOKINES CROSSTALK - all citokines at 72h
###########################################
# All cytokines at 72h Luminex
df_luminex_72h <- df_luminex %>%
  filter(Condition == "72h") 

cat("Total cytokines at 72h from Luminex:", nrow(df_luminex_72h), "\n")
Cy_list_expanded <- Cy_list %>%
  separate_rows(`CYT_GENE_SYMBOL_(Ligand)`, sep = ",\\s*")

# Join con Cy_list to get gene symbol
df_72h_full <- Cy_list_expanded %>%
  dplyr::select(Gene = `CYT_GENE_SYMBOL_(Ligand)`, 
                Luminex_Name = CYT_ALIAS) %>%
  distinct() %>%
  inner_join(df_luminex_72h, by = c("Luminex_Name" = "Cytokine")) %>%
  # ensembl_id
  mutate(ensembl_id = mapIds(org.Hs.eg.db,
                             keys = Gene,
                             column = "ENSEMBL",
                             keytype = "SYMBOL",
                             multiVals = "first"))

cat("Citochine with mapped gene symbol :", nrow(df_72h_full), "\n")

# abundance EC e LY a 72h
df_72h_full <- df_72h_full %>%
  left_join(df_ab_72h, by = "ensembl_id") %>%
  # Add LFC EC and LY at 72h
  left_join(
    as.data.frame(res_EC_72h) %>%
      rownames_to_column("ensembl_id") %>%
      dplyr::select(ensembl_id, Be_LFC = log2FoldChange),
    by = "ensembl_id"
  ) %>%
  left_join(
    as.data.frame(res_LY_72h) %>%
      rownames_to_column("ensembl_id") %>%
      dplyr::select(ensembl_id, Ly_LFC = log2FoldChange),
    by = "ensembl_id"
  ) %>%
  mutate(
    EC_Abundance     = replace_na(EC_Abundance, 0),
    LY_Abundance     = replace_na(LY_Abundance, 0),
    Be_LFC           = replace_na(Be_LFC, 0),
    Ly_LFC           = replace_na(Ly_LFC, 0),
    Expression_Ratio = log2((LY_Abundance + 1) / (EC_Abundance + 1)),
    LFC_Ratio        = abs(Ly_LFC) - abs(Be_LFC)
  )

cat("Citokines df final 72h:", nrow(df_72h_full), "\n")
cat("Range Expression_Ratio:", range(df_72h_full$Expression_Ratio, na.rm=TRUE), "\n")
cat("Range LFC_Ratio:", range(df_72h_full$LFC_Ratio, na.rm=TRUE), "\n")

#### ANNOTATION COMUNICATION
# this part of the plot tells me which cytokine is implicated with the crosstalk, and uses the results of the cellchat-nichenet analysis
# if you make changes at the code downhill you have to re-run this later

ligand_receptor_DE_LY <- read.csv("LigandReceptor_LY_to_EC.csv", header = TRUE)
ligand_receptor_DE_EC <- read.csv("LigandReceptor_EC_to_LY.csv", header = TRUE)

ligands_EC_in_luminex <- ligand_receptor_DE_EC %>%
  filter(Luminex == TRUE) %>%
  pull(ligand) %>% unique()

ligands_LY_in_luminex <- ligand_receptor_DE_LY %>%
  filter(Luminex == TRUE) %>%
  pull(ligand) %>% unique()

df_72h_annotated <- df_72h_full %>%
  mutate(
    communication = case_when(
      Gene %in% ligands_EC_in_luminex & 
        Gene %in% ligands_LY_in_luminex ~ "Both",
      Gene %in% ligands_EC_in_luminex   ~ "β Cell → LY",
      Gene %in% ligands_LY_in_luminex   ~ "LY → β Cell",
      TRUE                              ~ "Not involved"
    )
  )

# Verify
cat("\n--- Annotation comunication ---\n")
cat("EC→LY:", sum(df_72h_annotated$communication == "β Cell → LY"), "\n")
cat("LY→EC:", sum(df_72h_annotated$communication == "LY → β Cell"), "\n")
cat("Both:", sum(df_72h_annotated$communication == "Both"), "\n")
cat("Not involved:", sum(df_72h_annotated$communication == "Not involved"), "\n")

#################
# PLOT — colors
################

comm_colors <- c(
  "β Cell → LY"   = "#E89C4A",
  "LY → β Cell"   = "#1D9E75",
  "Both"          = "#8F6BC7",
  "Not involved"  = "grey80"
)

#  visible points
soglia_plot <- function(x) {
  filter(x, Mean_Value > 5 | abs(Expression_Ratio) > 2)
}
p_cytokine <- ggplot(df_72h_annotated,
                     aes(x = LFC_Ratio, y = Expression_Ratio)) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  
  annotate("text", x = Inf, y = Inf,
           label = "Lymphocytes Dominant production", hjust = 1.05, vjust = 1.5, color = "#1D9E75",
            size = 5, fontface = "bold") + 
  annotate("text", x = -Inf, y = -Inf,
           label = "\u03b2-Cell Dominant production", hjust = -0.1, vjust = -0.5, color = "#E89C4A",
            size = 5, fontface = "bold") +
  
  geom_point(
    data = function(x) soglia_plot(x) %>% filter(communication == "Not involved"),
    aes(size = Mean_Value),
    color = "grey80", alpha = 0.15, shape = 16
  ) +
  geom_point(
    data = function(x) soglia_plot(x) %>% filter(communication != "Not involved"),
    aes(size = Mean_Value, color = communication),
    alpha = 0.85, shape = 16
  ) +
  geom_text_repel(
    data = function(x) soglia_plot(x),
    aes(label = Gene, color = communication),
    size = 3,
    fontface = "bold",
    max.overlaps = 25,
    box.padding = 0.6,
    segment.size = 0.4,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = comm_colors,
    name = "Communication direction",
    breaks = c("β Cell → LY", "LY → β Cell", "Both", "Not involved")
  ) +
  scale_size_continuous(
    range = c(6, 16),
    name = "Protein concentration\n(pg/ml)",
    breaks = c(100, 2000, 4000, 8000),
    labels = c("1000", "2000", "4,000", "8,000")
  ) +
  
  guides(
    color = guide_legend(
      title = "Communication direction (Cellchat + Nichenet) ",
      override.aes = list(size = 5),
      order = 1
    ),
    size = guide_legend(
      title = "Protein concentration\n(pg/ml)",
      order = 2
    )
  ) +
  coord_cartesian(
    ylim = c(min(df_72h_annotated$Expression_Ratio, na.rm = TRUE) - 0.5,
             max(df_72h_annotated$Expression_Ratio, na.rm = TRUE) + 1.5),
    clip = "off"
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",      
    legend.spacing.y = unit(0.3, "cm"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10)
  ) + 
  labs(
    x = "Reactivity (|Lymphocytes| - |β Cells|)",
    y = "RNAseq: log2(Lymphocytes / β Cells)"
  )
ggsave(
  filename = "results/DE_RESULTS/cytokine_landscape_72h.png",
  plot = p_cytokine,
  width = 11, height = 8, dpi = 300
)

########################
##  CELLCHAT ANALYSIS ##
########################

### METADATA
# using the df from normalised count before
meta_final <- meta_joint %>%
  mutate(
    timepoint = str_extract(sample_id, "\\d{1,2}h"),
    timepoint = case_when(
      timepoint == "00h" ~ "0h",
      is.na(timepoint) ~ "0h",
      timepoint == "06h" ~ "6h",
      TRUE ~ timepoint
    ),
    condition = case_when(
      str_detect(sample_id, "_noEC") ~ "Mono",
      str_detect(sample_id, "Endo_0h") ~ "Mono",
      str_detect(sample_id, "_EC") ~ "CoCulture",
      str_detect(sample_id, "CD3_") & !str_detect(sample_id, "noEC") ~ "CoCulture",
      TRUE ~ "Mono"
    ),
    cell_group = paste(cell_type, condition, sep = "_")
  ) %>%
  mutate(timepoint = factor(timepoint, levels = c("0h", "6h", "24h", "72h"))) %>%
  as.data.frame()

rownames(meta_final) <- meta_final$sample_id

# I need it for nichenet
write.csv(meta_final, "results/CELLCHAT_RESULTS/meta_final.csv")

############################
### NORMALIATION FOR CELLCHAT 
#############################
#(different from previous plot)

## RAW COUNTS
counts_EC_raw <- counts(dds_EC, normalized = FALSE)
counts_LY_raw <- counts(dds_LY, normalized = FALSE)
common_genes <- intersect(
  rownames(counts_EC_raw),
  rownames(counts_LY_raw)
)
combined_raw <- cbind(
  counts_EC_raw[common_genes, ],
  counts_LY_raw[common_genes, ]
)

# ENSEMBL → SYMBOL
gene_map <- mapIds(
  org.Hs.eg.db,
  keys      = rownames(combined_raw),
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

df_raw <- as.data.frame(combined_raw)
df_raw$symbol <- gene_map

# AGGREGATION DUPLICATED GENES
df_aggregated <- df_raw %>%
  filter(!is.na(symbol), symbol != "") %>%
  group_by(symbol) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  as.data.frame()  

rownames(df_aggregated) <- df_aggregated$symbol
df_aggregated$symbol <- NULL
mat_raw_aggregated <- as.matrix(df_aggregated)

# Verify
cat("first 10 rownames:", rownames(mat_raw_aggregated)[1:10], "\n")

cat("Dimensions after aggregation:",
    dim(mat_raw_aggregated)[1],
    "genes x",
    dim(mat_raw_aggregated)[2],
    "samples\n")

#TMM NORMALIZATION (edgeR)
dge <- DGEList(counts = mat_raw_aggregated)
dge <- calcNormFactors(
  dge,
  method = "TMM"
)

#CPM + LOG TRANSFORMATION
cpm_norm <- cpm(
  dge,
  normalized.lib.sizes = TRUE,
  log = FALSE
)

data_cellchat <- log1p(cpm_norm)
cat("Expression range:", range(data_cellchat), "\n" )

# METADATA ANLIGNMENT
common_samples <- intersect(
  colnames(data_cellchat),
  rownames(meta_final)
)

data_cellchat <- data_cellchat[, common_samples]
meta_final <- meta_final[
  common_samples,
]

# ONLY COCULTURE
meta_coculture <- meta_final %>%
  filter(condition == "CoCulture")

data_coculture <- data_cellchat[
  ,
  rownames(meta_coculture)
]

# CHECK
cat("Raw EC counts range:", range(counts_EC_raw),"\n")
cat( "Raw LY counts range:", range(counts_LY_raw), "\n" )
cat( "CellChat matrix range:",  range(data_cellchat),  "\n" )
cat( "Number of genes:", nrow(data_cellchat), "\n" )
cat( "Number of samples:", ncol(data_cellchat), "\n" )


####### WE TAKE THE THREE SAMPLES AT 72 E TRE A 24
#3 donors
donatori_comuni <- c("CC24_CD4T05", "CC30_CD4T02", "CC30_CD4T04")

#################
### CELLCHAT LOOP 
#################

# Loop
object.list <- list()
tempi_validi <- c("24h", "72h")
CellChatDB <- CellChatDB.human

for (t in tempi_validi) {
  message(paste("\n--- Elaborazione Timepoint:", t, "---"))
  
  meta_sub <- meta_final %>%
    filter(timepoint == t, condition == "CoCulture") %>%
    filter(
      cell_type == "Beta" |
        (cell_type == "Lymphocyte" & 
           # only common donors with str_detect
           Reduce("|", lapply(donatori_comuni, 
                              function(d) str_detect(sample_id, d))))
    )
  
  cat("Samples Beta:", sum(meta_sub$cell_type == "Beta"),
      "| Samples Lymphocyte:", sum(meta_sub$cell_type == "Lymphocyte"), "\n")
  
  print(meta_sub$sample_id[
    meta_sub$cell_type == "Lymphocyte"
  ])
  
  if(length(unique(meta_sub$cell_type)) < 2) {
    message(paste("Salto", t, "- missing cell type"))
    next
  }
  
  data_sub <- data_cellchat[, rownames(meta_sub)]
  cellchat <- createCellChat(object = data_sub,
                             meta   = meta_sub,
                             group.by = "cell_type")
  cellchat@DB <- subsetDB(CellChatDB,
                          search = c("Secreted Signaling"
                                     ,"ECM-Receptor"
                                     ))
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat,
                                         thresh.p = 1,
                                         thresh.fc = 0.1)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat,
                                raw.use = TRUE,
                                type    = "triMean")
  
  if (!is.null(cellchat@net$prob) && sum(cellchat@net$prob) > 0) {
    cellchat <- filterCommunication(cellchat, min.cells = 1)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    object.list[[t]] <- cellchat
    message(paste("Succes for", t,
                  "- Interaction foung:", sum(cellchat@net$count)))
  } else {
    message(paste("No interaction found for", t))
  }
}

# check results from loop 
cat("Timepoint completed:", length(object.list), "\n")
cat("Names:", names(object.list), "\n")

############################
# ANALYSIS : SETUP AND MERGE
###########################

names(object.list) <- c("CoCulture_24h", "CoCulture_72h")

cellchat.24h <- object.list[["CoCulture_24h"]]
cellchat.72h <- object.list[["CoCulture_72h"]]

# Pathway analysis
pathways_24h <- cellchat.24h@netP$pathways
pathways_72h <- cellchat.72h@netP$pathways

pathway_comuni    <- intersect(pathways_24h, pathways_72h)
pathway_solo_24h  <- setdiff(pathways_24h, pathways_72h)
pathway_solo_72h  <- setdiff(pathways_72h, pathways_24h)

cat("Pathway comuni:", length(pathway_comuni), "\n")
cat("Solo a 24h:", paste(pathway_solo_24h, collapse=", "), "\n")
cat("Solo a 72h:", paste(pathway_solo_72h, collapse=", "), "\n")

# Save pathway
saveRDS(list(
  comuni   = pathway_comuni,
  solo_24h = pathway_solo_24h,
  solo_72h = pathway_solo_72h
), "results/CELLCHAT_RESULTS/pathway_timepoint_specific.rds")

# Merge
cellchat.merged <- mergeCellChat(
  object.list,
  add.names = names(object.list),
  cell.prefix = TRUE
)
saveRDS(cellchat.merged,
        "results/CELLCHAT_RESULTS/LY_EC_CELLCHAT_merged.rds")
saveRDS(object.list,
        "results/CELLCHAT_RESULTS/LY_EC_CELLCHAT_list.rds")

#######################
# EXTRACT INTERACTIONS
#######################

prepare_interactions <- function(cellchat_obj, timepoint_label) {
  subsetCommunication(cellchat_obj) %>%
    mutate(
      interaction_name = paste(ligand, "->", receptor),
      source    = as.character(source),
      target    = as.character(target),
      timepoint = timepoint_label
    ) %>%
    filter(
      pval < 0.05,
      annotation %in% c("Secreted Signaling" # justify witch db used
                        , "ECM-Receptor"
                        )
    )
}

df_24h <- prepare_interactions(cellchat.24h, "24h")
df_72h <- prepare_interactions(cellchat.72h, "72h")
df_all <- bind_rows(df_24h, df_72h)

df_all %>%
  group_by(timepoint, annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  print()

write.csv(df_all,
          "results/CELLCHAT_RESULTS/LY_EC_CELLCHAT.csv",
          row.names = FALSE)

###############################
# PLOT 1: Overview interactions
###############################
gg1 <- compareInteractions(cellchat.merged, 
                           show.legend = TRUE, group = c(1,2))
gg2 <- compareInteractions(cellchat.merged, 
                           show.legend = TRUE, group = c(1,2),
                           measure = "weight")
gg1 + gg2

#######################################################
# PLOT 2: Information flow — common pathway 24h vs 72h
#######################################################
df_ranknet <- rankNet(cellchat.merged,
                      mode        = "comparison",
                      stacked     = FALSE,
                      do.stat     = TRUE,
                      return.data = TRUE)$signaling.contribution

soglia_75 <- quantile(df_ranknet$contribution[df_ranknet$contribution > 0], 0.75)

df_ranknet_filtered <- df_ranknet %>%
  filter(contribution > 0,
         !is.infinite(contribution.relative.1)) %>%
  group_by(name) %>%
  filter(n() == 2) %>%                              # present in both 
  filter(max(contribution) > soglia_75) %>%         # top 25%
  ungroup() %>%
  mutate(group = case_when(
    group == "CoCulture_24h" ~ "24h",
    group == "CoCulture_72h" ~ "72h",
    TRUE ~ as.character(group)
  ))

pathway_order <- df_ranknet_filtered %>%
  group_by(name) %>%
  summarise(tot = sum(contribution), .groups = "drop") %>%
  arrange(tot) %>%
  pull(name)

# pvalue check
pval_col <- intersect(colnames(df_ranknet_filtered),
                      c("pvalues", "pvalue", "p.value", "pval"))[1]
df_ranknet_filtered %>%
  mutate(name   = factor(name, levels = pathway_order),
         is_sig = .data[[pval_col]] < 0.05) %>%
  ggplot(aes(x = name, y = contribution, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = ifelse(is_sig, "*", "")),
            position = position_dodge(width = 0.7),
            vjust = -0.3, size = 4, color = "black", fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("24h" = "steelblue", "72h" = "darkred")) +
  theme_minimal(base_size = 11) +
  theme(legend.position    = "top",
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_text(face = "bold", size = 8)) +
  labs(title    = "Information Flow: Signaling Pathways 24h vs 72h",
       subtitle = "* p < 0.05 | pathway presenti in entrambi i timepoint",
       x = "", y = "Information Flow", fill = "Timepoint")

###################################
### SINGL L-R PAIRS ACCORDING TO DB 
###################################

plot_interactions <- function(df, src, tgt, db, title_label, n = 10) {
  plot_df <- df %>%
    filter(source == src, target == tgt, annotation == db) %>%
    group_by(timepoint) %>%
    slice_max(prob, n = n) %>%
    ungroup()
  
  if (nrow(plot_df) == 0) {
    message("Nessuna interazione trovata per ", src, "→", tgt, " (", db, ")")
    return(NULL)
  }
  
  ggplot(plot_df, aes(x = reorder(interaction_name, prob),
                      y = prob, fill = timepoint)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~ timepoint, scales = "free_y") +
    scale_fill_manual(values = c("24h" = "steelblue", "72h" = "darkred")) +
    theme_minimal(base_size = 11) +
    labs(title    = title_label,
         subtitle = paste("Top", n, db, "| pval < 0.05"),
         x = "", y = "Communication Probability") +
    theme(axis.text.y     = element_text(size = 9, face = "bold"),
          strip.text      = element_text(size = 11, face = "bold"),
          legend.position = "none")
}

# Secreted Signaling
print(plot_interactions(df_all, "Beta", "LY", "Secreted Signaling",
                        "Secreted Signals: β-cells → Lymphocytes"))
print(plot_interactions(df_all, "LY", "Beta", "Secreted Signaling",
                        "Secreted Signals: Lymphocytes → β-cells"))

# ECM-Receptor
print(plot_interactions(df_all, "Beta", "Lymphocyte", "ECM-Receptor",
                        "ECM Signals: β-cells → Lymphocytes"))
print(plot_interactions(df_all, "Lymphocyte", "Beta", "ECM-Receptor",
                        "ECM Signals: Lymphocytes → β-cells"))



#### PATHWAY ANALYSIS####
########################################################
# PLOT 3: Pathway timepoint-specifici — 2 plot separated
########################################################

extract_pathway_strength <- function(cellchat_obj, timepoint_label) {
  prob  <- cellchat_obj@netP$prob
  ec_ly <- prob["Beta", "Lymphocyte", ]
  ly_ec <- prob["Lymphocyte", "Beta", ]
  bind_rows(
    data.frame(pathway   = names(ec_ly), strength = ec_ly,
               direction = "Beta→Lymphocyte",     timepoint = timepoint_label),
    data.frame(pathway   = names(ly_ec), strength = ly_ec,
               direction = "Lymphocyte→Beta",     timepoint = timepoint_label)
  ) %>% filter(strength > 0)
}

df_strength <- bind_rows(
  extract_pathway_strength(cellchat.24h, "24h"),
  extract_pathway_strength(cellchat.72h, "72h")
)

########################################
# PLOT 3A — Pathway detected only at 24h
########################################

df_24h_specific <- df_strength %>%
  filter(timepoint == "24h",
         pathway %in% pathway_solo_24h)

p3a <- ggplot(df_24h_specific,
              aes(x = reorder(pathway, strength), 
                  y = strength, 
                  fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  facet_wrap(~ direction) +
  scale_fill_manual(values = c("EC→LY" = "#E45756", 
                               "LY→EC" = "#4C78A8")) +
  theme_minimal(base_size = 11) +
  theme(legend.position  = "none",
        strip.text        = element_text(face = "bold", size = 11),
        axis.text.y       = element_text(face = "bold", size = 9),
        panel.grid.major.y = element_blank()) +
  labs(title    = "Pathways active exclusively at 24h",
       subtitle = "Split by communication direction",
       x = "", y = "Communication Probability")

####################################
# PLOT 3B — Pathway esclusivi a 72h
####################################
df_72h_specific <- df_strength %>%
  filter(timepoint == "72h",
         pathway %in% pathway_solo_72h)

p3b <- ggplot(df_72h_specific,
              aes(x = reorder(pathway, strength), 
                  y = strength, 
                  fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  facet_wrap(~ direction) +
  scale_fill_manual(values = c("EC→LY" = "#E45756", 
                               "LY→EC" = "#4C78A8")) +
  theme_minimal(base_size = 11) +
  theme(legend.position  = "none",
        strip.text        = element_text(face = "bold", size = 11),
        axis.text.y       = element_text(face = "bold", size = 9),
        panel.grid.major.y = element_blank()) +
  labs(title    = "Pathways active exclusively at 24h",
       subtitle = "Split by communication direction",
       x = "", y = "Communication Probability")



##################################################
# PLOT 3: 4 plot separati — timepoint x direzione
##################################################
# Helper function
plot_pathway <- function(df, tp, dir, fill_color, title_label) {
  
  df_sub <- df %>%
    filter(timepoint == tp, direction == dir)
  
  if (nrow(df_sub) == 0) {
    message("Nessun pathway per ", tp, " ", dir)
    return(NULL)
  }
  
  ggplot(df_sub,
         aes(x = reorder(pathway, strength), y = strength)) +
    geom_bar(stat = "identity", fill = fill_color, width = 0.6) +
    geom_text(aes(label = round(strength, 3)),
              hjust = -0.1, size = 3, color = "grey30") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y        = element_text(face = "bold", size = 10),
      panel.grid.major.y = element_blank(),
      plot.title         = element_text(face = "bold", size = 11)
    ) +
    labs(title = title_label,
         x = "", y = "Communication Probability")
}

# 4 plot
p_24h_EC_LY <- plot_pathway(df_24h_specific, "24h", "Beta→Lymphocyte",
                            "#E45756", "24h | EC → LY")
p_24h_LY_EC <- plot_pathway(df_24h_specific, "24h", "Lymphocyte→Beta",
                            "#4C78A8", "24h | LY → EC")
p_72h_EC_LY <- plot_pathway(df_72h_specific, "72h", "Beta→Lymphocyte",
                            "#E45756", "72h | EC → LY")
p_72h_LY_EC <- plot_pathway(df_72h_specific, "72h", "Lymphocyte→Beta",
                            "#4C78A8", "72h | LY → EC")

# Layout 2x2
  (p_72h_EC_LY | p_72h_LY_EC) +
  plot_annotation(
    title    = "Timepoint-specific signaling pathways",
    subtitle = "Pathways active exclusively at 72h",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey40")
    )
  )
