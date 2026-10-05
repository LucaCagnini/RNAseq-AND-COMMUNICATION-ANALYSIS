#######################
### NICHENET ANALYSIS #
#######################

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
library("patchwork")
library("readxl")
library("readr")
library("tidyr")


BiocManager::install(c("SingleCellExperiment", "Seurat"))
# install nichenetr da GitHub
install.packages(c("spatstat.explore", "spatstat.random"))
install.packages("remotes")
remotes::install_github("saeyslab/nichenetr")

library(nichenetr)

#########################
# 1. DOWNLOAD THE NETWORK
#########################

ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
lr_network           <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
weighted_networks    <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

#######################
# 2. PREPARE EPRESSION
#######################
# Nichenet just uses the expressed genes as background
#SEPARATED NORMALISATION 

dds_EC <- estimateSizeFactors(dds_EC)
dds_LY <- estimateSizeFactors(dds_LY)
counts_EC_norm <- counts(dds_EC, normalized = TRUE)
counts_LY_norm <- counts(dds_LY, normalized = TRUE)

#MAPPING SYMBOLS
map_symbols <- function(counts_matrix) {
  gene_map <- mapIds(org.Hs.eg.db,
                     keys      = rownames(counts_matrix),
                     column    = "SYMBOL",
                     keytype   = "ENSEMBL",
                     multiVals = "first")
  counts_sym <- counts_matrix[!is.na(gene_map), ]
  rownames(counts_sym) <- gene_map[!is.na(gene_map)]
  counts_sym[!duplicated(rownames(counts_sym)), ]
}

counts_EC_sym <- map_symbols(counts_EC_norm)
counts_LY_sym <- map_symbols(counts_LY_norm)

#SUBSET 72h COCULTURE

meta_final <- read.csv("results/CELLCHAT_RESULTS/meta_final.csv", header = TRUE)
rownames(meta_final) <- meta_final$sample_id

meta_72h <- meta_final %>%
  filter(timepoint == "72h", condition == "CoCulture")

idx_EC_72h <- meta_72h$sample_id[meta_72h$cell_type == "Beta"]
idx_LY_72h <- meta_72h$sample_id[meta_72h$cell_type == "Lymphocyte"]

expr_EC <- counts_EC_sym[, idx_EC_72h]
expr_LY <- counts_LY_sym[, idx_LY_72h]

cat("Campioni EC a 72h:", ncol(expr_EC), "\n")
cat("Campioni LY a 72h:", ncol(expr_LY), "\n")

#GENES EXPRESSED (soglia coerente con baseMean > 10)
expressed_genes_sender   <- rownames(expr_EC)[rowMeans(expr_EC) > 10]
expressed_genes_receiver <- rownames(expr_LY)[rowMeans(expr_LY) > 10]

cat("Genes expressed EC:", length(expressed_genes_sender), "\n")
cat("Genes expressed LY:", length(expressed_genes_receiver), "\n")

######################################
# GENESET OF INTEREST : LY AS RECIVERS
######################################
# I need to get ly/ec degs 
degs_LY <- readRDS("DDS_RESULTS.rds")
res_LY_72h <- results(degs_LY, contrast = c("design", "72h_EC", "72h_noEC"))
res_LY_72h_df <- as.data.frame(res_LY_72h) %>%
  rownames_to_column("ensembl") %>%
  mutate(symbol = mapIds(org.Hs.eg.db,
                         keys = ensembl,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")) %>%
  filter(!is.na(symbol), !duplicated(symbol))

# UP genes - activated upon coculture
geneset_oi_up <- res_LY_72h_df %>%
  filter(padj < 0.05, log2FoldChange > 1) %>%
  pull(symbol)

# DOWN genes - suppressed upon  coculture
geneset_oi_down <- res_LY_72h_df %>%
  filter(padj < 0.05, log2FoldChange < - 1) %>%
  pull(symbol)

# Background = all genes expressed in LY
background_genes <- expressed_genes_receiver

# Filtering for overlap with NicheNet
geneset_oi_up   <- geneset_oi_up[geneset_oi_up %in% rownames(ligand_target_matrix)]
geneset_oi_down <- geneset_oi_down[geneset_oi_down %in% rownames(ligand_target_matrix)]
background_genes <- background_genes[background_genes %in% rownames(ligand_target_matrix)]

cat("Geneset UP:", length(geneset_oi_up), "\n")
cat("Geneset DOWN:", length(geneset_oi_down), "\n")
cat("Background:", length(background_genes), "\n")
cat("Ratio background/geneset UP:", round(length(background_genes) / length(geneset_oi_up), 1), "\n")

##################
# POTENTIAL LIGAND
##################
ligands   <- lr_network %>% pull(from) %>% unique()
receptors <- lr_network %>% pull(to)   %>% unique()

expressed_receptors_LY <- intersect(receptors, expressed_genes_receiver)
expressed_ligands_EC  <- intersect(ligands,   expressed_genes_sender)

# Sender-focused: ligands EC with receptor in the  LY
potential_ligands_focused <- lr_network %>%
  filter(from %in% expressed_ligands_EC,
         to   %in% expressed_receptors_LY) %>%
  pull(from) %>% unique()

# Sender-agnostic: ligands EC with receptor in the  LY
potential_ligands_agnostic <- lr_network %>%
  filter(to %in% expressed_receptors_LY) %>%
  pull(from) %>% unique()

cat("Ligandi focused (EC→LY):", length(potential_ligands_focused), "\n")
cat("Ligandi agnostic:", length(potential_ligands_agnostic), "\n")

#######################
# NICHENET UP AND DOWN
#######################
# How I run nichenet: I repeated the same code using LY as reciver and EC as reciver
# I runned using geneset target up genes before and then down genes 

##########################
# LY: LIGAND ACTIVITY - UP
##########################
# Focused
ligand_activities_focused <- predict_ligand_activities(
  geneset = geneset_oi_up,
  background_expressed_genes = background_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands_focused
) %>% arrange(desc(aupr_corrected)) %>%
  mutate(approach = "focused")

# Agnostic
ligand_activities_agnostic <- predict_ligand_activities(
  geneset = geneset_oi_up,
  background_expressed_genes = background_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands_agnostic
) %>% arrange(desc(aupr_corrected)) %>%
  mutate(approach = "agnostic")

# Remove genes with pearson NA
ligand_activities_agnostic <- ligand_activities_agnostic %>%
  filter(!is.na(pearson))
top_ligands_focused <- ligand_activities_focused %>%
  filter(aupr_corrected > quantile(aupr_corrected, 0.90)) %>%   # taking top 10%
  pull(test_ligand)

top_ligands_agnostic <- ligand_activities_agnostic %>%
  filter(aupr_corrected > quantile(aupr_corrected, 0.90)) %>%
  pull(test_ligand)
# common ligands 
comuni <- intersect(top_ligands_focused, top_ligands_agnostic)
cat("Ligandi in comune:", length(comuni), "\n")
print(comuni)

############################
# LY: LIGAND ACTIVITY - DOWN
############################
# Focused
ligand_activities_focused_down <- predict_ligand_activities(
  geneset = geneset_oi_down,
  background_expressed_genes = background_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands_focused
) %>% arrange(desc(aupr_corrected)) %>%
  mutate(approach = "focused")

# Agnostic
ligand_activities_agnostic_down <- predict_ligand_activities(
  geneset = geneset_oi_down,
  background_expressed_genes = background_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands_agnostic
) %>% arrange(desc(aupr_corrected)) %>%
  mutate(approach = "agnostic")

# Remove genes with pearson NA
ligand_activities_agnostic_down <- ligand_activities_agnostic_down %>%
  filter(!is.na(pearson))

# Top ligands DOWN 
top_ligands_focused_down <- ligand_activities_focused_down %>%
  filter(aupr_corrected > quantile(aupr_corrected, 0.90)) %>%
  pull(test_ligand)

top_ligands_agnostic_down <- ligand_activities_agnostic_down %>%
  filter(aupr_corrected > quantile(aupr_corrected, 0.90)) %>%
  pull(test_ligand)

# common ligand DOWN 
comuni_down <- intersect(top_ligands_focused_down, top_ligands_agnostic_down)
cat("Ligandi DOWN in common:", length(comuni_down), "\n")
print(comuni_down)

#################################
# PLOT 1: Ligand Activity barplot
#################################
ligand_activities_focused %>%
  filter(test_ligand %in% comuni) %>%
  arrange(desc(aupr_corrected)) %>%
  ggplot(aes(x = reorder(test_ligand, aupr_corrected),
             y = aupr_corrected)) +
  geom_bar(stat = "identity", fill = "darkred") +
  geom_hline(yintercept = mean(ligand_activities_focused$aupr_corrected),
             linetype = "dashed", color = "grey50") +
  coord_flip() +
  theme_minimal() +
  labs(title = "NicheNet Ligand Activity: EC → LY (72h)",
       subtitle = "Ligands shared between focused and agnostic",
       x = "Ligand", y = "AUPR corrected",
       caption = "Dashed line = mean AUPR across all ligands")


#####################################
# PREPARE EXPRESSION (Target: EC)
#####################################
# EC are reciving 
# Samples 72h coculture
idx_EC_72h <- rownames(meta_72h)[meta_72h$cell_type == "Beta"]
idx_LY_72h <- rownames(meta_72h)[meta_72h$cell_type == "Lymphocyte"]

# Matrici separate
expr_EC_target <- counts_EC_sym[, idx_EC_72h, drop = FALSE]
expr_LY_sender <- counts_LY_sym[, idx_LY_72h, drop = FALSE]

# Expressed genes
expressed_genes_sender_LY <- rownames(expr_LY_sender)[rowMeans(expr_LY_sender) > 10]
expressed_genes_receiver_EC <- rownames(expr_EC_target)[rowMeans(expr_EC_target) > 10]

# need them later for pathway in t1d_association
write.csv(expressed_genes_receiver_EC,"NICHENET_RESULTS/background_EC")
write.csv(expressed_genes_receiver, "NICHENET_RESULTS/background_LY")

##################################
# GENESET OF INTEREST (DEGs in EC)
##################################
# Results DE of dds_EC per le 72h
res_EC <- DESeq(res_EC)dds_EC  <- readRDS(here::here("data", "DDS_EC"))
degs_EC <- readRDS(here::here("data", "DEGS_EC"))
res_EC_72h <- results(res_EC, contrast=c("design", "CD3_72h", "Endo"))
res_EC_72h_df <- as.data.frame(res_EC_72h) %>%
  rownames_to_column("ensembl") %>%
  mutate(symbol = mapIds(org.Hs.eg.db, keys = ensembl, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")) %>%
  filter(!is.na(symbol), !duplicated(symbol))

# Up genes in Beta Cells induced by lymphocytes
geneset_ec_up <- res_EC_72h_df %>%
  filter(padj < 0.05, log2FoldChange > 1.5) %>%
  pull(symbol)

background_genes_ec <- expressed_genes_receiver_EC

# Filtro per matrice NicheNet
geneset_ec_up <- intersect(geneset_ec_up, rownames(ligand_target_matrix))
background_genes_ec <- intersect(background_genes_ec, rownames(ligand_target_matrix))

cat("Geni target (EC UP):", length(geneset_ec_up), "\n")

###############################
# POTENTIAL LIGANDS (LY -> EC)
###############################
expressed_receptors_EC <- intersect(receptors, expressed_genes_receiver_EC)
expressed_ligands_LY   <- intersect(ligands,   expressed_genes_sender_LY)

# ligand produce by LY that have a receptor on EC
potential_ligands_LY_to_EC <- lr_network %>%
  filter(from %in% expressed_ligands_LY,
         to   %in% expressed_receptors_EC) %>%
  pull(from) %>% unique()

cat("Ligandi potenziali LY -> EC:", length(potential_ligands_LY_to_EC), "\n")

##########################
# EC: LIGAND ACTIVITY - UP
##########################
ligand_activities_LY_focused <- predict_ligand_activities(
  geneset = geneset_ec_up,
  background_expressed_genes = background_genes_ec,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands_LY_to_EC
) %>%
  arrange(desc(aupr_corrected)) %>%
  mutate(approach = "focused")

ligand_activities_LY_agnostic <- predict_ligand_activities(
  geneset = geneset_ec_up,
  background_expressed_genes = background_genes_ec,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = unique(lr_network$from)
) %>%
  arrange(desc(aupr_corrected)) %>%
  mutate(approach = "agnostic")

## top focused and agnostic

top_focused <- ligand_activities_LY_focused %>%
  filter(aupr_corrected > quantile(aupr_corrected, 0.90),
         #pearson > 0.10
         ) %>%   
  pull(test_ligand)

top_agnostic <- ligand_activities_LY_agnostic %>%
  filter(aupr_corrected > quantile(aupr_corrected, 0.90),
           #pearson > 0.10
         ) %>%
  pull(test_ligand)

ligands_robust <- intersect(top_focused, top_agnostic)
print(ligands_robust)


###########################
#EC: LIGAND ACTIVITY DOWN 
###########################

geneset_ec_down <- res_EC_72h_df %>%
  filter(padj < 0.05, log2FoldChange < -1.5) %>%
  pull(symbol) %>%
  intersect(rownames(ligand_target_matrix))

cat("Geneset EC DOWN:", length(geneset_ec_down), "\n")

ligand_activities_LY_focused_down <- predict_ligand_activities(
  geneset = geneset_ec_down,
  background_expressed_genes = background_genes_ec,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands_LY_to_EC
) %>%
  arrange(desc(aupr_corrected)) %>%
  mutate(approach = "focused")

ligand_activities_LY_agnostic_down <- predict_ligand_activities(
  geneset = geneset_ec_down,
  background_expressed_genes = background_genes_ec,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = unique(lr_network$from)
) %>%
  arrange(desc(aupr_corrected)) %>%
  mutate(approach = "agnostic") %>%
  filter(!is.na(pearson))

# Top ligandi DOWN
top_focused_down <- ligand_activities_LY_focused_down %>%
  filter(aupr_corrected > quantile(aupr_corrected, 0.90)
         #pearson > 0.10
         ) %>%
  pull(test_ligand)

top_agnostic_down <- ligand_activities_LY_agnostic_down %>%
  filter(aupr_corrected > quantile(aupr_corrected, 0.90)
         #pearson > 0.10
         ) %>%
  pull(test_ligand)

ligands_robust_down <- intersect(top_focused_down, top_agnostic_down)
cat("Ligandi LY->EC DOWN robusti:", length(ligands_robust_down), "\n")
print(ligands_robust_down)

#################################
# PLOT 2: Ligand Activity barplot
#################################

plot_common <- ligand_activities_LY_focused %>%
  filter(test_ligand %in% ligands_robust) %>%
  arrange(desc(aupr_corrected))

ggplot(plot_common,
       aes(x = reorder(test_ligand, aupr_corrected),
           y = aupr_corrected)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "NicheNet Ligand Activity: LY → EC (72h)",
       subtitle = "Ligands shared between focused and agnostic",
       x = "Ligand",
       y = "AUPR corrected")


################################################
# DATA INTEGRATION: NICHENET, CELLCHAT, LUMINEX
################################################

###########################################
###  LIGAND RECEPTORS EC -> LY ####
# STEP 1 — CellChat layer (72h, pval < 0.05)
############################################
df_all <- read.csv("results/CELLCHAT_RESULTS/LY_EC_CELLCHAT.csv")
cc_edges <- df_all %>%
  filter(timepoint == "72h",
         source == "Beta",
         target == "Lymphocyte",
         pval < 0.05,
         annotation %in% c(
           "Secreted Signaling", 
                           "ECM-Receptor"
                           )) %>%
  dplyr::select(ligand, receptor, annotation, prob) %>%
  distinct()
cat("EC → LY interactions (Secreted + ECM ):", nrow(cc_edges), "\n")
##############################################
# STEP 2 — NicheNet layer (UP + DOWN separati)
##############################################
# UP
nn_ligands_up <- intersect(top_ligands_focused, top_ligands_agnostic)
nn_ligands_up <- intersect(nn_ligands_up, rownames(ligand_target_matrix))

# DOWN
nn_ligands_down <- intersect(top_ligands_focused_down, top_ligands_agnostic_down)
nn_ligands_down <- intersect(nn_ligands_down, rownames(ligand_target_matrix))

# all ligands in NicheNet with directions
nn_all <- bind_rows(
  data.frame(ligand = nn_ligands_up,   nn_direction = "UP"),
  data.frame(ligand = nn_ligands_down, nn_direction = "DOWN")
) %>%
  # if a ligand is found in both, keep that info
  group_by(ligand) %>%
  summarise(nn_direction = paste(sort(unique(nn_direction)), collapse = " | "),
            .groups = "drop")

nn_ligands_all <- nn_all$ligand

cat("Ligandi NicheNet UP:   ", length(nn_ligands_up), "\n")
cat("Ligandi NicheNet DOWN: ", length(nn_ligands_down), "\n")
cat("Ligandi NicheNet total:", length(nn_ligands_all), "\n")

#########################
# STEP 3 — Luminex layer
#########################

get_expressed_cytokines_72h <- function(lum_data, cy_list, min_prot = 10) {
  lum_72h <- lum_data %>%
    filter(Condition == "72h", Mean_Value > min_prot)
  
  # Expand multiple genes sep by , 
  cy_mapping <- cy_list %>%
    dplyr::select(
      Gene         = `CYT_GENE_SYMBOL_(Ligand)`,
      Luminex_Name = CYT_ALIAS
    ) %>%
    separate_rows(Gene, sep = ",\\s*") %>%   # split su virgola + spazi opzionali
    distinct()
  
  result <- lum_72h %>%
    left_join(cy_mapping, by = c("Cytokine" = "Luminex_Name")) %>%
    dplyr::select(
      Luminex_Name = Cytokine,
      Gene,
      Condition,
      Mean_Value,
      SD_Value
    ) %>%
    arrange(desc(Mean_Value))
  
  return(result)
}
# Call that
luminex_signals <- get_expressed_cytokines_72h(lum_data  = df_luminex,
                                             cy_list   = Cy_list,
                                             min_prot  = 10)

#########################
# STEP 4 — NORMALISATION
#########################
cc_ligands  <- unique(cc_edges$ligand)
lum_ligands <- unique(luminex_signals$Gene)
all_ligands <- unique(c(cc_ligands, nn_ligands_all, lum_ligands)) 

evidence_table <- data.frame(ligand = all_ligands) %>%
  mutate(
    CellChat = ligand %in% cc_ligands,
    NicheNet = ligand %in% nn_ligands_all,
    Luminex  = ligand %in% lum_ligands,
    evidence_count = as.integer(CellChat) +
      as.integer(NicheNet) +
      as.integer(Luminex)
  ) %>%
  # Add NicheNet didection 
  left_join(nn_all, by = "ligand") %>%
  mutate(nn_direction = replace_na(nn_direction, "Not in NicheNet"))

###########################
# STEP 5 — Evidence matrix
##########################

# AUPR from UP and from DOWN regulated genes
aupr_up <- ligand_activities_focused %>%
  dplyr::select(test_ligand, aupr_corrected) %>%
  rename(aupr_up = aupr_corrected)

aupr_down <- ligand_activities_focused_down %>%
  dplyr::select(test_ligand, aupr_corrected) %>%
  rename(aupr_down = aupr_corrected)

evidence_final <- evidence_table %>%
  filter(evidence_count > 0) %>%
  left_join(aupr_up,  by = c("ligand" = "test_ligand")) %>%
  left_join(aupr_down, by = c("ligand" = "test_ligand")) %>%
  mutate(
    aupr_corrected = case_when(
      nn_direction == "UP"       ~ aupr_up,
      nn_direction == "DOWN"     ~ aupr_down,
      nn_direction == "UP | DOWN" ~ pmax(aupr_up, aupr_down, na.rm = TRUE),
      TRUE ~ 0
    )
  ) %>%
  left_join(
    df_all %>%
      filter(timepoint == "72h", source == "EC", target == "LY",
             annotation %in% c("Secreted Signaling", "ECM-Receptor"),
             pval < 0.05) %>%
      group_by(ligand) %>%
      dplyr::summarise(
        max_prob      = max(prob),
        cc_annotation = paste(sort(unique(annotation)), collapse = " | "),
        .groups = "drop"
      ),
    by = "ligand"
  ) %>%
  mutate(
    aupr_corrected = replace_na(aupr_corrected, 0),
    aupr_up        = replace_na(aupr_up, 0),
    aupr_down      = replace_na(aupr_down, 0),
    max_prob       = replace_na(max_prob, 0),
    cc_annotation  = replace_na(cc_annotation, "Not in CellChat"),
    ligand         = reorder(ligand, evidence_count),
    classification = case_when(
      evidence_count == 3 ~ "High confidence",
      evidence_count == 2 ~ "Medium confidence",
      evidence_count == 1 ~ "Low confidence",
      TRUE ~ "None"
    )
  ) %>%
  # Order in a readable way
  dplyr::select(ligand, CellChat, NicheNet, Luminex, 
                evidence_count, classification,
                nn_direction,           # UP / DOWN / UP | DOWN / Not in NicheNet
                aupr_corrected,         # AUPR rilevante
                aupr_up, aupr_down,     # ntrambi per trasparenza
                max_prob, cc_annotation)

#########################
# STEP 6 — Evidence scoring
#########################
cat("\nDistribuzione classificazione EC→LY:\n")
print(table(evidence_final$classification))

cat("\n=== HIGH CONFIDENCE ===\n")
print(evidence_final %>% filter(classification == "High confidence"))

cat("\n=== MEDIUM CONFIDENCE ===\n")
print(evidence_final %>% filter(classification == "Medium confidence"))

cat("\n=== LOW CONFIDENCE ===\n")
print(evidence_final %>% filter(classification == "Low confidence"))


###################################
###  LIGAND RECEPTORS LY -> EC ####

# STEP 1 — CellChat
####################

cc_df_LY <- df_all %>%
  filter(timepoint == "72h",
         source == "Lymphocyte",
         target == "Beta",
         pval < 0.05,
         annotation %in% c("Secreted Signaling", "ECM-Receptor")) %>%
  distinct(ligand, receptor, .keep_all = TRUE)
cc_ligands_LY <- unique(cc_df_LY$ligand)
cat("CellChat LY→EC ligands:", length(cc_ligands_LY), "\n")

#########################
# STEP 2 — NicheNet UP + DOWN
#########################

# UP
top_focused_up <- ligand_activities_LY_focused %>%
  filter(aupr_corrected > quantile(aupr_corrected, 0.90)) %>%
  pull(test_ligand)
top_agnostic_up <- ligand_activities_LY_agnostic %>%
  filter(aupr_corrected > quantile(aupr_corrected, 0.90)) %>%
  pull(test_ligand)
nn_ligands_LY_up <- intersect(top_focused_up, top_agnostic_up)

# DOWN
top_focused_down <- ligand_activities_LY_focused_down %>%
  filter(aupr_corrected > quantile(aupr_corrected, 0.90)) %>%
  pull(test_ligand)
top_agnostic_down <- ligand_activities_LY_agnostic_down %>%
  filter(aupr_corrected > quantile(aupr_corrected, 0.90)) %>%
  pull(test_ligand)
nn_ligands_LY_down <- intersect(top_focused_down, top_agnostic_down)

# Table direction NicheNet
nn_all_LY <- bind_rows(
  data.frame(ligand = nn_ligands_LY_up,   nn_direction = "UP"),
  data.frame(ligand = nn_ligands_LY_down, nn_direction = "DOWN")
) %>%
  group_by(ligand) %>%
  summarise(nn_direction = paste(sort(unique(nn_direction)), collapse = " | "),
            .groups = "drop")

nn_ligands_LY_all <- nn_all_LY$ligand

cat("NicheNet LY→EC UP:    ", length(nn_ligands_LY_up), "\n")
cat("NicheNet LY→EC DOWN:  ", length(nn_ligands_LY_down), "\n")
cat("NicheNet LY→EC total:", length(nn_ligands_LY_all), "\n")

#########################
# STEP 3 — Luminex
#########################
# same as before
luminex_signals_LY <- get_expressed_cytokines_72h(lum_data  = df_luminex,
                                               cy_list   = Cy_list,
                                               min_prot  = 10)
#########################
# STEP 4 — Universo ligandi
#########################
cc_ligands_LY     <- unique(as.character(cc_df_LY$ligand))  
lum_ligands_LY    <- unique(as.character(luminex_signals_LY$Gene))
nn_ligands_LY_all <- unique(as.character(nn_all_LY$ligand))

all_ligands_LY <- Reduce(union, list(
  cc_ligands_LY,
  nn_ligands_LY_all,
  lum_ligands_LY
)) %>%
  intersect(rownames(ligand_target_matrix))

cat("Final ligand universe LY→EC:", length(all_ligands_LY), "\n")

#########################
# STEP 5 — AUPR UP e DOWN
#########################
aupr_up_LY <- ligand_activities_LY_focused %>%
  dplyr::select(test_ligand, aupr_corrected) %>%
  rename(aupr_up = aupr_corrected)

aupr_down_LY <- ligand_activities_LY_focused_down %>%
  dplyr::select(test_ligand, aupr_corrected) %>%
  rename(aupr_down = aupr_corrected)

#########################
# STEP 6 — Evidence table
#########################

evidence_table_LY <- tibble(ligand = all_ligands_LY) %>%
  mutate(
    CellChat = ligand %in% cc_ligands_LY,
    NicheNet = ligand %in% nn_ligands_LY_all,
    Luminex  = ligand %in% lum_ligands_LY,
    evidence_count = rowSums(across(c(CellChat, NicheNet, Luminex))),
    classification = case_when(
      evidence_count == 3 ~ "High confidence",
      evidence_count == 2 ~ "Medium confidence",
      evidence_count == 1 ~ "Low confidence",
      TRUE ~ "None"
    )
  ) %>%
  left_join(nn_all_LY, by = "ligand") %>%
  mutate(nn_direction = replace_na(nn_direction, "Not in NicheNet")) %>%
  arrange(desc(evidence_count), ligand)

#########################
# STEP 7 — Evidence final
#########################
evidence_final_LY <- evidence_table_LY %>%
  filter(evidence_count > 0) %>%
  left_join(aupr_up_LY,   by = c("ligand" = "test_ligand")) %>%
  left_join(aupr_down_LY, by = c("ligand" = "test_ligand")) %>%
  mutate(
    aupr_corrected = case_when(
      nn_direction == "UP"        ~ aupr_up,
      nn_direction == "DOWN"      ~ aupr_down,
      nn_direction == "UP | DOWN" ~ pmax(aupr_up, aupr_down, na.rm = TRUE),
      TRUE ~ 0
    ),
    aupr_up   = replace_na(aupr_up, 0),
    aupr_down = replace_na(aupr_down, 0)
  ) %>%
  left_join(
    df_all %>%
      filter(timepoint == "72h", source == "LY", target == "EC",
             annotation %in% c("Secreted Signaling", "ECM-Receptor"),
             pval < 0.05) %>%
      group_by(ligand) %>%
      dplyr::summarise(
        max_prob      = max(prob),
        cc_annotation = paste(sort(unique(annotation)), collapse = " | "),
        .groups = "drop"
      ),
    by = "ligand"
  ) %>%
  mutate(
    aupr_corrected = replace_na(aupr_corrected, 0),
    max_prob       = replace_na(max_prob, 0),
    cc_annotation  = replace_na(cc_annotation, "Not in CellChat"),
    ligand         = reorder(ligand, evidence_count)
  ) %>%
  dplyr::select(
    ligand, CellChat, NicheNet, Luminex,
    evidence_count, classification,
    nn_direction,
    aupr_corrected, aupr_up, aupr_down,
    max_prob, cc_annotation
  ) %>%
  arrange(desc(evidence_count), desc(aupr_corrected))

#########################
# STEP 8 — OUTPUT
#########################
cat("\nDistribuzione classificazione LY→EC:\n")
print(table(evidence_final_LY$classification))

cat("\n=== HIGH CONFIDENCE ===\n")
print(evidence_final_LY %>% filter(classification == "High confidence"))

cat("\n=== MEDIUM CONFIDENCE ===\n")
print(evidence_final_LY %>% filter(classification == "Medium confidence"))

cat("\n=== LOW CONFIDENCE ===\n")
print(evidence_final_LY %>% filter(classification == "Low confidence"))

##################################################
# filtering out the one who are in the two df only because the portein is found in the media
evidence_final_LY <- evidence_final_LY %>%
  filter(!(Luminex == TRUE & CellChat == FALSE & NicheNet == FALSE))

evidence_final <- evidence_final %>%
  filter(!(Luminex == TRUE & CellChat == FALSE & NicheNet == FALSE))

## save the two dataframe 
write_csv(evidence_final, "results/NICHENET_RESULTS/evidence_table_EC_communication.csv")
write_csv(evidence_final_LY, "results/NICHENET_RESULTS/evidence_table_LY_communication.csv")

evidence_final<- read.csv("results/NICHENET_RESULTS/evidence_table_EC_communication.csv", header=TRUE)
evidence_final_LY<- read.csv("results/NICHENET_RESULTS/evidence_table_LY_communication.csv", header = TRUE)

#####################################
# LIGAND-RECEPTOR ANALYSIS — EC ↔ LY
# INTEGRATING RECEPTORS FOR NICHENET AND INFORMATION ABOUT DEGS  
#####################################
# STEP 0 — DEG symbols 
####################################

degs_EC_up   <- res_EC_72h_df %>% filter(padj < 0.05, log2FoldChange >  1.5) %>% pull(symbol)
degs_EC_down <- res_EC_72h_df %>% filter(padj < 0.05, log2FoldChange < -1.5) %>% pull(symbol)
degs_EC_any  <- c(degs_EC_up, degs_EC_down)

degs_LY_up   <- res_LY_72h_df %>% filter(padj < 0.05, log2FoldChange >  1) %>% pull(symbol)
degs_LY_down <- res_LY_72h_df %>% filter(padj < 0.05, log2FoldChange <  -1) %>% pull(symbol)
degs_LY_any  <- c(degs_LY_up, degs_LY_down)

cat("DEG EC:", length(degs_EC_any), "(UP:", length(degs_EC_up), "DOWN:", length(degs_EC_down), ")\n")
cat("DEG LY:", length(degs_LY_any), "(UP:", length(degs_LY_up), "DOWN:", length(degs_LY_down), ")\n")

##################################################################
# STEP1 - raw lists: ligand and receptors expressed for cell-type 
#         (baseMean > 10, in evidence final/network)
##################################################################

all_receptors <- unique(lr_network$to)
# Ligands: the ones I found expressed (NicheNet + CellChat)
ligands_expressed_EC <- res_EC_72h_df %>%
  filter(baseMean > 10, symbol %in% as.character(evidence_final$ligand)) %>%
  pull(symbol) %>% unique()

ligands_expressed_LY <- res_LY_72h_df %>%
  filter(baseMean > 10, symbol %in% as.character(evidence_final_LY$ligand)) %>%
  pull(symbol) %>% unique()

# Receptors: all the ones on the LR pairs
receptors_expressed_EC <- res_EC_72h_df %>%
  filter(baseMean > 10, symbol %in% all_receptors) %>%
  pull(symbol) %>% unique()

receptors_expressed_LY <- res_LY_72h_df %>%
  filter(baseMean > 10, symbol %in% all_receptors) %>%
  pull(symbol) %>% unique()

cat("\n--- Geni espressi (baseMean > 10) ---\n")
cat("Ligandi espressi EC:", length(ligands_expressed_EC), "\n")
cat("Ligandi espressi LY:", length(ligands_expressed_LY), "\n")
cat("Recettori espressi EC:", length(receptors_expressed_EC), "\n")
cat("Recettori espressi LY:", length(receptors_expressed_LY), "\n")

#########################################################################
# STEP 2 — Confidence lists
# I'm only taking high/medium confidence ligands
########################################################################

high_conf_EC   <- evidence_final    %>% filter(classification == "High confidence")   %>% pull(ligand) %>% as.character()
medium_conf_EC <- evidence_final    %>% filter(classification == "Medium confidence") %>% pull(ligand) %>% as.character()
high_conf_LY   <- evidence_final_LY %>% filter(classification == "High confidence")   %>% pull(ligand) %>% as.character()
medium_conf_LY <- evidence_final_LY %>% filter(classification == "Medium confidence") %>% pull(ligand) %>% as.character()

#########################################################################
# STEP 3 — Ligand sets final
#########################################################################

ligands_EC_final <- ligands_expressed_EC
ligands_LY_final <- ligands_expressed_LY

cat("\n--- Ligand sets finali ---\n")
cat("Ligandi EC finali:", length(ligands_EC_final), "\n")
cat("Ligandi LY finali:", length(ligands_LY_final), "\n")

#########################################################################
# STEP 4 — Explosions complexes receptors CellChat + conversione HGNC
#########################################################################

cellchat_to_hgnc <- c(
  "TGFbR1" = "TGFBR1",
  "TGFbR"  = "TGFBR1",
  "R2"     = "TGFBR2"
)

expand_receptor_complex <- function(receptor_string) {
  components <- strsplit(receptor_string, "_")[[1]]
  components <- dplyr::recode(components, !!!cellchat_to_hgnc)
  return(components)
}

# Cellchat couples with HGNC
cc_pairs_expanded <- df_all %>%
  filter(timepoint == "72h", pval < 0.05,
         annotation %in% c("Secreted Signaling", "ECM-Receptor")) %>%
  dplyr::select(source, target, ligand, receptor) %>%
  distinct() %>%
  rowwise() %>%
  mutate(receptor_hgnc = list(expand_receptor_complex(receptor))) %>%
  unnest(receptor_hgnc) %>%
  ungroup()

#########################################################################
# STEP 5 — Helper: inclusion_reason
#########################################################################

get_inclusion_reason <- function(symbol_val, degs_strict, high_conf, medium_conf, res_df) {
  is_deg <- symbol_val %in% degs_strict
  is_hc  <- symbol_val %in% high_conf
  is_mc  <- symbol_val %in% medium_conf
  bm     <- res_df %>% filter(symbol == symbol_val) %>% pull(baseMean)
  expressed <- length(bm) > 0 && bm[1] > 10
  case_when(
    is_deg & is_hc ~ "DE + High confidence",
    is_deg & is_mc ~ "DE + Medium confidence",
    is_deg         ~ "DE (padj < 0.05)",
    is_hc          ~ "High confidence + expressed",
    is_mc          ~ "Medium confidence + expressed",
    TRUE           ~ "expressed only"
  )
}

#########################################################################
# Helper: stato DEG
#########################################################################

deg_status <- function(res_df, symbols) {
  res_df %>%
    filter(symbol %in% symbols) %>%
    dplyr::select(symbol, log2FoldChange, padj, baseMean) %>%
    mutate(
      direction = case_when(
        log2FoldChange > 0 & padj < 0.05 ~ "UP",
        log2FoldChange < 0 & padj < 0.05 ~ "DOWN",
        TRUE ~ "not significant"
      )
    )
}

#########################################################################
# PARTE 1 — LY → EC
#########################################################################
# CellChat couples LY→EC with receptors expressed on EC
cc_pairs_LY_EC <- cc_pairs_expanded %>%
  filter(source == "Lymphocyte", target == "Beta",
         ligand       %in% ligands_LY_final,
         receptor_hgnc %in% receptors_expressed_EC) %>%
  dplyr::select(from = ligand, to = receptor_hgnc) %>%
  distinct() %>%
  mutate(source_db = "CellChat")

# NicheNet: ligands from Nichenet with receptors expressed on EC 
nn_ligands_LY <- evidence_final_LY %>%
  filter(NicheNet == TRUE) %>%
  pull(ligand) %>% as.character()

nn_pairs_LY_EC <- lr_network %>%
  filter(from %in% nn_ligands_LY,
         from %in% ligands_LY_final,
         to   %in% receptors_expressed_EC) %>%
  dplyr::select(from, to) %>%
  distinct() %>%
  mutate(source_db = "NicheNet")

# Rete ibrida LY→EC
lr_subset_LY <- bind_rows(cc_pairs_LY_EC, nn_pairs_LY_EC) %>%
  distinct(from, to, .keep_all = TRUE)

cat("\n--- LY → EC ---\n")
cat("Ligandi dopo filtro:", length(unique(lr_subset_LY$from)), "\n")
cat("  da CellChat:", length(unique(cc_pairs_LY_EC$from)), "\n")
cat("  da NicheNet:", length(unique(nn_pairs_LY_EC$from)), "\n")
cat("Recettori espressi nelle EC:", length(unique(lr_subset_LY$to)), "\n")

deg_ligands_LY <- deg_status(res_LY_72h_df, unique(lr_subset_LY$from)) %>%
  rename(ligand_log2FC    = log2FoldChange,
         ligand_padj      = padj,
         ligand_baseMean  = baseMean,
         ligand_direction = direction) %>%
  mutate(inclusion_reason = map_chr(symbol, ~get_inclusion_reason(
    .x, degs_LY_any, high_conf_LY, medium_conf_LY, res_LY_72h_df)))

deg_receptors_EC <- deg_status(res_EC_72h_df, unique(lr_subset_LY$to)) %>%
  rename(receptor_log2FC    = log2FoldChange,
         receptor_padj      = padj,
         receptor_baseMean  = baseMean,
         receptor_direction = direction)

ligand_receptor_LY_EC <- lr_subset_LY %>%
  left_join(deg_ligands_LY,   by = c("from" = "symbol")) %>%
  left_join(deg_receptors_EC, by = c("to"   = "symbol")) %>%
  left_join(
    evidence_final_LY %>%
      mutate(ligand = as.character(ligand)) %>%
      dplyr::select(ligand, classification, evidence_count,
                    CellChat, NicheNet, Luminex, nn_direction,
                    aupr_corrected, max_prob),
    by = c("from" = "ligand")
  ) %>%
  rename(ligand = from, receptor = to) %>%
  arrange(desc(evidence_count), desc(ligand_log2FC))

cat("Interazioni LY->EC totali:", nrow(ligand_receptor_LY_EC), "\n")
cat("  Ligandi High confidence:",   ligand_receptor_LY_EC %>% filter(classification == "High confidence")   %>% pull(ligand) %>% unique() %>% length(), "\n")
cat("  Ligandi Medium confidence:", ligand_receptor_LY_EC %>% filter(classification == "Medium confidence") %>% pull(ligand) %>% unique() %>% length(), "\n")
cat("  Ligandi DE:",                ligand_receptor_LY_EC %>% filter(ligand_direction != "not significant") %>% pull(ligand) %>% unique() %>% length(), "\n")

#########################################################################
# PARTE 2 — EC → LY
#########################################################################
# CellChat couples EC→LY with receptors expressed on LY

cc_pairs_EC_LY <- cc_pairs_expanded %>%
  filter(source == "Beta", target == "Lymphocyte",
         ligand        %in% ligands_EC_final,
         receptor_hgnc %in% receptors_expressed_LY) %>%
  dplyr::select(from = ligand, to = receptor_hgnc) %>%
  distinct() %>%
  mutate(source_db = "CellChat")

# NicheNet: ligands from Nichenet with receptors expressed on LY
nn_ligands_EC <- evidence_final %>%
  filter(NicheNet == TRUE) %>%
  pull(ligand) %>% as.character()

nn_pairs_EC_LY <- lr_network %>%
  filter(from %in% nn_ligands_EC,
         from %in% ligands_EC_final,
         to   %in% receptors_expressed_LY) %>%
  dplyr::select(from, to) %>%
  distinct() %>%
  mutate(source_db = "NicheNet")

# Rete ibrida EC→LY
lr_subset_EC <- bind_rows(cc_pairs_EC_LY, nn_pairs_EC_LY) %>%
  distinct(from, to, .keep_all = TRUE)

cat("\n--- EC → LY ---\n")
cat("Ligandi dopo filtro:", length(unique(lr_subset_EC$from)), "\n")
cat("  da CellChat:", length(unique(cc_pairs_EC_LY$from)), "\n")
cat("  da NicheNet:", length(unique(nn_pairs_EC_LY$from)), "\n")
cat("Recettori espressi nei LY:", length(unique(lr_subset_EC$to)), "\n")

deg_ligands_EC <- deg_status(res_EC_72h_df, unique(lr_subset_EC$from)) %>%
  rename(ligand_log2FC    = log2FoldChange,
         ligand_padj      = padj,
         ligand_baseMean  = baseMean,
         ligand_direction = direction) %>%
  mutate(inclusion_reason = map_chr(symbol, ~get_inclusion_reason(
    .x, degs_EC_any, high_conf_EC, medium_conf_EC, res_EC_72h_df)))

deg_receptors_LY <- deg_status(res_LY_72h_df, unique(lr_subset_EC$to)) %>%
  rename(receptor_log2FC    = log2FoldChange,
         receptor_padj      = padj,
         receptor_baseMean  = baseMean,
         receptor_direction = direction)

ligand_receptor_EC_LY <- lr_subset_EC %>%
  left_join(deg_ligands_EC,   by = c("from" = "symbol")) %>%
  left_join(deg_receptors_LY, by = c("to"   = "symbol")) %>%
  left_join(
    evidence_final %>%
      mutate(ligand = as.character(ligand)) %>%
      dplyr::select(ligand, classification, evidence_count,
                    CellChat, NicheNet, Luminex, nn_direction,
                    aupr_corrected, max_prob),
    by = c("from" = "ligand")
  ) %>%
  rename(ligand = from, receptor = to) %>%
  arrange(desc(evidence_count), desc(ligand_log2FC))

cat("Interactions EC->LY total:", nrow(ligand_receptor_EC_LY), "\n")
cat("  Ligandi DE strict:",         ligand_receptor_EC_LY %>% filter(inclusion_reason %in% c("DE (padj < 0.05)", "DE + High confidence", "DE + Medium confidence")) %>% pull(ligand) %>% unique() %>% length(), "\n")
cat("  Ligandi High confidence:",   ligand_receptor_EC_LY %>% filter(classification == "High confidence")  %>% pull(ligand) %>% unique() %>% length(), "\n")
cat("  Ligandi expressed only:",    ligand_receptor_EC_LY %>% filter(inclusion_reason == "expressed only") %>% pull(ligand) %>% unique() %>% length(), "\n")

#########################################################################
# EXPORT
#########################################################################

write.csv(ligand_receptor_LY_EC, "results/NICHENET_RESULTS/LigandReceptor_LY_to_EC.csv", row.names = FALSE)
write.csv(ligand_receptor_EC_LY, "results/NICHENET_RESULTS/LigandReceptor_EC_to_LY.csv", row.names = FALSE)


