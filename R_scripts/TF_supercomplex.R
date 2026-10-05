#########################################################################
# PIPELINE INTEGRATED: PROGENy + VIPER/DoRothEA
# LR pairs → pathway activity → TF activity → DEG targets
#########################################################################

library(progeny)
library(dorothea)
library(viper)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(DESeq2)
library(org.Hs.eg.db)


#helper functions DEG

Up_genes <- function(res, padj_cutoff = 0.05, lfc_cutoff = 0) {
  rownames(res)[which(res$padj < padj_cutoff &
                        res$log2FoldChange > lfc_cutoff)]
}

Dw_genes <- function(res, padj_cutoff = 0.05, lfc_cutoff = 0) {
  rownames(res)[which(res$padj < padj_cutoff &
                        res$log2FoldChange < -lfc_cutoff)]
}

load_dorothea_regulon <- function(confidence_levels = c("A", "B", "C")) {
  data(dorothea_hs, package = "dorothea", envir = environment())
  dorothea_hs %>%
    filter(confidence %in% confidence_levels) %>%
    df2regulon()
}

# 1.Data and normalisation

dds_LY <- readRDSdds_LY  <- readRDS(here::here("data", "FINAL_DDS.rds"))
degs_LY <- readRDS(here::here("data", "DDS_RESULTS.rds")
dds_LY <- DESeq(dds_LY)
vsd <- vst(dds_LY, blind = FALSE)
expr_mat <- assay(vsd)
write.csv(expr_mat, file = "cytosig_input.csv", row.names = TRUE)

# Ensembl → Symbol
counts_df <- as.data.frame(expr_mat) %>%
  rownames_to_column("ensembl") %>%
  mutate(symbol = mapIds(org.Hs.eg.db,
                         keys      = ensembl,
                         keytype   = "ENSEMBL",
                         column    = "SYMBOL",
                         multiVals = "first")) %>%
  filter(!is.na(symbol)) %>%
  group_by(symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

expr_mat <- counts_df %>%
  column_to_rownames("symbol") %>%
  as.matrix()

expr_mat <- expr_mat[apply(expr_mat, 1, var) > 0, ]
cat("Dimensions expression matrix:", dim(expr_mat), "\n")


### PROGENY 
run_progeny_paired <- function(expr_mat,
                               timepoint_pattern,
                               control_pattern = "noEC") {
  
  idx <- grep(timepoint_pattern, colnames(expr_mat))
  mat_tp <- expr_mat[, idx, drop = FALSE]
  
  pathway_scores <- t(
    progeny(
      mat_tp,
      scale = TRUE,
      organism = "Human",
      top = 500,
      perm = 1
    )
  )
  
  sample_info <- data.frame(
    sample = colnames(pathway_scores)
  ) %>%
    dplyr::mutate(
      condition = ifelse(grepl(control_pattern, sample), "noEC", "EC"),
      donor = stringr::str_remove(sample,
                                  paste0("_", timepoint_pattern, "_(EC|noEC)$"))
    )
  
  results <- lapply(rownames(pathway_scores), function(pathway) {
    
    df <- data.frame(
      donor = sample_info$donor,
      condition = sample_info$condition,
      score = pathway_scores[pathway, ]
    )
    
    paired_df <- df %>%
      tidyr::pivot_wider(names_from = condition,
                         values_from = score) %>%
      dplyr::filter(!is.na(EC) & !is.na(noEC)) %>%
      dplyr::mutate(delta = EC - noEC)
    
    data.frame(
      pathway = pathway,
      mean_delta = mean(paired_df$delta),
      sd_delta = sd(paired_df$delta),
      p_value = tryCatch(
        t.test(paired_df$EC, paired_df$noEC, paired = TRUE)$p.value,
        error = function(e) NA
      ),
      concordance = mean(paired_df$delta > 0)
    )
  })
  
  dplyr::bind_rows(results) %>%
    dplyr::mutate(
      padj = p.adjust(p_value, "BH"),
      abs_effect = abs(mean_delta)
    ) %>%
    dplyr::arrange(desc(abs_effect))
}

# running function for timepoints 

progeny_06h <- run_progeny_paired(
  expr_mat,
  timepoint_pattern = "06h"
)

progeny_24h <- run_progeny_paired(
  expr_mat,
  timepoint_pattern = "24h"
)

progeny_72h <- run_progeny_paired(
  expr_mat,
  timepoint_pattern = "72h"
)
cat("\n=== PROGENy 72h ===\n"); print(progeny_72h)
cat("\n=== PROGENy 24h ===\n"); print(progeny_24h)
cat("\n=== PROGENy 06h ===\n"); print(progeny_06h)


######################################
# PROGENy TEMPORAL TRAJECTORY 
######################################
library(forcats)

progeny_all <- bind_rows(
  progeny_06h %>% mutate(timepoint = "06h"),
  progeny_24h %>% mutate(timepoint = "24h"),
  progeny_72h %>% mutate(timepoint = "72h")
) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("06h", "24h", "72h")),
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

# --- Heatmap tile: tutti i pathway insieme ---
pathway_order <- progeny_all %>%
  filter(timepoint == "72h") %>%
  arrange(mean_delta) %>%
  pull(pathway)

progeny_all <- progeny_all %>%
  mutate(pathway = factor(pathway, levels = pathway_order))

ggplot(progeny_all, aes(x = timepoint, y = pathway, fill = mean_delta)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label), color = "black", size = 4, vjust = 0.75) +
  scale_fill_gradient2(
    low = "#0F6E56", mid = "white", high = "#BA7517",
    midpoint = 0, name = "Mean Δ\n(EC − noEC)"
  ) +
  labs(
    title = "PROGENy pathway activity over co-culture time",
    x = "Timepoint", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 11)
  )

# --- Line plot: trajectory keys pathways ---
key_pathways <- c("JAK-STAT", "TGFb", "MAPK", "PI3K", "NFkB", "EGFR")

plot_traj <- progeny_all %>%
  filter(pathway %in% key_pathways)

ggplot(plot_traj, aes(x = timepoint, y = mean_delta, group = pathway, color = pathway)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_delta - sd_delta, ymax = mean_delta + sd_delta),
                width = 0.1, alpha = 0.5) +
  labs(
    title = "Trajectory of key pathway activity (EC − noEC)",
    x = "Timepoint", y = "Mean Δ pathway score",
    color = "Pathway"
  ) +
  theme_classic(base_size = 13)



# TF activity 
# t-test as statistical framework
regulon <- load_dorothea_regulon(c("A", "B", "C"))

run_tf_activity <- function(expr_mat, timepoint_pattern,
                            timepoint_label,
                            control_pattern = "noEC") {
  
  idx <- grep(timepoint_pattern, colnames(expr_mat))
  mat_tp <- expr_mat[, idx, drop = FALSE]
  
  mat_tp <- mat_tp[matrixStats::rowSds(as.matrix(mat_tp)) > 0, ]
  
  tf_act <- viper(
    eset = mat_tp,
    regulon = regulon,
    nes = TRUE,
    method = "scale",
    minsize = 10,
    verbose = FALSE
  )
  
  df <- as.data.frame(tf_act) %>%
    tibble::rownames_to_column("TF") %>%
    tidyr::pivot_longer(-TF, names_to = "sample", values_to = "activity")
  
  sample_info <- data.frame(sample = colnames(mat_tp)) %>%
    dplyr::mutate(
      condition = ifelse(grepl(control_pattern, sample), "noEC", "EC"),
      donor = stringr::str_remove(sample,
                                  paste0("_", timepoint_pattern, "_(EC|noEC)$"))
    )
  
  df <- df %>%
    dplyr::left_join(sample_info, by = "sample")
  
  df %>%
    group_by(TF, donor, condition) %>%
    summarise(activity = mean(activity), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = condition, values_from = activity) %>%
    dplyr::filter(!is.na(EC) & !is.na(noEC)) %>%
    dplyr::mutate(delta = EC - noEC) %>%
    group_by(TF) %>%
    summarise(
      mean_delta = mean(delta),
      sd_delta = sd(delta),
      p_value = tryCatch(
        t.test(EC, noEC, paired = TRUE)$p.value,
        error = function(e) NA
      ),
      concordance = mean(delta > 0),
      .groups = "drop"
    ) %>%
    dplyr::mutate(timepoint = timepoint_label) %>%
    arrange(desc(abs(mean_delta)))
}

tf_72h <- run_tf_activity(expr_mat, "72h", "72h")
tf_24h <- run_tf_activity(expr_mat, "24h", "24h")
tf_06h <- run_tf_activity(expr_mat, "06h", "06h")

tf_all_timepoints <- dplyr::bind_rows(tf_06h, tf_24h, tf_72h)

## Absolute activity for trajectories
run_tf_absolute <- function(expr_mat, timepoint_pattern, timepoint_label, regulon) {
  
  idx    <- grep(timepoint_pattern, colnames(expr_mat))
  mat_tp <- expr_mat[, idx, drop = FALSE]
  mat_tp <- mat_tp[matrixStats::rowSds(as.matrix(mat_tp)) > 0, ]
  
  tf_act <- viper(
    eset    = mat_tp,
    regulon = regulon,
    nes     = TRUE,
    method  = "scale",
    minsize = 10,
    verbose = FALSE
  )
  
  as.data.frame(tf_act) %>%
    tibble::rownames_to_column("TF") %>%
    tidyr::pivot_longer(-TF, names_to = "sample", values_to = "activity") %>%
    dplyr::mutate(
      condition = ifelse(grepl("noEC", sample), "noEC", "EC"),
      donor     = stringr::str_remove(sample,
                                      paste0("_", timepoint_pattern, "_(EC|noEC)$")),
      timepoint = timepoint_label
    )
}

tf_absolute_06h <- run_tf_absolute(expr_mat, "06h", "06h", regulon)
tf_absolute_24h <- run_tf_absolute(expr_mat, "24h", "24h", regulon)
tf_absolute_72h <- run_tf_absolute(expr_mat, "72h", "72h", regulon)

tf_absolute_all <- dplyr::bind_rows(tf_absolute_06h, tf_absolute_24h, tf_absolute_72h)


cat("\n=== TOP TF attivati (tutti i timepoint, t-test paired) ===\n")
tf_activated <- tf_all_timepoints %>%
        filter(mean_delta > 0, p_value < 0.05) %>%
        dplyr::select(TF, mean_delta, p_value, concordance, timepoint) %>%
        arrange(timepoint, p_value)

cat("\n=== TOP TF inattivati (tutti i timepoint, t-test paired) ===\n")
tf_inactivated <- tf_all_timepoints %>%
        filter(mean_delta < 0, p_value < 0.05) %>%
        dplyr::select(TF, mean_delta, p_value, concordance, timepoint) %>%
        arrange(timepoint, p_value)

plot_down <- tf_trajectory_summary %>%
  filter(TF %in% c("BACH2","FOXO1","FOXO2","FOXO3","FOXO4","IRF9","STAT1","STAT2" )) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("06h","24h","72h"))
  )

ggplot(plot_down,
       aes(x = timepoint,
           y = mean_activity,
           group = interaction(TF, condition),
           color = TF,
           linetype = condition)) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  theme_classic(base_size = 13) +
  
  labs(
    title = "TF repressed in β-cell co-culture",
    x = "Timepoint",
    y = "VIPER NES",
    color = "Transcription factor",
    linetype = "Condition"
  )

##### CAUSAL CHAIN ####
ligand_receptor_DE_LY <- read.csv("results/NICHENET_RESULTS/LigandReceptor_EC_to_LY.csv", header = TRUE)

#### manually curated map of receptor associated to 
receptor_pathway_map <- tibble::tribble(
  ~receptor, ~pathway,    ~direction, ~evidence,
  
  # IL23 signaling
  "IL23R",    "JAK-STAT",  1, "IL23 canonical signaling",
  "IL12RB1",  "JAK-STAT",  1, "IL23 receptor complex",
  "IL12RB2",  "JAK-STAT",  1, "IL12 signaling",
  
  # ICOS axis
  "ICOS",     "PI3K",      1, "ICOS-AKT signaling",
  "ICOS",     "MAPK",      1, "ICOS signaling",
  
  # CD28 axis
  "CD28",     "PI3K",      1, "CD28 costimulation",
  "CD28",     "NFkB",      1, "CD28 costimulation",
  "CD28",     "MAPK",      1, "CD28 costimulation",
  
  # CTLA4 inhibitory signaling
  "CTLA4",    "PI3K",     -1, "CTLA4 inhibitory signaling",
  "CTLA4",    "NFkB",     -1, "CTLA4 inhibitory signaling",
  
  # IL7 signaling
  "IL7R",     "JAK-STAT",  1, "IL7 canonical signaling",
  "IL7R",     "PI3K",      1, "IL7 survival signaling",
  
  # Complement receptors
  "C3AR1",    "NFkB",      1, "Complement signaling",
  "C3AR1",    "MAPK",      1, "Complement signaling",
  
  "CR2",      "PI3K",      1, "B/T cell costimulation",
  
  # Siglec
  "SIGLEC6",  "PI3K",     -1, "Inhibitory ITIM receptor",
  "SIGLEC6",  "NFkB",     -1, "Inhibitory ITIM receptor"
)


compute_lr_activity <- function(lr_df){
  
  lr_df %>%
    mutate(
      high_confidence = CellChat & NicheNet,
      ligand_up       = ligand_log2FC > 0,
      receptor_up     = receptor_log2FC > 0
    )
}

lr_df <- compute_lr_activity(ligand_receptor_DE_LY)

lr_receptor_summary <- lr_df %>%
  group_by(receptor) %>%
  summarise(
    receptor_log2FC = max(receptor_log2FC, na.rm = TRUE),
    receptor_padj   = min(receptor_padj, na.rm = TRUE),
    
    receptor_up     = any(receptor_up),
    
    evidence_count  = max(evidence_count, na.rm = TRUE),
    
    n_ligands       = n(),
    
    high_confidence = any(high_confidence),
    
    .groups = "drop"
  )

# progeny score 
link_receptor_pathway <- function(
    receptor_summary,
    progeny_df,
    map_df){
  
  map_df %>%
    
    left_join(
      receptor_summary,
      by = "receptor"
    ) %>%
    
    left_join(
      progeny_df,
      by = "pathway"
    ) %>%
    
    mutate(
      
      pathway_direction =
        sign(mean_delta),
      
      expected_direction =
        ifelse(
          receptor_up,
          direction,
          -direction
        ),
      
      coherence =
        pathway_direction ==
        sign(expected_direction)
      
    ) %>%
    
    arrange(
      desc(high_confidence),
      desc(abs(mean_delta))
    )
}


receptor_pathway_72h <- link_receptor_pathway(
  receptor_summary = lr_receptor_summary,
  progeny_df       = progeny_72h,
  map_df           = receptor_pathway_map
)

## 
library(tibble)

pathway_tf_map <- tibble::tribble(
  ~pathway, ~TF,      ~expected, ~regulatory_logic, ~evidence,
  
  # === JAK-STAT / interferon module ===
  "JAK-STAT", "STAT1",    1, "activation",  "JAK-STAT canonical effector",
  "JAK-STAT", "STAT2",    1, "activation",  "ISGF3 complex, type I IFN",
  "JAK-STAT", "IRF9",     1, "activation",  "ISGF3 complex, type I IFN",
  "JAK-STAT", "IRF4",     1, "activation",  "Th17/Tfh effector priming",           # [T2]
  "JAK-STAT", "BATF",     1, "activation",  "Th17/Tfh differentiation driver",     # [T2]
  "JAK-STAT", "PRDM1",    1, "activation",  "Terminal effector differentiation",    # [T2]
  "JAK-STAT", "STAT6",    1, "activation",  "Th2 differentiation, IL-4 signaling", # [T2]
  
  # === MAPK / proliferative program ===
  "MAPK", "ETS2",  1, "activation", "ERK-ETS2 activation axis",
  "MAPK", "ELK1",  1, "activation", "ERK-ELK1 SRF cofactor",
  "MAPK", "FOS",   1, "activation", "AP-1 MAPK effector",
  "MAPK", "JUN",   1, "activation", "AP-1 MAPK effector",
  "MAPK", "MYC",   1, "activation", "MAPK-RAS-MYC proliferative axis",
  "MAPK", "MEF2C", 1, "activation", "TCR-mediated activation, cytokine production", # [T2]
  
  # === EGFR / growth signalling ===
  "EGFR", "ETS2", 1, "activation", "EGFR-ERK-ETS2 axis",
  "EGFR", "ELK1", 1, "activation", "EGFR-ERK-ELK1 axis",
  "EGFR", "MYC",  1, "activation", "EGFR-RAS-MYC axis",
  "EGFR", "JUN",  1, "activation", "EGFR-AP1 signaling",
  
  # === PI3K / metabolic and survival program ===
  "PI3K", "MYC",    1, "activation", "PI3K-mTOR-MYC axis",
  "PI3K", "SREBF2", 1, "activation", "PI3K-mTOR-SREBP lipogenesis",
  "PI3K", "FOXO1", -1, "repression", "PI3K-AKT inhibits FOXO1",
  "PI3K", "FOXO3", -1, "repression", "PI3K-AKT inhibits FOXO3",
  "PI3K", "FOXO4", -1, "repression", "PI3K-AKT inhibits FOXO4",
  
  # === NFkB / immune checkpoint balance ===
  "NFkB", "RELA",  1, "activation", "NFkB canonical subunit",
  "NFkB", "NFKB1", 1, "activation", "NFkB canonical subunit",
  "NFkB", "BACH2", -1, "repression", "NFkB represses BACH2",
  
  # === TGFb / tolerance axis ===
  "TGFb", "SMAD2", 1, "activation", "TGFb canonical effector",
  "TGFb", "SMAD3", 1, "activation", "TGFb canonical effector",
  "TGFb", "BACH2", 1, "activation", "TGFb maintains BACH2 expression",
  
  # === Hypoxia / stress metabolic response ===
  "Hypoxia", "HIF1A",  1, "activation", "HIF1A master hypoxia regulator",
  "Hypoxia", "MYC",    1, "activation", "HIF1A-MYC metabolic cooperation",
  "Hypoxia", "SREBF2", 1, "activation", "Hypoxia-lipid metabolism",
  
  # === WNT / T cell development ===
  "WNT", "TCF7", 1, "activation", "WNT-TCF1 T cell naive/stem-like identity",
  "WNT", "TCF3", 1, "activation", "E2A, T cell identity and development"  # [T2]
)




link_pathway_tf <- function(progeny_df, tf_df, pathway_tf_map) {
  
  # rimuovi TF wtih low importance
  tf_df_filtered <- tf_df %>%
    filter(abs(mean_delta) > 0.5, p_value < 0.05)
  
  pathway_tf_map %>%
    left_join(
      progeny_df %>%
        dplyr::select(pathway,
                      pathway_delta = mean_delta,
                      pathway_padj  = padj),
      by = "pathway"
    ) %>%
    left_join(
      tf_df_filtered %>%
        dplyr::select(TF,
                      tf_delta = mean_delta,
                      tf_pval  = p_value),
      by = "TF"
    ) %>%
    filter(!is.na(tf_delta)) %>%
    mutate(
	# for repression the direction is inverted
      expected_tf_direction = ifelse(
        regulatory_logic == "activation",
        sign(pathway_delta),
        -sign(pathway_delta)
      ),
      coherence = sign(tf_delta) == expected_tf_direction
    ) %>%
    arrange(desc(coherence), desc(abs(pathway_delta)))
}

pathway_tf_72h <- link_pathway_tf(
  progeny_df     = progeny_72h,
  tf_df          = tf_72h,
  pathway_tf_map = pathway_tf_map
)


## TF + pathway 
library(dplyr)
library(ggplot2)

#selected_tf built from casual map pathway -> TF already defined
selected_tf <- pathway_tf_map %>%
  dplyr::select(TF, pathway, regulatory_logic) %>%
  distinct()

tf_plot_summary <-
  tf_absolute_all %>%
  
  filter(TF %in% selected_tf$TF) %>%
  
  left_join(
    selected_tf,
    by = "TF"
  ) %>%
  
  group_by(pathway,
           TF,
           regulatory_logic,
           timepoint,
           condition) %>%
  
  summarise(
    mean_activity = mean(activity),
    sd_activity   = sd(activity),
    .groups = "drop"
  ) %>%
  
  mutate(
    timepoint = factor(
      timepoint,
      levels = c("06h","24h","72h")
    )
  )

ggplot(
  tf_plot_summary,
  aes(
    x = timepoint,
    y = mean_activity,
    colour = TF,
    group = interaction(TF, condition),
    linetype = condition
  )
) +
  
  geom_line(linewidth = 1) +
  
  geom_point(size = 2.8) +
  
  facet_wrap(
    ~ pathway,
    scales = "free_y"
  ) +
  
  scale_linetype_manual(
    values = c(
      EC = "solid",
      noEC = "dashed"
    )
  ) +
  
  theme_classic(base_size = 10) +
  
  labs(
    x = "Time point",
    y = "VIPER activity (NES)",
    colour = "Transcription factor",
    linetype = "Condition"
  )
# TEST 

library(dplyr)
library(ggplot2)

selected_tf <- pathway_tf_map %>%
  dplyr::select(TF, pathway, regulatory_logic) %>%
  distinct()

# Statistiche per TF x timepoint (dal t-test paired)
tf_sig <- tf_all_timepoints %>%
  dplyr::select(TF, timepoint, p_value, concordance)

tf_plot_summary <- tf_absolute_all %>%
  filter(TF %in% selected_tf$TF) %>%
  left_join(selected_tf, by = "TF") %>%
  group_by(pathway, TF, regulatory_logic, timepoint, condition) %>%
  summarise(
    mean_activity = mean(activity),
    sd_activity   = sd(activity),
    .groups = "drop"
  ) %>%
  mutate(timepoint = factor(timepoint, levels = c("06h","24h","72h"))) %>%
  left_join(tf_sig, by = c("TF", "timepoint")) %>%
  mutate(
    significant = !is.na(p_value) & p_value < 0.05,
    robust = significant & concordance %in% c(0, 1)
  )


## 
plot_pathway_tf <- function(pathway_name, data = tf_plot_summary) {
  
  df <- data %>% filter(pathway == pathway_name)
  
  ggplot(df, aes(x = timepoint, y = mean_activity,
                 group = interaction(TF, condition),
                 color = TF, linetype = condition)) +
    
    geom_line(linewidth = 1) +
    
    geom_point(aes(shape = robust), size = 3) +
    
    scale_linetype_manual(values = c(EC = "solid", noEC = "dashed")) +
    
    scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 1),
      name = "p < 0.05 &\nfull concordance"
    ) +
    
    theme_classic(base_size = 13) +
    
    labs(
      title = paste0(pathway_name, " — TF activity over time"),
      x = "Time point", y = "VIPER activity (NES)",
      colour = "Transcription factor", linetype = "Condition"
    )
}

pathways <- unique(tf_plot_summary$pathway)

for (p in pathways) {
  print(plot_pathway_tf(p))
}

# plottare il delta 

tf_delta_summary <- tf_all_timepoints %>%
  filter(TF %in% selected_tf$TF) %>%
  left_join(selected_tf, by = "TF") %>%
  mutate(
    timepoint = factor(timepoint, levels = c("06h","24h","72h")),
    robust = p_value < 0.05 & concordance %in% c(0, 1)
  )

plot_pathway_delta <- function(pathway_name, data = tf_delta_summary) {
  df <- data %>% filter(pathway == pathway_name)
  
  ggplot(df, aes(x = timepoint, y = mean_delta, group = TF, color = TF)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(linewidth = 1) +
    geom_point(aes(shape = robust), size = 3) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "p<0.05 &\nfull concordance") +
    theme_classic(base_size = 13) +
    labs(title = paste0(pathway_name, " — TF delta (EC − noEC)"),
         x = "Time point", y = "Δ VIPER activity", color = "TF")
}

for (p in pathways) {
  print(plot_pathway_delta(p))
}

library(forcats)

tf_delta_summary <- tf_all_timepoints %>%
  filter(TF %in% selected_tf$TF) %>%
  left_join(selected_tf, by = "TF") %>%
  mutate(
    timepoint = factor(timepoint, levels = c("06h","24h","72h")),
    sig_label = ifelse(p_value < 0.05 & concordance %in% c(0,1), "*", "")
  )

ggplot(tf_delta_summary, aes(x = timepoint, y = fct_rev(TF), fill = mean_delta)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 5, vjust = 0.75) +
  facet_grid(pathway ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#0F6E56", mid = "white", high = "#BA7517",
                       midpoint = 0, name = "Δ activity\n(EC − noEC)") +
  theme_minimal(base_size = 11) +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(title = "TF activity delta by pathway over time", x = "Time point", y = NULL)

## association 

library(patchwork)

plot_pathway_with_tf <- function(pathway_name) {
  
  p_top <- progeny_all %>%
    filter(pathway == pathway_name) %>%
    ggplot(aes(x = timepoint, y = mean_delta, group = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(linewidth = 1, color = "#1D9E75") +
    geom_point(size = 3, color = "#1D9E75") +
    geom_errorbar(aes(ymin = mean_delta - sd_delta, ymax = mean_delta + sd_delta), width = 0.1) +
    theme_classic(base_size = 11) +
    labs(title = pathway_name, y = "Pathway Δ", x = NULL)
  
  p_bottom <- tf_delta_summary %>%
    filter(pathway == pathway_name) %>%
    ggplot(aes(x = timepoint, y = fct_rev(TF), fill = mean_delta)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sig_label), size = 5, vjust = 0.75) +
    scale_fill_gradient2(low = "#0F6E56", mid = "white", high = "#BA7517", midpoint = 0) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank()) +
    labs(y = NULL, x = "Time point", fill = "Δ TF\nactivity")
  
  p_top / p_bottom + plot_layout(heights = c(1, 2))
}

plot_pathway_with_tf("MAPK")

