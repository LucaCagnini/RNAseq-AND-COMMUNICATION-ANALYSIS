# **RNA-seq Data and Cell–Cell Communication Analysis**

## **Introduction**

R scripts for the analysis of bulk RNA-seq data from two cell types grown in a co-culture system, with a focus on differential expression, transcription factor (TF) and pathway activity, and inferred cell–cell communication.
This repository contains code only. Input data and analysis results are not included (see Data).


## **Table of Contents**

- [Overview](#overview)
- [Repository structure](#repository-structure)
- [Analysis workflow](#analysis-workflow)
- [Requirements](#requirements)
- [How to run](#how-to-run)
- [Notes](#notes)
- [Author](#author)
- [License](#license)

## **Overview**

Two cell types are co-cultured without direct contact and profiled by bulk RNA-seq at multiple timepoints:

- EC: pancreatic β-cell line 
- LY: primary activated human lymphocytes

The analysis answers three questions:

- How does each cell type change its transcriptional program in co-culture? (differential expression)
- Which TFs and pathways drive these changes? (TF and pathway activity)
- Which ligand–receptor interactions could mediate the crosstalk between the two cell types? (communication inference)


## **Repository structure**

```
RNAseq_analysis/
├── R_scripts/
│   ├── functions_dds.R        # helper functions shared by the other scripts
│   ├── DE_Final.R             # differential expression analysis
│   ├── TF_supercomplex.R      # TF and pathway activity
│   ├── CellChat_Final.R       # ligand–receptor inference with CellChat
│   └── Nichenet_Final.R       # ligand activity inference with NicheNet
├── data/                      # input data (not tracked by git)
├── results/                   # outputs (not tracked by git)
│   ├── DE_RESULTS/
│   ├── TF_RESULTS/
│   ├── CELLCHAT_RESULTS/
│   └── NICHENET_RESULTS/
├── .here                      # marks the project root for the `here` package
├── .gitignore
└── README.md
```

## **Analysis workflow**

| Step | Script | What it does |
|------|--------|--------------|
| 0 | `functions_dds.R` | Helper functions, loaded with `source()` by the other scripts. |
| 1 | `DE_Final.R` | Differential expression with DESeq2 using a paired design by donor across timepoints. Produces DEG tables per cell type and timepoint. |
| 2 | `TF_supercomplex.R` | Transcription factor activity (VIPER with DoRothEA regulons) and pathway activity (PROGENy), plus enrichment analyses. |
| 3 | `CellChat_Final.R` | Ligand–receptor inference with CellChat, restricted to secreted signaling and ECM–receptor interactions (cell–cell contact is excluded because of the transwell setup). |
| 4 | `Nichenet_Final.R` | Ligand activity prediction with NicheNet, adapted for bulk RNA-seq, run separately on up- and down-regulated gene sets to preserve directionality. |

Run `DE_Final.R` first: the downstream scripts use its outputs. Steps 2–4 are independent of each other.

<!-- TODO: check the order and dependencies between scripts and correct this table if needed. -->

## **Requirements**

- R (developed with R 4.5.x)
- Main packages: `DESeq2`, `edgeR`, `CellChat`, `nichenetr`, `progeny`, `viper`, `dorothea`, `clusterProfiler`, `ReactomePA`, `tidyverse`, `ggplot2`, `patchwork`, `ComplexHeatmap`, `pheatmap`, `ggrepel`, `here`, `readxl`, `conflicted`

<!-- TODO: add exact package versions, e.g. with renv::snapshot() (renv.lock) or by pasting sessionInfo(). -->

Several packages mask each other's functions (e.g. `plyr` and `dplyr`), so the scripts use explicit `package::function()` calls where needed.

## **Data**

No data are distributed with this repository. To run the scripts you need:

| Input | Description |
|-------|-------------|
| `FINAL_DDS.rds` | DESeq2 object for the lymphocyte (LY) samples. |
| `DDS_RESULTS.rds` | DESeq2 results for LY. |
| `DDS_EC` | DESeq2 object for the β-cell (EC) samples. |
| `Cytokines_Concentration_Clean.txt` | Secreted protein concentrations (multiplex immunoassay). |
| `Table_genes_receptors.xlsx` | Table of genes and receptors used to link the analyses to the secretome panel. |
| NicheNet prior files | Ligand–target matrix, ligand–receptor network and weighted networks, downloaded from the [NicheNet Zenodo repository](https://zenodo.org/). |

Place the input files in `data/` at the project root. Paths are built with `here::here()`, so the scripts work from any working directory as long as the `.here` file is present.

<!-- TODO: describe how raw reads were processed upstream (e.g. nf-core/rnaseq) and where a collaborator can obtain the DESeq2 objects. -->



## **Notes**

- Only the RNA-seq, TF/pathway and communication analyses are included.
- Interpretation of the results is in the accompanying thesis.

## **Author**

<!-- TODO: name, affiliation and contact. -->

## License

<!-- TODO: choose a license before making the repository public, or state "All rights reserved". -->
