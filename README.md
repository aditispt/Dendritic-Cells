# Decoding Functional and Developmental Trajectories of Tissue-Resident Uterine Dendritic Cells

[📄 Paper (iScience, 2025)](https://www.cell.com/iscience/fulltext/S2589-0042(25)01612-8) •
[📊 GEO: GSE288249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288249)

Single-cell dendritic cell (uDC) workflow (R/Seurat) for preprocessing, clustering, stage-aware expression dynamics, CellChat receptor–ligand analysis, co-expression ratios, and GO enrichment.

## Citation
Singh A, Mor G, *et al.* **Decoding Functional and Developmental Trajectories of Uterine Dendritic Cells.** *iScience.* 2025. 

# Dendritic-Cells
ScRNAseq of dendritic cells

Single-cell dendritic cell (uDC) workflow for preprocessing, clustering, stage-aware expression dynamics, CellChat receptor–ligand analysis, co-expression ratios, and GO enrichment.

Contents
1. Overview
2. Input Data
3. Pipeline Steps
4. Quick Start
5. Outputs
6. Repo Structure
7. Parameters & Notes
8. Troubleshooting
9. Reproducibility
10. Citations


1. Overview
This repository hosts an R/Seurat pipeline to:
Load count and stage metadata
Build a Seurat object, perform QC, normalize, cluster, and run UMAP
Visualize marker genes by cluster
Analyze AXL and other markers across menstrual stages
Run CellChat signaling analysis
Compute co-expression ratios and GO enrichment per cluster
Tested with Seurat v5-style calls (compatible with Seurat v4+).

2. Counts and stage metadata are publicly available on GEO:
GEO Accession: GSE288249

Expected files:
DC_count.txt — gene × cell matrix (rows = genes, columns = cells). Column names must match row names in DC_stage.txt.
DC_stage.txt — rows = cells, includes column x with stage labels:

The script creates:
newclass (1–10) mapping the above stages
One-hot columns for each stage (e.g., Decb, Decp, …)
A Stages row in the combined matrix for convenience

3. Pipeline Steps
Load packages: Seurat, tidyverse, ggplot2, ggpubr, Matrix, dplyr, patchwork, openxlsx, CellChat, clusterProfiler, org.Hs.eg.db, enrichplot, pheatmap, ggrepel.

Read data and align:
stopifnot(all(rownames(DC_stage) == colnames(DC_count)))

Create Seurat object (min.cells = 3, min.features = 200)
QC: mitochondrial % (percent.mt), violin & scatter plots
Normalize → HVGs → Scale → PCA; elbow plot & loadings
Neighbors / Clusters / UMAP: dims 1:12, resolution 0.5
Markers: violin plots for curated lists (cDC1, cDC2, pDC, DC3, moDC, progenitors, subtype panels)

Expression dynamics:
AXL across newclass stages
Combined plot with secondary y-axis for cell counts
Export AXL_gene_expression.xlsx
Function plot_gene_temporal() for multiple markers
CellChat: human, Secreted Signaling; object creation, overexpression, probabilities, pathways, aggregation, netVisual_bubble()
Co-expression ratios: clusters 0/2/3 (AXL with CD1C, PCNA, XCR1)
GO enrichment: per cluster (BP ontology, BH correction; dotplots)

4. Install dependencies
install.packages(c(
  "Seurat","tidyverse","ggpubr","Matrix","dplyr","patchwork",
  "ggplot2","openxlsx","pheatmap","ggrepel","tidyr"
))

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler","org.Hs.eg.db","enrichplot"))

#CellChat:
#install.packages("CellChat")
#or:
#remotes::install_github("sqjin/CellChat")

5. Outputs
QC & DR: PCA, elbow, UMAP, violin/scatter plots
Marker panels: multi-gene violins per curated lists
AXL results:
AXL_gene_expression.xlsx
Combined_AXL_expression_and_cell_count.png
Temporal marker plots: via plot_gene_temporal()
CellChat: bubble plot of pathways/interactions
Co-expression ratios: line/point plot by cluster
GO enrichment: per-cluster dotplots (top terms)

6. Repo Structure
├─ README.md
├─ scripts/
│  └─ dc_pipeline.R
├─ data/                    # optional local copies of DC_count.txt / DC_stage.txt
├─ outputs/                 # optional figures/tables
└─ renv/ renv.lock          # if using renv

7. Parameters & Notes
Dims for PCA/UMAP/Neighbors: 1:12 (adjust if elbow suggests otherwise)
Clustering resolution: 0.5
Stage order: menstrual → proliferative → secretory → decidua
CellChat groups: converts seurat_clusters to 1-indexed integers
Gene names: assumes human symbol case; adjust mitochondrial pattern (^MT-) if needed

8. Troubleshooting
Row/column mismatch: check colnames(DC_count) == rownames(DC_stage)
Low AXL expression: code filters AXL > 0; remove filter to see all cells
Secondary y-axis scaling: adjust scaling_factor if curves are flat
CellChat errors:
Ensure CellChatDB.human is loaded
Data should be log-normalized
Very small clusters may be removed by min.cells = 10; lower if needed

9. Reproducibility
install.packages("renv")
renv::init()
renv::snapshot()
sessionInfo()

10. Citations
Please cite:
Seurat: Hao et al., Cell 2021; Stuart et al., Cell 2019
CellChat: Jin et al., Nat Commun 2021
clusterProfiler / enrichplot / org.Hs.eg.db: Yu et al., OMICS 2012; Wu et al., Innovation 2021
openxlsx: Alexander Walker (CRAN)
ggplot2 / tidyverse: Wickham et al.

