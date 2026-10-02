# =============================================================================
# Single-cell RNA-seq analysis: dissociated human liver
# -----------------------------------------------------------------------------
# Part of: "Evaluating the Reliability of Cell-Type Biomarkers in Spatial
#           Transcriptomics"
#
# This script processes dissociated human liver single-cell RNA-seq data
# (Liver Cell Atlas; Guilliams et al., 2022) to:
#   1. Perform quality control, normalisation, clustering and annotation
#   2. Curate candidate cell-type markers from the Human Protein Atlas
#   3. Score each marker by fold-change, specificity (pct.1/pct.2) and
#      classification performance (AUC)
#   4. Produce feature plots, violin plots and dot plots supporting the
#      final validated marker set
#
# Platform: Seurat (v5) on the BlueBEAR HPC service, University of Birmingham
# =============================================================================


# ---- Libraries --------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(clustree)       # visualise clustering across resolutions
library(purrr)
library(Polychrome)     # colour palettes for many clusters

# presto speeds up differential expression (FindMarkers / FindAllMarkers)
install.packages("devtools")
devtools::install_github("immunogenomics/presto")

# Use cairo for raster graphics (required on the HPC for reliable plotting)
options(bitmapType = "cairo")


# =============================================================================
# 1. LOAD DATA
# =============================================================================

# Clear any existing objects from a previous run
rm(liv, liv.data)

# Load the raw 10x count matrix (gene.column = 1 for this reference build)
liv.data <- Read10X(
  data.dir = "/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/Data/ABI000_Liver/livercellatlasdownload/rawData_human/countTable_human",
  gene.column = 1
)

# Create the Seurat object; keep genes in >=3 cells and cells with >=200 genes
liv <- CreateSeuratObject(counts = liv.data, project = "liver",
                          min.cells = 3, min.features = 200)
liv

# Load the Human Protein Atlas candidate marker set (gene + cell type)
downloaded_markers_filepath <- "/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/cell_markers_HPA.csv"
downloaded_markers <- read.csv(downloaded_markers_filepath)


# =============================================================================
# 2. QUALITY CONTROL
# =============================================================================

# Percentage of mitochondrial reads per cell (high = low-quality/dying cells)
liv[["percent.mt"]] <- PercentageFeatureSet(liv, pattern = "^MT-")
summary(liv$percent.mt)

# QC violin plots: genes/cell, counts/cell, mitochondrial %
vln <- VlnPlot(liv, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               ncol = 3, pt.size = 0)
vln
ggsave("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/plots/liver plots/QC_vln_plot.png",
       plot = vln, width = 12, height = 4, dpi = 150)

# QC scatter plots to inspect relationships between metrics
plot1 <- FeatureScatter(liv, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = FALSE)
plot2 <- FeatureScatter(liv, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE)
scatters <- plot1 + plot2
ggsave("QC_scatter_liv.png", plot = scatters, width = 12, height = 5, dpi = 150)

# Apply QC filter: drop low-gene cells and high-mitochondrial cells
liv <- subset(liv, subset = nFeature_RNA > 200 & percent.mt < 20)

# Report how many cells remain after QC
n_cells_after_QC <- ncol(liv)
cat("Number of cells passing QC:", n_cells_after_QC, "\n")

summary(liv$nFeature_RNA)
summary(liv$percent.mt)

# Compare cell counts before and after QC
cat("Cells before QC:",
    ncol(CreateSeuratObject(counts = liv.data, project = "liver",
                            min.cells = 3, min.features = 200)), "\n")
cat("Cells after QC:", ncol(liv), "\n")


# =============================================================================
# 3. NORMALISATION, FEATURE SELECTION, SCALING, PCA
# =============================================================================

# Log-normalise the counts
liv <- NormalizeData(liv, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify the 2,000 most variable genes
liv <- FindVariableFeatures(liv, selection.method = "vst", nfeatures = 2000)

# Label the 10 most variable genes on the variable-feature plot
top10 <- head(VariableFeatures(liv), 10)
plot1 <- VariableFeaturePlot(liv)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
varfeat <- plot1 + plot2
ggsave("VariableFeatures_liv.png", plot = varfeat, width = 12, height = 5, dpi = 150)

# Scale the variable features (only variable features are needed for PCA)
all.genes <- rownames(liv)
liv <- ScaleData(liv, features = VariableFeatures(liv))

# Run PCA on the variable features
liv <- RunPCA(liv, features = VariableFeatures(object = liv))
print(liv[["pca"]], dims = 1:5, nfeatures = 5)

# Inspect PCA loadings and structure
VizDimLoadings(liv, dims = 1:2, reduction = "pca")
DimPlot(liv, reduction = "pca", raster = FALSE) + NoLegend()
DimHeatmap(liv, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(liv, dims = 1:15, cells = 500, balanced = TRUE)

# Elbow plot: used to decide how many PCs to retain
ElbowPlot(liv)


# =============================================================================
# 4. CLUSTERING AND UMAP
# =============================================================================

# Build the shared nearest-neighbour graph (first 15 PCs)
liv <- FindNeighbors(liv, dims = 1:15)

# Cluster across a range of resolutions to choose a sensible value
liv <- FindClusters(liv, resolution = seq(0.1, 1.5, 0.1), algorithm = 1)

# clustree: shows how clusters split/merge as resolution increases
clustree_plot <- clustree(liv)
ggsave("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/plots/liver plots/clustree.png",
       plot = clustree_plot, height = 13, width = 8)

# Compute the UMAP embedding (same PCs as FindNeighbors)
liv <- RunUMAP(liv, dims = 1:15)

# Compare several candidate resolutions side by side
resolutions_to_explore <- c("RNA_snn_res.0.3", "RNA_snn_res.0.4",
                            "RNA_snn_res.0.5", "RNA_snn_res.0.6")
(umap_multiple_resolution <- DimPlot(liv, reduction = "umap",
                                     group.by = sort(resolutions_to_explore),
                                     ncol = 2, label = TRUE))
ggsave("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/plots/liver plots/umap_resolution_comparison.png",
       plot = umap_multiple_resolution, width = 12, height = 10)

# Set the chosen clustering resolution (0.3)
Idents(liv) <- "RNA_snn_res.0.3"
head(Idents(liv), 5)


# =============================================================================
# 5. MARKER EXPLORATION BY LINEAGE (pre-filtering)
# =============================================================================
# Inspect HPA candidate markers visually, grouped by cell type, before applying
# any quantitative thresholds.

plotdir <- "/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/plots/liver plots/HPA_marker_feature_plots"

# Split HPA markers into lineage groups and keep only genes present in the data
myeloid_markers <- subset(downloaded_markers, cell.type == "Myeloid")$gene
intersect(myeloid_markers, rownames(liv))

plasma_markers <- subset(downloaded_markers, cell.type == "Plasma cells")$gene
intersect(plasma_markers, rownames(liv))

b_cell_markers <- subset(downloaded_markers, cell.type == "B-cells")$gene

t_cell_markers <- subset(downloaded_markers, cell.type == "T-cells")$gene
intersect(t_cell_markers, rownames(liv))

endo_markers <- subset(downloaded_markers, cell.type == "Endothelial cells")$gene

# Reload the full HPA marker list and tag its source
all_lineage_markers <- read.csv("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/cell_markers_HPA.csv")
all_lineage_markers$Type <- "HPA"

# Identify the lineages present
lineages_with_markers <- unique(all_lineage_markers$cell.type)
lineages_with_markers <- setdiff(lineages_with_markers, "")

# Base output directory for the per-lineage marker plots
currLevel_resolution_analysis_outdir <- "/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/plots"

# Loop over each lineage, saving violin + feature plots of its markers
for(marker_type in unique(all_lineage_markers$Type)){

  resolution_markerplots_outdir <-
    file.path(currLevel_resolution_analysis_outdir,
              paste0(marker_type, "LineageMarkers_Plots"))
  if(!dir.exists(resolution_markerplots_outdir)){
    dir.create(resolution_markerplots_outdir, recursive = TRUE)
  }

  for(i in 1:length(lineages_with_markers)){
    lineage_name <- lineages_with_markers[i]

    # Markers for this lineage, present in the dataset
    lineage_markers <- subset(all_lineage_markers,
                              cell.type == lineage_name & Type == marker_type)$gene
    lineage_markers <- na.omit(lineage_markers)
    lineage_markers <- intersect(lineage_markers, rownames(liv))

    # Skip lineages with no detectable markers (write a warning placeholder)
    n_lineage_markers <- length(lineage_markers)
    if(n_lineage_markers < 1) {
      file.create(file.path(resolution_markerplots_outdir,
                            paste("WARNING - No", marker_type,
                                  "Markers For", lineage_name,
                                  "Are Present In Dataset")))
      next
    }

    # Lay out the grid based on marker count
    n_col <- ceiling(sqrt(n_lineage_markers)) + 1
    n_row <- ceiling(n_lineage_markers / n_col)

    message("##############################")
    message(paste0("Plotting ", lineage_name, " ", marker_type, " markers"))
    message("##############################")

    # Violin plot of this lineage's markers
    lineage_check_violin <- VlnPlot(liv, lineage_markers, alpha = 0, ncol = n_col)
    ggsave(filename = file.path(resolution_markerplots_outdir,
                                paste0(marker_type, "Markers_", lineage_name, "_Vln.png")),
           plot = lineage_check_violin,
           height = n_row * 3 + 1, width = n_col * 3 + 2, limitsize = FALSE)

    # Feature plot of this lineage's markers
    lineage_check_feature <- FeaturePlot(liv, lineage_markers, ncol = n_col,
                                         reduction = "umap")
    ggsave(filename = file.path(resolution_markerplots_outdir,
                                paste0(marker_type, "Markers_", lineage_name, "_Ftr.png")),
           plot = lineage_check_feature,
           height = n_row * 3 + 1, width = n_col * 3.5 + 2, limitsize = FALSE)
  }
}

# Collect the lineage marker groups into one table and save
marker_list <- list(
  Myeloid     = myeloid_markers,
  Plasma      = plasma_markers,
  B_cells     = b_cell_markers,
  T_cells     = t_cell_markers,
  Endothelial = endo_markers
)
liver_hpa_markers <- bind_rows(lapply(names(marker_list), function(ct){
  data.frame(cell_type = ct, gene = marker_list[[ct]])
}))
write.csv(liver_hpa_markers, "liver_HPA_markers_present.csv", row.names = FALSE)

# Save paged feature + violin plots for every marker, grouped by lineage
outdir <- "/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/plots/ALL_MARKERS_per_cell_type_feature_plots/"

for(ct in names(marker_list)){
  genes <- intersect(marker_list[[ct]], rownames(liv))
  pages <- split(genes, ceiling(seq_along(genes) / 9))   # 9 genes per page
  for(i in seq_along(pages)){
    FeaturePlot(liv, features = pages[[i]], ncol = 3, raster = FALSE)
    ggsave(file.path(outdir, paste0("feature_", ct, "_", i, ".png")),
           width = 14, height = 14, limitsize = FALSE)
  }
}

for(ct in names(marker_list)){
  genes <- intersect(marker_list[[ct]], rownames(liv))
  pages <- split(genes, ceiling(seq_along(genes) / 6))   # 6 genes per page
  for(i in seq_along(pages)){
    p <- VlnPlot(liv, features = pages[[i]], ncol = 3, pt.size = 0)
    ggsave(file.path(outdir, paste0("vln_", ct, "_", i, ".png")),
           plot = p, width = 14, height = 14, limitsize = FALSE)
  }
}


# =============================================================================
# 6. DIFFERENTIAL EXPRESSION: MARKER SCORING (log2FC, pct, AUC)
# =============================================================================

# Positive markers for every cluster (Wilcoxon, default thresholds)
liv.markers <- FindAllMarkers(liv, only.pos = TRUE)
liv.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)
View(liv.markers)

# Save the processed object (checkpoint) and confirm it reloads
saveRDS(liv, "liver_runthrough_seurat.rds")
file.info("liver_runthrough_seurat.rds")$size
test <- readRDS("liver_runthrough_seurat.rds"); rm(test)

# ---- AUC (ROC) scoring of the HPA markers -----------------------------------
# Re-run per cluster using the ROC test to obtain AUC (classification power)
liv <- readRDS("liver_runthrough_seurat.rds")
Idents(liv) <- "RNA_snn_res.0.3"

# HPA marker genes actually present in the data
hpa_present <- intersect(unique(downloaded_markers$gene), rownames(liv))

clusters <- as.character(levels(Idents(liv)))
roc_list <- vector("list", length(clusters))
names(roc_list) <- clusters

# For each cluster, score every HPA marker by AUC (no thresholds applied)
for (i in seq_along(clusters)) {
  cluster_id <- clusters[i]
  m <- FindMarkers(liv, ident.1 = cluster_id, features = hpa_present,
                   test.use = "roc", only.pos = TRUE,
                   logfc.threshold = 0, min.pct = 0)
  m$gene <- rownames(m)
  m$cluster <- cluster_id
  roc_list[[i]] <- m
}
liv.markers.roc <- dplyr::bind_rows(roc_list)
table(liv.markers.roc$cluster)

# ---- Combined table: AUC + log2FC + pct.1 + pct.2 per HPA marker ------------
# Wilcoxon stats (log2FC, pct) for every HPA marker, no filtering
wil <- FindAllMarkers(liv, only.pos = TRUE, features = hpa_present,
                      logfc.threshold = 0, min.pct = 0)

# Join AUC with log2FC/pct and attach each gene's HPA cell type
auc_table <- liv.markers.roc %>%
  select(gene, cluster, myAUC, power) %>%
  full_join(wil %>% select(gene, cluster, avg_log2FC, pct.1, pct.2, p_val_adj),
            by = c("gene", "cluster")) %>%
  left_join(downloaded_markers[, c("gene", "cell.type")], by = "gene") %>%
  rename(HPA_cell_type = cell.type) %>%
  mutate(pct.dif = round(pct.1 - pct.2, 3)) %>%
  arrange(HPA_cell_type, gene, desc(myAUC))
write.csv(auc_table, file.path(base_out, "HPA_markers_AUC_full.csv"), row.names = FALSE)


# =============================================================================
# 7. AUC THRESHOLD SENSITIVITY
# =============================================================================
# Explore how many markers are retained as the AUC cutoff is made stricter.

total_markers <- nrow(liv.markers.roc)
auc_thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9)

# Retained vs removed at each AUC threshold
auc_summary <- lapply(auc_thresholds, function(th){
  retained <- sum(liv.markers.roc$myAUC >= th)
  removed  <- total_markers - retained
  data.frame(AUC_cutoff = th, Retained = retained, Removed = removed,
             Pct_retained = round(100 * retained / total_markers, 1))
}) %>% bind_rows()
print(auc_summary)
write.csv(auc_summary, file.path(base_out, "AUC_markers_retained_vs_removed.csv"),
          row.names = FALSE)

# Per-cluster retention at each threshold
auc_by_cluster <- lapply(auc_thresholds, function(th){
  liv.markers.roc %>% filter(myAUC >= th) %>% count(cluster) %>%
    mutate(AUC_cutoff = th)
}) %>% bind_rows()
write.csv(auc_by_cluster, file.path(base_out, "AUC_markers_retained_per_cluster.csv"),
          row.names = FALSE)

# Bar plot of markers retained at each AUC cutoff
auc_retention_plot <- ggplot(auc_summary, aes(x = factor(AUC_cutoff), y = Retained)) +
  geom_col() +
  geom_text(aes(label = Retained), vjust = -0.3) +
  labs(x = "AUC threshold", y = "Markers retained",
       title = "HPA markers retained at each AUC cutoff") +
  theme_minimal()
ggsave(file.path(base_out, "plots", "AUC_markers_retained.png"),
       auc_retention_plot, width = 8, height = 6)


# =============================================================================
# 8. COMBINED FEATURE & VIOLIN PLOTS FOR AUC-PASSING MARKERS
# =============================================================================
# For each cluster, plot the markers passing an AUC cutoff, annotated with AUC.

library(patchwork)
auc_cutoff <- 0.6

# Keep markers passing the cutoff and present in the data
auc_pass <- liv.markers.roc %>%
  filter(myAUC >= auc_cutoff) %>%
  filter(gene %in% rownames(liv))

# Output folders
dir.create(file.path(base_out, "plots", "AUC_markers_by_cluster", "FeaturePlots"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_out, "plots", "AUC_markers_by_cluster", "ViolinPlots"),
           recursive = TRUE, showWarnings = FALSE)

for(cl in sort(unique(auc_pass$cluster))){

  cluster_markers <- auc_pass %>% filter(cluster == cl) %>% arrange(desc(myAUC))
  genes <- cluster_markers$gene
  if(length(genes) == 0) next

  # Feature plots, one per marker, titled with its AUC
  feature_list <- lapply(genes, function(g){
    auc <- round(cluster_markers$myAUC[cluster_markers$gene == g], 2)
    FeaturePlot(liv, features = g, reduction = "umap",
                cols = c("lightgrey", "red")) +
      ggtitle(paste0(g, "\nAUC = ", auc))
  })
  feature_page <- wrap_plots(feature_list) +
    plot_annotation(title = paste0("Cluster ", cl, " - Markers Passing AUC >= ", auc_cutoff),
                    subtitle = paste(length(genes), "marker genes"))
  ggsave(file.path(base_out, "plots", "AUC_markers_by_cluster", "FeaturePlots",
                   paste0("Cluster_", cl, "_FeaturePlots.png")),
         feature_page, width = 12, height = 8, dpi = 300)

  # Violin plots, one per marker, titled with its AUC
  violin_list <- lapply(genes, function(g){
    auc <- round(cluster_markers$myAUC[cluster_markers$gene == g], 2)
    VlnPlot(liv, features = g, pt.size = 0) +
      ggtitle(paste0(g, "\nAUC = ", auc))
  })
  violin_page <- wrap_plots(violin_list) +
    plot_annotation(title = paste0("Cluster ", cl, " - Markers Passing AUC >= ", auc_cutoff),
                    subtitle = paste(length(genes), "marker genes"))
  ggsave(file.path(base_out, "plots", "AUC_markers_by_cluster", "ViolinPlots",
                   paste0("Cluster_", cl, "_ViolinPlots.png")),
         violin_page, width = 12, height = 8, dpi = 300)
}


# =============================================================================
# 9. SELECTED MARKER FEATURE / VIOLIN PLOTS (for figures)
# =============================================================================
# Curated per-lineage plots of the chosen marker genes, for presentation.

# ---- Feature plots ----
T_cell_plot <- FeaturePlot(liv, features = c("CD3D","CD3G","CD8A","KLRB1","CD6","THEMIS","CTLA4","PDCD1"),
                           ncol = 4, raster = FALSE)
ggsave(file.path(plotdir, "feature_Tcells.png"), plot = T_cell_plot, width = 16, height = 8, limitsize = FALSE)

B_cell_plot <- FeaturePlot(liv, features = c("MS4A1","CD79B","BANK1","BLK","FCRL1","CR2","PAX5","FCER2"),
                           ncol = 4, raster = FALSE)
ggsave(file.path(plotdir, "feature_Bcells.png"), plot = B_cell_plot, width = 16, height = 8, limitsize = FALSE)

Endo_plot <- FeaturePlot(liv, features = c("CDH5","VWF","CLEC14A","KDR","FLT1","ERG","CD34"),
                         ncol = 4, raster = FALSE)
ggsave(file.path(plotdir, "feature_Endothelial.png"), plot = Endo_plot, width = 16, height = 8, limitsize = FALSE)

Myeloid_plot <- FeaturePlot(liv, features = c("CD14","FCGR3A","MARCO","C1QB","VSIG4","CD300E","MS4A7","TREM2"),
                            ncol = 4, raster = FALSE)
ggsave(file.path(plotdir, "feature_Myeloid.png"), plot = Myeloid_plot, width = 16, height = 8, limitsize = FALSE)

Plasma_plot <- FeaturePlot(liv, features = c("CD79A","MZB1","TNFRSF17","DERL3","FKBP11","CD27"),
                           ncol = 4, raster = FALSE)
ggsave(file.path(plotdir, "feature_Plasma.png"), plot = Plasma_plot, width = 16, height = 8, limitsize = FALSE)

# Save the full marker table
write.csv(liv.markers,
          "/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/liv_markers_all.csv",
          row.names = FALSE)

# ---- Violin plots (log counts) justifying the chosen markers ----
Idents(liv) <- "RNA_snn_res.0.3"

vln_plot <- VlnPlot(liv, features = c("CD3D","CD3G","CD8A","CD6","KLRB1"),
                    layer = "counts", log = TRUE, pt.size = 0, raster = FALSE)
ggsave("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/plots/reason for using specific markers/vln_Tcells.png",
       plot = vln_plot, width = 16, height = 8, limitsize = FALSE)

vln_plot <- VlnPlot(liv, features = c("MS4A1","CD79B","BANK1","BLK","FCRL1","CR2","PAX5","FCER2"),
                    layer = "counts", log = TRUE, pt.size = 0, raster = FALSE)
ggsave("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/plots/reason for using specific markers/vln_Bcells.png",
       plot = vln_plot, width = 16, height = 8, limitsize = FALSE)

vln_plot <- VlnPlot(liv, features = c("CD79A","MZB1","TNFRSF17","DERL3","CD27"),
                    layer = "counts", log = TRUE, pt.size = 0, raster = FALSE)
ggsave("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/plots/reason for using specific markers/vln_plasmacells.png",
       plot = vln_plot, width = 16, height = 8, limitsize = FALSE)

vln_plot <- VlnPlot(liv, features = c("CD14","MARCO","C1QB","VSIG4","MS4A7","TREM2"),
                    layer = "counts", log = TRUE, pt.size = 0, raster = FALSE)
ggsave("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/plots/reason for using specific markers/vln_myeloidcells.png",
       plot = vln_plot, width = 16, height = 8, limitsize = FALSE)

vln_plot <- VlnPlot(liv, features = c("CDH5","CLEC14A","KDR","FLT1","ERG"),
                    layer = "counts", log = TRUE, pt.size = 0, raster = FALSE)
ggsave("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/plots/reason for using specific markers/vln_endothelial_cells.png",
       plot = vln_plot, width = 16, height = 8, limitsize = FALSE)


# =============================================================================
# 10. MANUAL CLUSTER ANNOTATION
# =============================================================================

Idents(liv) <- "RNA_snn_res.0.3"

# B cells and plasma cells co-cluster in cluster 7 -> subcluster to separate them
liv <- FindSubCluster(liv, cluster = "7", graph.name = "RNA_snn",
                      subcluster.name = "sub_cluster", resolution = 0.015)
Idents(liv) <- "sub_cluster"
table(Idents(liv))

# Visual confirmation of the B (MS4A1) vs plasma (MZB1) split
feature_plot <- FeaturePlot(liv, c("MS4A1","MZB1"))
ggsave("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/plots/b_cell_distiguishing_markers.png",
       plot = feature_plot, width = 16, height = 8, limitsize = FALSE)

# Map each (sub)cluster to a fine cell type.
# Only the five lineages of interest are named; all other clusters (hepatocytes
# and other non-target cells) are labelled "Unassigned" (deliberate scoping).
current_cluster_ids <- c("0","1","2","3","4","5","6","7_0","7_1","8","9",
                         "10","11","12","13","14","15","16","17","18","19")
new_cluster_ids <- c(
  "T-cells","Monocytes (myeloid)","Dendritic cells (myeloid)","Unassigned",
  "Kupffer cells/Macrophages (myeloid)","Unassigned","Unassigned",
  "B cells","Plasma cells","Unassigned","Endothelial cells_1","Unassigned",
  "Endothelial cells_2","Macrophages/Monocytes (myeloid)","Unassigned",
  "Unassigned","Proliferating T cells","Proliferating myeloid","Unassigned",
  "Unassigned","Unassigned"
)
Idents(liv) <- plyr::mapvalues(as.character(Idents(liv)),
                               from = current_cluster_ids, to = new_cluster_ids)
liv$celltype <- Idents(liv)

# Annotated UMAP (fine cell types)
clusters <- levels(Idents(liv))
plot2 <- DimPlot(liv, reduction = "umap", label = TRUE, label.size = 3,
                 repel = TRUE, raster = FALSE)
ggsave("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/plots/liv_manual_annot_UMAP_plot.png",
       plot = plot2, width = 12, height = 7)

# Zoomed UMAP of the myeloid compartment only
myeloid_plot <- DimPlot(
  subset(liv, celltype_merged == "Myeloid"),
  reduction = "umap",
  cols = setNames(Polychrome::palette36.colors(36), unique(liv$celltype_merged)),
  label = TRUE, repel = TRUE, raster = FALSE
) + NoLegend()
ggsave("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/plots/liv_myeloid_UMAP_plot.png",
       plot = myeloid_plot, width = 12, height = 7)
saveRDS(liv, file = "liver_runthrough_seurat.rds")


# =============================================================================
# 11. VALIDATED-MARKER DOT PLOTS BY BROAD LINEAGE
# =============================================================================

# Load the final validated marker set (gene + cell_type)
keepers <- read.csv("sc_liver_markers_final.csv")

# Rebuild the subclustering and fine annotation (needed after a fresh session)
Idents(liv) <- "RNA_snn_res.0.3"
liv <- FindSubCluster(liv, cluster = "7", graph.name = "RNA_snn",
                      subcluster.name = "sub_cluster", resolution = 0.015)
Idents(liv) <- "sub_cluster"

current_cluster_ids <- c("0","1","2","3","4","5","6","7_0","7_1","8","9",
                         "10","11","12","13","14","15","16","17","18","19")

# Collapse fine cell types into the five broad lineages (+ Unassigned)
fine_to_broad <- c(
  "T-cells"="T cells", "Proliferating T cells"="T cells",
  "B cells"="B cells", "Plasma cells"="Plasma cells",
  "Monocytes (myeloid)"="Myeloid", "Dendritic cells (myeloid)"="Myeloid",
  "Kupffer cells/Macrophages (myeloid)"="Myeloid",
  "Macrophages/Monocytes (myeloid)"="Myeloid", "Proliferating myeloid"="Myeloid",
  "Endothelial cells_1"="Endothelial", "Endothelial cells_2"="Endothelial",
  "Unassigned"="Unassigned"
)
liv$celltype_merged <- unname(fine_to_broad[as.character(liv$celltype)])

# Remove the two markers excluded for cross-lineage (LSEC) expression -> 42 markers
keepers <- keepers %>% filter(!gene %in% c("MS4A6A","LGMN"))
nrow(keepers)

# Fix lineage order and set broad annotation as the identity
lin_order <- c("T cells","B cells","Plasma cells","Myeloid","Endothelial")
liv$celltype_merged <- factor(liv$celltype_merged, levels = c(lin_order,"Unassigned"))
Idents(liv) <- "celltype_merged"

# Drop Unassigned so only the five lineages are plotted
liv_lin <- subset(liv, celltype_merged %in% lin_order)
liv_lin$celltype_merged <- factor(liv_lin$celltype_merged, levels = lin_order)
Idents(liv_lin) <- "celltype_merged"

# Helper: return the markers for given lineages, ordered by lineage
gene_by_lin <- function(types) {
  keepers %>% filter(cell_type %in% types) %>%
    mutate(cell_type = factor(cell_type, levels = lin_order)) %>%
    arrange(cell_type) %>% pull(gene) %>% intersect(rownames(liv_lin))
}

# Dot plot A: endothelial + myeloid markers
dpA <- DotPlot(liv_lin, features = gene_by_lin(c("Endothelial","Myeloid"))) +
  RotatedAxis() + theme(axis.text.x = element_text(size = 9)) +
  ggtitle("Endothelial and myeloid markers")
ggsave("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/plots/SC_dotplot_endo_myeloid.png",
       dpA, width = 15, height = 5, dpi = 300)

# Dot plot B: plasma + T + B markers
dpB <- DotPlot(liv_lin, features = gene_by_lin(c("Plasma cells","T cells","B cells"))) +
  RotatedAxis() + theme(axis.text.x = element_text(size = 9)) +
  ggtitle("Plasma, T and B cell markers")
ggsave("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/plots/SC_dotplot_lymphoid.png",
       dpB, width = 9, height = 5, dpi = 300)


# =============================================================================
# 12. ZOOMED UMAPs (B/plasma and myeloid subpopulations)
# =============================================================================

# Report subcluster centres (for positioning annotation labels)
emb <- Embeddings(liv, "umap")
for (sc in c("7_0","7_1")) {
  cells <- rownames(liv@meta.data)[liv$sub_cluster == sc]
  cat(sc, "centre: x =", round(median(emb[cells,1]),1),
      " y =", round(median(emb[cells,2]),1), "\n")
}

# ---- B cell / plasma cell zoom (subcluster 7_0 vs 7_1) ----
library(ggplot2)
bplasma_zoom <- DimPlot(
  subset(liv, sub_cluster %in% c("7_0","7_1")),
  reduction = "umap", group.by = "zoom_label",
  label = FALSE, pt.size = 0.6, raster = FALSE,
  cols = c("B cells" = "cadetblue3", "Plasma cells" = "azure4")
) +
  ggtitle("B-cell and Plasma cell subclustering") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  annotate("segment", x = -3, y = 1.5, xend = -1.7, yend = -1.0,
           arrow = arrow(length = unit(0.25,"cm")), linewidth = 0.6) +
  annotate("text", x = -3, y = 1.9, label = "B cells", fontface = "bold", size = 5) +
  annotate("segment", x = 2, y = -3.5, xend = 0.3, yend = -1.6,
           arrow = arrow(length = unit(0.25,"cm")), linewidth = 0.6) +
  annotate("text", x = 2, y = -3.9, label = "Plasma cells", fontface = "bold", size = 5)
ggsave("UMAP_zoom_bplasma.png", bplasma_zoom, width = 6, height = 6, dpi = 300)

# ---- Myeloid subpopulation zoom ----
myeloid_types <- c("Monocytes (myeloid)","Dendritic cells (myeloid)",
                   "Kupffer cells/Macrophages (myeloid)",
                   "Macrophages/Monocytes (myeloid)",
                   "Proliferating myeloid")

# Pastel palette for the myeloid subtypes
pastel <- c(
  "Monocytes (myeloid)"                 = "#F4A9A0",
  "Dendritic cells (myeloid)"           = "#A9C9E8",
  "Kupffer cells/Macrophages (myeloid)" = "#B5D8B1",
  "Macrophages/Monocytes (myeloid)"     = "#CDB4DB",
  "Proliferating myeloid"               = "#F7C59F"
)

# Sanity check: palette names must match the cell-type labels exactly
setdiff(unique(as.character(subset(liv, celltype %in% myeloid_types)$celltype)), names(pastel))

myeloid_zoom <- DimPlot(
  subset(liv, celltype %in% myeloid_types),
  reduction = "umap", group.by = "celltype",
  label = FALSE, pt.size = 0.5, raster = FALSE
) +
  scale_colour_manual(values = pastel) +
  ggtitle("Myeloid subpopulations") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank())

# Manual annotation labels with arrows pointing to each subpopulation
myeloid_zoom <- myeloid_zoom +
  annotate("text", x = -8, y = -2.5, label = "Monocytes", fontface="bold", size=4.5, hjust=0) +
  annotate("segment", x = -5.5, y = -3.5, xend = -4.5, yend = -5.5,
           arrow = arrow(length = unit(0.25,"cm")), linewidth = 0.6) +
  annotate("text", x = 4, y = -9, label = "Dendritic cells", fontface="bold", size=4.5, hjust=0) +
  annotate("segment", x = 3.3, y = -9, xend = -0.3, yend = -9,
           arrow = arrow(length = unit(0.25,"cm")), linewidth = 0.6) +
  annotate("text", x = 4.8, y = -5.5, label = "Proliferating\nmyeloid", fontface="bold", size=4.5, hjust=0) +
  annotate("segment", x = 4.2, y = -6, xend = 1.8, yend = -6.5,
           arrow = arrow(length = unit(0.25,"cm")), linewidth = 0.6) +
  annotate("text", x = 4.5, y = 4, label = "Macrophages /\nMonocytes", fontface="bold", size=4.5, hjust=0) +
  annotate("segment", x = 6, y = 2.5, xend = 9, yend = -0.5,
           arrow = arrow(length = unit(0.25,"cm")), linewidth = 0.6) +
  annotate("text", x = -8.5, y = -16.5, label = "Kupffer cells /\nMacrophages", fontface="bold", size=4.5, hjust=0) +
  annotate("segment", x = -4, y = -15.5, xend = -3, yend = -15,
           arrow = arrow(length = unit(0.25,"cm")), linewidth = 0.6)
ggsave("UMAP_zoom_myeloid.png", myeloid_zoom, width = 8, height = 7, dpi = 300)

# Combine both zoomed panels into one figure (A/B)
library(patchwork)
combined <- myeloid_zoom + bplasma_zoom +
  plot_annotation(
    tag_levels = "A",
    title = "Single-cell RNA-seq: Myeloid and Lymphoid Subpopulations",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )
ggsave("UMAP_zoom_combined.png", combined, width = 14, height = 6, dpi = 300)


# =============================================================================
# 13. CELL-TYPE COMPOSITION BAR PLOT
# =============================================================================

library(ggplot2); library(dplyr); library(scales)

# Broad-lineage cell counts, ordered by frequency
sc_counts <- as.data.frame(table(liv$celltype_merged)) %>%
  setNames(c("cell_type","count")) %>%
  arrange(desc(count))
sc_counts$cell_type <- factor(sc_counts$cell_type, levels = sc_counts$cell_type)

p_sc <- ggplot(sc_counts, aes(cell_type, count, fill = cell_type)) +
  geom_col() +
  scale_y_continuous(breaks = seq(0, 125000, 25000), labels = comma) +
  labs(x = NULL, y = "Cell count",
       title = "Single-cell annotation frequency for each cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave("SC_cellcount_bar.png", p_sc, width = 6, height = 5, dpi = 300)

# =============================================================================
# END
# =============================================================================
