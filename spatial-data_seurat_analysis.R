# Spatial transcriptomics analysis: CosMx human liver (1,000-plex)
# Part of: "Evaluating the Reliability of Cell-Type Biomarkers in Spatial
#           Transcriptomics"
#
# This script validates the single-cell-derived marker set (from
# 01_singlecell_analysis.R) against CosMx spatial molecular imaging of human
# liver (Bruker Spatial Biology, 1,000-plex panel). It:
#   1. Assesses QC of the provided spatial object
#   2. Groups the supplied (fine) cell-type annotation into five broad lineages
#   3. Tests whether each marker is upregulated in its expected lineage spatially
#   4. Performs an independent clustering + annotation of the same data
#   5. Produces annotated UMAPs, isolation panels, dot plots and bar charts
#
# Platform: Seurat (v5) on the BlueBEAR HPC service, University of Birmingham


# Libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(clustree)       # clustering across resolutions
library(purrr)
library(Polychrome)
library(scales)         # comma-formatted axis labels

options(bitmapType = "cairo")   # reliable raster output on the HPC

# Output directory for all spatial plots
sp_base <- "/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/spatial"
dir.create(file.path(sp_base, "plots"), recursive = TRUE, showWarnings = FALSE)

# Load the provided CosMx spatial object (Bruker, Slide 3, normal liver)
sp <- readRDS("/rds/projects/g/gilberts-spatial-biology-image-analysis/data/CosMx/Data/Default/01_ConstructedSeuratObjects_Slide/seuratObject_Slide3_Bruker_Liver.RDS")


# 1. QUALITY CONTROL
# The object arrives pre-QC'd by the facility; here QC metrics are reported
# (rather than re-derived) to confirm data quality.

cells_before_QC <- ncol(sp)
cat("Number of spatial cells:", cells_before_QC, "\n")

# Mitochondrial transcript percentage per cell
sp[["percent.mt"]] <- PercentageFeatureSet(sp, pattern = "^MT-")

# QC summary statistics
cat("\nNumber of detected genes per cell:\n");   print(summary(sp$nFeature_RNA))
cat("\nNumber of transcripts per cell:\n");       print(summary(sp$nCount_RNA))
cat("\nMitochondrial transcript percentage:\n");  print(summary(sp$percent.mt))

# QC violin plots
cosmx_qc_vln <- VlnPlot(sp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        ncol = 3, pt.size = 0)
cosmx_qc_vln
ggsave(file.path(sp_base, "plots", "SPATIAL_QC_violin_plot.png"),
       cosmx_qc_vln, width = 12, height = 4, dpi = 300)

# QC scatter plots
qc_scatter1 <- FeatureScatter(sp, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = FALSE)
qc_scatter2 <- FeatureScatter(sp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE)
cosmx_qc_scatter <- qc_scatter1 + qc_scatter2
ggsave(file.path(sp_base, "plots", "SPATIAL_QC_scatter_plot.png"),
       cosmx_qc_scatter, width = 12, height = 5, dpi = 300)


# 2. MARKER AVAILABILITY & SUPPLIED ANNOTATION

# Load the HPA candidate marker set
downloaded_markers_filepath <- "/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/cell_markers_HPA.csv"
downloaded_markers <- read.csv(downloaded_markers_filepath)

# How many HPA markers are present on the CosMx panel?
hpa_present <- intersect(unique(downloaded_markers$gene), rownames(sp))
length(hpa_present)

# The supplied (reference) fine cell-type annotation
table(sp$cellType)


# 3. UMAP & TISSUE MAP (supplied annotation)

Idents(sp) <- "cellType"

# UMAP coloured by the supplied reference annotation
umap_plot <- DimPlot(sp,
                     reduction = "approximateUMAP_bySlide",  # provided UMAP embedding
                     group.by  = "cellType",
                     label = TRUE, label.size = 3, repel = TRUE, raster = FALSE) +
  NoLegend()
ggsave(file.path(sp_base, "plots", "SPATIAL_UMAP_cellType.png"),
       umap_plot, width = 12, height = 9)

# Spatial tissue map: cells plotted at their physical coordinates
df <- sp@meta.data[, c("x_slide_mm", "y_slide_mm", "cellType")]
tissue_plot <- ggplot(df, aes(x_slide_mm, y_slide_mm, colour = cellType)) +
  geom_point(size = 0.1) +
  coord_fixed() +                 # preserve true tissue proportions
  theme_void() +
  guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave(file.path(sp_base, "plots", "SPATIAL_tissue_cellType.png"),
       tissue_plot, width = 14, height = 12, dpi = 200)


# 4. MARKER VALIDATION AGAINST THE SUPPLIED ANNOTATION
# Group the fine CosMx annotation into five broad lineages, then test whether
# each single-cell marker is upregulated in its expected lineage spatially.

# Final single-cell marker list (gene + cell_type)
markers <- read.csv("sc_liver_markers_final.csv")

# Mapping of fine CosMx cell types -> five broad lineages
# (e.g. inflammatory + non-inflammatory macrophages -> Myeloid)
mapping <- read.csv("cosmx_celltype_mapping.csv")

# Normalise the RNA assay for plotting
DefaultAssay(sp) <- "RNA"
sp <- NormalizeData(sp)

# Translate each cell's fine type into its broad lineage group
sp$generic_group <- unname(setNames(mapping$generic_group, mapping$cosmx_cellType)[as.character(sp$cellType)])
sp$generic_group[is.na(sp$generic_group)] <- "Unassigned"   # unmapped -> Unassigned

# Keep only the five lineages of interest
keep_groups <- c("T cells","B cells","Plasma cells","Myeloid","Endothelial")
sp <- subset(sp, subset = generic_group %in% keep_groups)

Idents(sp) <- "generic_group"
table(sp$generic_group)

# Dot plot of markers across broad lineages
# (marker_order / avg / avg_mat / markers_on are produced in the validation
#  workflow; ensure they are defined before running these lines)
DotPlot(sp, features = marker_order) + RotatedAxis()
ggsave(file.path(sp_base, "plots/SPATIAL_marker_validation_dotplot.png"),
       width = 22, height = 7, limitsize = FALSE)

# Which markers are / are not on the panel
table(markers$cell_type[markers$gene %in% rownames(avg)])               # on panel
markers[!markers$gene %in% rownames(avg), c("gene","cell_type")]        # not on panel

# For each marker, find which lineage it is highest in on the spatial data
markers_on$top_group <- apply(avg_mat[markers_on$gene, , drop = FALSE], 1,
                              function(x) colnames(avg_mat)[which.max(x)])

# A marker "validates" if its peak lineage matches its expected cell type
markers_on$hits_expected <- markers_on$top_group == markers_on$cell_type
table(markers_on$cell_type, markers_on$hits_expected)
write.csv(markers_on, file.path(sp_base, "SPATIAL_marker_validation_verdict.csv"),
          row.names = FALSE)

# Markers that did NOT validate spatially
markers_on[markers_on$hits_expected == FALSE, c("gene","cell_type","top_group","myAUC")]

# Inspect the two failed markers' expression across lineages
round(avg[c("SPP1","TNFRSF17"), ], 2)


# 5. INDEPENDENT CLUSTERING & ANNOTATION (own analysis)
# Re-process the full spatial object from scratch and annotate it independently,
# rather than relying solely on the supplied reference annotation.

# Reload the full object (all 332,877 cells)
sp_full <- readRDS("/rds/projects/g/gilberts-spatial-biology-image-analysis/data/CosMx/Data/Default/01_ConstructedSeuratObjects_Slide/seuratObject_Slide3_Bruker_Liver.RDS")
ncol(sp_full)   # 332877

DefaultAssay(sp_full) <- "RNA"

# Standard processing pipeline
sp_full <- NormalizeData(sp_full)
sp_full <- FindVariableFeatures(sp_full)
sp_full <- ScaleData(sp_full)
sp_full <- RunPCA(sp_full)
ElbowPlot(sp_full, ndims = 30)

# Build graph and cluster (dims must match RunUMAP below)
sp_full <- FindNeighbors(sp_full, dims = 1:20)
sp_full <- FindClusters(sp_full, resolution = 0.3, algorithm = 1)

# clustree across resolutions
clustree_plot <- clustree(sp_full)
ggsave(file.path(sp_base, "plots", "SPATIAL_res3_clustree.png"),
       clustree_plot, height = 12, width = 8)

# UMAP (same dims as FindNeighbors)
sp_full <- RunUMAP(sp_full, dims = 1:20, min.dist = 0.1, n.neighbors = 30)

# Pin the chosen resolution (0.3)
sp_full$seurat_clusters <- sp_full$RNA_snn_res.0.3
Idents(sp_full) <- "seurat_clusters"
table(sp_full$seurat_clusters)

# UMAP of the chosen clustering
p <- DimPlot(sp_full, reduction = "umap", label = TRUE, raster = TRUE) + NoLegend()
ggsave(file.path(sp_base, "plots", "SPATIAL_FINAL_UMAP_clusters.png"),
       p, width = 10, height = 7, dpi = 300)


# 6. MODULE-SCORE ANNOTATION
# Score each cluster against the HPA lineage marker sets and assign the
# best-matching lineage.

downloaded_markers <- read.csv("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/liver/cell_markers_HPA.csv")
Idents(sp_full) <- "seurat_clusters"

# 1. HPA markers grouped by cell type, keep only genes on the CosMx panel
marker_list <- split(downloaded_markers$gene, downloaded_markers$cell.type)
marker_list <- lapply(marker_list, function(g) intersect(unique(g), rownames(sp_full)))
marker_list <- marker_list[sapply(marker_list, length) > 0]

# 2. Score every cell for each lineage's marker set (ctrl=20 for the ~1k panel)
sp_full <- AddModuleScore(sp_full, features = marker_list, name = "ctscore", ctrl = 20)
score_cols <- paste0("ctscore", seq_along(marker_list))

# 3. Mean score per cluster, z-scored so lineages compete fairly;
#    assign the top-scoring lineage, or "Unassigned/Other" if none is clear
cluster_scores <- sp_full@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean), .groups = "drop")
z <- scale(cluster_scores[score_cols])
cluster_scores$HPA_call <- ifelse(apply(z, 1, max) < 0.5, "Unassigned/Other",
                                  names(marker_list)[apply(z, 1, which.max)])
print(cluster_scores)

# NOTE: module scoring alone is unreliable here (hepatocytes have no HPA signal),
# so cluster labels were assigned by cross-referencing the scores with direct
# marker inspection (cross-tab against the reference cellType).
sp_full$seurat_clusters <- sp_full$RNA_snn_res.0.3
Idents(sp_full) <- "seurat_clusters"

# Manual cluster-to-lineage assignment (from scores + marker inspection).
# Only the lineages of interest are named; all others -> "Unassigned"
# (predominantly hepatocytes; deliberately outside the scope of this study).
new_ids <- c(
  "0" = "Unassigned", "1" = "Unassigned", "2" = "Unassigned",
  "3" = "T cells",    "4" = "Unassigned", "5" = "Myeloid",
  "6" = "Endothelial","7" = "Unassigned", "8" = "Endothelial",
  "9" = "T cells"
)
sp_full$celltype <- unname(new_ids[as.character(sp_full$seurat_clusters)])
sp_full$celltype[is.na(sp_full$celltype)] <- "Unassigned"
sp_full$celltype_final <- as.character(sp_full$celltype)

# B and plasma cells could not be resolved by clustering on the 1k panel, so
# they are placed using the supplied REFERENCE annotation (mixed-method step).
b_ref      <- colnames(sp_full)[sp_full$cellType == "Mature.B.cells"]
plasma_ref <- colnames(sp_full)[sp_full$cellType == "Antibody.secreting.B.cells"]
sp_full$celltype_final[b_ref]      <- "B cells"
sp_full$celltype_final[plasma_ref] <- "Plasma cells"

Idents(sp_full) <- "celltype_final"
table(sp_full$celltype_final)

# Save the annotated object (checkpoint)
saveRDS(sp_full, file.path(sp_base, "SPATIAL_annotated_full.rds"))

# Annotated UMAP
p <- DimPlot(sp_full, reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE)
ggsave(file.path(sp_base, "plots", "SPATIAL_annotated_UMAP_final.png"), p, width = 12, height = 8)


# 7. CELL-TYPE ISOLATION PANELS
# One UMAP per lineage, highlighting that lineage's cells, to show spatial
# distribution and how well each resolves.

# (reload the annotated object if continuing in a fresh session)
readRDS("/rds/projects/g/gilberts-spatial-biology-image-analysis/ABI000_SingleCellMarkerTesting/2026_GAdams/spatial/SPATIAL_annotated_full.rds")

# Shared lineage colour palette
lineage_cols <- c(
  "T cells"      = "#4E79A7",
  "B cells"      = "#762A83",
  "Plasma cells" = "darkslategray",
  "Myeloid"      = "#F28E2B",
  "Endothelial"  = "#E15759",
  "Unassigned"   = "darkgrey"
)

celltypes <- setdiff(names(lineage_cols), "Unassigned")   # skip Unassigned

# Helper: highlight one lineage on the UMAP
make_highlight <- function(obj, ct, reduction = "umap") {
  cells <- colnames(obj)[obj$celltype_final == ct]
  DimPlot(obj, reduction = reduction, cells.highlight = cells,
          cols = "lightgrey", cols.highlight = lineage_cols[[ct]],
          sizes.highlight = 0.8, raster = TRUE) +
    ggtitle(ct) + NoLegend() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank())
}

# Combine all lineage panels into one figure
panels <- lapply(celltypes, function(ct) make_highlight(sp_full, ct))
combined <- wrap_plots(panels, ncol = 3) +
  plot_annotation(title = "Spatial distribution of annotated cell types",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)))
ggsave("SPATIAL_celltype_isolation.png", combined, width = 12, height = 8, dpi = 300)


# 8. CELL-TYPE COMPOSITION BAR PLOT

sp_counts <- as.data.frame(table(sp_full$celltype_final)) %>%
  setNames(c("cell_type","count")) %>%
  arrange(desc(count))
sp_counts$cell_type <- factor(sp_counts$cell_type, levels = sp_counts$cell_type)

p_sp <- ggplot(sp_counts, aes(cell_type, count, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = lineage_cols) +
  scale_y_continuous(breaks = seq(0, 300000, 25000), labels = comma) +
  labs(x = NULL, y = "Cell count",
       title = "Spatial annotation frequency for each cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave("SPATIAL_cellcount_bar.png", p_sp, width = 6, height = 5, dpi = 300)


# 9. IMPROVED ANNOTATED UMAP (ordered, coloured)
# Final publication UMAP: fixed lineage order, rare types drawn on top.

lin_levels <- c("T cells","B cells","Plasma cells","Myeloid","Endothelial","Unassigned")
sp_full$celltype_final <- factor(sp_full$celltype_final, levels = lin_levels)
Idents(sp_full) <- "celltype_final"

p <- DimPlot(sp_full, reduction = "umap", group.by = "celltype_final",
             cols = lineage_cols, pt.size = 0.4, raster = TRUE, raster.dpi = c(600,600),
             order = c("B cells","Plasma cells","T cells","Myeloid","Endothelial","Unassigned")) +
  ggtitle("Spatial liver — annotated cell types") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")
ggsave(file.path(sp_base, "plots", "SPATIAL_annotated_UMAP_improved.png"),
       p, width = 11, height = 8, dpi = 300)


# 10. FINAL LINEAGE MARKER DOT PLOT (validated markers, ordered by lineage)
# A clean, self-contained dot plot of the panel markers across broad lineages.

# Reload and re-group (standalone block)
sp <- readRDS("/rds/projects/g/gilberts-spatial-biology-image-analysis/data/CosMx/Data/Default/01_ConstructedSeuratObjects_Slide/seuratObject_Slide3_Bruker_Liver.RDS")
DefaultAssay(sp) <- "RNA"
sp <- NormalizeData(sp)

mapping <- read.csv("cosmx_celltype_mapping.csv")
sp$generic_group <- unname(setNames(mapping$generic_group, mapping$cosmx_cellType)[as.character(sp$cellType)])
sp$generic_group[is.na(sp$generic_group)] <- "Unassigned"

keep_groups <- c("T cells","B cells","Plasma cells","Myeloid","Endothelial")
sp <- subset(sp, subset = generic_group %in% keep_groups)
Idents(sp) <- "generic_group"
table(Idents(sp))

# Panel markers, grouped by lineage
features <- c(
  # Endothelial
  "ADGRL4","CD34","CDH5","CLEC14A","FLT1","KDR","PECAM1","TIE1",
  # Myeloid
  "C1QB","C1QC","CD14","CD5L","FCGR3A","MARCO",
  # Plasma cells
  "FKBP11","IGHA1","IGHG1","IGHG2","IGKC","JCHAIN","MZB1",
  # T cells
  "CD3D","CD3G","CD8A","IL17A","PDCD1",
  # B cells
  "CD27","MS4A1"
)

# Report which are present / missing on the panel
features_present <- intersect(features, rownames(sp))
features_missing <- setdiff(features, rownames(sp))
cat("Genes present:", length(features_present), "\n"); print(features_present)
cat("\nGenes missing from spatial panel:", length(features_missing), "\n"); print(features_missing)

# Keep lineage order, present genes only
feature_order <- intersect(features, features_present)

# Dot plot
spatial_marker_dotplot <- DotPlot(sp, features = feature_order, group.by = "generic_group") +
  RotatedAxis() +
  ggtitle("CosMx Spatial Marker Validation DotPlot") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))
spatial_marker_dotplot
ggsave(filename = file.path(sp_base, "plots", "SPATIAL_lineage_marker_dotplot.png"),
       plot = spatial_marker_dotplot, width = 16, height = 8, dpi = 300)

# END
