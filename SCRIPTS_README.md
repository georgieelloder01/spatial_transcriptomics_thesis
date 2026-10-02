# Analysis scripts

This repository contains two R (Seurat) scripts that make up the analysis pipeline. They are designed to be run in order: the single-cell script produces the validated marker set that the spatial script then tests.

| Script | Purpose |
|--------|---------|
| `singlecell-data_seurat_analysis.R` | Processes the dissociated single-cell liver data and derives the validated cell-type marker set. |
| `spatial-data_seurat_analysis.R` | Validates that marker set against CosMx spatial liver data and independently annotates the spatial data. |

---

## `singlecell-data_seurat_analysis.R`

Processes dissociated human liver single-cell RNA-seq data (Liver Cell Atlas) and curates a validated set of cell-type markers.

**Sections:**
1. **Load data** — read the 10x count matrix and the Human Protein Atlas candidate markers.
2. **Quality control** — filter cells on gene count and mitochondrial percentage.
3. **Normalisation, feature selection, scaling, PCA** — standard Seurat pre-processing.
4. **Clustering and UMAP** — cluster across resolutions (clustree), embed with UMAP, and select resolution 0.3.
5. **Marker exploration by lineage** — feature and violin plots of candidate markers before filtering.
6. **Differential expression / marker scoring** — score each marker by fold-change, specificity (pct.1/pct.2) and classification performance (AUC).
7. **AUC threshold sensitivity** — how many markers are retained as the AUC cutoff is made stricter.
8. **Combined feature & violin plots** — per-cluster plots of markers passing the AUC cutoff.
9. **Selected marker plots** — curated per-lineage feature/violin plots.
10. **Manual cluster annotation** — assign fine cell types; subcluster to separate B cells from plasma cells.
11. **Validated-marker dot plots** — markers collapsed to five broad lineages (two cross-lineage markers excluded, giving the final set of 42).
12. **Zoomed UMAPs** — annotated myeloid and B/plasma subpopulations.
13. **Composition bar plot** — cell counts per lineage.

**Key output:** `sc_liver_markers_final.csv` — the 42 validated single-cell markers.

---

## `spatial-data_seurat_analysis.R`

Validates the single-cell marker set against CosMx spatial molecular imaging of human liver (1,000-plex panel) and performs an independent annotation of the spatial data.

**Sections:**
1. **Quality control** — report QC metrics of the provided spatial object.
2. **Marker availability** — how many markers are present on the CosMx panel.
3. **UMAP & tissue map** — visualise the supplied annotation and the physical tissue layout.
4. **Marker validation** — group the supplied annotation into five broad lineages and test whether each marker peaks in its expected lineage; record the verdict.
5. **Independent clustering & annotation** — re-process and cluster the full object independently.
6. **Module-score annotation** — assign lineages by scoring clusters against the marker sets, cross-checked with marker inspection; rare B and plasma cells placed from the reference annotation.
7. **Isolation panels** — one UMAP per lineage showing spatial distribution.
8. **Composition bar plot** — cell counts per lineage.
9. **Improved annotated UMAP** — final annotated figure.
10. **Lineage marker dot plot** — panel markers across broad lineages, reporting which genes are present vs missing.

**Key outputs:** `SPATIAL_marker_validation_verdict.csv` (per-marker spatial verdict) and `SPATIAL_annotated_full.rds` (annotated spatial object).

---

## Requirements

- R with Seurat (v5), dplyr, ggplot2, patchwork, clustree, purrr, Polychrome, scales

## Notes

- File paths are specific to the BlueBEAR HPC environment and will need updating to run elsewhere.
- Raw data are not included; the single-cell and Human Protein Atlas data are publicly available from their original sources.
- Some validation objects (`avg`, `avg_mat`, `markers_on`, `marker_order`) are produced within the validation workflow and must be defined before the dependent plotting lines are run.
