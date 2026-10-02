Evaluating the Reliability of Cell-Type Biomarkers in Spatial Transcriptomics

MSc Bioinformatics thesis project assessing whether cell-type marker genes derived from dissociated single-cell RNA sequencing retain their specificity when applied to spatial transcriptomics data.

Overview

Cell-type identity in transcriptomic data is assigned using marker genes, most of which are derived from dissociated single-cell RNA sequencing. These markers are increasingly applied to spatial transcriptomics, where cells remain in intact tissue and only a restricted, pre-selected gene panel is measured. This project tests whether such markers transfer reliably across platforms, and distinguishes genuine marker failure from platform-driven limitations.

Candidate markers for five lineages — T cells, B cells, plasma cells, myeloid cells and endothelial cells — were curated from the Human Protein Atlas, filtered in human liver single-cell data, and validated against CosMx spatial molecular imaging of human liver. A second, higher-plex (6,000-gene) spatial dataset was used to separate the causes of marker failure.

Key findings
Of 42 markers retained in the dissociated data, 19 were present on the 1,000-plex CosMx panel; 17 of these validated in both modalities (~90% of testable markers).
Panel coverage, not marker quality, was the dominant constraint on transferability.
Two markers failed spatially for distinct reasons: TNFRSF17 (low transcript abundance and reduced detection sensitivity) and SPP1 (a resolution mismatch between the myeloid subpopulation it marks and the broad lineage level achievable spatially).
A higher-plex (6,000-gene) panel recovered rare populations (B and plasma cells) that were unresolvable on the smaller panel, confirming that coverage was the limiting factor.

Methods
Marker curation: candidate markers curated from the Human Protein Atlas and filtered in single-cell data using fold-change, expression specificity (pct.1/pct.2), classification performance (AUC), and visual inspection of feature and violin plots.
Data processing: normalisation, clustering and cell-type annotation of single-cell and spatial data using the Seurat framework.
Cross-platform validation: testing whether each marker was specifically upregulated in its expected lineage spatially.
Failure diagnosis: retesting failed markers on a second, higher-plex dataset to distinguish coverage, resolution, and intrinsic-biology causes.
Platform comparison: genome-wide comparison of detection (pct.1), background (pct.2) and fold-change between platforms.

Data
Single-cell: Liver Cell Atlas (Guilliams et al., 2022), human liver.
Spatial (primary): CosMx 1,000-plex human liver (Bruker Spatial Biology).
Spatial (secondary): CosMx 6,000-plex human liver (in-house).
Markers: Human Protein Atlas Single Cell Type Atlas (Karlsson et al., 2021).

Tools
R / Seurat
High-performance computing (BlueBEAR, University of Birmingham)
Repository structure
├── scripts/        # analysis scripts (processing, annotation, validation)
├── data/           # input marker sets and comparison tables
├── plots/          # figures (UMAPs, dot plots, correlation plots)
└── results/        # output tables (marker verdicts, failed-marker stats)
