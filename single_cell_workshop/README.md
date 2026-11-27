# Single-cell workshop

Open single_cell_workshop.html for the exercises, adapted from Mert Cihan's
workshop using real GSE248995 brain data. One route, about 60 minutes after
setup.

## Start

1. Fully extract this folder beside AVA_training.
2. Open single_cell_workshop.Rproj.
3. Open student_workbook.R. From Exercise 2 on, let an AI assistant write the
   code from the prompt in the HTML; the supplied code is under "Show the code".

start_here.R finds Rlib inside the project, beside it, one folder further up,
or in Documents/alien-lecture.
It checks Seurat, ggplot2, and patchwork without installing packages.

## Files

- data/brain_cells_start.rds: 1,344 cells from 15 samples; original counts retained.
- data/sample_manifest.csv: original and teaching sample sizes.
- data/marker_reference.csv: short marker panel for naming clusters.
- checkpoints/brain_cells_qc.rds: resume after QC (Exercise 2).
- checkpoints/brain_cells_umap.rds: resume after clustering and UMAP (Exercise 4).
- output/cluster_markers.csv: fallback if marker computation stalls (Exercise 5).

The original count object contains 10,834 cells. The teaching subset was sampled
reproducibly before clustering. Counts and gene identifiers are unchanged.

## Solutions and editing

The instructor bundle additionally contains single_cell_workshop_solutions.Rmd,
single_cell_workshop_solutions.html, and reference_analysis.R.
The solution Rmd runs the analysis when knitted and shows all plots and answer
sketches. Its example names go to output/example_cell_types.csv; students save
theirs as output/my_cell_types.csv.

Student and solution Rmds are separate documents. Use Knit to rebuild the
corresponding HTML, and update both sources when changing shared text.

Tested locally with R 4.4.3 and Seurat 5.4.0. The classroom package bundle is
not included in this archive.
