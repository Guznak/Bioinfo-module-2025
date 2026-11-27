# Cronus mission: single cells from the crew biopsies
# Open cronus_mission.Rproj first, then run:  source("start_single_cells.R")
# Checks the packages and builds the Seurat object `cells` from raw counts.
# It never installs or updates anything.

library_candidates <- c("Rlib", file.path("..", "Rlib"), file.path("..", "..", "Rlib"),
                        "~/alien-lecture/Rlib", "~/alien-lecture/Rlib/Rlib")  # project, parent, two up, course folder
library_candidates <- library_candidates[dir.exists(library_candidates)]
if (length(library_candidates)) {
  .libPaths(c(normalizePath(library_candidates), .libPaths()))
}
needed <- c("Seurat", "tidyverse", "patchwork")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing course packages: ", paste(missing, collapse = ", "),
                         ". Ask the instructor to complete package setup.")
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
})

sc_dir <- file.path("data", "single_cell")
counts <- readRDS(file.path(sc_dir, "cronus_cell_counts.rds"))
cell_info <- read_csv(file.path(sc_dir, "cronus_cell_metadata.csv"), show_col_types = FALSE)
donors <- read_csv(file.path(sc_dir, "cronus_donors.csv"), show_col_types = FALSE)

# Join the donor table to every cell by donor ID, then keep the cell order.
cell_info <- cell_info %>%
  left_join(donors, by = "donor_id") %>%
  column_to_rownames("cell_id")
stopifnot(identical(rownames(cell_info), colnames(counts)),
          !anyNA(cell_info$group))

cells <- CreateSeuratObject(counts = counts, meta.data = cell_info,
                            project = "Cronus_T3_biopsies")
cells[["percent.mt"]] <- PercentageFeatureSet(cells, pattern = "^MT-")
dir.create("output", showWarnings = FALSE)
rm(counts, cell_info)
message("Ready: ", nrow(cells), " genes; ", ncol(cells), " cells; ",
        length(unique(cells$donor_id)), " crew biopsies. Starting from raw counts.")
