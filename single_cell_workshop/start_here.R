library_candidates <- c("Rlib", file.path("..", "Rlib"), file.path("..", "..", "Rlib"),
                        "~/alien-lecture/Rlib", "~/alien-lecture/Rlib/Rlib")  # project, parent, two up, course folder
library_candidates <- library_candidates[dir.exists(library_candidates)]
if (length(library_candidates)) {
  .libPaths(c(normalizePath(library_candidates), .libPaths()))
}
needed <- c("Seurat", "ggplot2", "patchwork")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing course packages: ", paste(missing, collapse = ", "),
                         ". Ask the instructor to complete package setup.")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
sc <- readRDS("data/brain_cells_start.rds")
stopifnot(inherits(sc, "Seurat"), identical(colnames(sc), rownames(sc[[]])))
dir.create("output", showWarnings = FALSE)
message("Ready: ", nrow(sc), " genes; ", ncol(sc), " cells; ",
        length(unique(sc$sample_id)), " samples. Starting from counts.")
