# Open AVA_training.Rproj first. This checks; it never installs automatically.
library_candidates <- c("Rlib", file.path("..", "Rlib"), file.path("..", "..", "Rlib"),
                        "~/alien-lecture/Rlib", "~/alien-lecture/Rlib/Rlib")  # project, parent, two up, course folder
library_candidates <- library_candidates[dir.exists(library_candidates)]
if (length(library_candidates)) {
  .libPaths(c(normalizePath(library_candidates), .libPaths()))
}
required <- c("tidyverse", "limma", "knitr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing course packages: ", paste(missing, collapse = ", "),
       ". Ask your instructor to complete setup before class.")
}
stopifnot(file.exists("data/bulk/GSE64912.csv"),
          file.exists("data/bulk/GSE68799.csv"),
          file.exists("data/ml/training_archive.csv"))
message("Core setup is ready. Optional single-cell package available: ",
        requireNamespace("Seurat", quietly = TRUE))
