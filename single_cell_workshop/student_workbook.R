# Single-cell workshop: your working script.
# Open single_cell_workshop.Rproj first. From Exercise 2 on, paste the code
# your AI assistant writes from the prompt in single_cell_workshop.html (or the
# code under "Show the code") below each heading, run it, and write your
# answers as comments.

# Course package library (Rlib in Documents/alien-lecture). Run this first.
# if you named the folder differently replace the name ~/ is generally the location of you Dcouments folder
setwd("~/Alien lecture rna-seq")
Rlib <- "~/alien-lecture/Rlib"
if (dir.exists(file.path(Rlib, "Rlib"))) Rlib <- file.path(Rlib, "Rlib")  # unpacked one level deeper
.libPaths(Rlib)

# if loading of libraries fails, restart your R Session and try again (Session -> Restart R)
# keep in mind that all the objects you loaded will then be gone and the scripts will have to run again
library(tidyverse)
library(Seurat)
library(ggplot2)
library(patchwork)



# Exercise 1: inspect the object (source("start_here.R") loads sc)

# Exercise 2: QC plots, compare the two rules, apply the working rule
# Stuck? sc <- readRDS("checkpoints/brain_cells_qc.rds")

# Exercise 3: normalise, variable genes, scale, PCA, elbow plot

# Exercise 4: neighbours, clusters (0.4 and 1.0), UMAP, colour by sample
# Stuck? sc <- readRDS("checkpoints/brain_cells_umap.rds")

# Exercise 5: dot plot, top-5 genes per cluster (output/top5_markers.csv),
# then name every cluster ("Unassigned" is allowed)

# Exercise 6: which samples do the clusters come from?

# Exercise 7: save the figure, then the four statements for Ava
# QC decision:
# Resolution effect:
# One cell type and its markers:
# One cluster you could not name, and the evidence you would need:
