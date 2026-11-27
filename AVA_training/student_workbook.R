# AVA / Recovery training
# Open AVA_training.Rproj. Save your work in this script.
# The HTML contains explanations, examples, AI prompts, and a "Done when"
# check for each exercise. Worked solutions are in a separate teaching version.
# Work in order. Run small blocks with Ctrl+Enter.


# Course package library (Rlib in Documents/alien-lecture). Run this first.
setwd("~/alien lecture")
Rlib <- "~/alien-lecture/Rlib"
if (dir.exists(file.path(Rlib, "Rlib"))) Rlib <- file.path(Rlib, "Rlib")  # unpacked one level deeper
.libPaths(Rlib)

# if loading of libraries fails, restart your R Session and try again (Session -> Restart R)
# keep in mind that all the objects you loaded will then be gone and the scripts will have to run again
library(tidyverse)
library(limma)
out.dir <- "output"
in.dir <- file.path("data", "bulk")
dir.create(out.dir, showWarnings = FALSE)

setwd("~/alien lecture")
Rlib <- "~/alien-lecture/Rlib"
if (dir.exists(file.path(Rlib, "Rlib"))) Rlib <- file.path(Rlib, "Rlib")  # unpacked one level deeper
.libPaths(Rlib)


# For each exercise: question -> code -> check -> interpretation.
# Paste only code you can connect to the objects already created here.
# Stuck? source("catch_up/after_exercise_3.R") or
#        source("catch_up/after_exercise_9.R") recreate the objects you need.

# Exercise 0: read, predict, then run (mini table)

# Exercise 1: working directory and input files

# Exercise 2: import and inspect the archive

# Exercise 3: metadata_combined and expr (genes x samples)

# Exercise 4: select a study without losing sample alignment

# Exercise 5: PCA (prcomp needs samples x genes)

# Exercise 6: expression and PCA plots; save a figure

# Exercise 7: variance explained

# Exercise 8: exploratory limma comparison (create Condition first)

# Exercise 9: batch-aware model
# Exercise 9b: run the supplied two-volcano comparison; predict first
# Exercise 9c: join gene names; plot one gene; check the AI's description
# Exercise 9d: four top-20 lists (logFC or t, with or without study); compare
# Before/after PCA: run the supplied code

# Exercise 10: train a logistic regression on independent archive people

# Exercise 11: evaluate held-out predictions against a baseline

# Exercise 12: twelve programs instead of three; which programs carry it?

# M06 and settlement: run the supplied code

# Exercise 13: treatment thresholds 0.3, 0.5, 0.7

# Final transmission: write your evidence and uncertainty here.
# Single-cell workshop: open the separately supplied single_cell_workshop folder.
# Work in the student_workbook.R there.
