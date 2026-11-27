# Catch-up after Exercise 9
# Use this if your limma models are stuck, so you can continue with the
# volcano comparison (9b), the gene names (9c), the gene lists (9d) and the
# before/after PCA. Open AVA_training.Rproj first, then run:
#
#   source("catch_up/after_exercise_9.R")
#
# It first repeats catch_up/after_exercise_3.R. Afterwards you also have:
# metadata_combined$Condition, metadata_combined$Batch, design, design_batch,
# de_unadjusted (Exercise 8) and de_batch (Exercise 9).
# This is one possible solution. Read it and compare it with your own attempt.

source(file.path("catch_up", "after_exercise_3.R"))

# Exercise 8: Tumor minus Normal, ignoring study.
metadata_combined$Condition <- factor(metadata_combined$class,
                                      levels = c("Normal", "Tumor"))
design <- model.matrix(~ Condition, data = metadata_combined)
stopifnot(identical(rownames(design), colnames(expr)))
fit <- eBayes(lmFit(expr, design))
de_unadjusted <- topTable(fit, coef = "ConditionTumor",
                          number = Inf, adjust.method = "BH")

# Exercise 9: the same comparison, with study as a term in the model.
metadata_combined$Batch <- factor(metadata_combined$dataset)
design_batch <- model.matrix(~ Condition + Batch, data = metadata_combined)
stopifnot(identical(rownames(design_batch), colnames(expr)),
          qr(design_batch)$rank == ncol(design_batch))
fit_batch <- eBayes(lmFit(expr, design_batch))
de_batch <- topTable(fit_batch, coef = "ConditionTumor",
                     number = Inf, adjust.method = "BH")

message("Catch-up complete: de_unadjusted and de_batch each have ",
        nrow(de_batch), " genes.")
