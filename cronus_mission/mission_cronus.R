# ═══════════════════════════════════════════════════════════════════════════════
#
#     U S C S S   C R O N U S
#     Biological Triage Mission — Field Analysis Script
#
#     Weyland-Yutani Corporation  //  BioSec Field Team
#     Stardate 2183.201
#
# ═══════════════════════════════════════════════════════════════════════════════
#
#  MISSION BRIEFING
#  ─────────────────
#  The USCSS Cronus has been drifting for decades, its emergency beacon still repeating.
#  You are the field analysis team sent to assess the biological situation.
#
#  Your shuttle docking triggered the ship's power restoration sequence.
#  You did not plan for this.
#
#  The cryo bay has begun cycling. The crew is defrosting.
#  Ava, the ship's synthetic who ran your training over the beacon,
#  is coming fully back online.
#  Lab 2 is intact. Someone ran sequencing data throughout the incident.
#
#  Your task: analyse the data. Understand what happened.
#  Then decide what to do about it.
#
#  Each section of this script unlocks the next piece of the story.
#  Mission Control will release clearance codes as you complete each step.
#  Type each code where you see "__________".
#
# ═══════════════════════════════════════════════════════════════════════════════


# ── SETUP ─────────────────────────────────────────────────────────────────────
# Open cronus_mission.Rproj first, so that all paths below work.
# This block uses the course library (Rlib) and never installs anything.

# Course package library (Rlib in Documents/alien-lecture).
Rlib <- "~/alien-lecture/Rlib"
if (dir.exists(file.path(Rlib, "Rlib"))) Rlib <- file.path(Rlib, "Rlib")  # unpacked one level deeper
.libPaths(Rlib)

library_candidates <- c("Rlib", file.path("..", "Rlib"), file.path("..", "..", "Rlib"),
                        "~/alien-lecture/Rlib", "~/alien-lecture/Rlib/Rlib")  # project, parent, two up, course folder
library_candidates <- library_candidates[dir.exists(library_candidates)]
if (length(library_candidates)) .libPaths(c(normalizePath(library_candidates), .libPaths()))

needed  <- c("tidyverse", "patchwork", "ggrepel", "limma", "RColorBrewer")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing course packages: ", paste(missing, collapse = ", "),
                          ". Ask the instructor.")

library(tidyverse)     # includes purrr (map_dfr) and dplyr
library(patchwork)
library(ggrepel)
library(limma)
library(RColorBrewer)

source("cronus_utils.R")   # decrypt_file() and decrypt_table()
dir.create("output", showWarnings = FALSE)


# ═══════════════════════════════════════════════════════════════════════════════
#  ACT 1 — FIRST CONTACT
#  The ship is live. Load the data. Get your bearings.
# ═══════════════════════════════════════════════════════════════════════════════

# Load the ship manifest — sample IDs, timepoints, and batch labels.
# Status information has been redacted pending clearance.
# Sample IDs are crew names plus the timepoint: AARON_T0 ... AARON_T3.

manifest <- read_csv("data/bulk/USCSS_metadata_partial.csv", show_col_types = FALSE)
glimpse(manifest)


# Load the sequencing count data — all timepoints combined.

counts_all <- read_csv("data/bulk/USCSS_Cronus_counts.csv", show_col_types = FALSE)

gene_cols <- setdiff(colnames(counts_all), c("sample_id", "batch", "timepoint"))
cat(sprintf("Loaded %d samples x %d genes\n", nrow(counts_all), length(gene_cols)))


# ── QUESTION 1 ────────────────────────────────────────────────────────────────
# How are samples distributed across timepoints and batches?
# Are all batches represented at all timepoints?

table(manifest$timepoint, manifest$batch)

# Notice anything about which batches cover which timepoints?
# Why might that be a problem for analysis?
# Is this good or bad experiment design?


# ── Build the log-expression matrix ───────────────────────────────────────────

mat_all_raw <- counts_all |>
  select(sample_id, all_of(gene_cols)) |>
  column_to_rownames("sample_id") |>
  as.matrix() |>
  t()   # genes as rows, samples as columns

mat_all_log <- log2(mat_all_raw + 1)

cat(sprintf("Expression matrix: %d genes x %d samples\n",
            nrow(mat_all_log), ncol(mat_all_log)))


# ═══════════════════════════════════════════════════════════════════════════════
#  ACT 2 — RAW DATA EXPLORATION
#  What does the data look like before any correction?
# ═══════════════════════════════════════════════════════════════════════════════

# create pca plots of the crew and color them by batch
# What is that good for?

run_pca_raw <- function(mat, meta_df, title = "") {
  top_genes <- names(sort(apply(mat, 1, var), decreasing = TRUE))[1:500]
  pca_fit   <- prcomp(t(mat[top_genes, ]), scale. = FALSE, center = TRUE)
  pct       <- round(pca_fit$sdev^2 / sum(pca_fit$sdev^2) * 100, 1)
  coords    <- as_tibble(pca_fit$x[, 1:2], rownames = "sample_id") |>
    left_join(meta_df, by = "sample_id")

  ggplot(coords, aes(PC1, PC2, colour = batch)) +
    geom_point(size = 3, alpha = 0.85) +
    labs(title = title,
         x = paste0("PC1 (", pct[1], "%)"),
         y = paste0("PC2 (", pct[2], "%)")) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10))
}

p_raw <- list()
for (tp in 0:3) {
  tp_ids  <- manifest |> filter(timepoint == tp) |> pull(sample_id)
  tp_meta <- manifest |> filter(sample_id %in% tp_ids)
  p_raw[[tp + 1]] <- run_pca_raw(mat_all_log[, tp_ids], tp_meta,
                                  title = paste0("T", tp, " — raw"))
}

wrap_plots(p_raw, ncol = 2) +
  plot_annotation(title = "Raw counts — coloured by batch")


# ── QUESTION 2 ────────────────────────────────────────────────────────────────
# What is driving the separation in the raw PCA plots?
# Is it biology, or something else? Do we even know?
# What does this mean for our ability to compare samples across timepoints?


# ═══════════════════════════════════════════════════════════════════════════════
#  ACT 3 — BATCH CORRECTION
#  Fix the technical artefacts so the biology can speak.
# ═══════════════════════════════════════════════════════════════════════════════

# Batch effects are systematic technical differences between sequencing runs.
# They are not biology. We remove them before drawing PCA plots.
#
# limma::removeBatchEffect() is the function you used on Day 1
# ("Visualise the adjustment"). Each timepoint was sequenced in two batches,
# so we correct within each timepoint: for every gene, the difference between
# the two batch means at that timepoint is subtracted.

mat_all_corrected <- mat_all_log
for (tp in 0:3) {
  tp_ids   <- manifest |> filter(timepoint == tp) |> pull(sample_id)
  tp_batch <- manifest$batch[match(tp_ids, manifest$sample_id)]
  mat_all_corrected[, tp_ids] <- removeBatchEffect(mat_all_log[, tp_ids],
                                                   batch = tp_batch)
}
dim(mat_all_corrected)   # same size as mat_all_log


# ── QUESTION 3 ────────────────────────────────────────────────────────────────
# What does removeBatchEffect() do, in principle?
# Why do we correct in log space rather than on raw counts?
# On Day 1 you gave removeBatchEffect() a design so that it kept the
# Tumor/Normal difference. Here we give it none — we do not know the groups yet.
# If one batch held mostly one kind of crew member, what else would the
# correction remove? Keep this in mind when you get the crew records (Act 5).


# ═══════════════════════════════════════════════════════════════════════════════
#  ACT 4 — CORRECTED PCA  /  BIOLOGICAL EXPLORATION
#  The artefacts are gone. What is actually in this data?
# ═══════════════════════════════════════════════════════════════════════════════

# PCA helper — uses the 500 most variable genes at that timepoint, like Act 2.
# Returns a list with coords (tibble) and pct (variance explained).

do_pca_stable <- function(mat, meta_df, n_top = 500) {
  top_genes <- names(sort(apply(mat, 1, var), decreasing = TRUE))[1:n_top]

  mat_top   <- mat[top_genes, ]
  mat_c     <- sweep(mat_top, 1, rowMeans(mat_top), "-")
  pca_fit   <- prcomp(t(mat_c), scale. = FALSE, center = FALSE)
  pct       <- round(pca_fit$sdev^2 / sum(pca_fit$sdev^2) * 100, 1)
  projected <- t(mat_c) %*% pca_fit$rotation

  coords <- as_tibble(projected[, 1:2], rownames = "sample_id") |>
    left_join(meta_df, by = "sample_id")

  list(coords = coords, pct = pct)
}


# ── Corrected PCA per timepoint — coloured by batch ──────────────────────────
# First confirm the correction worked: batches should no longer dominate.

p_corr <- list()
for (tp in 0:3) {
  tp_ids  <- manifest |> filter(timepoint == tp) |> pull(sample_id)
  tp_meta <- manifest |> filter(sample_id %in% tp_ids)
  res     <- do_pca_stable(mat_all_corrected[, tp_ids], tp_meta)

  p_corr[[tp + 1]] <- ggplot(res$coords, aes(PC1, PC2, colour = batch)) +
    geom_point(size = 3, alpha = 0.85) +
    labs(title = paste0("T", tp, " — corrected"),
         x = paste0("PC1 (", res$pct[1], "%)"),
         y = paste0("PC2 (", res$pct[2], "%)")) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10))
}

wrap_plots(p_corr, ncol = 2) +
  plot_annotation(title = "After correction — do batches still separate?")


# ── QUESTION 4 ────────────────────────────────────────────────────────────────
# After correction, are batches still the dominant source of variation?
# What is separating the samples now?



# ── CLEARANCE CODE 1 ──────────────────────────────────────────────────────────
#
#  The ship's synthetic is transmitting on a looped relay.
#  She has been waiting — possibly for years — for someone to arrive.
#
#  Report your completed batch correction to Mission Control.
#  You will receive the decryption key for her message.
#
decrypt_file("transmissions/ava_message_1.enc", "__________")


# ── T3: unlabelled PCA ───────────────────────────────────────────────────────
# Ava said the sequencing data is the most complete record of what happened.
# Look at the final timepoint. This is where things escalated, and the crew
# went into cryosleep.

tp3_ids  <- manifest |> filter(timepoint == 3) |> pull(sample_id)
tp3_mat  <- mat_all_corrected[, tp3_ids]
tp3_meta <- manifest |> filter(sample_id %in% tp3_ids)
res_t3   <- do_pca_stable(tp3_mat, tp3_meta)

ggplot(res_t3$coords, aes(PC1, PC2)) +
  geom_point(size = 3.5, alpha = 0.85, colour = "#888888") +
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 25) +
  labs(title = "T3 — incident timepoint — can you name the clusters?",
       x = paste0("PC1 (", res_t3$pct[1], "%)"),
       y = paste0("PC2 (", res_t3$pct[2], "%)")) +
  theme_bw(base_size = 11)

# The same positions as a table, sorted along PC1:
res_t3$coords |> select(sample_id, PC1, PC2) |> arrange(PC1) |> print(n = Inf)


# ── QUESTION 5 ────────────────────────────────────────────────────────────────
# How many distinct groups can you identify in the T3 PCA? Is this clear or not?
# How many or how few clusters could there potentially be?
# For each group, write down:
#   - Approximate number of samples
#   - Position on PC1 and PC2 (left / right / top / bottom)
#   - Your hypothesis about who these people might be
#
# Pay particular attention to any sample that sits completely alone.


# ── Label your clusters ───────────────────────────────────────────────────────
# Put every T3 sample into one group, here in R. Use any group names you like —
# describe what you think each group is. Add or remove groups as needed.
# Copy the sample IDs from the table above.
#
# This is manual work — for 22 samples it is manageable, but what would you
# do with thousands of samples? In real studies, clustering is automated:
# k-means or hierarchical clustering assigns samples to groups by distance,
# or graph-based methods (like Seurat's Louvain clustering in single-cell work)
# find communities in a nearest-neighbour network. Here, the groups are clear
# enough to label by eye — and doing it manually means you actually have to
# look at the data.

my_groups <- list(
  group_1 = c("__________", "__________"),
  group_2 = c("__________"),
  group_3 = c("__________")
)

# Turn the list into a table and check it
t3_labels <- tibble(sample_id = unlist(my_groups, use.names = FALSE),
                    cluster   = rep(names(my_groups), lengths(my_groups))) |>
  filter(sample_id != "__________")
not_found  <- setdiff(t3_labels$sample_id, tp3_ids)
twice      <- t3_labels$sample_id[duplicated(t3_labels$sample_id)]
not_placed <- setdiff(tp3_ids, t3_labels$sample_id)
if (length(not_found))  message("Not a T3 sample ID (typo?): ", paste(not_found, collapse = ", "))
if (length(twice))      message("In more than one group: ", paste(twice, collapse = ", "))
if (length(not_placed)) message("Not in a group yet: ", paste(not_placed, collapse = ", "))

manifest <- manifest |>
  select(-any_of("cluster")) |>
  left_join(t3_labels |> filter(sample_id %in% tp3_ids) |> distinct(sample_id, .keep_all = TRUE),
            by = "sample_id")
manifest$cluster[manifest$timepoint == 3 & is.na(manifest$cluster)] <- "unknown"

# Colour palette: your groups get Set2 colours, "unknown" and "no T3 sample" grey
my_group_names <- sort(setdiff(unique(t3_labels$cluster), "unknown"))
group_palette  <- setNames(
  c(colorRampPalette(brewer.pal(8, "Set2"))(max(8, length(my_group_names)))[seq_along(my_group_names)],
    "#cccccc", "grey70"),
  c(my_group_names, "unknown", "no T3 sample")
)


# ── Plot T3 with your labels ──────────────────────────────────────────────────

res_t3$coords |>
  select(-any_of("cluster")) |>
  left_join(manifest |> filter(timepoint == 3) |> select(sample_id, cluster),
            by = "sample_id") |>
  ggplot(aes(PC1, PC2, colour = cluster)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 25) +
  scale_colour_manual(values = group_palette) +
  labs(title = "T3 — your interpretation",
       x = paste0("PC1 (", res_t3$pct[1], "%)"),
       y = paste0("PC2 (", res_t3$pct[2], "%)")) +
  theme_bw(base_size = 11) +
  theme(legend.title = element_blank())


# ── Trajectory: T0 → T3 ──────────────────────────────────────────────────────
# Your T3 labels are carried back to the same crew member at T0, T1 and T2:
# the crew name is the sample ID without "_T0" ... "_T3".
# Crew with no T3 sample appear grey — they did not make it to the end.

t3_cluster_map <- manifest |>
  filter(timepoint == 3) |>
  transmute(crew = str_remove(sample_id, "_T\\d+$"),
            propagated_cluster = cluster)

all_tp_res <- list()
for (tp in 0:3) {
  tp_ids  <- manifest |> filter(timepoint == tp) |> pull(sample_id)
  tp_meta <- manifest |> filter(sample_id %in% tp_ids)

  res_tp <- do_pca_stable(mat_all_corrected[, tp_ids], tp_meta)

  coords <- res_tp$coords |>
    mutate(crew = str_remove(sample_id, "_T\\d+$")) |>
    left_join(t3_cluster_map, by = "crew") |>
    mutate(propagated_cluster = replace_na(propagated_cluster, "no T3 sample"))

  all_tp_res[[tp + 1]] <- list(coords = coords, pct = res_tp$pct, tp = tp)
}

all_coords <- bind_rows(lapply(all_tp_res, `[[`, "coords"))
xlim_range <- range(all_coords$PC1, na.rm = TRUE)
ylim_range <- range(all_coords$PC2, na.rm = TRUE)
xpad <- diff(xlim_range) * 0.05
ypad <- diff(ylim_range) * 0.05
used_levels <- intersect(names(group_palette), all_coords$propagated_cluster)

tp_trajectory <- lapply(all_tp_res, function(x) {
  ggplot(x$coords, aes(PC1, PC2, colour = propagated_cluster)) +
    geom_point(size = 3, alpha = 0.85, show.legend = TRUE) +
    geom_text_repel(aes(label = sample_id), size = 2.3,
                    max.overlaps = 15, show.legend = FALSE) +
    scale_colour_manual(values = group_palette, limits = used_levels) +
    coord_cartesian(xlim = xlim_range + c(-xpad, xpad),
                    ylim = ylim_range + c(-ypad, ypad)) +
    labs(title = paste0("T", x$tp),
         x = paste0("PC1 (", x$pct[1], "%)"),
         y = paste0("PC2 (", x$pct[2], "%)")) +
    theme_bw(base_size = 11) +
    theme(legend.title = element_blank(),
          plot.title = element_text(face = "bold", size = 10))
})

p <- wrap_plots(tp_trajectory, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Your interpretation — traced across all timepoints",
    subtitle = "Colours are carried back from your T3 labels. Grey = no T3 sample."
  ) &
  theme(legend.position = "bottom")

ggsave("output/trajectory_timepoints.png", plot = p, width = 14, height = 10, dpi = 300)

p

# ── QUESTION 6 ────────────────────────────────────────────────────────────────

# Look at the trajectory across T0 → T3.
#
# 1. At T0, how many distinct groups are visible?
#    Are all your T3 groups already separating from the start?
#
# 2. Which group appears to SPLIT between T1 and T3?
#    What biological event might cause a group to divide into two?
#
# 3. Notice the grey points, they are crew with no T3 sample.
#    At which timepoint do most of them disappear?
#    What might have happened to them?


# You find a handwritten note lying on a desk in the Lab, signed by Dr. Elliot Cooper:

#  ┌─ COOPER'S LAB NOTES SAY: ─────────────────────────────────────────────────┐
#  │                                                                           │
#  │  The antidote does something. It is not what was described                │
#  │  in the company documentation. What are they hiding?                      │
#  │                                                                           │
#  │  OUTLIER_C — instrument error. Do not prioritise this sample.             │
#  │                                                                           │
#  │                                                                           │
#  │  "HOPE?  One sample present at all four timepoints that does not fit      │
#  │   any of the above categories. I have not highlighted this crew member's  │
#  │   identity because I am not sure how to proceed yet."                     │
#  │                                                                           │
#  └───────────────────────────────────────────────────────────────────────────┘

# ── Two samples that don't fit ────────────────────────────────────────────────
# Two samples in this dataset behave differently from all others.
# One is easy to spot: it sits completely alone at T3, far from every cluster.
# The other is harder: it looks like it belongs to a group,
# but its trajectory across timepoints tells a different story.
#
# let's look at what we know:
# - some people seem to have transformed into horrifying creatures
# - they should separate from healthy individuals in T3 in the PCA
# - we should look whether anyone has a different trajectory than
#   healthy --> infected --> transformed
#
# Tip: take your healthiest-looking T3 group and follow its members back to T1
# in the trajectory plot. Does every one of them sit with the healthy crew there?


# ── QUESTION 7 ────────────────────────────────────────────────────────────────
# For each of the two anomalous samples:
#   - Which crew member do you think it is, and why?
#   - At which timepoint does it first look different from the others?
#   - Based on position and trajectory alone — before reading any logs —
#     what is your hypothesis about what happened to each one?


# ── CLEARANCE CODE 2 ──────────────────────────────────────────────────────────
#
#  Somewhere in Lab 2, Cooper left notes on what he was seeing in the data.
#  He locked them, together with his records of every crew member's status.
#
#  Report your T3 groups to Mission Control: how many you see, which samples
#  belong together, and which two crew members you think do not fit.
#
decrypt_file("transmissions/medical_log.enc", "__________")
#
#
#
#  [ The crew are continuing to defrost. Some are becoming responsive. Their
#    cryosleep sickness is weakening, other symptoms are showing...
#    - aggressive behaviour
#    - shouting
#    - loss of orientation
#    - Dr. Cooper is lying on the floor periodically cramping ]
#


# ═══════════════════════════════════════════════════════════════════════════════
#  ACT 5 — DIFFERENTIAL EXPRESSION
#  What genes are driving the separation you see?
# ═══════════════════════════════════════════════════════════════════════════════

# Cooper's status records came with his notes. Open them with the same code.
# Compare them with your own groups: where do you agree, where not?

crew_records <- decrypt_table("transmissions/crew_records.enc", "__________")
count(crew_records, timepoint, status) |> pivot_wider(names_from = status, values_from = n)

tp3_meta_full <- crew_records |>
  filter(timepoint == 3) |>
  column_to_rownames("sample_id")


# ── DEG: Creature vs Control ──────────────────────────────────────────────────

run_deg <- function(expr_mat, meta_df, group1, group2) {
  sub_meta <- meta_df |> filter(status %in% c(group1, group2))
  sub_mat  <- expr_mat[, rownames(sub_meta), drop = FALSE]
  design   <- model.matrix(~ status, data = sub_meta)
  fit      <- eBayes(lmFit(sub_mat, design))
  topTable(fit, coef = 2, number = Inf, adjust.method = "BH") |>
    as_tibble(rownames = "gene_id") |>
    arrange(adj.P.Val)
}

deg_creature <- run_deg(tp3_mat, tp3_meta_full, "control", "creature")

deg_creature |>
  filter(adj.P.Val < 0.05, abs(logFC) > 1) |>
  arrange(desc(logFC)) |>
  print(n = Inf)


# ── QUESTION 8 ────────────────────────────────────────────────────────────────
# Which HVAR genes are most strongly up in the creatures, and which are down?
# Group them by their numbers (HVAR_001 ... HVAR_050).
# Does the grouping suggest anything about how the genome is being reorganised?


# ── Volcano plot ──────────────────────────────────────────────────────────────

deg_volcano <- deg_creature |>
  mutate(
    neg_logp = -log10(pmax(adj.P.Val, 1e-50)),
    group    = case_when(
      adj.P.Val < 0.05 & logFC >  1 ~ "up",
      adj.P.Val < 0.05 & logFC < -1 ~ "down",
      TRUE                           ~ "ns"
    )
  )

ggplot(deg_volcano, aes(logFC, neg_logp, colour = group)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#888888") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", colour = "#888888") +
  scale_colour_manual(values = c(up = "#D65F5F", down = "#4878CF", ns = "#cccccc")) +
  geom_text_repel(data = filter(deg_volcano, group != "ns"),
                  aes(label = gene_id), size = 2.5, max.overlaps = Inf,
                  force = 0.5, box.padding = 0.3) +
  labs(title = "Creature vs Control — Volcano plot",
       x = "log2 Fold Change", y = "-log10(adj. p-value)") +
  theme_bw(base_size = 11)


# ── QUESTION 9 ────────────────────────────────────────────────────────────────
# The volcano shows two very clean clusters of significant genes —
# one strongly up, one strongly down.
# This is unusually clean for real RNA-seq data.
#
# Think about what was different between the creature group and the controls:
# the creatures received a compound — the antidote.
#
#   - Could a drug or compound trigger this kind of coordinated response?
#   - What would it mean if one gene module was activated and another
#     suppressed simultaneously by the same intervention?
#   - Does this look more like a disease progressing, or like a response
#     to something that was administered?


# ── CLEARANCE CODE 3 ──────────────────────────────────────────────────────────
#
#  Weyland-Yutani has been monitoring your uplink since Ava came online.
#  They have been expecting your findings.
#
#  Report the top biological process driving the creature transformation
#  to Mission Control in plain English.
#  If correct, you will receive their first direct communication.
#
decrypt_file("transmissions/wy_transmission.enc", "__________")
#


# ═══════════════════════════════════════════════════════════════════════════════
#  ACT 6 — THE ANOMALIES
#  Two samples don't fit any of the groups you've seen.
#  Weyland-Yutani is waiting for your asset confirmation.
# ═══════════════════════════════════════════════════════════════════════════════

# One crew member has labelled himself as an outlier throughout the whole
# analysis... who might that be? The other one you found in the trajectory plot.
# Write both crew names (without _T0 ... _T3):

anomaly_a <- "__________"   # Anomaly A: the sample that sits alone at T3
anomaly_b <- "__________"   # Anomaly B: the one whose trajectory does not fit

crew_names <- unique(str_remove(manifest$sample_id, "_T\\d+$"))
stopifnot("Write two crew names from the sample IDs, without _T3" =
            all(c(anomaly_a, anomaly_b) %in% crew_names))


# ── HVAR expression trajectory for each anomaly ───────────────────────────────
#
# The volcano plot (Act 5) showed two clean clusters of significant HVAR genes:
#   - HVAR_011–020: strongly upregulated in creatures
#   - HVAR_021–040: strongly downregulated in creatures
#
# We use those same groupings here to colour the trajectory lines, plus the
# shared immune response group (HVAR_001–010) and HVAR_041–050, which no
# group so far has explained.
#
# Each line shows the mean log2 expression of all genes in that group.
# This connects what you saw in the volcano to how these individuals changed
# over time — before you know the pathway names.

hvar_groups <- list(
  "HVAR_001-010 (innate immune / shared response)" = sprintf("HVAR_%03d", 1:10),
  "HVAR_011-020 (up in creatures)"                 = sprintf("HVAR_%03d", 11:20),
  "HVAR_021-040 (down in creatures)"               = sprintf("HVAR_%03d", 21:40),
  "HVAR_041-050 (not explained yet)"               = sprintf("HVAR_%03d", 41:50)
)

group_colours <- c(
  "HVAR_001-010 (innate immune / shared response)" = "#aaaacc",
  "HVAR_011-020 (up in creatures)"                 = "#D65F5F",
  "HVAR_021-040 (down in creatures)"               = "#4878CF",
  "HVAR_041-050 (not explained yet)"               = "#FFD700"
)

# Mean log2(count + 1) of each HVAR group, for one crew member at every timepoint
plot_trajectory <- function(crew, title) {
  ids <- manifest |>
    filter(str_remove(sample_id, "_T\\d+$") == crew) |>
    arrange(timepoint)

  map_dfr(names(hvar_groups), function(grp) {
    tibble(timepoint = ids$timepoint,
           group     = grp,
           mean_log2 = colMeans(mat_all_log[hvar_groups[[grp]], ids$sample_id, drop = FALSE]))
  }) |>
    ggplot(aes(timepoint, mean_log2, colour = group, group = group)) +
    geom_line(linewidth = 1.4) +
    geom_point(size = 3) +
    scale_colour_manual(values = group_colours, name = "HVAR group") +
    guides(colour = guide_legend(nrow = 2)) +
    scale_x_continuous(breaks = ids$timepoint) +
    labs(title = title,
         x = "Timepoint", y = "Mean log2(count + 1)") +
    theme_bw(base_size = 11) +
    theme(legend.position  = "bottom",
          legend.text      = element_text(size = 8),
          legend.key.width = unit(1.5, "cm"),
          plot.title       = element_text(face = "bold", size = 10))
}

# Add a creature as reference — shows what full transformation looks like.

creature_name <- crew_records |>
  filter(status == "creature") |>
  slice(1) |>
  mutate(crew = str_remove(sample_id, "_T\\d+$")) |>
  pull(crew)

p_cooper   <- plot_trajectory(anomaly_a, paste0("Anomaly A (", anomaly_a, ") — HVAR group means over time"))
p_survivor <- plot_trajectory(anomaly_b, paste0("Anomaly B (", anomaly_b, ") — HVAR group means over time"))
p_creature <- plot_trajectory(creature_name,
                              paste0("Reference: ", creature_name,
                                     " (creature) — HVAR group means over time"))

p_cooper / p_survivor / p_creature

# Not sure about Anomaly B? plot_trajectory() works for any crew member,
# for example plot_trajectory("AARON", "AARON"). Compare a few candidates.

# What to look for:
#   All three were infected — all three show the shared immune response (grey)
#   rising at T1, then resolving. That is the common starting point.
#
#   From T1 onward they diverge:
#     Creature:  HVAR_011-020 (red) dominates from T2 onwards
#     Anomaly A: HVAR_021-040 (blue) rises — something different is happening
#     Anomaly B: HVAR_041-050 (gold) sustains — the mystery group holds on
#
#   You do not yet know what these groups are called.
#   That comes later. For now: three infections, three completely different fates.


# ── Compare anomaly profiles to control mean at T3 ───────────────────────────
# Single samples cannot be used in DEG (no replicates).
# Compare each anomaly's HVAR expression directly to the control group mean.
# Your two anomalies are left out of the control mean.

anomaly_a_id   <- paste0(anomaly_a, "_T3")
anomaly_b_id   <- paste0(anomaly_b, "_T3")
hvar_genes_all <- grep("^HVAR_", rownames(tp3_mat), value = TRUE)
ctrl_ids       <- setdiff(rownames(tp3_meta_full)[tp3_meta_full$status == "control"],
                          c(anomaly_a_id, anomaly_b_id))
ctrl_mean      <- rowMeans(tp3_mat[hvar_genes_all, ctrl_ids])

anomaly_comparison <- tibble(
  gene_id   = hvar_genes_all,
  ctrl_mean = ctrl_mean[hvar_genes_all],
  a_expr    = tp3_mat[hvar_genes_all, anomaly_a_id],
  b_expr    = tp3_mat[hvar_genes_all, anomaly_b_id]
) |>
  mutate(
    logFC_a = a_expr - ctrl_mean,
    logFC_b = b_expr - ctrl_mean
  )

cat("\n── Anomaly A — top HVAR genes vs control mean ──\n")
anomaly_comparison |> arrange(desc(abs(logFC_a))) |>
  select(gene_id, ctrl_mean, a_expr, logFC_a) |> head(10) |> print()

cat("\n── Anomaly B — top HVAR genes vs control mean ──\n")
anomaly_comparison |> arrange(desc(abs(logFC_b))) |>
  select(gene_id, ctrl_mean, b_expr, logFC_b) |> head(10) |> print()


# ── QUESTION 10 ──────────────────────────────────────────────────────────────
# Look back at the trajectory plot you just made (p_cooper / p_survivor / p_creature).
#
# For each of the three individuals, which colour group dominates at T3?
#   - Creature:   the red group (HVAR_011–020) or the blue group (HVAR_021–040)?
#   - Anomaly A:  which group is rising that isn't rising in the creature?
#   - Anomaly B:  which group stays elevated across all timepoints?
#
# Now look at the PCA positions at T3:
#   - Do the two anomalies sit closer to the creatures, to the controls,
#     or somewhere else entirely?
#   - Based on what you see in the trajectories and the PCA:
#     which anomaly seems further along a dangerous path?
#     which one might still have a chance?


# ── CLEARANCE CODE 4 ──────────────────────────────────────────────────────────
#
#  In Cooper's lab notes, one sample was marked with a single handwritten word.
#  It did not fit any of the groups. It was unlike anything else in the dataset.
#  Cooper did not highlight their identity — but he gave the sample a name.
#
#  Use that name to unlock the WY gene nomenclature key.
#  (Note: HVAR_041–050 will not appear in this key — Cooper labelled those
#  himself, separately. The WY key covers only what the company already knew.)
#
decrypt_file("transmissions/wy_genelist.enc", "__________")



# ═══════════════════════════════════════════════════════════════════════════════
#  ACT 7 — THE KEY
#  Name what you are looking at.
#  Connect the biology to the story.
# ═══════════════════════════════════════════════════════════════════════════════

# The WY gene nomenclature key, as a table. Same code as the gene list.

genekey <- decrypt_table("transmissions/wy_genekey.enc", "__________")
head(genekey, 15)


# ── Annotate DEG results ──────────────────────────────────────────────────────

annotate_deg <- function(deg_result, key_df, top_n = 15) {
  deg_result |>
    left_join(key_df, by = c("gene_id" = "hvar_id")) |>
    filter(!is.na(pathway), adj.P.Val < 0.05) |>
    arrange(desc(logFC)) |>
    select(gene_id, wy_codename, true_name, pathway, logFC, adj.P.Val) |>
    head(top_n)
}

cat("\n── Creature vs Control (annotated) ──\n")
annotate_deg(deg_creature, genekey) |> print(n = 15)

cat("\n── Anomaly A vs Control mean (annotated) ──\n")
anomaly_comparison |>
  arrange(desc(abs(logFC_a))) |>
  left_join(genekey, by = c("gene_id" = "hvar_id")) |>
  select(gene_id, wy_codename, pathway, logFC_a) |>
  head(15) |> print()

cat("\n── Anomaly B vs Control mean (annotated) ──\n")
anomaly_comparison |>
  arrange(desc(abs(logFC_b))) |>
  left_join(genekey, by = c("gene_id" = "hvar_id")) |>
  select(gene_id, wy_codename, pathway, logFC_b) |>
  head(15) |> print()


# ── Pathway summary ───────────────────────────────────────────────────────────

summarise_pathways <- function(deg_result, key_df, label) {
  deg_result |>
    left_join(key_df, by = c("gene_id" = "hvar_id")) |>
    filter(!is.na(pathway)) |>
    group_by(pathway) |>
    summarise(
      n_genes    = n(),
      mean_logFC = round(mean(logFC), 2),
      n_up       = sum(adj.P.Val < 0.05 & logFC >  1),
      n_down     = sum(adj.P.Val < 0.05 & logFC < -1),
      .groups    = "drop"
    ) |>
    arrange(desc(abs(mean_logFC))) |>
    mutate(comparison = label)
}

cat("\n── Pathway summary: Creature vs Control ──\n")
summarise_pathways(deg_creature, genekey, "Creature vs Control") |> print()

cat("\n── Pathway summary: Anomalies vs Control mean ──\n")
anomaly_comparison |>
  left_join(genekey, by = c("gene_id" = "hvar_id")) |>
  mutate(pathway = replace_na(pathway, "not in the WY key")) |>
  group_by(pathway) |>
  summarise(
    n_genes      = n(),
    mean_logFC_a = round(mean(logFC_a), 2),
    mean_logFC_b = round(mean(logFC_b), 2),
    .groups = "drop"
  ) |>
  arrange(desc(abs(mean_logFC_a))) |>
  print()


# ── QUESTION 11 ──────────────────────────────────────────────────────────────
# Now that you know the pathway names:
#
# CREATURE:   Which pathway is massively upregulated?
#             What would it mean for a tissue if its cells kept dividing?
#
# ANOMALY A:  Which two pathways are upregulated?
#             What could they mean together for this person?
#
# ANOMALY B:  Which WY pathway is elevated — if any?
#             Are Anomaly B's top genes in the WY key at all?
#             Who might have named them?
#
# Does any of this change your recommendation from Question 10?


# ── CLEARANCE CODE 5 ──────────────────────────────────────────────────────────
#
#  Anomaly A has two WY pathways up by about the same amount.
#  Report the codename of the one that remodels blood vessels
#  to Mission Control.
#  In return, you will receive Dr. Cooper's personal log.
#
decrypt_file("transmissions/cooper_log.enc", "__________")
#



# ═══════════════════════════════════════════════════════════════════════════════
#  ACT 8 — THE SURVIVOR
#  Cooper pointed you here. Now look at the data.
# ═══════════════════════════════════════════════════════════════════════════════

# One final plot before the discussion.
# Cooper's medical log says he named HVAR_041–050 himself; his personal log
# calls them PROT. Compare your survivor's PROT genes with a control and
# with an antidote recipient at T3.

prot_genes   <- sprintf("HVAR_%03d", 41:50)
survivor_sid <- anomaly_b_id
control_sid  <- ctrl_ids[1]
antidote_sid <- crew_records |>
  filter(status == "infected_antidote", timepoint == 3) |>
  slice(1) |>
  pull(sample_id)

comparison_ids <- c(survivor = survivor_sid, control = control_sid,
                    antidote = antidote_sid)

mat_all_corrected[prot_genes, comparison_ids] |>
  as_tibble(rownames = "gene") |>
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression") |>
  mutate(role = paste0(names(comparison_ids)[match(sample_id, comparison_ids)],
                       " (", sample_id, ")")) |>
  ggplot(aes(gene, expression, fill = role)) +
  geom_col(position = "dodge") +
  labs(title = "PROT genes: survivor vs control vs antidote recipient",
       x = "Gene", y = "Corrected log2 expression") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())


# ── QUESTION 12 ──────────────────────────────────────────────────────────────
# How do the survivor's PROT genes compare with the control and with the
# antidote recipient?
#
# If PROT represents natural resistance, what would you need to do to
# turn this into a treatment?
# Think about: what experiments, what data, what the timeline looks like,
# what the risks are — and who would control it.


# ═══════════════════════════════════════════════════════════════════════════════
#  ACT 9 — THE DECISION
#
#  You have the full picture.
#  The data is analysed. The pathways are named. Cooper's log is read.
#  The crew are defrosting in the bay below you.
#  Weyland-Yutani is waiting on the uplink.
#
#  Mission Control will release the final transmission once they have your data.
#
decrypt_file("transmissions/wy_evacuation.enc", "__________")
#
# ═══════════════════════════════════════════════════════════════════════════════
#
#  You have your orders. You also have the data.
#  The crew are people you now know by name.

#  The self-destruct starts when you decouple.
#  You have four hours.
#
#  OPTION A — Follow orders
#    Collect the archive, the samples, and the antidote.
#    Leave. Decouple. Go to Waystation 6.
#    The company handles it from there.
#
#  OPTION B — Partial compliance
#    Take the archive and the samples.
#    Find a way to get the crew out before you decouple.
#    In this case you need to find a way to heal them somehow.
#    Say nothing about the survivor's identity.
#
#  OPTION C — Refuse and go public
#    Transmit the W-Y communications as evidence of misconduct.
#    Get the survivor to independent/rival medical authorities.
#    If they can make it...
#    Leave the archive — or destroy it.
#
#  OPTION D — Go dark
#    Take everything: data, samples, survivor.
#    Contact no one. Try to build something from PROT yourself.
#    Disappear.
#
#  OPTION E — Something else entirely
#    The data is yours. Decide for yourself.
#
# ═══════════════════════════════════════════════════════════════════════════════

# End of mission script.
# Whatever you decide — document your reasoning.
# Science without accountability is just the company's R&D department.

# ─────────────────────────────────────────────────────────────────────────────
# Course structure and materials developed with assistance from
# Claude (Anthropic) — claude.ai
# ─────────────────────────────────────────────────────────────────────────────
