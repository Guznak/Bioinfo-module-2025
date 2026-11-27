# Catch-up after Exercise 3
# Use this if your own import or matrix code is stuck, so you can continue with
# Exercise 4. Open AVA_training.Rproj first, then run:
#
#   source("catch_up/after_exercise_3.R")
#
# Afterwards you have: combined_data, metadata_combined, expr, in.dir, out.dir.
# This is one possible solution. Read it and compare it with your own attempt:
# which step was different?

source("check_setup.R")
library(tidyverse)
library(limma)

out.dir <- "output"
in.dir <- file.path("data", "bulk")
dir.create(out.dir, showWarnings = FALSE)

# Exercise 2: read both studies and keep the study name in a column.
archive_types <- cols(.default = col_double(),
                      ID = col_character(), class = col_character())
csv_files <- c("GSE64912.csv", "GSE68799.csv")
data_list <- map(file.path(in.dir, csv_files),
                 ~ read_csv(.x, col_types = archive_types, num_threads = 1))
stopifnot(identical(names(data_list[[1]]), names(data_list[[2]])))
names(data_list) <- str_remove(csv_files, "\\.csv$")
data_list <- imap(data_list, ~ mutate(.x, dataset = .y))
combined_data <- bind_rows(data_list)
rm(data_list)

# Exercise 3: sample information (samples x annotations) and
# expression (genes x samples), linked by the sample IDs.
metadata_combined <- combined_data %>%
  select(ID, class, dataset) %>%
  column_to_rownames("ID") %>%
  as.data.frame()
expr <- combined_data %>%
  select(-class, -dataset) %>%
  column_to_rownames("ID") %>%
  as.matrix() %>%
  t()

stopifnot(!anyDuplicated(rownames(metadata_combined)),
          is.numeric(expr), all(is.finite(expr)),
          identical(colnames(expr), rownames(metadata_combined)))
message("Catch-up complete: expr has ", nrow(expr), " genes x ", ncol(expr),
        " samples; metadata_combined has ", nrow(metadata_combined), " samples.")
