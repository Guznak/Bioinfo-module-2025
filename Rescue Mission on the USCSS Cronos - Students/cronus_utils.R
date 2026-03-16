# ═══════════════════════════════════════════════════════════════════════════════
#  USCSS CRONUS — UTILITY FUNCTIONS
#  Load this file at the top of your mission script:
#    source("cronus_utils.R")
# ═══════════════════════════════════════════════════════════════════════════════


# ── decrypt_file() ─────────────────────────────────────────────────────────────
#
#  Decrypts a .enc file using a clearance code (password).
#  The encryption is a simple XOR cipher — the same key that locked the file
#  unlocks it. You can read this function and understand exactly how it works.
#
#  Usage:
#    decrypt_file("ava_message_1.enc", "BATCHCLEAR")
#
decrypt_file <- function(filename, password) {

  # Read the hex-encoded encrypted content
  hex_str <- readLines(filename, warn = FALSE)
  hex_str <- paste(hex_str, collapse = "")

  # Convert hex string to raw bytes
  hex_pairs  <- substring(hex_str, seq(1, nchar(hex_str)-1, 2),
                                   seq(2, nchar(hex_str),   2))
  raw_cipher <- as.raw(strtoi(hex_pairs, base = 16L))

  # Build the repeating key from the password (uppercase, trimmed)
  pw      <- toupper(trimws(password))
  key_raw <- charToRaw(pw)
  key_rep <- key_raw[((seq_along(raw_cipher) - 1) %% length(key_raw)) + 1]

  # XOR to decrypt
  decrypted_raw <- as.raw(bitwXor(as.integer(raw_cipher),
                                  as.integer(key_rep)))

  # Convert back to text and print
  result <- rawToChar(decrypted_raw)

  if (grepl("[\x01-\x08\x0e-\x1f]", result)) {
    cat("\n[ INCORRECT CLEARANCE CODE — ACCESS DENIED ]\n\n")
    return(invisible(NULL))
  }

  cat(result)
  return(invisible(result))
}


# ── remove_batch_effect() ──────────────────────────────────────────────────────
#
#  Removes batch effects from a log-expression matrix.
#  Wraps limma::removeBatchEffect() with validation.
#
#  Arguments:
#    expr_matrix  : genes × samples matrix (log-scale)
#    batch        : character vector of batch labels (length = ncol(expr_matrix))
#    verbose      : print progress messages (default TRUE)
#
#  Returns:
#    A corrected genes × samples matrix, same dimensions as input.
#
remove_batch_effect <- function(expr_matrix, batch, verbose = TRUE) {

  # ── Input validation (same checks limma would run) ──────────────────────────
  if (!is.matrix(expr_matrix))
    stop("expr_matrix must be a matrix.")

  if (ncol(expr_matrix) != length(batch))
    stop(sprintf(
      "Length of 'batch' (%d) must equal number of columns in expr_matrix (%d).",
      length(batch), ncol(expr_matrix)
    ))

  if (anyNA(expr_matrix))
    stop("expr_matrix contains NA values. Please handle missing data first.")

  if (length(unique(batch)) < 2)
    stop("At least two distinct batch labels are required.")

  batch <- as.factor(batch)

  if (verbose) {
    cat(sprintf(
      "Removing batch effects: %d genes x %d samples, %d batches (%s)\n",
      nrow(expr_matrix), ncol(expr_matrix),
      nlevels(batch), paste(levels(batch), collapse = ", ")
    ))
  }

  # ── Load pre-computed corrected matrix ──────────────────────────────────────
  #  (In a real analysis this would call limma::removeBatchEffect directly.
  #   For this course dataset, we return a pre-engineered matrix that
  #   guarantees clean PCA geometry for teaching purposes.)

  corrected_path <- "USCSS_corrected_precomputed.csv"

  if (!file.exists(corrected_path))
    stop(paste(
      "Corrected matrix file not found:", corrected_path,
      "\nMake sure all course data files are in your working directory."
    ))

  corrected <- read.csv(corrected_path, row.names = 1, check.names = FALSE)
  corrected <- as.matrix(corrected)

  # Return only columns matching the input samples, in input order
  sample_ids <- colnames(expr_matrix)
  missing    <- setdiff(sample_ids, colnames(corrected))
  if (length(missing) > 0)
    stop(paste("These sample IDs were not found in the corrected matrix:",
               paste(missing, collapse = ", ")))

  corrected <- corrected[, sample_ids, drop = FALSE]

  if (verbose)
    cat("Batch correction complete.\n")

  return(corrected)
}
