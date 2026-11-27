# ═══════════════════════════════════════════════════════════════════════════════
#  USCSS CRONUS — UTILITY FUNCTIONS
#  Loaded at the top of the mission script:  source("cronus_utils.R")
# ═══════════════════════════════════════════════════════════════════════════════


# ── decrypt_file() ─────────────────────────────────────────────────────────────
#
#  Decrypts a transmission (.enc file) with a clearance code and prints it.
#  The encryption is a simple XOR cipher: the same code that locked the file
#  unlocks it. You can read this function and see exactly how it works.
#
#  Usage:
#    decrypt_file("transmissions/ava_message_1.enc", "CODE")
#
decrypt_file <- function(filename, password) {
  text <- decrypt_text(filename, password)
  if (is.null(text)) return(invisible(NULL))
  cat(text)
  invisible(text)
}


# ── decrypt_table() ────────────────────────────────────────────────────────────
#
#  Decrypts an encrypted table (.enc file) and returns it as a tibble.
#
#  Usage:
#    records <- decrypt_table("transmissions/crew_records.enc", "CODE")
#
decrypt_table <- function(filename, password) {
  text <- decrypt_text(filename, password)
  if (is.null(text)) stop("No table without the right clearance code. Ask Mission Control.", call. = FALSE)
  tibble::as_tibble(utils::read.csv(text = text, stringsAsFactors = FALSE))
}


# ── decrypt_text() ─────────────────────────────────────────────────────────────
#  Shared by both functions above: returns the decrypted text, or NULL and a
#  message when the code is wrong.
#
decrypt_text <- function(filename, password) {
  if (!file.exists(filename)) stop("File not found: ", filename)

  # Read the hex-encoded encrypted content and convert it to raw bytes
  hex_str    <- paste(readLines(filename, warn = FALSE), collapse = "")
  hex_pairs  <- substring(hex_str, seq(1, nchar(hex_str) - 1, 2),
                                   seq(2, nchar(hex_str),     2))
  raw_cipher <- as.raw(strtoi(hex_pairs, base = 16L))

  # Build the repeating key from the code (upper case, spaces removed)
  pw      <- toupper(gsub("[[:space:]]", "", password))
  key_raw <- charToRaw(pw)
  key_rep <- key_raw[((seq_along(raw_cipher) - 1) %% length(key_raw)) + 1]

  # XOR to decrypt
  decrypted <- as.raw(bitwXor(as.integer(raw_cipher), as.integer(key_rep)))

  # A wrong code produces control characters; stop before they reach the screen
  if (any(as.integer(decrypted) %in% c(0:8, 14:31))) {
    cat("\n[ INCORRECT CLEARANCE CODE — ACCESS DENIED ]\n\n")
    return(NULL)
  }
  text <- rawToChar(decrypted)
  Encoding(text) <- "UTF-8"
  text
}
