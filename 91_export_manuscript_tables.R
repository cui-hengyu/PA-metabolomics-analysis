# Export the six requested manuscript tables from the completed 794-feature run.
# No refitting, new randomization, or change to any model specification.
if (!exists(".integration_code_root", inherits = TRUE)) {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  stopifnot(length(a) == 1L)
  .integration_code_root <- dirname(normalizePath(sub("^--file=", "", a), winslash = "/"))
}
local({
run_root <- .integration_code_root
nested <- file.path(run_root, "results", "analysis", "nested")
total <- file.path(nested, "total_PA")
nonocc <- file.path(nested, "nonoccupational_PA", "forced_15_seed257_complete_analysis")
table_dir <- file.path(run_root, "table")
dir.create(table_dir, showWarnings = FALSE)
read_csv <- function(path) read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                                   fileEncoding = "UTF-8")
write_csv_bom <- function(x, name) {
  path <- file.path(table_dir, name)
  con <- file(path, open = "wb")
  on.exit(close(con))
  writeBin(as.raw(c(239, 187, 191)), con)
  lines <- capture.output(write.table(x, file = "", sep = ",", row.names = FALSE,
    col.names = TRUE, na = "", qmethod = "double"))
  writeLines(enc2utf8(lines), con, useBytes = TRUE)
}

screen <- read_csv(file.path(nested, "1_total_PA_signature", "02_Spearman_FDR_lt_0.05_full_nested.csv"))
total_coef <- read_csv(file.path(total, "01_signature_metabolite_coefficients.csv"))
nonocc_coef <- read_csv(file.path(nonocc, "02_signature_metabolites_and_coefficients.csv"))
# The original non-occupational export uses empty strings for missing HMDB IDs;
# use the same blank CSV representation as the other tables.
nonocc_coef$HMDB_ID[which(nonocc_coef$HMDB_ID == "")] <- NA_character_
crc176 <- read_csv(file.path(total, "screened_CRC", "01_screened_metabolites_CRC_full_model.csv"))
stopifnot(nrow(screen) == 176L, nrow(total_coef) == 18L, nrow(nonocc_coef) == 16L,
          nrow(crc176) == 176L, all(screen$p_fdr < .05),
          setequal(screen$COMPOUND.Name, crc176$Metabolite),
          all(total_coef$Metabolite %in% screen$COMPOUND.Name))
excluded <- c("3-amino-2-piperidone", "4-hydroxyphenylacetylglutamine",
              "N,N-dimethyl-pro-pro", "O-acetylhomoserine", "alpha-tocopherol")
stopifnot(!any(screen$COMPOUND.Name %in% excluded),
          !any(total_coef$Metabolite %in% excluded), !any(nonocc_coef$Metabolite %in% excluded))

# Table 3: exactly the unadjusted and fully adjusted signature models, not the
# mutually adjusted model. Ranges and the score SD refer to the 288-person sample.
t3 <- read_csv(file.path(total, "05_signature_CRC_tertiles_and_perSD.csv"))
stopifnot(identical(t3$Model, c("Unadjusted", "Full model")))
table3 <- data.frame(Variable = "Metabolic signature", Model = c("Univariate", "Full model"),
  T1 = t3$T1, T2 = t3$T2, T3 = t3$T3,
  Per_1SD_OR_95CI = t3[["Per 1-SD"]], P = t3$P,
  T1_range = t3$T1_range, T2_range = t3$T2_range, T3_range = t3$T3_range,
  Signature_SD = t3$Signature_SD, N = 256L, Matched_pairs = 128L,
  SD_and_tertile_reference = t3$Range_and_SD_reference)
write_csv_bom(table3, "Table_3.csv")

# Table S2: retain full numerical precision, with nominal P also supplied.
screen <- screen[order(screen$p_fdr, screen$P.value, screen$COMPOUND.Name), ]
s2 <- data.frame(Metabolite = screen$COMPOUND.Name,
  Super_pathway = screen$SUPER.META.PATHWAY, HMDB_ID = screen$HMDB.ID,
  Spearman_r = screen$Estimate, Q_value = screen$p_fdr, P_value = screen$P.value)
write_csv_bom(s2, "Table_S2.csv")

coef_columns <- c("Metabolite", "Super_pathway", "HMDB_ID", "LASSO_coefficient")
write_csv_bom(total_coef[, coef_columns], "Table_S3.csv")

# Table S4 is deliberately absent: use the 176-metabolite CRC analysis instead.
# Table S5: two model rows, paired exposure columns, and unrounded P values.
mutual <- read_csv(file.path(total, "TableS5_PA_signature_CRC_model1_model2_numeric.csv"))
stopifnot(nrow(mutual) == 4L, all(mutual$N == 256L), all(mutual$Matched_pairs == 128L))
pa <- mutual[match(paste(c("Model 1", "Model 2"), "Physical activity"),
                   paste(mutual$Model, mutual$Exposure)), ]
sig <- mutual[match(paste(c("Model 1", "Model 2"), "Metabolic signature"),
                    paste(mutual$Model, mutual$Exposure)), ]
stopifnot(!anyNA(pa$P), !anyNA(sig$P), abs(sig$P[1] - t3$P[2]) < 1e-12)
fmt_or <- function(z) sprintf("%.2f (%.2f-%.2f)", z$OR, z$CI_lower, z$CI_upper)
s5 <- data.frame(Model = c("Model 1", "Model 2"),
  Total_PA_per_1SD_OR_95CI = fmt_or(pa), Total_PA_P = pa$P,
  Signature_per_1SD_OR_95CI = fmt_or(sig), Signature_P = sig$P,
  Total_PA_OR = pa$OR, Total_PA_CI_lower = pa$CI_lower, Total_PA_CI_upper = pa$CI_upper,
  Signature_OR = sig$OR, Signature_CI_lower = sig$CI_lower, Signature_CI_upper = sig$CI_upper,
  Total_PA_SD_MET_h_day = pa$Exposure_SD, Signature_SD = sig$Exposure_SD,
  N = pa$N, Matched_pairs = pa$Matched_pairs)
write_csv_bom(s5, "Table_S5.csv")
write_csv_bom(nonocc_coef[, coef_columns], "Table_S6.csv")

t7 <- read_csv(file.path(nonocc, "10_signature_CRC_conditional_logistic.csv"))
stopifnot(identical(t7$Model, c("Unadjusted", "Full")))
s7 <- data.frame(Variable = "Metabolic signature", Model = c("Univariate", "Full model"),
  T1 = t7$T1, T2 = t7$T2, T3 = t7$T3, Per_1SD_OR_95CI = t7$Per_1SD_OR_95CI,
  P = t7$P, T1_range = t7$T1_range, T2_range = t7$T2_range, T3_range = t7$T3_range,
  Signature_SD = t7$Signature_SD, N = t7$N, Matched_pairs = t7$Matched_pairs,
  SD_and_tertile_reference = "Full nested derivation sample (N=288)")
write_csv_bom(s7, "Table_S7.csv")

# Re-open the delivered CSVs, checking row counts, full-precision values, and keys.
tables <- list(Table_3.csv = table3, Table_S2.csv = s2, Table_S3.csv = total_coef[, coef_columns],
               Table_S5.csv = s5, Table_S6.csv = nonocc_coef[, coef_columns], Table_S7.csv = s7)
for (f in names(tables)) {
  got <- read.csv(file.path(table_dir, f), check.names = FALSE, stringsAsFactors = FALSE,
                  fileEncoding = "UTF-8-BOM", na.strings = "")
  expected <- tables[[f]]
  stopifnot(identical(names(got), names(expected)), nrow(got) == nrow(expected))
  for (col in names(expected)) {
    if (is.numeric(expected[[col]])) {
      stopifnot(isTRUE(all.equal(got[[col]], expected[[col]], tolerance = 1e-12)))
    } else {
      if (!identical(got[[col]], expected[[col]])) {
        stop(sprintf("CSV round-trip mismatch: %s / %s", f, col))
      }
    }
  }
}
cat("Paper CSV tables exported and re-read successfully:\n")
print(data.frame(Table = names(tables), Rows = vapply(tables, nrow, integer(1))), row.names = FALSE)
inventory <- unlist(lapply(c("results", "table"), function(subdir)
  list.files(file.path(run_root, subdir), recursive = TRUE, full.names = TRUE)))
inventory <- inventory[!dir.exists(inventory)]
write.csv(data.frame(Relative_path = substring(inventory, nchar(run_root) + 2L),
  Bytes = file.info(inventory)$size), file.path(run_root, "RESULT_FILES.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")
})
