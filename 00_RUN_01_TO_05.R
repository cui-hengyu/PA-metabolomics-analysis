# One entry point: exact original 01 -> 02 -> 03 -> 04 -> 05 order.
options(stringsAsFactors = FALSE, warn = 1)
a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(a) == 1L)
.pipeline_root <- dirname(normalizePath(sub("^--file=", "", a), winslash = "/"))
.pipeline_input_root <- dirname(.pipeline_root)
.libPaths(c(file.path(.pipeline_input_root, "_library"), .libPaths()))
# Register compatible metrics for pre-device grid layout; Cairo PDF exports
# still use the installed Arial font itself.
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))
if (!"Arial" %in% names(pdfFonts())) pdfFonts(Arial = pdfFonts("Helvetica")[[1L]])
# Grid may request a measurement device before ggsave; avoid an extra Rplots.pdf.
options(device = function(...) grDevices::pdf(file = NULL, ...))
.pipeline_script_names <- c("01_data_preprocessing.R", "02_full_cohort_analysis.R",
  "03_nested_sample_screening.R", "04_total_pa_metabolite_score.R", "05_leisure_pa_metabolite_score.R")
.pipeline_protected_inputs <- c(file.path(.pipeline_input_root, .pipeline_script_names),
  file.path(.pipeline_input_root, "data", c("cohort.xlsx", "metabolomics.xlsx")))
.pipeline_before_md5 <- tools::md5sum(.pipeline_protected_inputs)
dir.create(file.path(.pipeline_root, "results", "provenance"), recursive = TRUE, showWarnings = FALSE)
for (.pipeline_file in c(.pipeline_script_names, "90_summary_and_figures.R", "91_export_manuscript_tables.R")) invisible(parse(file = file.path(.pipeline_root, .pipeline_file), encoding = "UTF-8"))
.pipeline_phases <- list()
for (.pipeline_i in seq_along(.pipeline_script_names)) {
  .pipeline_started <- Sys.time()
  cat(sprintf("\n[PIPELINE %d/5] %s -- %s\n", .pipeline_i, .pipeline_script_names[.pipeline_i], format(.pipeline_started)))
  flush.console()
  source(file.path(.pipeline_root, .pipeline_script_names[.pipeline_i]), local = .GlobalEnv, encoding = "UTF-8")
  .pipeline_phases[[.pipeline_i]] <- data.frame(Script = .pipeline_script_names[.pipeline_i], Started = format(.pipeline_started),
    Finished = format(Sys.time()), Status = "Completed")
  write.csv(do.call(rbind, .pipeline_phases), file.path(.pipeline_root, "results", "provenance", "completed_stages.csv"),
    row.names = FALSE, fileEncoding = "UTF-8")
}
cat("\n[FINAL] Exporting named figures and a compact results summary.\n")
source(file.path(.pipeline_root, "90_summary_and_figures.R"), local = .GlobalEnv, encoding = "UTF-8")
cat("\n[TABLES] Exporting manuscript Table 3 and Tables S2, S3, S5, S6, S7 (no S4).\n")
source(file.path(.pipeline_root, "91_export_manuscript_tables.R"), local = .GlobalEnv, encoding = "UTF-8")
.pipeline_after_md5 <- tools::md5sum(.pipeline_protected_inputs)
stopifnot(identical(.pipeline_before_md5, .pipeline_after_md5))
write.csv(data.frame(File = names(.pipeline_before_md5), MD5_before = unname(.pipeline_before_md5), MD5_after = unname(.pipeline_after_md5)),
  file.path(.pipeline_root, "results", "provenance", "original_input_checksums.csv"), row.names = FALSE, fileEncoding = "UTF-8")
capture.output(sessionInfo(), file = file.path(.pipeline_root, "results", "sessionInfo.txt"))
cat("\nCOMPLETE: all five requested stages finished; original scripts and Excel inputs unchanged.\n")
