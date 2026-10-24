# Re-draw only the two Manhattan plots from the completed Spearman CSVs.
# The plotting implementation remains in 03_nested_sample_screening.R, so a
# later full pipeline run produces the same uniformly sized inset tables.
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(args) == 1L)
run_root <- dirname(normalizePath(sub("^--file=", "", args), winslash = "/"))
.libPaths(c(file.path(dirname(run_root), "_library"), .libPaths()))
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))
exprs <- parse(file.path(run_root, "03_nested_sample_screening.R"), encoding = "UTF-8")
is_stage <- vapply(exprs, function(x)
  is.call(x) && identical(x[[1L]], as.name(".integration_run_stage")), logical(1))
# Initialize paths and the stage runner without running sample preparation.
# source() path discovery can use this script's --file argument (same directory).
for (i in seq_len(which(is_stage)[1L]-1L)) eval(exprs[[i]], envir = .GlobalEnv)
write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
target <- which(vapply(exprs, function(x) is.call(x) &&
  identical(x[[1L]], as.name(".integration_run_stage")) &&
  identical(x[[2L]], "manhattan_plots.R"), logical(1)))
stopifnot(length(target) == 1L)
eval(exprs[[target]], envir = .GlobalEnv)
jobs_to_copy <- data.frame(
  source = c(file.path(run_root,"results/analysis/nested/1_total_PA_signature/03_Spearman_FDR_Manhattan"),
             file.path(run_root,"results/analysis/nested/nonoccupational_PA/03_Spearman_nominal_P_Manhattan")),
  final = c("Total_PA_Manhattan", "Nonoccupational_PA_Manhattan"))
for (i in seq_len(nrow(jobs_to_copy))) {
  stopifnot(file.copy(paste0(jobs_to_copy$source[i], ".pdf"),
    file.path(run_root, "figures", paste0(jobs_to_copy$final[i], ".pdf")), overwrite = TRUE),
    file.copy(paste0(jobs_to_copy$source[i], ".png"),
    file.path(run_root, "figures", paste0(jobs_to_copy$final[i], "_600dpi.png")), overwrite = TRUE))
}
cat("Both Manhattan inset tables use one common font size; PDF/600-dpi PNG copies updated.\n")
