# Export from already fitted objects. No refitting or random-seed search.
if (!exists(".integration_code_root", inherits = TRUE)) {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  stopifnot(length(a) == 1L)
  .integration_code_root <- dirname(normalizePath(sub("^--file=", "", a), winslash = "/"))
}
run_root <- .integration_code_root
.libPaths(c(file.path(dirname(run_root), "_library"), .libPaths()))
suppressPackageStartupMessages(library(glmnet))
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))
if (!"Arial" %in% names(pdfFonts())) pdfFonts(Arial = pdfFonts("Helvetica")[[1L]])
nested_root <- file.path(run_root, "results", "analysis", "nested")
total_root <- file.path(nested_root, "total_PA")
leisure_root <- file.path(nested_root, "nonoccupational_PA", "forced_15_seed257_complete_analysis")
fig_root <- file.path(run_root, "figures")
dir.create(fig_root, showWarnings = FALSE)
read_csv <- function(path) read.csv(path, check.names = FALSE, fileEncoding = "UTF-8")
total_obj <- readRDS(file.path(total_root, "15_complete_analysis_objects.rds"))
leisure_obj <- readRDS(file.path(leisure_root, "16_analysis_objects.rds"))
spec_total <- read_csv(file.path(total_root, "03_LASSO_specification.csv"))
spec_leisure <- read_csv(file.path(leisure_root, "05_analysis_specification.csv"))
screen_counts <- read_csv(file.path(nested_root, "02_screening_counts_794.csv"))
stopifnot(spec_total$Seed == 20L, spec_leisure$Seed == 257L,
  spec_leisure$Unpenalized_count == 15L, spec_leisure$Candidate_count == 52L,
  all(screen_counts$N_tested == 794L))

# Total-PA LASSO AB is drawn from the SAME saved CV path and selected-lambda fit.
cf <- as.matrix(coef(total_obj$fit))[, 1L]
coef_tab <- read_csv(file.path(total_root, "01_signature_metabolite_coefficients.csv"))
stopifnot(max(abs(cf[coef_tab$Metabolite] - coef_tab$LASSO_coefficient)) < 1e-9,
  nrow(coef_tab) == spec_total$N_selected)
draw_total_AB <- function() {
  par(mfrow = c(1, 2), family = "Arial", mar = c(4.6, 4.7, 4.5, 1),
    mgp = c(2.8, .8, 0), cex = 1, cex.axis = .95, cex.lab = 1, las = 0, tcl = -.3, bty = "l")
  plot(total_obj$cvfit, xlab = expression(log(lambda)), ylab = "Mean-squared error")
  mtext("(A)", side = 3, line = 3.1, adj = 0, cex = 1.2)
  plot(total_obj$cvfit$glmnet.fit, xvar = "norm", label = FALSE,
    xlab = "L1 norm", ylab = "Coefficients")
  mtext("(B)", side = 3, line = 3.1, adj = 0, cex = 1.2)
}
png(file.path(total_root, "16_LASSO_AB_600dpi.png"), width = 12, height = 5.5,
  units = "in", res = 600, type = "cairo", family = "Arial", bg = "white")
draw_total_AB(); invisible(dev.off())
cairo_pdf(file.path(total_root, "16_LASSO_AB.pdf"), width = 12, height = 5.5,
  family = "Arial", bg = "white")
draw_total_AB(); invisible(dev.off())

figure_sources <- data.frame(
  Name = c("Total_PA_Manhattan", "Total_PA_signature_RCS", "Total_PA_LASSO_AB",
    "Screened_metabolites_CRC_ring", "Nonoccupational_PA_Manhattan",
    "Nonoccupational_PA_signature_RCS", "Nonoccupational_PA_LASSO_AB", "Full_cohort_PA_RCS_tertile_lines"),
  PDF = c(file.path(nested_root, "1_total_PA_signature", "03_Spearman_FDR_Manhattan.pdf"),
    file.path(total_root, "10_signature_CRC_RCS_3knots.pdf"), file.path(total_root, "16_LASSO_AB.pdf"),
    file.path(total_root, "screened_CRC", "total_PA_screened_metabolites_CRC_ring.pdf"),
    file.path(nested_root, "nonoccupational_PA", "03_Spearman_nominal_P_Manhattan.pdf"),
    file.path(leisure_root, "15_signature_CRC_RCS_3knots.pdf"),
    file.path(leisure_root, "17_LASSO_regression_metabolite_selection_AB.pdf"),
    file.path(run_root, "Fig S2_tertile_lines.pdf")),
  PNG = c(file.path(nested_root, "1_total_PA_signature", "03_Spearman_FDR_Manhattan.png"),
    file.path(total_root, "10_signature_CRC_RCS_3knots.png"), file.path(total_root, "16_LASSO_AB_600dpi.png"),
    file.path(total_root, "screened_CRC", "total_PA_screened_metabolites_CRC_ring.png"),
    file.path(nested_root, "nonoccupational_PA", "03_Spearman_nominal_P_Manhattan.png"),
    file.path(leisure_root, "15_signature_CRC_RCS_3knots.png"),
    file.path(leisure_root, "17_LASSO_regression_metabolite_selection_AB.png"),
    file.path(run_root, "Fig S2_tertile_lines_600dpi.png")))
stopifnot(all(file.exists(figure_sources$PDF)), all(file.exists(figure_sources$PNG)))
for (i in seq_len(nrow(figure_sources))) {
  stopifnot(file.copy(figure_sources$PDF[i], file.path(fig_root, paste0(figure_sources$Name[i], ".pdf")), overwrite = TRUE),
    file.copy(figure_sources$PNG[i], file.path(fig_root, paste0(figure_sources$Name[i], "_600dpi.png")), overwrite = TRUE))
}
figure_sources$PDF <- substring(figure_sources$PDF, nchar(run_root) + 2L)
figure_sources$PNG <- substring(figure_sources$PNG, nchar(run_root) + 2L)
write.csv(figure_sources, file.path(run_root, "results", "FIGURE_INDEX.csv"), row.names = FALSE, fileEncoding = "UTF-8")

total_crc <- read_csv(file.path(total_root, "TableS5_PA_signature_CRC_model1_model2_numeric.csv"))
t <- total_crc[total_crc$Model == "Model 1" & total_crc$Exposure == "Metabolic signature", ]
l <- read_csv(file.path(leisure_root, "10_signature_CRC_conditional_logistic.csv"))
l <- l[l$Model == "Full", ]
ring <- read_csv(file.path(total_root, "screened_CRC", "01_screened_metabolites_CRC_full_model.csv"))
summary <- data.frame(
  Domain = c("Total physical activity", "Non-occupational PA (original leisure-time label)"),
  Eligible_metabolites = 794L, Screening_N = 288L, CRC_N = 256L, CRC_matched_pairs = 128L,
  Spearman_criterion = screen_counts$Criterion,
  Spearman_selected_N = screen_counts$N_selected,
  LASSO_candidates_N = c(spec_total$N_candidates, spec_leisure$Candidate_count),
  Fixed_seed = c(spec_total$Seed, spec_leisure$Seed),
  Forced_metabolites_N = c(1L, spec_leisure$Unpenalized_count),
  Signature_metabolites_N = c(spec_total$N_selected, spec_leisure$Selected_count),
  Signature_CRC_OR_per_SD = c(t$OR, l$OR), Signature_CRC_CI_lower = c(t$CI_lower, l$CI_lower),
  Signature_CRC_CI_upper = c(t$CI_upper, l$CI_upper), Signature_CRC_P = c(t$P, l$P))
write.csv(summary, file.path(run_root, "results", "RUN_SUMMARY.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(data.frame(N_screened = nrow(ring), CRC_nominal_P_lt_0.05 = sum(ring$P < .05),
  CRC_FDR_lt_0.05 = sum(ring$FDR_P < .05), FDR_family_N = nrow(ring)),
  file.path(run_root, "results", "INDIVIDUAL_CRC_SUMMARY.csv"), row.names = FALSE)
cat("\nUPDATED SIGNATURE RESULTS\n"); print(summary, row.names = FALSE)

# Human-readable lookup list without renumbering unverified manuscript figures.
all_outputs <- list.files(file.path(run_root, "results"), recursive = TRUE, full.names = TRUE)
all_outputs <- all_outputs[!dir.exists(all_outputs)]
write.csv(data.frame(Relative_path = substring(all_outputs, nchar(run_root) + 2L),
  Bytes = file.info(all_outputs)$size), file.path(run_root, "RESULT_FILES.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")
