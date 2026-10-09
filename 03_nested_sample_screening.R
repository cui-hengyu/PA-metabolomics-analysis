# Nested sample and metabolite screening.
# Historical scores are retained; final scores are generated in scripts 04 and 06.

# Resolve the package directory from source() or Rscript.
.integration_code_root <- local({
  frames <- sys.frames()
  sources <- vapply(frames, function(x) {
    value <- x$ofile
    if (is.null(value) || length(value) != 1L) "" else as.character(value)
  }, character(1))
  sources <- sources[nzchar(sources)]
  if (length(sources)) {
    path <- tail(sources, 1L)
  } else {
    argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(argument) != 1L) stop("Run with Rscript or source() the saved file.")
    path <- sub("^--file=", "", argument)
  }
  dirname(normalizePath(path, winslash = "/", mustWork = TRUE))
})
.integration_data <- file.path(.integration_code_root, "data")
.integration_results <- file.path(.integration_code_root, "results")
.integration_analysis <- file.path(.integration_results, "analysis")
.integration_nested <- file.path(.integration_analysis,
  "13_full_nested_derived_signatures_nominalP0.05_lag1_CRC")
.integration_total <- file.path(.integration_nested, "seed20_GPA_only_181met_signature_complete_analysis")
.integration_leisure <- file.path(.integration_nested, "3_leisure_PA_signature")
dir.create(.integration_results, recursive = TRUE, showWarnings = FALSE)
if (.Platform$OS.type == "windows" && getRversion() >= "4.3.0" &&
    getRversion() < "4.4.0" && dir.exists(file.path(.integration_code_root, "_library"))) {
  .libPaths(c(file.path(.integration_code_root, "_library"), .libPaths()))
}

# Run a stage and restore the caller directory.
.integration_run_stage <- function(stage, working_directory, code,
                                   envir = parent.frame()) {
  dir.create(working_directory, recursive = TRUE, showWarnings = FALSE)
  message("[Stage] ", stage)
  previous_directory <- getwd()
  on.exit(setwd(previous_directory), add = TRUE)
  setwd(working_directory)
  eval(substitute(code), envir = envir)
}

# Build missing prerequisites within this package.
.integration_require_output <- function(script, output) {
  if (!file.exists(output)) {
    message("[Prerequisite] ", script)
    source(file.path(.integration_code_root, script), local = .GlobalEnv, encoding = "UTF-8")
  }
  if (!file.exists(output)) stop("Required output was not generated: ", output)
  invisible(output)
}

# 1. Nested baseline and PA analyses.
.integration_run_stage("code_baseline_and_PA.R",
  file.path(.integration_analysis, "01_baseline_and_PA"), {
# Exclude CRC cases diagnosed within one year and their matched controls.
# Screen at FDR <0.05; glycine is not forced.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(ggplot2)
  library(rms)
  library(tableone)
  library(gridExtra)
})

if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

project_dir <- .integration_code_root
output_root <- .integration_analysis

dirs <- c(
  baseline_pa = file.path(output_root, "01_baseline_and_PA"),
  total = file.path(output_root, "02_total_PA_metabolomics"),
  nonoccupational = file.path(output_root, "03_nonoccupational_PA_metabolomics"),
  rings = file.path(output_root, "04_ring_plots"),
  mediation = file.path(output_root, "05_mediation_and_mutual_adjustment"),
  summary = file.path(output_root, "06_summary")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_est <- function(or, lo, hi) {
  if (any(!is.finite(c(or, lo, hi)))) return("")
  sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
}
extract_or <- function(fit, term) {
  tab <- summary(fit)$coefficients
  if (!term %in% rownames(tab)) return(c(OR = NA, Lower = NA, Upper = NA, P = NA))
  b <- tab[term, "coef"]
  se <- tab[term, "se(coef)"]
  c(OR = exp(b), Lower = exp(b - 1.96 * se), Upper = exp(b + 1.96 * se),
    P = tab[term, "Pr(>|z|)"])
}
fit_clogit <- function(formula, data) {
  suppressWarnings(clogit(formula, data = data, method = "exact", na.action = na.fail))
}

# Data preparation from raw Excel.
source(file.path(.integration_code_root, "01_data_preprocessing.R"), local = TRUE, encoding = "UTF-8")
stopifnot(exists("data1"), exists("metarowdata_0.2_imputed1"),
          exists("scaled_metarowdata"), exists("metarowdata2"))

smoking_levels <- c(
  "Never", "Former, <20 pack-years", "Former, >=20 pack-years",
  "Current, <20 pack-years", "Current, >=20 pack-years", "Unknown"
)
data1$smoking5 <- with(data1, ifelse(
  B1 == 2, "Never",
  ifelse(B3 == 3 & packyear1 < 20, "Former, <20 pack-years",
  ifelse(B3 == 3 & packyear1 >= 20, "Former, >=20 pack-years",
  ifelse(B3 %in% c(1, 2) & packyear1 < 20, "Current, <20 pack-years",
  ifelse(B3 %in% c(1, 2) & packyear1 >= 20,
         "Current, >=20 pack-years", "Unknown"))))
))
data1$smoking5[is.na(data1$smoking5)] <- "Unknown"  # Original missing-value code: 999.
data1$smoking5 <- factor(data1$smoking5, levels = smoking_levels)

pa_days_or_zero <- function(x) ifelse(is.na(x), 0, x)
pa_fields <- c("E2", "E33", "E5", "E66", "E8", "E99", "E11", "E122", "E14", "E155")
if (!all(pa_fields %in% names(data1))) {
  stop("Missing PA variables: ", paste(setdiff(pa_fields, names(data1)), collapse = ", "))
}
data1$occupational_PA_MET_h_day <-
  (8 * pa_days_or_zero(data1$E2) * data1$E33 +
     4 * pa_days_or_zero(data1$E5) * data1$E66) / 60 / 7
data1$leisure_PA_MET_h_day <-
  (8 * pa_days_or_zero(data1$E11) * data1$E122 +
     4 * pa_days_or_zero(data1$E14) * data1$E155) / 60 / 7
data1$transport_PA_MET_h_day <-
  4 * pa_days_or_zero(data1$E8) * data1$E99 / 60 / 7
data1$nonoccupational_PA_MET_h_day <-
  data1$leisure_PA_MET_h_day + data1$transport_PA_MET_h_day

dat <- metarowdata_0.2_imputed1
idx <- match(dat$ID, data1$ID)
stopifnot(!anyNA(idx), !anyDuplicated(data1$ID))
update_vars <- c(
  "AGE", "A2", "BMI", "education", "smoking5", "drink", "DMTIME", "H1H2",
  "redmeat", "fruit", "veg", "occupational_PA_MET_h_day",
  "nonoccupational_PA_MET_h_day"
)
for (v in update_vars) dat[[v]] <- data1[[v]][idx]
dat$total_PA_MET_h_day <- as.numeric(dat$MET)

metabolite_original <- colnames(scaled_metarowdata)
if (length(metabolite_original) != 799L) {
  warning("Current metabolite count is ", length(metabolite_original), ", not 799.")
}
metabolite_clean <- gsub("\\.", " ", metabolite_original)
met_idx <- match(metabolite_original, names(dat))
if (anyNA(met_idx)) stop("Metabolite columns could not be matched to the nested dataset.")
names(dat)[met_idx] <- metabolite_clean

dat$Group <- as.integer(as.character(dat$Group))
dat$Match <- factor(dat$Match)
dat$education <- factor(dat$education)
dat$smoking5 <- droplevels(factor(dat$smoking5, levels = smoking_levels))
dat$drink <- factor(dat$drink)
dat$H1H2 <- factor(dat$H1H2)
dat$A2 <- factor(dat$A2)

case_within_one_year <- dat$Group == 1 & as.numeric(dat$endtime5) <= 1
removed_matches <- unique(as.character(dat$Match[case_within_one_year]))
lag1_dat <- droplevels(dat[!(as.character(dat$Match) %in% removed_matches), , drop = FALSE])
pair_tab <- table(lag1_dat$Match, lag1_dat$Group)
if (ncol(pair_tab) != 2L || any(pair_tab != 1L)) {
  stop("Lag-1 sample is not composed exclusively of intact 1:1 matched sets.")
}

base_covars <- c(
  "AGE", "BMI", "education", "smoking5", "drink", "DMTIME", "H1H2",
  "redmeat", "fruit", "veg"
)
pa_exposures <- c(
  total = "total_PA_MET_h_day",
  occupational = "occupational_PA_MET_h_day",
  nonoccupational = "nonoccupational_PA_MET_h_day"
)
pa_labels <- c(
  total = "Physical activity",
  occupational = "Occupational physical activity",
  nonoccupational = "Non-occupational physical activity"
)
required <- c("ID", "Group", "Match", "endtime5", "A2", base_covars,
              unname(pa_exposures), metabolite_clean)
missing_count <- colSums(is.na(lag1_dat[, required, drop = FALSE]))
if (any(missing_count > 0)) {
  stop("Missing values remain in analysis variables: ",
       paste(names(missing_count)[missing_count > 0], collapse = ", "))
}

sample_flow <- data.frame(
  Stage = c("Original nested case-control sample",
            "Matched sets removed because the case was diagnosed within 1 year",
            "Primary lag-1 analysis sample"),
  Participants = c(nrow(dat), 2 * length(removed_matches), nrow(lag1_dat)),
  CRC_cases = c(sum(dat$Group == 1), length(removed_matches), sum(lag1_dat$Group == 1)),
  Controls = c(sum(dat$Group == 0), length(removed_matches), sum(lag1_dat$Group == 0)),
  Matched_sets = c(nlevels(dat$Match), length(removed_matches), nlevels(lag1_dat$Match))
)
write_csv(sample_flow, file.path(dirs["summary"], "00_sample_flow.csv"))
write_csv(data.frame(Removed_match = removed_matches),
          file.path(dirs["summary"], "01_removed_matched_sets.csv"))

# Save code provenance in each analysis directory.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- file.path(.integration_code_root, "03_nested_sample_screening.R")
  code_names <- c(
    baseline_pa = "code_baseline_and_PA.R",
    total = "code_total_PA_metabolomics.R",
    nonoccupational = "code_nonoccupational_PA_metabolomics.R",
    rings = "code_ring_plots.R",
    mediation = "code_mediation_and_mutual_adjustment.R"
  )
  for (nm in names(code_names)) {
    destination <- file.path(dirs[nm], code_names[nm])
    same_file <- identical(
      normalizePath(script_path, winslash = "/", mustWork = TRUE),
      normalizePath(destination, winslash = "/", mustWork = FALSE)
    )
    if (!same_file) file.copy(script_path, destination, overwrite = TRUE)
  }
}

# Baseline table.
baseline_vars <- c(
  "A2", "AGE", "BMI", "education", "smoking5", "drink", "DMTIME", "H1H2",
  "redmeat", "fruit", "veg", unname(pa_exposures)
)
baseline_factor_vars <- c("A2", "education", "smoking5", "drink", "H1H2")
# Report continuous variables as mean (SD).
baseline_continuous_vars <- c(
  "AGE", "BMI", "DMTIME", "redmeat", "fruit", "veg",
  "total_PA_MET_h_day", "occupational_PA_MET_h_day",
  "nonoccupational_PA_MET_h_day"
)
stopifnot(setequal(baseline_continuous_vars,
                   setdiff(baseline_vars, baseline_factor_vars)))
# Retain all five smoking categories, including empty levels.
baseline_dat <- lag1_dat
baseline_dat$smoking5 <- factor(as.character(baseline_dat$smoking5), levels = smoking_levels)
baseline_obj <- CreateTableOne(
  vars = baseline_vars, strata = "Group", data = baseline_dat,
  factorVars = baseline_factor_vars, addOverall = TRUE
)
baseline_print <- print(
  baseline_obj, exact = NULL,
  quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = TRUE,
  smd = FALSE
)
baseline_table <- data.frame(Variable = rownames(baseline_print), baseline_print,
                             row.names = NULL, check.names = FALSE)
# Tests: Wilcoxon rank-sum for continuous variables; Pearson chi-square for categories.
format_baseline_p <- function(p) {
  if (!is.finite(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}
for (v in baseline_continuous_vars) {
  target_row <- which(baseline_table$Variable == paste0(v, " (mean (SD))"))
  if (length(target_row) == 1L) {
    p_wilcox <- wilcox.test(
      baseline_dat[[v]] ~ baseline_dat$Group,
      exact = FALSE
    )$p.value
    baseline_table$p[target_row] <- format_baseline_p(p_wilcox)
    baseline_table$test[target_row] <- "nonnorm"
  }
}
# Restore empty smoking levels omitted by tableone.
smoke_start <- which(baseline_table$Variable == "smoking5 (%)")
if (length(smoke_start) == 1L) {
  later_nonblank <- which(seq_len(nrow(baseline_table)) > smoke_start & baseline_table$Variable != "")
  smoke_end <- if (length(later_nonblank)) min(later_nonblank) - 1L else nrow(baseline_table)
  old_p <- baseline_table$p[smoke_start]
  old_test <- baseline_table$test[smoke_start]
  fmt_n_pct <- function(x, denom) sprintf("%d (%.1f)", x, 100 * x / denom)
  smoke_rows <- do.call(rbind, lapply(seq_along(smoking_levels), function(i) {
    lev <- smoking_levels[i]
    n_all <- sum(as.character(baseline_dat$smoking5) == lev)
    n0 <- sum(as.character(baseline_dat$smoking5[baseline_dat$Group == 0]) == lev)
    n1 <- sum(as.character(baseline_dat$smoking5[baseline_dat$Group == 1]) == lev)
    data.frame(
      Variable = if (i == 1L) "smoking5 (%)" else "",
      level = lev,
      Overall = fmt_n_pct(n_all, nrow(baseline_dat)),
      `0` = fmt_n_pct(n0, sum(baseline_dat$Group == 0)),
      `1` = fmt_n_pct(n1, sum(baseline_dat$Group == 1)),
      p = if (i == 1L) old_p else "",
      test = if (i == 1L) old_test else "",
      check.names = FALSE
    )
  }))
  before <- if (smoke_start > 1L) baseline_table[seq_len(smoke_start - 1L), , drop = FALSE] else baseline_table[0, ]
  after <- if (smoke_end < nrow(baseline_table)) baseline_table[(smoke_end + 1L):nrow(baseline_table), , drop = FALSE] else baseline_table[0, ]
  baseline_table <- rbind(before, smoke_rows, after)
}
names(baseline_table)[names(baseline_table) == "0"] <- "Controls"
names(baseline_table)[names(baseline_table) == "1"] <- "CRC cases"
write_csv(baseline_table, file.path(dirs["baseline_pa"], "01_baseline_table_by_CRC.csv"))
baseline_missing <- data.frame(
  Variable = baseline_vars,
  Missing_N = vapply(lag1_dat[, baseline_vars, drop = FALSE], function(x) sum(is.na(x)), numeric(1)),
  Missing_percent = vapply(lag1_dat[, baseline_vars, drop = FALSE], function(x) mean(is.na(x)) * 100, numeric(1))
)
write_csv(baseline_missing, file.path(dirs["baseline_pa"], "02_baseline_missingness.csv"))

# Allow baseline-only execution.
if (identical(Sys.getenv("JBPCD_BASELINE_ONLY"), "1")) {
  message("Baseline table regenerated with Wilcoxon rank-sum and chi-square tests.")
  quit(save = "no", status = 0)
}

# PA associations: tertiles, per SD, and splines.
make_tertiles <- function(x) {
  cuts <- unname(quantile(x, c(1 / 3, 2 / 3), type = 7, na.rm = TRUE))
  if (cuts[1] >= cuts[2]) stop("Two distinct tertile cutpoints could not be defined.")
  g <- cut(x, breaks = c(-Inf, cuts, Inf), labels = c("T1", "T2", "T3"),
           include.lowest = TRUE, right = TRUE)
  ranges <- c(
    sprintf("[%.2f, %.2f]", min(x), cuts[1]),
    sprintf("(%.2f, %.2f]", cuts[1], cuts[2]),
    sprintf("(%.2f, %.2f]", cuts[2], max(x))
  )
  list(group = g, cuts = cuts, ranges = ranges)
}

pa_groups <- lapply(pa_exposures, function(v) make_tertiles(lag1_dat[[v]]))
for (nm in names(pa_exposures)) lag1_dat[[paste0(nm, "_tertile")]] <- pa_groups[[nm]]$group

fit_pa_row <- function(nm, adjusted) {
  exposure <- pa_exposures[[nm]]
  group_var <- paste0(nm, "_tertile")
  lag1_dat[[paste0(nm, "_z")]] <- as.numeric(scale(lag1_dat[[exposure]]))
  zvar <- paste0(nm, "_z")
  covars <- if (adjusted) base_covars else character()
  f_sd <- as.formula(paste("Group ~", paste(c(zvar, covars, "strata(Match)"), collapse = " + ")))
  f_t <- as.formula(paste("Group ~", paste(c(group_var, covars, "strata(Match)"), collapse = " + ")))
  m_sd <- fit_clogit(f_sd, lag1_dat)
  m_t <- fit_clogit(f_t, lag1_dat)
  e_sd <- extract_or(m_sd, zvar)
  e_t2 <- extract_or(m_t, paste0(group_var, "T2"))
  e_t3 <- extract_or(m_t, paste0(group_var, "T3"))
  data.frame(
    Variable = pa_labels[[nm]], Model = if (adjusted) "Fully adjusted" else "Unadjusted",
    T1 = "1.00 (ref.)", T2 = fmt_est(e_t2["OR"], e_t2["Lower"], e_t2["Upper"]),
    T3 = fmt_est(e_t3["OR"], e_t3["Lower"], e_t3["Upper"]),
    `Per 1-SD` = fmt_est(e_sd["OR"], e_sd["Lower"], e_sd["Upper"]),
    `P for per 1-SD` = unname(e_sd["P"]), SD_MET_h_day = sd(lag1_dat[[exposure]]),
    N = nrow(lag1_dat), Matched_sets = nlevels(lag1_dat$Match), check.names = FALSE
  )
}
pa_results <- do.call(rbind, lapply(names(pa_exposures), function(nm) {
  rbind(fit_pa_row(nm, FALSE), fit_pa_row(nm, TRUE))
}))
rownames(pa_results) <- NULL
write_csv(pa_results, file.path(dirs["baseline_pa"], "03_PA_CRC_tertiles_and_perSD.csv"))

pa_cutpoints <- do.call(rbind, lapply(names(pa_exposures), function(nm) {
  z <- pa_groups[[nm]]
  data.frame(
    Variable = pa_labels[[nm]], Unit = "MET-h/day", T1 = z$ranges[1], T2 = z$ranges[2],
    T3 = z$ranges[3], Cut1 = z$cuts[1], Cut2 = z$cuts[2],
    SD = sd(lag1_dat[[pa_exposures[[nm]]]]), Minimum = min(lag1_dat[[pa_exposures[[nm]]]]),
    Maximum = max(lag1_dat[[pa_exposures[[nm]]]])
  )
}))
rownames(pa_cutpoints) <- NULL
write_csv(pa_cutpoints, file.path(dirs["baseline_pa"], "04_PA_tertile_cutpoints_and_SD.csv"))

run_rcs <- function(data, exposure, label, output_prefix, out_dir,
                    linear_p = NA_real_, unit = "MET-h/day") {
  rcs_data <- droplevels(data[, c("Group", exposure, base_covars), drop = FALSE])
  dd_name <- paste0("dd_", gsub("[^A-Za-z0-9]", "_", output_prefix))
  dd <- datadist(rcs_data)
  dd$limits["Adjust to", exposure] <- min(rcs_data[[exposure]])
  assign(dd_name, dd, envir = .GlobalEnv)
  old_dd <- getOption("datadist")
  options(datadist = dd_name)
  on.exit({options(datadist = old_dd); rm(list = dd_name, envir = .GlobalEnv)}, add = TRUE)
  form <- as.formula(paste(
    "Group ~ rms::rcs(", exposure, ", 3) +", paste(base_covars, collapse = " + ")
  ))
  fit <- suppressWarnings(rms::lrm(form, data = rcs_data, x = TRUE, y = TRUE))
  av <- anova(fit)
  pcol <- grep("^P$", colnames(av), value = TRUE)[1]
  erow <- grep(paste0("^", exposure, "$"), trimws(rownames(av)))[1]
  nrow_idx <- grep("Nonlinear", trimws(rownames(av)), ignore.case = TRUE)[1]
  p_overall <- if (length(erow) && is.finite(erow)) as.numeric(av[erow, pcol]) else NA_real_
  p_nonlinear <- if (length(nrow_idx) && is.finite(nrow_idx)) as.numeric(av[nrow_idx, pcol]) else NA_real_
  xpred <- seq(min(rcs_data[[exposure]]), max(rcs_data[[exposure]]), length.out = 300)
  args <- c(list(fit), setNames(list(xpred), exposure), list(fun = exp, ref.zero = TRUE))
  pred <- as.data.frame(do.call(rms::Predict, args))
  pred_out <- data.frame(Exposure_value = pred[[exposure]], OR = pred$yhat,
                         Lower = pred$lower, Upper = pred$upper)
  p <- ggplot(pred_out, aes(Exposure_value, OR)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#5B9BD5", alpha = 0.25) +
    geom_line(color = "#2C78B8", linewidth = 0.9) +
    geom_hline(yintercept = 1, linetype = 2, linewidth = 0.35) +
    theme_classic(base_family = "Arial", base_size = 11) +
    labs(x = paste0(label, " (", unit, ")"),
         y = "Odds ratio") +
    annotate("text", x = quantile(xpred, 0.25), y = Inf,
             label = paste0("P for linear: ", fmt_p(linear_p)), vjust = 1.2, size = 3.4, family = "Arial") +
    annotate("text", x = quantile(xpred, 0.25), y = Inf,
             label = paste0("P for non-linear: ", fmt_p(p_nonlinear)), vjust = 3.0, size = 3.4, family = "Arial")
  ggsave(file.path(out_dir, paste0(output_prefix, "_RCS.pdf")), p, width = 6, height = 5.5,
         device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(output_prefix, "_RCS.png")), p, width = 6, height = 5.5, dpi = 600)
  write_csv(pred_out, file.path(out_dir, paste0(output_prefix, "_RCS_predictions.csv")))
  stat <- data.frame(
    Test = c("Linear association (per 1-SD)", "Overall spline association", "Non-linear association"),
    P_value = c(linear_p, p_overall, p_nonlinear)
  )
  write_csv(stat, file.path(out_dir, paste0(output_prefix, "_RCS_statistics.csv")))
  list(fit = fit, predictions = pred_out, statistics = stat)
}

pa_rcs <- lapply(names(pa_exposures), function(nm) {
  linear_p <- pa_results$`P for per 1-SD`[
    pa_results$Variable == pa_labels[[nm]] & pa_results$Model == "Fully adjusted"
  ][1]
  run_rcs(lag1_dat, pa_exposures[[nm]], pa_labels[[nm]], paste0(nm, "_PA"),
          dirs["baseline_pa"], linear_p = linear_p, unit = "MET-h/day")
})
names(pa_rcs) <- names(pa_exposures)

# Shared metabolomics functions.
annotation <- metarowdata2
annotation$COMPOUND.Name <- gsub("\\.", " ", annotation$COMPOUND.Name)
annotation <- annotation[!duplicated(annotation$COMPOUND.Name), , drop = FALSE]

# Total-PA screening: FDR <0.05; two forced metabolites have zero penalty.
spearman_fdr_threshold <- 0.05
total_pa_forced_metabolites <- c(
  "1-linoleoyl-GPA (18:2)*",
  "creatine"
)

spearman_screen <- function(exposure, out_dir) {
  x <- lag1_dat[[exposure]]
  out <- data.frame(
    COMPOUND.Name = metabolite_clean,
    Estimate = NA_real_, P.value = NA_real_, stringsAsFactors = FALSE
  )
  for (i in seq_along(metabolite_clean)) {
    ct <- suppressWarnings(cor.test(x, lag1_dat[[metabolite_clean[i]]],
                                    method = "spearman", exact = FALSE))
    out$Estimate[i] <- unname(ct$estimate)
    out$P.value[i] <- ct$p.value
  }
  out$p_fdr <- p.adjust(out$P.value, method = "fdr")
  ann_idx <- match(out$COMPOUND.Name, annotation$COMPOUND.Name)
  for (v in setdiff(names(annotation), "COMPOUND.Name")) out[[v]] <- annotation[[v]][ann_idx]
  out$SUPER.META.PATHWAY[is.na(out$SUPER.META.PATHWAY) | out$SUPER.META.PATHWAY == ""] <-
    "Partially Characterized Molecules"
  selected <- out[out$p_fdr < spearman_fdr_threshold, , drop = FALSE]
  write_csv(out, file.path(out_dir, "01_Spearman_all_metabolites.csv"))
  stale_screen <- file.path(out_dir, "02_Spearman_FDR_lt_0.10.csv")
  if (file.exists(stale_screen)) file.remove(stale_screen)
  write_csv(selected, file.path(out_dir, "02_Spearman_FDR_lt_0.05.csv"))
  target <- out[out$COMPOUND.Name == "1-linoleoyl-GPA (18:2)*", , drop = FALSE]
  write_csv(target, file.path(out_dir, "03_target_1-linoleoyl-GPA_screening.csv"))
  forced_screening <- out[match(total_pa_forced_metabolites, out$COMPOUND.Name), , drop = FALSE]
  write_csv(forced_screening, file.path(out_dir, "03b_two_forced_metabolites_screening.csv"))
  list(all = out, selected = selected)
}

# Manhattan plot: pathway groups with positions randomized using a fixed seed.
# Inset: total and FDR-significant counts by pathway.
plot_total_pa_spearman_manhattan <- function(screen_df, out_dir,
                                             random_seed = 20260911L) {
  pathway_order <- c(
    "Lipid", "Amino Acid", "Xenobiotics", "Peptide", "Nucleotide",
    "Cofactors and Vitamins", "Carbohydrate", "Energy",
    "Partially Characterized Molecules"
  )
  pathway_colors <- c(
    "Lipid" = "#A54B37",
    "Amino Acid" = "#FA3100",
    "Xenobiotics" = "#FDAE61",
    "Peptide" = "#E8D300",
    "Nucleotide" = "#4DAF4A",
    "Cofactors and Vitamins" = "#AFCFEA",
    "Carbohydrate" = "#6A3D9A",
    "Energy" = "#37437A",
    "Partially Characterized Molecules" = "#999999"
  )
  cutoff <- spearman_fdr_threshold
  d <- screen_df
  d$SUPER.META.PATHWAY[is.na(d$SUPER.META.PATHWAY) |
                         d$SUPER.META.PATHWAY == ""] <-
    "Partially Characterized Molecules"
  unexpected <- setdiff(unique(d$SUPER.META.PATHWAY), pathway_order)
  if (length(unexpected) > 0L) {
    stop("Unexpected SUPER.META.PATHWAY values: ", paste(unexpected, collapse = ", "))
  }
  d$SUPER.META.PATHWAY <- factor(d$SUPER.META.PATHWAY,
                                 levels = pathway_order, ordered = TRUE)

  set.seed(random_seed)
  d$random_key_within_pathway <- runif(nrow(d))
  d <- d[order(d$SUPER.META.PATHWAY, d$random_key_within_pathway), , drop = FALSE]
  d$plot_pos <- seq_len(nrow(d))
  d$within_pathway_position <- ave(
    d$plot_pos, d$SUPER.META.PATHWAY,
    FUN = function(z) seq_along(z)
  )
  d$minus_log10_FDR_P <- -log10(pmax(d$p_fdr, .Machine$double.xmin))

  pathway_split <- split(d$plot_pos, d$SUPER.META.PATHWAY, drop = TRUE)
  pathway_summary <- data.frame(
    SUPER.META.PATHWAY = names(pathway_split),
    count_all = vapply(pathway_split, length, integer(1)),
    start = vapply(pathway_split, min, numeric(1)),
    end = vapply(pathway_split, max, numeric(1)),
    stringsAsFactors = FALSE
  )
  pathway_summary$mid <- (pathway_summary$start + pathway_summary$end) / 2
  pathway_summary$count_FDR_lt_0.05 <- vapply(
    pathway_summary$SUPER.META.PATHWAY,
    function(pathway) sum(d$p_fdr[as.character(d$SUPER.META.PATHWAY) == pathway] < cutoff,
                           na.rm = TRUE),
    integer(1)
  )
  pathway_summary <- rbind(
    pathway_summary,
    data.frame(
      SUPER.META.PATHWAY = "Total", count_all = nrow(d), start = NA_real_,
      end = NA_real_, mid = NA_real_,
      count_FDR_lt_0.05 = sum(d$p_fdr < cutoff, na.rm = TRUE)
    )
  )

  write_csv(d, file.path(out_dir, "01a_Spearman_manhattan_plot_data.csv"))
  write_csv(pathway_summary,
            file.path(out_dir, "01b_Spearman_manhattan_pathway_summary.csv"))

  pathway_display <- c(
    "Lipid" = "Lipid",
    "Amino Acid" = "Amino acid",
    "Xenobiotics" = "Xenobiotics",
    "Peptide" = "Peptide",
    "Nucleotide" = "Nucleotide",
    "Cofactors and Vitamins" = "Cofactors and vitamins",
    "Carbohydrate" = "Carbohydrate",
    "Energy" = "Energy",
    "Partially Characterized Molecules" = "Unclassified metabolites"
  )
  table_body <- pathway_summary[
    match(pathway_order, pathway_summary$SUPER.META.PATHWAY), , drop = FALSE
  ]
  table_body$display_name <- unname(pathway_display[table_body$SUPER.META.PATHWAY])
  table_body <- rbind(
    table_body,
    transform(
      pathway_summary[pathway_summary$SUPER.META.PATHWAY == "Total", , drop = FALSE],
      display_name = "Total"
    )
  )

  # Use plot coordinates for editable inset text.
  table_left <- 602
  table_right <- 797
  table_top <- 3.96
  table_bottom <- 2.70
  table_row_h <- (table_top - table_bottom) / (nrow(table_body) + 1)
  table_col_edges <- c(table_left, 700, 746, table_right)
  table_hlines <- data.frame(
    x = table_left,
    xend = table_right,
    y = seq(table_bottom, table_top, by = table_row_h),
    yend = seq(table_bottom, table_top, by = table_row_h)
  )
  table_vlines <- data.frame(
    x = table_col_edges,
    xend = table_col_edges,
    y = table_bottom,
    yend = table_top
  )
  table_body$y <- table_top - table_row_h * (seq_len(nrow(table_body)) + 0.5)
  table_body$point_color <- unname(pathway_colors[table_body$SUPER.META.PATHWAY])
  table_body$point_color[table_body$SUPER.META.PATHWAY == "Total"] <- NA_character_

  p <- ggplot(d, aes(x = plot_pos, y = minus_log10_FDR_P,
                     color = SUPER.META.PATHWAY)) +
    geom_hline(yintercept = -log10(cutoff), linetype = 2,
               color = "black", linewidth = 0.45) +
    geom_point(size = 2.4, alpha = 0.90) +
    scale_color_manual(values = pathway_colors, drop = FALSE) +
    scale_x_continuous(
      limits = c(0, 800), breaks = NULL, labels = NULL,
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0, 4), breaks = 0:4, expand = c(0, 0)) +
    labs(
      x = "Identified metabolites",
      y = "-log10(FDR P)",
      color = NULL
    ) +
    annotate(
      "rect", xmin = table_left, xmax = table_right,
      ymin = table_bottom, ymax = table_top,
      fill = "white", color = NA
    ) +
    geom_segment(
      data = table_hlines,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE, color = "grey75", linewidth = 0.35
    ) +
    geom_segment(
      data = table_vlines,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE, color = "grey75", linewidth = 0.35
    ) +
    annotate(
      "text", x = table_left + 10,
      y = table_top - table_row_h / 2,
      label = "Super pathway", hjust = 0,
      family = "Arial", size = 3.05, color = "black"
    ) +
    annotate(
      "text", x = mean(table_col_edges[2:3]),
      y = table_top - table_row_h / 2,
      label = "No. of\nmetabolites", lineheight = 0.86,
      family = "Arial", size = 2.70, color = "black"
    ) +
    annotate(
      "text", x = mean(table_col_edges[3:4]),
      y = table_top - table_row_h / 2,
      label = "No. statistically\nsignificant", lineheight = 0.86,
      family = "Arial", size = 2.55, color = "black"
    ) +
    geom_point(
      data = table_body[table_body$SUPER.META.PATHWAY != "Total", , drop = FALSE],
      aes(x = table_left + 6, y = y, color = SUPER.META.PATHWAY),
      inherit.aes = FALSE, size = 3.2, show.legend = FALSE
    ) +
    geom_text(
      data = table_body,
      aes(x = table_left + 11, y = y, label = display_name),
      inherit.aes = FALSE, hjust = 0,
      family = "Arial", size = 2.85, color = "black"
    ) +
    geom_text(
      data = table_body,
      aes(x = mean(table_col_edges[2:3]), y = y, label = count_all),
      inherit.aes = FALSE,
      family = "Arial", size = 2.85, color = "black"
    ) +
    geom_text(
      data = table_body,
      aes(x = mean(table_col_edges[3:4]), y = y,
          label = count_FDR_lt_0.05),
      inherit.aes = FALSE,
      family = "Arial", size = 2.85, color = "black"
    ) +
    guides(color = "none") +
    theme_bw(base_family = "Arial", base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.65),
      axis.text = element_text(family = "Arial", color = "black"),
      axis.title = element_text(family = "Arial", color = "black"),
      legend.position = "none",
      plot.margin = margin(8, 14, 8, 8)
    )

  pdf_out <- file.path(out_dir, "total_PA_Spearman_manhattan_with_pathway_summary.pdf")
  png_out <- file.path(out_dir, "total_PA_Spearman_manhattan_with_pathway_summary.png")
  ggsave(pdf_out, p,
         width = 12, height = 8, units = "in", device = cairo_pdf)
  ggsave(png_out, p,
         width = 12, height = 8, units = "in", dpi = 600, bg = "white")
  # Keep the new PDF if the original is locked.
  suppressWarnings(try(file.copy(
    png_out, file.path(out_dir, "total_PA_Spearman_manhattan.png"), overwrite = TRUE
  ), silent = TRUE))
  suppressWarnings(try(file.copy(
    pdf_out, file.path(out_dir, "total_PA_Spearman_manhattan.pdf"), overwrite = TRUE
  ), silent = TRUE))
  invisible(list(plot = p, data = d, summary = pathway_summary))
}

metabolites_crc <- function(candidates, out_dir, file_prefix = "PA_related") {
  if (length(candidates) == 0L) {
    out <- data.frame(Metabolite = character(), OR = numeric(), CI_lower = numeric(),
                      CI_upper = numeric(), OR_CI = character(), P = numeric(), FDR_P = numeric())
    write_csv(out, file.path(out_dir, paste0("04_", file_prefix, "_metabolites_CRC.csv")))
    return(out)
  }
  out <- data.frame(Metabolite = candidates, OR = NA_real_, CI_lower = NA_real_,
                    CI_upper = NA_real_, OR_CI = "", P = NA_real_)
  for (i in seq_along(candidates)) {
    met <- candidates[i]
    form <- as.formula(paste0(
      "Group ~ scale(`", met, "`) + ", paste(base_covars, collapse = " + "),
      " + strata(Match)"
    ))
    fit <- fit_clogit(form, lag1_dat)
    term <- rownames(summary(fit)$coefficients)[1]
    e <- extract_or(fit, term)
    out$OR[i] <- e["OR"]; out$CI_lower[i] <- e["Lower"]; out$CI_upper[i] <- e["Upper"]
    out$P[i] <- e["P"]; out$OR_CI[i] <- fmt_est(e["OR"], e["Lower"], e["Upper"])
  }
  out$FDR_P <- p.adjust(out$P, method = "fdr")
  ann_idx <- match(out$Metabolite, annotation$COMPOUND.Name)
  out$SUPER.META.PATHWAY <- annotation$SUPER.META.PATHWAY[ann_idx]
  out$SUPER.META.PATHWAY[is.na(out$SUPER.META.PATHWAY) | out$SUPER.META.PATHWAY == ""] <-
    "Partially Characterized Molecules"
  out <- out[order(out$SUPER.META.PATHWAY, out$P, out$OR), ]
  rownames(out) <- NULL
  write_csv(out, file.path(out_dir, paste0("04_", file_prefix, "_metabolites_CRC.csv")))
  write_csv(out[out$P < 0.05, , drop = FALSE],
            file.path(out_dir, paste0("05_", file_prefix, "_CRC_raw_P_lt_0.05.csv")))
  write_csv(out[out$FDR_P < 0.05, , drop = FALSE],
            file.path(out_dir, paste0("06_", file_prefix, "_CRC_FDR_lt_0.05.csv")))
  out
}

plot_ring <- function(crc_results, file_stub) {
  if (nrow(crc_results) == 0L) return(invisible(NULL))
  class_order <- c("Lipid", "Amino Acid", "Xenobiotics", "Peptide", "Nucleotide",
                   "Cofactors and Vitamins", "Carbohydrate", "Energy",
                   "Partially Characterized Molecules")
  class_labels <- c(
    "Lipid" = "Lipid", "Amino Acid" = "Amino acid",
    "Xenobiotics" = "Xenobiotics", "Peptide" = "Peptide",
    "Nucleotide" = "Nucleotide",
    "Cofactors and Vitamins" = "Cofactors and\nvitamins",
    "Carbohydrate" = "Carbohydrate", "Energy" = "Energy",
    "Partially Characterized Molecules" = "Partially characterized\nmolecules"
  )
  class_colors <- c(
    "Lipid" = "#5B8DB8", "Amino Acid" = "#C87587",
    "Xenobiotics" = "#D6A269", "Peptide" = "#B8A46E",
    "Nucleotide" = "#76A98C", "Cofactors and Vitamins" = "#299B99",
    "Carbohydrate" = "#9B8DB8", "Energy" = "#7286A7",
    "Partially Characterized Molecules" = "#A8ADB2"
  )
  association_levels <- c(
    "Positive association (OR > 1)",
    "Negative association (OR < 1)",
    "Non-significant (P >= 0.05)"
  )
  association_colors <- c(
    "Positive association (OR > 1)" = "#D75B4E",
    "Negative association (OR < 1)" = "#3F7FB3",
    "Non-significant (P >= 0.05)" = "#C8CDD2"
  )

  d <- crc_results
  d$Class <- factor(d$SUPER.META.PATHWAY, levels = class_order)
  d <- d[order(d$Class, d$P, d$OR), , drop = FALSE]
  d$class_number <- as.integer(d$Class)
  gap_slots <- 2L
  d$index <- seq_len(nrow(d)) + (d$class_number - 1L) * gap_slots
  total_slots <- max(d$index) + gap_slots
  d$logP <- -log10(pmax(d$P, .Machine$double.xmin))
  d$Association <- ifelse(
    d$P < 0.05 & d$OR > 1, "Positive association (OR > 1)",
    ifelse(d$P < 0.05 & d$OR < 1, "Negative association (OR < 1)",
           "Non-significant (P >= 0.05)")
  )
  d$Association <- factor(d$Association, levels = association_levels)

  # Ring layers: class, -log10(P) by OR direction, and names with P <0.05.
  class_ring_inner <- 2.46
  class_ring_outer <- 2.72
  bar_base <- 2.82
  d$bar_top <- bar_base + d$logP
  p05_height <- bar_base - log10(0.05)

  sector_split <- split(d$index, d$Class, drop = TRUE)
  sectors <- data.frame(
    Class = factor(names(sector_split), levels = class_order),
    xmin = vapply(sector_split, min, numeric(1)) - 0.44,
    xmax = vapply(sector_split, max, numeric(1)) + 0.44,
    mid = vapply(sector_split, function(z) mean(range(z)), numeric(1)),
    row.names = NULL
  )

  sig <- d[d$P < 0.05, , drop = FALSE]
  if (nrow(sig) > 0L) {
    sig <- sig[order(sig$index), , drop = FALSE]
    circular_gap <- c(Inf, diff(sig$index))
    sig$label_cluster <- cumsum(circular_gap > 3)
    sig$label_index <- sig$index
    sig$label_radius <- sig$bar_top + 0.16
    for (cluster_id in unique(sig$label_cluster)) {
      ii <- which(sig$label_cluster == cluster_id)
      cluster_size <- length(ii)
      if (cluster_size == 2L) {
        ordered_ii <- ii[order(sig$label_radius[ii], sig$index[ii])]
        sig$label_radius[ordered_ii[2]] <- max(
          sig$label_radius[ordered_ii[2]], sig$label_radius[ordered_ii[1]] + 1.35
        )
      } else if (cluster_size == 3L) {
        ordered_ii <- ii[order(sig$label_radius[ii], sig$index[ii])]
        for (j in 2:cluster_size) {
          previous <- ordered_ii[j - 1L]
          current <- ordered_ii[j]
          sig$label_radius[current] <- max(
            sig$label_radius[current], sig$label_radius[previous] + 0.60
          )
        }
      } else if (cluster_size >= 4L) {
        angular_order <- ii[order(sig$index[ii])]
        outer_track <- angular_order[seq(1L, cluster_size, by = 2L)]
        sig$label_radius[outer_track] <- max(sig$label_radius[outer_track])
      }
    }
    raw_angle <- 90 - 360 * (sig$index - 0.5) / total_slots
    sig$label_angle <- ifelse(raw_angle < -90, raw_angle + 180, raw_angle)
    sig$label_hjust <- ifelse(raw_angle < -90, 1, 0)
  }

  max_logp <- max(c(d$logP, -log10(0.05)), na.rm = TRUE)
  outer_limit <- max(c(d$bar_top, if (nrow(sig) > 0L) sig$label_radius else NA_real_),
                     na.rm = TRUE) + 0.80
  grid_values <- seq_len(max(1L, ceiling(max_logp)))
  grid_values <- grid_values[bar_base + grid_values <= outer_limit]
  grid_data <- data.frame(y = bar_base + grid_values)
  grid_label_data <- data.frame(
    x = 0.5 + total_slots / 4,
    y = bar_base + grid_values,
    label = as.character(grid_values)
  )

  p <- ggplot() +
    geom_hline(
      data = grid_data, aes(yintercept = y),
      linewidth = 0.28, color = "#E3E6E8"
    ) +
    geom_rect(
      data = sectors,
      aes(xmin = xmin, xmax = xmax, ymin = class_ring_inner,
          ymax = class_ring_outer, fill = Class),
      color = "white", linewidth = 0.45, inherit.aes = FALSE
    ) +
    geom_hline(yintercept = bar_base, linewidth = 0.30, color = "#B8BDC2") +
    geom_segment(
      data = d,
      aes(x = index, xend = index, y = bar_base, yend = bar_top,
          color = Association),
      linewidth = 1.55, lineend = "butt"
    ) +
    geom_hline(
      yintercept = p05_height, linetype = 2,
      linewidth = 0.55, color = "#D75B4E"
    ) +
    geom_label(
      data = grid_label_data,
      aes(x = x, y = y, label = label),
      family = "Arial", size = 3.0, fontface = "plain",
      color = "#686D72", fill = "white", linewidth = 0,
      label.padding = grid::unit(0.06, "lines"),
      label.r = grid::unit(0, "lines"),
      inherit.aes = FALSE
    ) +
    {if (nrow(sig) > 0L) geom_text(
      data = sig,
      aes(index, label_radius, label = Metabolite,
          angle = label_angle, hjust = label_hjust),
      size = 2.65, family = "Arial", color = "#202020",
      inherit.aes = FALSE
    )} +
    coord_polar(start = 0, clip = "off") +
    scale_x_continuous(limits = c(0.5, total_slots + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, outer_limit), expand = c(0, 0)) +
    scale_color_manual(values = association_colors, drop = FALSE, guide = "none") +
    scale_fill_manual(values = class_colors, labels = class_labels, drop = TRUE) +
    guides(
      fill = guide_legend(order = 1, ncol = 2, byrow = TRUE,
                          override.aes = list(color = NA))
    ) +
    labs(
      fill = "Metabolite class",
      caption = "Bars show -log10(P); dashed circle indicates P = 0.05."
    ) +
    theme_void(base_family = "Arial") +
    theme(
      legend.position = c(0.50, 0.50),
      legend.justification = c(0.50, 0.50),
      legend.box = "vertical",
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 11.5, face = "bold", hjust = 0.5,
                                  family = "Arial"),
      legend.text = element_text(size = 9.5, lineheight = 0.95,
                                 family = "Arial"),
      legend.spacing.y = grid::unit(0.03, "cm"),
      legend.key.height = grid::unit(0.40, "cm"),
      legend.key.width = grid::unit(0.47, "cm"),
      legend.margin = margin(1, 1, 1, 1),
      plot.caption = element_text(size = 8.0, hjust = 0.5, color = "#4B4B4B",
                                  family = "Arial", margin = margin(t = 4)),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(28, 34, 16, 34)
    )

  direction_key <- grid::grobTree(
    grid::textGrob("Effect direction", x = 0.79, y = 0.94, just = "left",
                   gp = grid::gpar(fontfamily = "Arial", fontsize = 14,
                                   fontface = "bold")),
    grid::segmentsGrob(x0 = 0.79, x1 = 0.83, y0 = 0.90, y1 = 0.90,
                       gp = grid::gpar(
                         col = association_colors["Negative association (OR < 1)"],
                         lwd = 5)),
    grid::textGrob("OR < 1", x = 0.85, y = 0.90, just = "left",
                   gp = grid::gpar(fontfamily = "Arial", fontsize = 12)),
    grid::segmentsGrob(x0 = 0.79, x1 = 0.83, y0 = 0.865, y1 = 0.865,
                       gp = grid::gpar(
                         col = association_colors["Positive association (OR > 1)"],
                         lwd = 5)),
    grid::textGrob("OR > 1", x = 0.85, y = 0.865, just = "left",
                   gp = grid::gpar(fontfamily = "Arial", fontsize = 12))
  )
  p <- ggplotGrob(p)
  p <- gtable::gtable_add_grob(
    p, direction_key, t = 1, l = 1,
    b = nrow(p), r = ncol(p), clip = "off"
  )

  pdf_path <- file.path(dirs["rings"], paste0(file_stub, ".pdf"))
  png_path <- file.path(dirs["rings"], paste0(file_stub, ".png"))
  ggsave(pdf_path, p, width = 14, height = 13.2,
         limitsize = FALSE, device = cairo_pdf, bg = "white")
  ggsave(png_path, p, width = 14, height = 13.2,
         dpi = 600, limitsize = FALSE, bg = "white")
  write_csv(d, file.path(dirs["rings"], paste0(file_stub, "_plot_data.csv")))
  write_csv(sectors, file.path(dirs["rings"], paste0(file_stub, "_sector_positions.csv")))
  file.copy(pdf_path,
            file.path(dirs["rings"], paste0(file_stub, "_revised.pdf")),
            overwrite = TRUE)
  invisible(p)
}

run_lasso_signature <- function(exposure, exposure_label, screen, crc_results, out_dir, prefix,
                                forced_metabolites = character()) {
  stale_signature_files <- c(
    "08_single_metabolite_CRC_result.csv", "08_LASSO_selected_metabolites.csv",
    "09_LASSO_cross_validation_curve.pdf", "09_LASSO_cross_validation_curve.png",
    "10_LASSO_coefficient_paths.pdf", "10_LASSO_coefficient_paths.png",
    "11_LASSO_nonpenalized_covariates.csv", "12_LASSO_lambda.csv",
    "13_signature_scores.csv", "14_signature_cutpoints_and_SD.csv",
    "15_signature_CRC_tertiles_and_perSD.csv", "16_signature_metabolites_CRC.csv",
    "17_LASSO_and_CRC_models.rds",
    paste0(prefix, "_signature_RCS.pdf"), paste0(prefix, "_signature_RCS.png"),
    paste0(prefix, "_signature_RCS_predictions.csv"),
    paste0(prefix, "_signature_RCS_tests.csv"),
    paste0(prefix, "_signature_RCS_knots.csv")
  )
  stale_signature_paths <- file.path(out_dir, stale_signature_files)
  invisible(file.remove(stale_signature_paths[file.exists(stale_signature_paths)]))

  candidates <- screen$selected$COMPOUND.Name
  missing_forced <- setdiff(forced_metabolites, candidates)
  if (length(missing_forced) > 0L) {
    stop("Forced metabolites are absent from the PA-related candidate set: ",
         paste(missing_forced, collapse = ", "))
  }
  status <- data.frame(
    Exposure = exposure_label,
    N_PA_related_metabolites = length(candidates),
    Spearman_FDR_threshold = spearman_fdr_threshold,
    Forced_metabolites = if (length(forced_metabolites))
      paste(forced_metabolites, collapse = "; ") else "None",
    N_forced_metabolites = length(forced_metabolites),
    All_forced_in_PA_related_set = length(missing_forced) == 0L,
    N_LASSO_selected_metabolites = NA_integer_,
    All_forced_retained = NA,
    Signature_constructed = FALSE,
    Reason = ""
  )
  if (length(candidates) <= 1L) {
    status$N_LASSO_selected_metabolites <- length(candidates)
    status$Reason <- if (length(candidates) == 0L) "No PA-related metabolite after FDR < 0.05 screening" else
      "Only one PA-related metabolite; single-metabolite CRC analysis retained"
    write_csv(status, file.path(out_dir, "07_signature_status.csv"))
    if (length(candidates) == 1L) {
      write_csv(crc_results[crc_results$Metabolite == candidates, , drop = FALSE],
                file.path(out_dir, "08_single_metabolite_CRC_result.csv"))
    }
    return(list(status = status, signature = NULL, selected = candidates,
                lasso = NULL, crc = NULL, rcs = NULL, constituent_crc = crc_results))
  }

  x_met <- as.matrix(lag1_dat[, candidates, drop = FALSE])
  x_unpen <- model.matrix(~ AGE + BMI + A2 + smoking5, data = lag1_dat)[, -1, drop = FALSE]
  x <- cbind(x_met, x_unpen)
  metabolite_penalty <- ifelse(colnames(x_met) %in% forced_metabolites, 0, 1)
  penalty <- c(metabolite_penalty, rep(0, ncol(x_unpen)))
  y <- lag1_dat[[exposure]]
  set.seed(222)
  cvfit <- cv.glmnet(x, y, alpha = 1, family = "gaussian", penalty.factor = penalty,
                     nfolds = 10)
  final_fit <- glmnet(x, y, alpha = 1, family = "gaussian", lambda = cvfit$lambda.min,
                      penalty.factor = penalty)
  all_coef <- as.matrix(coef(final_fit))
  selected <- intersect(rownames(all_coef)[all_coef[, 1] != 0], candidates)
  status$N_LASSO_selected_metabolites <- length(selected)
  status$All_forced_retained <- all(forced_metabolites %in% selected)
  if (length(selected) <= 1L) {
    status$Reason <- "LASSO retained one or no metabolite; no multimetabolite signature constructed"
    write_csv(status, file.path(out_dir, "07_signature_status.csv"))
    selected_crc <- crc_results[crc_results$Metabolite %in% selected, , drop = FALSE]
    write_csv(selected_crc, file.path(out_dir, "08_single_metabolite_CRC_result.csv"))
    return(list(status = status, signature = NULL, selected = selected,
                lasso = cvfit, crc = NULL, rcs = NULL, constituent_crc = selected_crc))
  }

  status$Signature_constructed <- TRUE
  status$Reason <- "Multimetabolite signature constructed"
  coef_table <- data.frame(Metabolite = selected,
                           LASSO_coefficient = as.numeric(all_coef[selected, 1]),
                           Penalty_factor = metabolite_penalty[match(selected, colnames(x_met))])
  ann_idx <- match(coef_table$Metabolite, annotation$COMPOUND.Name)
  coef_table$HMDB.ID <- annotation$HMDB.ID[ann_idx]
  coef_table$KEGG.ID <- annotation$KEGG.ID[ann_idx]
  coef_table$SUPER.META.PATHWAY <- annotation$SUPER.META.PATHWAY[ann_idx]
  coef_table$SUB.META.PATHWAY <- annotation$SUB.META.PATHWAY[ann_idx]
  coef_table <- coef_table[, c(
    "Metabolite", "HMDB.ID", "KEGG.ID", "SUPER.META.PATHWAY", "SUB.META.PATHWAY",
    "LASSO_coefficient", "Penalty_factor"
  )]
  nonpen_table <- data.frame(Feature = colnames(x_unpen),
                             Coefficient = as.numeric(all_coef[colnames(x_unpen), 1]),
                             Penalty_factor = 0)
  signature <- as.numeric(as.matrix(lag1_dat[, selected, drop = FALSE]) %*% all_coef[selected, 1])
  score_name <- paste0(prefix, "_signature_score")
  score_dat <- lag1_dat
  score_dat[[score_name]] <- signature
  score_dat[[paste0(prefix, "_signature_z")]] <- as.numeric(scale(signature))
  zname <- paste0(prefix, "_signature_z")
  tert <- make_tertiles(signature)
  tname <- paste0(prefix, "_signature_tertile")
  score_dat[[tname]] <- tert$group

  f_sd0 <- as.formula(paste("Group ~", zname, "+ strata(Match)"))
  f_sd1 <- as.formula(paste("Group ~", paste(c(zname, base_covars, "strata(Match)"), collapse = " + ")))
  f_t0 <- as.formula(paste("Group ~", tname, "+ strata(Match)"))
  f_t1 <- as.formula(paste("Group ~", paste(c(tname, base_covars, "strata(Match)"), collapse = " + ")))
  m_sd0 <- fit_clogit(f_sd0, score_dat); m_sd1 <- fit_clogit(f_sd1, score_dat)
  m_t0 <- fit_clogit(f_t0, score_dat); m_t1 <- fit_clogit(f_t1, score_dat)
  build_sig_row <- function(label, msd, mt) {
    esd <- extract_or(msd, zname)
    et2 <- extract_or(mt, paste0(tname, "T2")); et3 <- extract_or(mt, paste0(tname, "T3"))
    data.frame(
      Variable = exposure_label, Model = label, T1 = "1.00 (ref.)",
      T2 = fmt_est(et2["OR"], et2["Lower"], et2["Upper"]),
      T3 = fmt_est(et3["OR"], et3["Lower"], et3["Upper"]),
      `Per 1-SD` = fmt_est(esd["OR"], esd["Lower"], esd["Upper"]),
      `P for per 1-SD` = esd["P"], SD = sd(signature), N = nrow(score_dat),
      check.names = FALSE
    )
  }
  sig_crc <- rbind(build_sig_row("Unadjusted", m_sd0, m_t0),
                   build_sig_row("Fully adjusted", m_sd1, m_t1))

  constituent_crc <- metabolites_crc(selected, out_dir, "signature")
  names(score_dat)[names(score_dat) == score_name] <- "signature_score"
  sig_linear_p <- sig_crc$`P for per 1-SD`[sig_crc$Model == "Fully adjusted"][1]
  sig_label <- if (prefix == "total_PA") {
    "Physical activity metabolic signature"
  } else {
    "Non-occupational physical activity metabolic signature"
  }
  sig_rcs <- run_rcs(score_dat, "signature_score", sig_label,
                     paste0(prefix, "_signature"), out_dir,
                     linear_p = sig_linear_p, unit = "score")

  pdf(file.path(out_dir, "09_LASSO_cross_validation_curve.pdf"), width = 7, height = 6)
  plot(cvfit); dev.off()
  png(file.path(out_dir, "09_LASSO_cross_validation_curve.png"), width = 4200, height = 3600, res = 600)
  plot(cvfit); dev.off()
  pdf(file.path(out_dir, "10_LASSO_coefficient_paths.pdf"), width = 7, height = 6)
  plot(cvfit$glmnet.fit); dev.off()
  png(file.path(out_dir, "10_LASSO_coefficient_paths.png"), width = 4200, height = 3600, res = 600)
  plot(cvfit$glmnet.fit); dev.off()

  write_csv(status, file.path(out_dir, "07_signature_status.csv"))
  write_csv(coef_table, file.path(out_dir, "08_LASSO_selected_metabolites.csv"))
  write_csv(nonpen_table, file.path(out_dir, "11_LASSO_nonpenalized_covariates.csv"))
  write_csv(data.frame(Lambda_min = cvfit$lambda.min, Lambda_1se = cvfit$lambda.1se),
            file.path(out_dir, "12_LASSO_lambda.csv"))
  write_csv(data.frame(ID = score_dat$ID, Group = score_dat$Group, Match = score_dat$Match,
                       Signature_score = signature), file.path(out_dir, "13_signature_scores.csv"))
  write_csv(data.frame(Variable = exposure_label, T1 = tert$ranges[1], T2 = tert$ranges[2],
                       T3 = tert$ranges[3], Cut1 = tert$cuts[1], Cut2 = tert$cuts[2],
                       SD = sd(signature)), file.path(out_dir, "14_signature_cutpoints_and_SD.csv"))
  write_csv(sig_crc, file.path(out_dir, "15_signature_CRC_tertiles_and_perSD.csv"))
  write_csv(constituent_crc, file.path(out_dir, "16_signature_metabolites_CRC.csv"))
  saveRDS(list(cvfit = cvfit, final_fit = final_fit, selected = selected,
               signature = signature, signature_crc = list(m_sd0, m_sd1, m_t0, m_t1)),
          file.path(out_dir, "17_LASSO_and_CRC_models.rds"))
  list(status = status, signature = signature, score_data = score_dat, selected = selected,
       lasso = cvfit, crc = sig_crc, rcs = sig_rcs, constituent_crc = constituent_crc,
       score_name = "signature_score", tertile = tert)
}

# Total and non-occupational PA metabolomics.
if (identical(Sys.getenv("JBPCD_TOTAL_PA_MANHATTAN_ONLY"), "1")) {
  existing_total_screen <- read.csv(
    file.path(dirs["total"], "01_Spearman_all_metabolites.csv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  total_manhattan <- plot_total_pa_spearman_manhattan(
    existing_total_screen, dirs["total"], random_seed = 20260911L
  )
  message("Total-PA Spearman Manhattan plot regenerated from the existing screen.")
  quit(save = "no", status = 0)
}
total_screen <- spearman_screen(pa_exposures["total"], dirs["total"])
total_manhattan <- plot_total_pa_spearman_manhattan(
  total_screen$all, dirs["total"], random_seed = 20260911L
)
total_crc <- metabolites_crc(total_screen$selected$COMPOUND.Name, dirs["total"])
plot_ring(total_crc, "total_PA_related_metabolites_CRC_ring")
total_signature <- run_lasso_signature(
  pa_exposures["total"], "Physical activity metabolic signature",
  total_screen, total_crc, dirs["total"], "total_PA",
  forced_metabolites = total_pa_forced_metabolites
)

nonocc_screen <- spearman_screen(pa_exposures["nonoccupational"], dirs["nonoccupational"])
nonocc_crc <- metabolites_crc(nonocc_screen$selected$COMPOUND.Name, dirs["nonoccupational"])
plot_ring(nonocc_crc, "nonoccupational_PA_related_metabolites_CRC_ring")
nonocc_signature <- run_lasso_signature(
  pa_exposures["nonoccupational"], "Non-occupational physical activity metabolic signature",
  nonocc_screen, nonocc_crc, dirs["nonoccupational"], "nonoccupational_PA",
  forced_metabolites = intersect(
    "1-linoleoyl-GPA (18:2)*", nonocc_screen$selected$COMPOUND.Name
  )
)

signature_exposure_correlation <- function(sig_obj, exposure, out_file) {
  if (is.null(sig_obj$signature)) {
    if (file.exists(out_file)) file.remove(out_file)
    return(NULL)
  }
  x <- lag1_dat[[exposure]]
  y <- sig_obj$signature
  keep <- complete.cases(x, y)
  pearson_test <- cor.test(x[keep], y[keep], method = "pearson")
  spearman_test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  pearson_r <- unname(pearson_test$estimate)
  out <- data.frame(
    N = sum(keep),
    Pearson_r = pearson_r,
    Pearson_R_squared = pearson_r^2,
    Pearson_P = pearson_test$p.value,
    Spearman_rho = unname(spearman_test$estimate),
    Spearman_P = spearman_test$p.value,
    R_squared_definition = "Pearson_r squared",
    stringsAsFactors = FALSE
  )
  write_csv(out, out_file)
  out
}

total_signature_pa_correlation <- signature_exposure_correlation(
  total_signature, pa_exposures["total"],
  file.path(dirs["total"], "18_signature_total_PA_correlation.csv")
)

# Mediation and mutually adjusted models.
bootstrap_spearman <- function(x, y, nboot = 1000L, seed = 123L) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  set.seed(seed)
  b <- replicate(nboot, {
    ii <- sample(seq_along(x), replace = TRUE)
    cor(x[ii], y[ii], method = "spearman", use = "complete.obs")
  })
  c(Rho = unname(ct$estimate), Lower = unname(quantile(b, 0.025, na.rm = TRUE)),
    Upper = unname(quantile(b, 0.975, na.rm = TRUE)), P = ct$p.value)
}

mediation_one <- function(exposure, mediator, mediator_label, data) {
  d <- data
  d$exposure_z <- as.numeric(scale(d[[exposure]]))
  d$mediator_z <- as.numeric(scale(mediator))
  total_fit <- fit_clogit(as.formula(paste(
    "Group ~", paste(c("exposure_z", base_covars, "strata(Match)"), collapse = " + ")
  )), d)
  direct_fit <- fit_clogit(as.formula(paste(
    "Group ~", paste(c("exposure_z", "mediator_z", base_covars, "strata(Match)"), collapse = " + ")
  )), d)
  mediator_fit <- fit_clogit(as.formula(paste(
    "Group ~", paste(c("mediator_z", base_covars, "strata(Match)"), collapse = " + ")
  )), d)
  et <- extract_or(total_fit, "exposure_z")
  ed <- extract_or(direct_fit, "exposure_z")
  em <- extract_or(mediator_fit, "mediator_z")
  rho <- bootstrap_spearman(d[[exposure]], mediator)
  prop <- (coef(total_fit)["exposure_z"] - coef(direct_fit)["exposure_z"]) /
    coef(total_fit)["exposure_z"] * 100
  data.frame(
    Mediator = mediator_label,
    Exposure_mediator_Spearman = sprintf("%.2f (%.2f-%.2f)", rho["Rho"], rho["Lower"], rho["Upper"]),
    Exposure_mediator_P = rho["P"],
    Total_effect_exposure = fmt_est(et["OR"], et["Lower"], et["Upper"]),
    Total_effect_P = et["P"],
    Direct_effect_exposure = fmt_est(ed["OR"], ed["Lower"], ed["Upper"]),
    Direct_effect_P = ed["P"],
    Mediator_CRC = fmt_est(em["OR"], em["Lower"], em["Upper"]),
    Mediator_CRC_P = em["P"],
    Mediation_proportion_percent = unname(prop),
    stringsAsFactors = FALSE
  )
}

mutual_adjustment <- function(exposure, exposure_label, sig_obj, prefix) {
  output_file <- file.path(dirs["mediation"], paste0(prefix, "_mutually_adjusted_CRC.csv"))
  if (file.exists(output_file)) file.remove(output_file)
  if (is.null(sig_obj$signature)) return(NULL)
  d <- lag1_dat
  d$signature_score <- sig_obj$signature
  d$exposure_z <- as.numeric(scale(d[[exposure]]))
  d$signature_z <- as.numeric(scale(d$signature_score))
  d$exposure_tertile <- make_tertiles(d[[exposure]])$group
  d$signature_tertile <- make_tertiles(d$signature_score)$group
  sig_sd <- fit_clogit(as.formula(paste(
    "Group ~", paste(c("signature_z", "exposure_z", base_covars, "strata(Match)"), collapse = " + ")
  )), d)
  sig_t <- fit_clogit(as.formula(paste(
    "Group ~", paste(c("signature_tertile", "exposure_z", base_covars, "strata(Match)"), collapse = " + ")
  )), d)
  pa_sd <- fit_clogit(as.formula(paste(
    "Group ~", paste(c("exposure_z", "signature_z", base_covars, "strata(Match)"), collapse = " + ")
  )), d)
  pa_t <- fit_clogit(as.formula(paste(
    "Group ~", paste(c("exposure_tertile", "signature_z", base_covars, "strata(Match)"), collapse = " + ")
  )), d)
  make_row <- function(variable, msd, mt, sdterm, tprefix) {
    esd <- extract_or(msd, sdterm)
    e2 <- extract_or(mt, paste0(tprefix, "T2")); e3 <- extract_or(mt, paste0(tprefix, "T3"))
    data.frame(
      Variable = variable, Model = "Mutually adjusted full model", T1 = "1.00 (ref.)",
      T2 = fmt_est(e2["OR"], e2["Lower"], e2["Upper"]),
      T3 = fmt_est(e3["OR"], e3["Lower"], e3["Upper"]),
      `Per 1-SD` = fmt_est(esd["OR"], esd["Lower"], esd["Upper"]),
      `P for per 1-SD` = esd["P"], check.names = FALSE
    )
  }
  out <- rbind(
    make_row(exposure_label, pa_sd, pa_t, "exposure_z", "exposure_tertile"),
    make_row(paste0(exposure_label, " metabolic signature"), sig_sd, sig_t,
             "signature_z", "signature_tertile")
  )
  write_csv(out, output_file)
  out
}

run_mediation_set <- function(exposure, exposure_label, sig_obj, prefix) {
  output_file <- file.path(dirs["mediation"], paste0(prefix, "_mediation_results.csv"))
  if (file.exists(output_file)) file.remove(output_file)
  if (is.null(sig_obj$signature)) return(NULL)
  out <- mediation_one(exposure, sig_obj$signature,
                       paste0(exposure_label, " metabolic signature"), lag1_dat)
  sig_crc <- sig_obj$constituent_crc
  significant_mets <- sig_crc$Metabolite[sig_crc$P < 0.05]
  if (length(significant_mets)) {
    for (met in significant_mets) {
      out <- rbind(out, mediation_one(exposure, lag1_dat[[met]], met, lag1_dat))
    }
  }
  write_csv(out, output_file)
  out
}

total_mediation <- run_mediation_set(pa_exposures["total"], "Physical activity",
                                     total_signature, "total_PA")
nonocc_mediation <- run_mediation_set(pa_exposures["nonoccupational"],
                                      "Non-occupational physical activity",
                                      nonocc_signature, "nonoccupational_PA")
total_mutual <- mutual_adjustment(pa_exposures["total"], "Physical activity",
                                  total_signature, "total_PA")
nonocc_mutual <- mutual_adjustment(pa_exposures["nonoccupational"],
                                   "Non-occupational physical activity",
                                   nonocc_signature, "nonoccupational_PA")

# Combined summary.
pipeline_summary <- data.frame(
  Analysis = c("Lag-1 matched sample", "Total PA-related metabolites",
               "Total PA LASSO-selected metabolites", "Total PA signature constructed",
               "Non-occupational PA-related metabolites",
               "Non-occupational PA LASSO-selected metabolites",
               "Non-occupational PA signature constructed"),
  Value = c(
    nrow(lag1_dat), nrow(total_screen$selected), length(total_signature$selected),
    total_signature$status$Signature_constructed,
    nrow(nonocc_screen$selected), length(nonocc_signature$selected),
    nonocc_signature$status$Signature_constructed
  ),
  stringsAsFactors = FALSE
)
if (!is.null(total_signature_pa_correlation)) {
  pipeline_summary <- rbind(
    pipeline_summary,
    data.frame(
      Analysis = c("Total PA signature Pearson r", "Total PA signature Pearson R-squared",
                   "Total PA signature Spearman rho"),
      Value = c(total_signature_pa_correlation$Pearson_r,
                total_signature_pa_correlation$Pearson_R_squared,
                total_signature_pa_correlation$Spearman_rho),
      stringsAsFactors = FALSE
    )
  )
}
write_csv(pipeline_summary, file.path(dirs["summary"], "02_pipeline_summary.csv"))
saveRDS(
  list(lag1_dat = lag1_dat, sample_flow = sample_flow, baseline = baseline_table,
       pa_results = pa_results, pa_cutpoints = pa_cutpoints,
       total_screen = total_screen, total_crc = total_crc, total_signature = total_signature,
       nonocc_screen = nonocc_screen, nonocc_crc = nonocc_crc,
       nonocc_signature = nonocc_signature, total_mediation = total_mediation,
       nonocc_mediation = nonocc_mediation, total_mutual = total_mutual,
       nonocc_mutual = nonocc_mutual,
       total_signature_pa_correlation = total_signature_pa_correlation),
  file.path(dirs["summary"], "03_complete_analysis_objects.rds")
)

# Supplementary lag-1 MSEA.
if (identical(Sys.getenv("RUN_TOP3_MSEA_AFTER_MAIN", unset = "1"), "1")) {
  local({

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(patchwork)
})

if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

BiocParallel::register(BiocParallel::SerialParam())
set.seed(20260911)

analysis_root <- .integration_analysis
output_dir <- file.path(analysis_root, "07_MSEA_top3_pathways")
smpdb_dir <- file.path(.integration_data, "pathways")

objects_file <- file.path(
  analysis_root, "06_summary", "03_complete_analysis_objects.rds"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, filename) {
  write.csv(x, file.path(output_dir, filename), row.names = FALSE,
            fileEncoding = "UTF-8")
}

target_pathways <- data.frame(
  SMPDB_ID = c("SMP0000039", "SMP0000020", "SMP0000004"),
  Pathway = c(
    "Glycerolipid Metabolism",
    "Arginine and Proline Metabolism",
    "Glycine and Serine Metabolism"
  ),
  stringsAsFactors = FALSE
)

superseded_files <- c(
  "08_MSEA_total_PA_exploratory_minGSSize2.csv",
  "09_MSEA_CRC_exploratory_minGSSize2.csv"
)
old_pathway_plots <- list.files(
  output_dir, pattern = "^(PA|CRC)_SMP[0-9]+_.*[.]png$", full.names = TRUE
)
invisible(unlink(c(file.path(output_dir, superseded_files), old_pathway_plots)))

make_rank_table <- function(df, statistic_col) {
  keep <- !is.na(df$HMDB.ID) & nzchar(df$HMDB.ID) &
    !is.na(df[[statistic_col]]) & is.finite(df[[statistic_col]])
  out <- df[keep, c("COMPOUND.Name", "HMDB.ID", statistic_col), drop = FALSE]
  out <- out[out$COMPOUND.Name != "N-delta-acetylornithine", , drop = FALSE]
  out <- out[!out$HMDB.ID %in% c("HMDB0240429", "HMDB0010383"), , drop = FALSE]
  if (anyDuplicated(out$HMDB.ID)) {
    stop("Duplicate HMDB IDs remain after applying the prespecified exclusions.")
  }
  names(out)[3] <- "Rank_statistic"
  out <- out[order(out$Rank_statistic, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

run_msea <- function(rank_table, min_size) {
  gene_list <- rank_table$Rank_statistic
  names(gene_list) <- rank_table$HMDB.ID
  gene_list <- sort(gene_list, decreasing = TRUE)
  set.seed(20260911)
  clusterProfiler::GSEA(
    geneList = gene_list,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    minGSSize = min_size,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "none",
    eps = 0,
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea"
  )
}

tidy_gsea <- function(gsea_object, analysis_label, min_size) {
  out <- as.data.frame(gsea_object)
  if (nrow(out) == 0L) {
    return(data.frame(
      Analysis = character(), minGSSize = integer(), SMPDB_ID = character(),
      Pathway = character(), Set_size = integer(), Enrichment_score = numeric(),
      NES = numeric(), P = numeric(), Rank_at_max = integer(),
      Leading_edge = character(), Core_enrichment = character()
    ))
  }
  out <- data.frame(
    Analysis = analysis_label,
    minGSSize = min_size,
    SMPDB_ID = out$ID,
    Pathway = out$Description,
    Set_size = out$setSize,
    Enrichment_score = out$enrichmentScore,
    NES = out$NES,
    P = out$pvalue,
    Rank_at_max = out$rank,
    Leading_edge = out$leading_edge,
    Core_enrichment = out$core_enrichment,
    stringsAsFactors = FALSE
  )
  out$Pathway_order <- match(out$SMPDB_ID, target_pathways$SMPDB_ID)
  out <- out[order(out$Pathway_order), setdiff(names(out), "Pathway_order"),
             drop = FALSE]
  rownames(out) <- NULL
  out
}

if (!file.exists(objects_file)) {
  stop("Required complete-analysis object was not found: ", objects_file)
}
analysis_objects <- readRDS(objects_file)
lag1_dat <- analysis_objects$lag1_dat
pa_screen <- analysis_objects$total_screen$all

if (nrow(lag1_dat) != 256L || nlevels(factor(lag1_dat$Match)) != 128L) {
  stop("The saved analysis object is not the expected 128-pair lag-1 sample.")
}
if (nrow(pa_screen) != 799L) {
  stop("The saved total-PA Spearman screen does not contain 799 metabolites.")
}

base_covars <- c(
  "AGE", "BMI", "education", "smoking5", "drink", "DMTIME", "H1H2",
  "redmeat", "fruit", "veg"
)
required <- c("Group", "Match", base_covars, pa_screen$COMPOUND.Name)
missing_names <- setdiff(required, names(lag1_dat))
if (length(missing_names) > 0L) {
  stop("Variables missing from lag1_dat: ", paste(missing_names, collapse = ", "))
}
if (anyNA(lag1_dat[, required, drop = FALSE])) {
  stop("Missing values remain in the MSEA analysis variables.")
}

pathway_members <- do.call(rbind, lapply(seq_len(nrow(target_pathways)), function(i) {
  pathway_file <- file.path(
    smpdb_dir, paste0(target_pathways$SMPDB_ID[i], "_metabolites.csv")
  )
  if (!file.exists(pathway_file)) stop("SMPDB pathway file not found: ", pathway_file)
  raw <- read.csv(pathway_file, check.names = FALSE, fileEncoding = "UTF-8")
  required_cols <- c("SMPDB ID", "Pathway Name", "Metabolite Name", "HMDB ID")
  if (!all(required_cols %in% names(raw))) {
    stop("Unexpected SMPDB file structure: ", pathway_file)
  }
  raw <- raw[!is.na(raw[["HMDB ID"]]) & nzchar(raw[["HMDB ID"]]), required_cols]
  raw <- raw[!duplicated(raw[["HMDB ID"]]), , drop = FALSE]
  names(raw) <- c("SMPDB_ID", "Pathway", "SMPDB_metabolite", "HMDB_ID")
  raw
}))

term2gene <- unique(pathway_members[, c("SMPDB_ID", "HMDB_ID")])
names(term2gene) <- c("term", "gene")
term2name <- target_pathways
names(term2name) <- c("term", "name")

pa_rank_input <- pa_screen[, c("COMPOUND.Name", "HMDB.ID", "Estimate")]
pa_rank <- make_rank_table(pa_rank_input, "Estimate")
write_csv(pa_rank, "01_total_PA_ranked_metabolites.csv")

crc_all <- data.frame(
  COMPOUND.Name = pa_screen$COMPOUND.Name,
  logOR_per_SD = NA_real_,
  OR_per_SD = NA_real_,
  CI_lower = NA_real_,
  CI_upper = NA_real_,
  P = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(crc_all))) {
  metabolite <- crc_all$COMPOUND.Name[i]
  model_formula <- as.formula(paste0(
    "Group ~ scale(`", metabolite, "`) + ",
    paste(base_covars, collapse = " + "), " + strata(Match)"
  ))
  fit <- suppressWarnings(clogit(
    model_formula, data = lag1_dat, method = "exact", na.action = na.fail
  ))
  coefficient_table <- summary(fit)$coefficients
  beta <- coefficient_table[1, "coef"]
  se <- coefficient_table[1, "se(coef)"]
  crc_all$logOR_per_SD[i] <- beta
  crc_all$OR_per_SD[i] <- exp(beta)
  crc_all$CI_lower[i] <- exp(beta - 1.96 * se)
  crc_all$CI_upper[i] <- exp(beta + 1.96 * se)
  crc_all$P[i] <- coefficient_table[1, "Pr(>|z|)"]
}

annotation_index <- match(crc_all$COMPOUND.Name, pa_screen$COMPOUND.Name)
crc_all$HMDB.ID <- pa_screen$HMDB.ID[annotation_index]
crc_all$SUPER.META.PATHWAY <- pa_screen$SUPER.META.PATHWAY[annotation_index]
crc_all$FDR_P <- p.adjust(crc_all$P, method = "BH")
write_csv(crc_all, "02_all_metabolites_adjusted_CRC_associations.csv")

crc_rank_input <- crc_all[, c("COMPOUND.Name", "HMDB.ID", "logOR_per_SD")]
crc_rank <- make_rank_table(crc_rank_input, "logOR_per_SD")
write_csv(crc_rank, "03_CRC_ranked_metabolites.csv")

measured_annotation <- unique(pa_rank[, c("COMPOUND.Name", "HMDB.ID")])
names(measured_annotation) <- c("Measured_metabolite", "HMDB_ID")
matched_members <- merge(pathway_members, measured_annotation, by = "HMDB_ID")
matched_members$Pathway_order <- match(matched_members$SMPDB_ID,
                                       target_pathways$SMPDB_ID)
matched_members <- matched_members[order(matched_members$Pathway_order,
                                         matched_members$HMDB_ID),
                                   setdiff(names(matched_members), "Pathway_order")]
write_csv(matched_members, "04_pathway_members_measured_in_JBPCD.csv")

overlap <- merge(
  aggregate(HMDB_ID ~ SMPDB_ID + Pathway, pathway_members,
            function(x) length(unique(x))),
  aggregate(HMDB_ID ~ SMPDB_ID + Pathway, matched_members,
            function(x) length(unique(x))),
  by = c("SMPDB_ID", "Pathway"), all.x = TRUE,
  suffixes = c("_SMPDB", "_measured")
)
names(overlap)[3:4] <- c("SMPDB_members", "Measured_unique_HMDB_members")
overlap$Measured_unique_HMDB_members[is.na(overlap$Measured_unique_HMDB_members)] <- 0L
overlap$Eligible_at_minGSSize_5 <- overlap$Measured_unique_HMDB_members >= 5L
overlap$Pathway_order <- match(overlap$SMPDB_ID, target_pathways$SMPDB_ID)
overlap <- overlap[order(overlap$Pathway_order), setdiff(names(overlap), "Pathway_order")]
write_csv(overlap, "05_pathway_overlap_and_eligibility.csv")

pa_gsea_min5 <- run_msea(pa_rank, min_size = 5)
crc_gsea_min5 <- run_msea(crc_rank, min_size = 5)
pa_res_min5 <- tidy_gsea(pa_gsea_min5, "Total PA (Spearman rho)", 5)
crc_res_min5 <- tidy_gsea(crc_gsea_min5, "CRC (adjusted log OR)", 5)
write_csv(pa_res_min5, "06_MSEA_total_PA_minGSSize5.csv")
write_csv(crc_res_min5, "07_MSEA_CRC_minGSSize5.csv")

summary_long <- rbind(pa_res_min5, crc_res_min5)
summary_long <- merge(
  summary_long,
  overlap[, c("SMPDB_ID", "Measured_unique_HMDB_members",
              "Eligible_at_minGSSize_5")],
  by = "SMPDB_ID", all.x = TRUE
)
summary_long$Direction <- ifelse(
  summary_long$P < 0.05,
  ifelse(
    summary_long$NES > 0,
    ifelse(grepl("Total PA", summary_long$Analysis),
           "Significant positive enrichment with total PA",
           "Significant enrichment toward higher CRC risk"),
    ifelse(grepl("Total PA", summary_long$Analysis),
           "Significant inverse enrichment with total PA",
           "Significant enrichment toward lower CRC risk")
  ),
  ifelse(summary_long$NES > 0,
         "Not statistically significant; positive NES",
         "Not statistically significant; negative NES")
)
summary_long$Pathway_order <- match(summary_long$SMPDB_ID,
                                    target_pathways$SMPDB_ID)
summary_long$Analysis_order <- match(summary_long$Analysis,
  c("Total PA (Spearman rho)", "CRC (adjusted log OR)"))
summary_long <- summary_long[order(summary_long$Pathway_order,
                                   summary_long$Analysis_order),
                             setdiff(names(summary_long),
                                     c("Pathway_order", "Analysis_order"))]
write_csv(summary_long, "10_MSEA_top3_pathways_summary.csv")

save_gsea_set <- function(gsea_object, result_table, analysis_stub) {
  pathway_plots <- vector("list", nrow(target_pathways))
  for (i in seq_len(nrow(target_pathways))) {
    pathway_id <- target_pathways$SMPDB_ID[i]
    if (!pathway_id %in% result_table$SMPDB_ID) next
    pathway_name <- target_pathways$Pathway[i]
    pathway_p <- result_table$P[match(pathway_id, result_table$SMPDB_ID)]
    pathway_nes <- result_table$NES[match(pathway_id, result_table$SMPDB_ID)]
    plot_parts <- enrichplot::gseaplot2(
      gsea_object,
      geneSetID = pathway_id,
      title = pathway_name,
      base_size = 12,
      pvalue_table = FALSE
    )
    plot_parts[[1]] <- plot_parts[[1]] +
      labs(title = paste0("(", LETTERS[i], ") ", pathway_name)) +
      annotate(
        "text", x = Inf, y = Inf,
        label = sprintf("NES: %.2f\nP: %.3f", pathway_nes, pathway_p),
        hjust = 1.08, vjust = 1.10, size = 3.4, family = "Arial",
        lineheight = 1.05
      ) +
      theme(
        text = element_text(family = "Arial"),
        plot.title = element_text(family = "Arial", face = "plain", size = 11,
                                  hjust = 0),
        plot.margin = margin(6, 12, 0, 6)
      )
    plot_parts[[2]] <- plot_parts[[2]] + theme(text = element_text(family = "Arial"))
    plot_parts[[3]] <- plot_parts[[3]] + theme(text = element_text(family = "Arial"))
    full_plot <- patchwork::wrap_plots(
      plot_parts, ncol = 1, heights = c(1.45, 0.45, 1)
    )
    filename <- paste0(
      ifelse(analysis_stub == "Total PA", "PA", "CRC"), "_",
      pathway_id, "_", gsub("[^A-Za-z0-9]+", "_", pathway_name), ".png"
    )
    png(file.path(output_dir, filename), width = 5100, height = 4320,
        res = 600, bg = "white", family = "Arial")
    print(full_plot)
    dev.off()
    pathway_plots[[i]] <- full_plot
  }
  combined_plot <- patchwork::wrap_plots(pathway_plots, nrow = 1)
  combined_stub <- ifelse(analysis_stub == "Total PA", "total_PA", "CRC")
  ggsave(
    file.path(output_dir, paste0("MSEA_", combined_stub, "_three_pathways.pdf")),
    combined_plot, width = 18, height = 6.5, units = "in", device = cairo_pdf
  )
  ggsave(
    file.path(output_dir, paste0("MSEA_", combined_stub, "_three_pathways.png")),
    combined_plot, width = 18, height = 6.5, units = "in", dpi = 600,
    bg = "white"
  )
  invisible(pathway_plots)
}

pa_pathway_plots <- save_gsea_set(pa_gsea_min5, pa_res_min5, "Total PA")
crc_pathway_plots <- save_gsea_set(crc_gsea_min5, crc_res_min5, "CRC")

plot_dat <- summary_long
plot_dat$Analysis_label <- factor(
  plot_dat$Analysis,
  levels = c("Total PA (Spearman rho)", "CRC (adjusted log OR)"),
  labels = c("Total physical activity", "Colorectal cancer risk")
)
plot_dat$Pathway <- factor(
  plot_dat$Pathway,
  levels = rev(target_pathways$Pathway)
)
plot_dat$Significance <- ifelse(plot_dat$P < 0.05, "P < 0.05", "P >= 0.05")

summary_plot <- ggplot(plot_dat, aes(y = Pathway, x = NES)) +
  geom_vline(xintercept = 0, color = "grey45", linetype = 2, linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = NES, yend = Pathway),
               color = "grey70", linewidth = 0.7) +
  geom_point(aes(fill = NES, shape = Significance),
             size = 4.2, color = "black", stroke = 0.45) +
  facet_wrap(~ Analysis_label, nrow = 1) +
  scale_fill_gradient2(low = "#3B73B9", mid = "white", high = "#C84B45",
                       midpoint = 0, name = "NES") +
  scale_shape_manual(values = c("P < 0.05" = 21, "P >= 0.05" = 22),
                     name = NULL) +
  labs(
    x = "Normalized enrichment score (NES)", y = NULL,
    caption = paste0(
      "Targeted MSEA based on all measured metabolites with mapped HMDB IDs ",
      "(minGSSize = 5); significance is based on nominal P values."
    )
  ) +
  theme_bw(base_family = "Arial", base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(family = "Arial", face = "bold"),
    axis.text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(family = "Arial", color = "black"),
    legend.text = element_text(family = "Arial"),
    legend.title = element_text(family = "Arial"),
    plot.caption = element_text(family = "Arial", size = 9, hjust = 0),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(file.path(output_dir, "MSEA_top3_pathways_summary.pdf"), summary_plot,
       width = 11, height = 5.5, units = "in", device = cairo_pdf)
ggsave(file.path(output_dir, "MSEA_top3_pathways_summary.png"), summary_plot,
       width = 11, height = 5.5, units = "in", dpi = 600, bg = "white")

saveRDS(
  list(
    target_pathways = target_pathways,
    pathway_members = pathway_members,
    overlap = overlap,
    pa_rank = pa_rank,
    crc_rank = crc_rank,
    pa_gsea_min5 = pa_gsea_min5,
    crc_gsea_min5 = crc_gsea_min5
  ),
  file.path(output_dir, "11_MSEA_analysis_objects.rds")
)

message("MSEA completed. Results saved to: ", output_dir)
  })
}

# PA histograms from this run's matched sample.
local({
suppressPackageStartupMessages({library(ggplot2); library(gridExtra)})
root <- .integration_analysis
out_dir <- file.path(root, "01_baseline_and_PA")
d <- readRDS(file.path(root, "06_summary", "03_complete_analysis_objects.rds"))$lag1_dat
stopifnot(nrow(d) == 256L, sum(d$Group == 1) == 128L, sum(d$Group == 0) == 128L)
variables <- c(total_PA="total_PA_MET_h_day", occupational_PA="occupational_PA_MET_h_day",
               nonoccupational_PA="nonoccupational_PA_MET_h_day")
labels <- c(total_PA="Total physical activity", occupational_PA="Occupational physical activity",
            nonoccupational_PA="Non-occupational physical activity")
if (.Platform$OS.type == "windows") windowsFonts(Arial=windowsFont("Arial"))

bin_width <- 2
plots <- list(); summaries <- list(); bins <- list()
for (nm in names(variables)) {
  x <- d[[variables[[nm]]]]
  observed <- x[is.finite(x)]
  stopifnot(all(observed >= 0))
  summaries[[nm]] <- data.frame(Variable=labels[[nm]], Unit="MET-h/day", N=length(x),
    Missing=sum(is.na(x)), Zero=sum(observed == 0), Nonzero=sum(observed > 0),
    Zero_percent=mean(observed == 0)*100, Mean=mean(observed), SD=sd(observed),
    Median=median(observed), Q1=unname(quantile(observed,.25)),
    Q3=unname(quantile(observed,.75)), Minimum=min(observed), Maximum=max(observed))
  breaks <- seq(0, (floor(max(observed)/bin_width)+1)*bin_width, by=bin_width)
  h <- hist(observed, breaks=breaks, right=FALSE, plot=FALSE)
  bins[[nm]] <- data.frame(Variable=labels[[nm]], Lower=head(h$breaks,-1),
    Upper=tail(h$breaks,-1), Count=h$counts)
  caption <- sprintf("n = %d; zero PA: %d (%.1f%%)", length(observed),
                     sum(observed == 0), mean(observed == 0)*100)
  plots[[nm]] <- ggplot(data.frame(PA=observed),aes(PA))+
    geom_histogram(breaks=breaks,closed="left",fill="#5B9BD5",color="white",linewidth=.25)+
    scale_x_continuous(limits=c(0,100),breaks=seq(0,100,20),expand=expansion(mult=c(0,.01)))+
    scale_y_continuous(expand=expansion(mult=c(0,.08)))+
    theme_classic(base_family="Arial",base_size=12)+
    labs(title=labels[[nm]],subtitle=caption,x="Physical activity (MET-h/day)",y="Number of participants")+
    theme(plot.title=element_text(size=13),plot.subtitle=element_text(size=10),
          plot.margin=margin(8,12,8,8))
  stopifnot(sum(h$counts) == length(observed), max(observed) <= 100)
  ggsave(file.path(out_dir,paste0("05_",nm,"_histogram.pdf")),plots[[nm]],width=7,height=5,device=cairo_pdf)
  ggsave(file.path(out_dir,paste0("05_",nm,"_histogram.png")),plots[[nm]],width=7,height=5,dpi=600)
}
summary_table <- do.call(rbind,summaries); rownames(summary_table) <- NULL
write.csv(summary_table,file.path(out_dir,"05_PA_distribution_summary.csv"),row.names=FALSE,fileEncoding="UTF-8")
write.csv(do.call(rbind,bins),file.path(out_dir,"05_PA_histogram_bin_counts.csv"),row.names=FALSE,fileEncoding="UTF-8")
combined <- arrangeGrob(grobs=plots,ncol=3)
ggsave(file.path(out_dir,"05_PA_distribution_histograms.pdf"),combined,width=15,height=5,device=cairo_pdf)
ggsave(file.path(out_dir,"05_PA_distribution_histograms.png"),combined,width=15,height=5,dpi=600)
print(summary_table)
})

cat("Completed lag-1 nested case-control pipeline.\n")
print(sample_flow)
print(pa_results)
print(pipeline_summary)
})

# 2. Full-sample derivation and lag-1 CRC analysis.
.integration_run_stage("00_full_nested_derivation_lag1_CRC_pipeline.R",
  .integration_nested, {
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(ggplot2)
})

if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))

output_root <- .integration_nested
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
fmt_est <- function(or, lo, hi) sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
extract_or <- function(fit, term) {
  z <- summary(fit)$coefficients
  if (!term %in% rownames(z)) return(c(OR=NA, Lower=NA, Upper=NA, P=NA))
  b <- z[term, "coef"]; se <- z[term, "se(coef)"]
  c(OR=exp(b), Lower=exp(b-1.96*se), Upper=exp(b+1.96*se), P=z[term,"Pr(>|z|)"])
}
fit_clogit <- function(formula, data) suppressWarnings(
  clogit(formula, data=data, method="exact", na.action=na.fail)
)

# Rebuild cohort and metabolomics data from raw Excel.
source(file.path(.integration_code_root, "01_data_preprocessing.R"),
       local=TRUE, encoding="UTF-8")
stopifnot(exists("data1"), exists("metarowdata_0.2_imputed1"),
          exists("scaled_metarowdata"), exists("metarowdata2"))

smoking_levels <- c("Never", "Former, <20 pack-years", "Former, >=20 pack-years",
                    "Current, <20 pack-years", "Current, >=20 pack-years", "Unknown")
data1$smoking5 <- with(data1, ifelse(
  B1 == 2, "Never",
  ifelse(B3 == 3 & packyear1 < 20, "Former, <20 pack-years",
  ifelse(B3 == 3 & packyear1 >= 20, "Former, >=20 pack-years",
  ifelse(B3 %in% c(1,2) & packyear1 < 20, "Current, <20 pack-years",
  ifelse(B3 %in% c(1,2) & packyear1 >= 20, "Current, >=20 pack-years", "Unknown"))))))
data1$smoking5[is.na(data1$smoking5)] <- "Unknown"
data1$smoking5 <- factor(data1$smoking5, levels=smoking_levels)

z0 <- function(x) ifelse(is.na(x), 0, x)
data1$occupational_PA_MET_h_day <-
  (8*z0(data1$E2)*data1$E33 + 4*z0(data1$E5)*data1$E66)/60/7
data1$leisure_only_PA_MET_h_day <-
  (8*z0(data1$E11)*data1$E122 + 4*z0(data1$E14)*data1$E155)/60/7
data1$transport_PA_MET_h_day <- 4*z0(data1$E8)*data1$E99/60/7
data1$nonoccupational_PA_MET_h_day <-
  data1$leisure_only_PA_MET_h_day + data1$transport_PA_MET_h_day

dat <- metarowdata_0.2_imputed1
idx <- match(dat$ID, data1$ID)
stopifnot(!anyNA(idx), !anyDuplicated(data1$ID))
update_vars <- c("AGE","A2","BMI","education","smoking5","drink","DMTIME","H1H2",
                 "redmeat","fruit","veg","occupational_PA_MET_h_day",
                 "nonoccupational_PA_MET_h_day")
for (v in update_vars) dat[[v]] <- data1[[v]][idx]
dat$total_PA_MET_h_day <- as.numeric(dat$MET)

metabolite_original <- colnames(scaled_metarowdata)
metabolites <- gsub("\\.", " ", metabolite_original)
mi <- match(metabolite_original, names(dat))
stopifnot(!anyNA(mi))
names(dat)[mi] <- metabolites

dat$Group <- as.integer(as.character(dat$Group))
dat$Match <- factor(dat$Match)
dat$education <- factor(dat$education)
dat$smoking5 <- droplevels(factor(dat$smoking5, levels=smoking_levels))
dat$drink <- factor(dat$drink)
dat$H1H2 <- factor(dat$H1H2)
dat$A2 <- factor(dat$A2)

# Derivation: all intact matched sets, including first-year cases.
derivation_dat <- droplevels(dat)
early_case <- derivation_dat$Group == 1 & as.numeric(derivation_dat$endtime5) <= 1
removed_matches <- unique(as.character(derivation_dat$Match[early_case]))
# CRC analysis: lag-1 matched sample.
lag1_dat <- droplevels(derivation_dat[!(as.character(derivation_dat$Match) %in% removed_matches),])

check_pairs <- function(d, expected_pairs) {
  pt <- table(d$Match, d$Group)
  stopifnot(nrow(pt)==expected_pairs, ncol(pt)==2L, all(pt==1L), nrow(d)==2L*expected_pairs)
}
check_pairs(derivation_dat, 144L)
check_pairs(lag1_dat, 128L)

base_covars <- c("AGE","BMI","education","smoking5","drink","DMTIME","H1H2",
                 "redmeat","fruit","veg")
pa <- c(total="total_PA_MET_h_day",
        occupational="occupational_PA_MET_h_day",
        leisure="nonoccupational_PA_MET_h_day")
pa_labels <- c(total="Total physical activity",
               occupational="Occupational physical activity",
               leisure="Leisure-time physical activity")

required <- c("ID","Group","Match","endtime5","A2",base_covars,unname(pa),metabolites)
stopifnot(!anyNA(derivation_dat[,required]), !anyNA(lag1_dat[,required]))

flow <- data.frame(
  Dataset=c("Signature derivation: full nested case-control sample",
            "Excluded matched sets: case diagnosed within 1 year",
            "CRC association: lag-1 nested case-control sample"),
  Participants=c(nrow(derivation_dat),2*length(removed_matches),nrow(lag1_dat)),
  CRC_cases=c(sum(derivation_dat$Group==1),length(removed_matches),sum(lag1_dat$Group==1)),
  Controls=c(sum(derivation_dat$Group==0),length(removed_matches),sum(lag1_dat$Group==0)),
  Matched_pairs=c(nlevels(derivation_dat$Match),length(removed_matches),nlevels(lag1_dat$Match))
)
write_csv(flow, file.path(output_root,"00_sample_flow_and_roles.csv"))
write_csv(data.frame(Removed_match=removed_matches), file.path(output_root,"01_removed_year1_matches.csv"))

annotation <- metarowdata2
annotation$COMPOUND.Name <- gsub("\\.", " ", annotation$COMPOUND.Name)
annotation <- annotation[!duplicated(annotation$COMPOUND.Name),]

screen_spearman <- function(exposure, out_dir, criterion=c("FDR","nominal_P")) {
  criterion <- match.arg(criterion)
  x <- derivation_dat[[exposure]]
  out <- data.frame(COMPOUND.Name=metabolites, Estimate=NA_real_, P.value=NA_real_)
  for (i in seq_along(metabolites)) {
    ct <- suppressWarnings(cor.test(x, derivation_dat[[metabolites[i]]],
                                    method="spearman", exact=FALSE))
    out$Estimate[i] <- unname(ct$estimate); out$P.value[i] <- ct$p.value
  }
  out$p_fdr <- p.adjust(out$P.value, method="BH")
  ai <- match(out$COMPOUND.Name, annotation$COMPOUND.Name)
  for (v in setdiff(names(annotation),"COMPOUND.Name")) out[[v]] <- annotation[[v]][ai]
  out$SUPER.META.PATHWAY[is.na(out$SUPER.META.PATHWAY)|out$SUPER.META.PATHWAY==""] <-
    "Partially Characterized Molecules"
  selected <- if (criterion == "FDR") out[out$p_fdr < 0.05,] else out[out$P.value < 0.05,]
  write_csv(out, file.path(out_dir,"01_Spearman_all_799_metabolites_full_nested.csv"))
  selected_name <- if (criterion == "FDR") "02_Spearman_FDR_lt_0.05_full_nested.csv" else "02_Spearman_nominal_P_lt_0.05_full_nested.csv"
  write_csv(selected, file.path(out_dir,selected_name))
  list(all=out, selected=selected)
}

plot_screen <- function(s, out_dir, label, seed=20260920L) {
  d <- s$all
  lev <- c("Lipid","Amino Acid","Xenobiotics","Peptide","Nucleotide",
           "Cofactors and Vitamins","Carbohydrate","Energy",
           "Partially Characterized Molecules")
  pal <- c("#A54B37","#FA3100","#FDAE61","#E8D300","#4DAF4A",
           "#AFCFEA","#6A3D9A","#37437A","#999999"); names(pal)<-lev
  set.seed(seed)
  d$key <- runif(nrow(d)); d$Class <- factor(d$SUPER.META.PATHWAY,levels=lev)
  d <- d[order(d$Class,d$key),]; d$Index <- seq_len(nrow(d)); d$logP <- -log10(d$P.value)
  p <- ggplot(d,aes(Index,logP,color=Class))+
    geom_hline(yintercept=-log10(.05),linetype=2,color="#444444")+
    geom_point(size=1.7,alpha=.9)+scale_color_manual(values=pal,drop=FALSE)+
    labs(x="Identified metabolites",y=expression(-log[10](italic(P))),title=label,color="Metabolite class")+
    theme_bw(base_family="Arial",base_size=12)+theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position="right")
  ggsave(file.path(out_dir,"03_Spearman_nominal_P_Manhattan.png"),p,width=11,height=6,dpi=400,bg="white")
  ggsave(file.path(out_dir,"03_Spearman_nominal_P_Manhattan.pdf"),p,width=11,height=6,device=cairo_pdf)
}

score_crc <- function(score, train_mean, train_sd, cuts, label) {
  d <- lag1_dat
  d$signature_score <- score
  d$signature_z <- (score-train_mean)/train_sd
  d$signature_tertile <- cut(score, breaks=c(-Inf,cuts,Inf), labels=c("T1","T2","T3"), include.lowest=TRUE)
  f0 <- fit_clogit(Group~signature_z+strata(Match),d)
  f1 <- fit_clogit(as.formula(paste("Group ~ signature_z +",paste(base_covars,collapse=" + "),"+ strata(Match)")),d)
  t0 <- fit_clogit(Group~signature_tertile+strata(Match),d)
  t1 <- fit_clogit(as.formula(paste("Group ~ signature_tertile +",paste(base_covars,collapse=" + "),"+ strata(Match)")),d)
  make_row <- function(name,fs,ft) {
    e<-extract_or(fs,"signature_z"); e2<-extract_or(ft,"signature_tertileT2"); e3<-extract_or(ft,"signature_tertileT3")
    data.frame(Variable=label,Model=name,T1="1.00 (ref.)",
      T2=fmt_est(e2["OR"],e2["Lower"],e2["Upper"]),
      T3=fmt_est(e3["OR"],e3["Lower"],e3["Upper"]),
      OR=e["OR"],CI_lower=e["Lower"],CI_upper=e["Upper"],OR_95CI=fmt_est(e["OR"],e["Lower"],e["Upper"]),P=e["P"],
      SD_reference="Full nested derivation sample",N=nrow(d),Matched_pairs=nlevels(d$Match),check.names=FALSE)
  }
  list(table=rbind(make_row("Unadjusted",f0,t0),make_row("Full model",f1,t1)),
       models=list(perSD_unadjusted=f0,perSD_full=f1,tertile_unadjusted=t0,tertile_full=t1), data=d)
}

individual_crc <- function(mets) {
  if (!length(mets)) return(data.frame())
  out <- lapply(mets,function(m) {
    d<-lag1_dat; d$met_z<-as.numeric(scale(d[[m]]))
    f<-fit_clogit(as.formula(paste("Group ~ met_z +",paste(base_covars,collapse=" + "),"+ strata(Match)")),d)
    e<-extract_or(f,"met_z")
    data.frame(Metabolite=m,OR=e["OR"],CI_lower=e["Lower"],CI_upper=e["Upper"],OR_95CI=fmt_est(e["OR"],e["Lower"],e["Upper"]),P=e["P"])
  })
  out<-do.call(rbind,out); out$FDR_P<-p.adjust(out$P,"BH"); out
}

run_one <- function(key) {
  exposure <- pa[[key]]; label <- pa_labels[[key]]
  out_dir <- file.path(output_root,paste0(match(key,names(pa)),"_",key,"_PA_signature"))
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  criterion <- if (key %in% c("total","occupational")) "FDR" else "nominal_P"
  s <- screen_spearman(exposure,out_dir,criterion)
  plot_screen(s,out_dir,paste0(label,": Spearman screen in full nested sample"),20260920L+match(key,names(pa)))
  candidates <- s$selected$COMPOUND.Name
  target_candidates <- if (key == "total") {
    c("1-linoleoyl-GPA (18:2)*", "1-linoleoyl-GPA (18:2)", "creatine", "glycine")
  } else {
    c("1-linoleoyl-GPA (18:2)*", "1-linoleoyl-GPA (18:2)")
  }
  forced_metabolites <- intersect(target_candidates, candidates)
  status <- data.frame(Exposure=label,Derivation_N=nrow(derivation_dat),Derivation_pairs=nlevels(derivation_dat$Match),
                       Screening=if (criterion=="FDR") "Spearman BH-FDR P < 0.05" else "Spearman nominal P < 0.05",N_candidates=length(candidates),
                       Forced_metabolites=if(length(forced_metabolites)) paste(forced_metabolites,collapse="; ") else "None (target did not pass screening)",
                       Conditional_forcing_rule=paste0(if (key == "total") "Force 1-linoleoyl-GPA (18:2), creatine, and glycine individually only when each passes the " else "Force 1-linoleoyl-GPA (18:2) only if it passes the ",if (criterion=="FDR") "BH-FDR P < 0.05" else "nominal P < 0.05"," screen"),
                       Unpenalized_covariates="AGE; BMI; A2; smoking5",
                       Validation_N=nrow(lag1_dat),Validation_pairs=nlevels(lag1_dat$Match))
  if (length(candidates)<2L) stop(label,": fewer than two FDR-significant candidates; cannot construct a multimetabolite signature.")
  x_met <- as.matrix(derivation_dat[,candidates,drop=FALSE])
  x_unpen <- model.matrix(~AGE+BMI+A2+smoking5,data=derivation_dat)[,-1,drop=FALSE]
  metabolite_penalty <- ifelse(colnames(x_met) %in% forced_metabolites,0,1)
  x <- cbind(x_met,x_unpen); penalty <- c(metabolite_penalty,rep(0,ncol(x_unpen)))
  set.seed(222)
  cv <- cv.glmnet(x,derivation_dat[[exposure]],alpha=1,family="gaussian",penalty.factor=penalty,nfolds=10)
  fit <- glmnet(x,derivation_dat[[exposure]],alpha=1,family="gaussian",lambda=cv$lambda.min,penalty.factor=penalty)
  cf <- as.matrix(coef(fit)); selected <- intersect(rownames(cf)[cf[,1]!=0],candidates)
  if (length(selected)<2L) stop(label,": LASSO retained fewer than two metabolites.")
  weights <- as.numeric(cf[selected,1]); names(weights)<-selected
  coef_tab <- data.frame(Metabolite=selected,LASSO_coefficient=weights,
                         Penalty_factor=metabolite_penalty[match(selected,colnames(x_met))])
  ai<-match(selected,annotation$COMPOUND.Name)
  coef_tab$HMDB.ID<-annotation$HMDB.ID[ai]
  coef_tab$SUPER.META.PATHWAY<-annotation$SUPER.META.PATHWAY[ai]
  coef_tab$SUB.META.PATHWAY<-annotation$SUB.META.PATHWAY[ai]
  coef_tab<-coef_tab[,c("Metabolite","HMDB.ID","LASSO_coefficient","Penalty_factor",
                        "SUPER.META.PATHWAY","SUB.META.PATHWAY")]
  train_score <- as.numeric(as.matrix(derivation_dat[,selected,drop=FALSE])%*%weights)
  valid_score <- as.numeric(as.matrix(lag1_dat[,selected,drop=FALSE])%*%weights)
  train_mean<-mean(train_score); train_sd<-sd(train_score)
  cuts<-as.numeric(quantile(train_score,c(1/3,2/3),type=2,names=FALSE))
  crc<-score_crc(valid_score,train_mean,train_sd,cuts,paste0(label,"-related metabolic signature"))
  train_min <- min(train_score); train_max <- max(train_score)
  valid_min <- min(valid_score); valid_max <- max(valid_score)
  range_labels <- c(
    sprintf("[%.2f, %.2f]", train_min, cuts[1]),
    sprintf("(%.2f, %.2f]", cuts[1], cuts[2]),
    sprintf("(%.2f, %.2f]", cuts[2], train_max)
  )
  crc$table$T1_score_range <- range_labels[1]
  crc$table$T2_score_range <- range_labels[2]
  crc$table$T3_score_range <- range_labels[3]
  crc$table$Signature_SD <- train_sd
  crc$table$Range_reference <- "Full nested derivation sample"
  ind<-individual_crc(selected)
  # Score performance: Pearson r and R-squared; screening: Spearman correlation.
  rho_train<-cor.test(train_score,derivation_dat[[exposure]],method="pearson")
  rho_valid<-cor.test(valid_score,lag1_dat[[exposure]],method="pearson")
  corr<-data.frame(Sample=c("Full nested derivation","Lag-1 CRC analysis"),N=c(nrow(derivation_dat),nrow(lag1_dat)),
                   Pearson_r=c(unname(rho_train$estimate),unname(rho_valid$estimate)),
                   P=c(rho_train$p.value,rho_valid$p.value),
                   Pearson_R_squared=c(unname(rho_train$estimate)^2,unname(rho_valid$estimate)^2),
                   Variance_explained_percent=100*c(unname(rho_train$estimate)^2,unname(rho_valid$estimate)^2))
  status$N_LASSO_selected<-length(selected)
  status$Target_forced_and_retained<-if(length(forced_metabolites)) all(forced_metabolites %in% selected) else NA
  status$Signature_constructed<-TRUE
  target_screen <- s$all[s$all$COMPOUND.Name %in% unique(target_candidates),,drop=FALSE]
  target_screen$Passed_domain_screen <- if (criterion=="FDR") target_screen$p_fdr < 0.05 else target_screen$P.value < 0.05
  target_screen$Forced_in_LASSO <- target_screen$COMPOUND.Name %in% forced_metabolites
  write_csv(target_screen,file.path(out_dir,"03b_prespecified_forced_metabolites_screening.csv"))
  write_csv(status,file.path(out_dir,"04_signature_status.csv"))
  write_csv(coef_tab,file.path(out_dir,"05_LASSO_selected_metabolites_and_weights.csv"))
  write_csv(data.frame(Feature=colnames(x_unpen),Coefficient=as.numeric(cf[colnames(x_unpen),1]),Penalty_factor=0),file.path(out_dir,"06_LASSO_nonpenalized_covariates.csv"))
  write_csv(data.frame(lambda_min=cv$lambda.min,lambda_1se=cv$lambda.1se),file.path(out_dir,"07_LASSO_lambda.csv"))
  write_csv(data.frame(Sample="Full nested derivation",ID=derivation_dat$ID,Group=derivation_dat$Group,Match=derivation_dat$Match,Signature_score=train_score),file.path(out_dir,"08a_signature_scores_full_nested.csv"))
  write_csv(data.frame(Sample="Lag-1 CRC analysis",ID=lag1_dat$ID,Group=lag1_dat$Group,Match=lag1_dat$Match,Signature_score=valid_score),file.path(out_dir,"08b_signature_scores_lag1.csv"))
  write_csv(data.frame(
    Training_mean=train_mean,Training_SD=train_sd,
    Training_min=train_min,Training_max=train_max,
    Lag1_min=valid_min,Lag1_max=valid_max,
    Tertile_cut1=cuts[1],Tertile_cut2=cuts[2],
    T1_score_range=range_labels[1],T2_score_range=range_labels[2],
    T3_score_range=range_labels[3],
    Range_reference="Full nested derivation sample"
  ),file.path(out_dir,"09_training_scaling_and_tertile_parameters.csv"))
  write_csv(corr,file.path(out_dir,"10_signature_PA_Pearson_r_and_r2.csv"))
  write_csv(crc$table,file.path(out_dir,"11_signature_CRC_lag1_conditional_logistic.csv"))
  write_csv(ind,file.path(out_dir,"12_signature_metabolites_CRC_lag1.csv"))
  pdf(file.path(out_dir,"13_LASSO_cross_validation.pdf"),7,6); plot(cv); dev.off()
  png(file.path(out_dir,"13_LASSO_cross_validation.png"),2100,1800,res=300); plot(cv); dev.off()
  pdf(file.path(out_dir,"14_LASSO_coefficient_paths.pdf"),7,6); plot(cv$glmnet.fit); dev.off()
  png(file.path(out_dir,"14_LASSO_coefficient_paths.png"),2100,1800,res=300); plot(cv$glmnet.fit); dev.off()
  saveRDS(list(screen=s,cvfit=cv,final_fit=fit,weights=weights,train_score=train_score,
               validation_score=valid_score,crc=crc,individual_crc=ind,status=status),
          file.path(out_dir,"15_analysis_objects.rds"))
  list(key=key,label=label,out_dir=out_dir,status=status,weights=weights,crc=crc$table,corr=corr)
}

results <- lapply(names(pa),run_one); names(results)<-names(pa)
summary_status <- do.call(rbind,lapply(results,function(x)x$status))
summary_crc <- do.call(rbind,lapply(results,function(x)x$crc[x$crc$Model=="Full model",]))
write_csv(summary_status,file.path(output_root,"02_three_signature_construction_summary.csv"))
write_csv(summary_crc,file.path(output_root,"03_three_signatures_CRC_fully_adjusted_summary.csv"))

# Forest plot of fully adjusted CRC associations.
pd <- summary_crc; pd$Short <- c("Total PA","Occupational PA","Leisure-time PA")
p <- ggplot(pd,aes(x=OR,y=factor(Short,levels=rev(Short))))+
  geom_vline(xintercept=1,linetype=2,color="#777777")+
  geom_errorbarh(aes(xmin=CI_lower,xmax=CI_upper),height=.12,color="#3B6F8F",linewidth=.8)+
  geom_point(shape=15,size=3,color="#3B6F8F")+
  scale_x_log10()+labs(x="OR per 1-SD increase in metabolic signature",y=NULL,
    title="Metabolic signatures and colorectal cancer risk")+
  theme_bw(base_family="Arial",base_size=12)+theme(panel.grid.minor=element_blank(),panel.grid.major.y=element_blank())
ggsave(file.path(output_root,"04_three_signatures_CRC_forest.png"),p,width=8,height=4.5,dpi=400,bg="white")
ggsave(file.path(output_root,"04_three_signatures_CRC_forest.pdf"),p,width=8,height=4.5,device=cairo_pdf)

script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(length(script_arg)==1L) {
  script_source <- file.path(.integration_code_root, "03_nested_sample_screening.R")
  script_target <- normalizePath(file.path(output_root,"00_full_nested_derivation_lag1_CRC_pipeline.R"),winslash="/",mustWork=FALSE)
  if (!identical(tolower(script_source),tolower(script_target))) file.copy(script_source,script_target,overwrite=TRUE)
}
saveRDS(list(derivation_dat=derivation_dat,lag1_dat=lag1_dat,results=results,flow=flow),file.path(output_root,"05_complete_analysis_objects.rds"))
writeLines(c(
  "Three PA metabolic signatures derived in the complete nested case-control sample and tested for CRC in the lag-1 sample.",
  "Derivation: 288 participants (144 intact matched pairs).",
  "CRC association: 256 participants (128 intact matched pairs) after removal of year-1 cases and matched controls.",
  "Screening: Spearman correlation across 799 metabolites; BH-FDR P < 0.05 for total and occupational PA, and nominal P < 0.05 for leisure-time PA.",
  "LASSO: Gaussian model of PA; AGE, BMI, sex (A2), and smoking5 dummy variables unpenalized.",
  "Conditional forcing: for total PA, 1-linoleoyl-GPA (18:2), creatine, and glycine are each unpenalized only when that metabolite passes BH-FDR P < 0.05. For occupational PA, only 1-linoleoyl-GPA (18:2) is conditionally forced after passing BH-FDR P < 0.05. For leisure-time PA, only 1-linoleoyl-GPA (18:2) is conditionally forced after passing nominal P < 0.05.",
  "Application: metabolite weights, score mean/SD, and tertile cutpoints fixed from the 288-person derivation sample.",
  paste0("CRC full model: conditional logistic regression with strata(Match), adjusted for ",paste(base_covars,collapse=", "),"."),
  "Note: the lag-1 sample is a subset of the derivation sample and therefore is not an independent external validation."
),file.path(output_root,"README.txt"),useBytes=TRUE)

message("Completed: ",output_root)
})

# 3. Manhattan plots with inset tables.
.integration_run_stage("06_rerender_manhattan_with_topright_table.R",
  .integration_nested, {
library(ggplot2)
library(readr)

root <- .integration_analysis
threshold <- 0.05

pathway_order <- c("Lipid","Amino Acid","Xenobiotics","Peptide","Nucleotide",
                   "Cofactors and Vitamins","Carbohydrate","Energy",
                   "Partially Characterized Molecules")
pathway_colors <- c("Lipid"="#A54B37","Amino Acid"="#FA3100","Xenobiotics"="#FDAE61",
                    "Peptide"="#E8D300","Nucleotide"="#4DAF4A",
                    "Cofactors and Vitamins"="#AFCFEA","Carbohydrate"="#6A3D9A",
                    "Energy"="#37437A","Partially Characterized Molecules"="#999999")
pathway_display <- c("Lipid"="Lipid","Amino Acid"="Amino acid","Xenobiotics"="Xenobiotics",
                     "Peptide"="Peptide","Nucleotide"="Nucleotide",
                     "Cofactors and Vitamins"="Cofactors and vitamins",
                     "Carbohydrate"="Carbohydrate","Energy"="Energy",
                     "Partially Characterized Molecules"="Unclassified metabolites")

draw_one <- function(input_csv, output_pdf, output_png, data_csv, summary_csv, seed,
                     criterion=c("nominal_P","FDR")) {
  criterion <- match.arg(criterion)
  d <- read.csv(input_csv, check.names=FALSE, stringsAsFactors=FALSE)
  d$SUPER.META.PATHWAY[is.na(d$SUPER.META.PATHWAY)|d$SUPER.META.PATHWAY==""] <-
    "Partially Characterized Molecules"
  stopifnot(all(d$SUPER.META.PATHWAY %in% pathway_order))
  d$SUPER.META.PATHWAY <- factor(d$SUPER.META.PATHWAY,levels=pathway_order,ordered=TRUE)
  set.seed(seed)
  d$random_key_within_pathway <- runif(nrow(d))
  d <- d[order(d$SUPER.META.PATHWAY,d$random_key_within_pathway),,drop=FALSE]
  d$plot_pos <- seq_len(nrow(d))
  value_col <- if (criterion == "FDR") "p_fdr" else "P.value"
  stopifnot(value_col %in% names(d))
  d$plot_value <- d[[value_col]]
  d$minus_log10_value <- -log10(pmax(d$plot_value,.Machine$double.xmin))

  summary <- do.call(rbind,lapply(pathway_order,function(z) {
    w <- as.character(d$SUPER.META.PATHWAY)==z
    data.frame(SUPER.META.PATHWAY=z,count_all=sum(w),
               count_significant=sum(w & d$plot_value<threshold,na.rm=TRUE))
  }))
  summary <- rbind(summary,data.frame(SUPER.META.PATHWAY="Total",count_all=nrow(d),
                                      count_significant=sum(d$plot_value<threshold,na.rm=TRUE)))
  write_csv(d,data_csv)
  write_csv(summary,summary_csv)

  table_body <- summary
  table_body$display_name <- c(unname(pathway_display[pathway_order]),"Total")
  y_max <- max(4,ceiling(max(d$minus_log10_value,na.rm=TRUE)+0.25))
  table_left <- 602; table_right <- 797
  table_top <- y_max*0.99; table_bottom <- y_max*0.64
  table_row_h <- (table_top-table_bottom)/(nrow(table_body)+1)
  table_col_edges <- c(table_left,700,746,table_right)
  table_hlines <- data.frame(x=table_left,xend=table_right,
    y=seq(table_bottom,table_top,length.out=nrow(table_body)+2),
    yend=seq(table_bottom,table_top,length.out=nrow(table_body)+2))
  table_vlines <- data.frame(x=table_col_edges,xend=table_col_edges,
                             y=table_bottom,yend=table_top)
  table_body$y <- table_top-table_row_h*(seq_len(nrow(table_body))+0.5)

  y_axis_label <- if (criterion == "FDR") expression(-log[10](FDR~italic(P))) else expression(-log[10](italic(P)))
  p <- ggplot(d,aes(plot_pos,minus_log10_value,color=SUPER.META.PATHWAY))+
    geom_hline(yintercept=-log10(threshold),linetype=2,color="black",linewidth=.45)+
    geom_point(size=2.4,alpha=.90)+
    scale_color_manual(values=pathway_colors,drop=FALSE)+
    scale_x_continuous(limits=c(0,800),breaks=NULL,labels=NULL,expand=c(0,0))+
    scale_y_continuous(limits=c(0,y_max),breaks=0:y_max,expand=c(0,0))+
    labs(x="Identified metabolites",y=y_axis_label,color=NULL)+
    annotate("rect",xmin=table_left,xmax=table_right,ymin=table_bottom,ymax=table_top,
             fill="white",color=NA)+
    geom_segment(data=table_hlines,aes(x=x,xend=xend,y=y,yend=yend),inherit.aes=FALSE,
                 color="grey75",linewidth=.35)+
    geom_segment(data=table_vlines,aes(x=x,xend=xend,y=y,yend=yend),inherit.aes=FALSE,
                 color="grey75",linewidth=.35)+
    annotate("text",x=table_left+10,y=table_top-table_row_h/2,label="Super pathway",
             hjust=0,family="Arial",size=3.05,color="black")+
    annotate("text",x=mean(table_col_edges[2:3]),y=table_top-table_row_h/2,
             label="No. of\nmetabolites",lineheight=.86,family="Arial",size=2.70,color="black")+
    annotate("text",x=mean(table_col_edges[3:4]),y=table_top-table_row_h/2,
             label="No. statistically\nsignificant",lineheight=.86,family="Arial",size=2.55,color="black")+
    geom_point(data=table_body[table_body$SUPER.META.PATHWAY!="Total",],
               aes(x=table_left+6,y=y,color=SUPER.META.PATHWAY),inherit.aes=FALSE,
               size=3.2,show.legend=FALSE)+
    geom_text(data=table_body,aes(x=table_left+11,y=y,label=display_name),inherit.aes=FALSE,
              hjust=0,family="Arial",size=2.85,color="black")+
    geom_text(data=table_body,aes(x=mean(table_col_edges[2:3]),y=y,label=count_all),
              inherit.aes=FALSE,family="Arial",size=2.85,color="black")+
    geom_text(data=table_body,aes(x=mean(table_col_edges[3:4]),y=y,label=count_significant),
              inherit.aes=FALSE,family="Arial",size=2.85,color="black")+
    guides(color="none")+
    theme_bw(base_family="Arial",base_size=13)+
    theme(panel.grid=element_blank(),panel.border=element_rect(color="black",fill=NA,linewidth=.65),
          axis.text=element_text(family="Arial",color="black"),
          axis.title=element_text(family="Arial",color="black"),legend.position="none",
          plot.margin=margin(8,14,8,8))
  ggsave(output_pdf,p,width=12,height=8,units="in",device=cairo_pdf)
  ggsave(output_png,p,width=12,height=8,units="in",dpi=600,bg="white")
}

jobs <- data.frame(
  rel=c(
    "13_full_nested_derived_signatures_nominalP0.05_lag1_CRC/1_total_PA_signature",
    "13_full_nested_derived_signatures_nominalP0.05_lag1_CRC/2_occupational_PA_signature",
    "13_full_nested_derived_signatures_nominalP0.05_lag1_CRC/3_leisure_PA_signature",
    "15_full_nested_leisure_signature_Pearson_P0.05_lag1_CRC/3_leisure_PA_signature"),
  method=c("Spearman","Spearman","Spearman","Pearson"),
  seed=c(20260921L,20260922L,20260923L,20260923L),
  criterion=c("FDR","FDR","nominal_P","nominal_P"),stringsAsFactors=FALSE)

# Rebuild the comparison Pearson screen from the current derivation sample.
pearson_dir <- file.path(root, "15_full_nested_leisure_signature_Pearson_P0.05_lag1_CRC",
                         "3_leisure_PA_signature")
dir.create(pearson_dir, recursive = TRUE, showWarnings = FALSE)
pearson_dat <- readRDS(file.path(.integration_nested, "05_complete_analysis_objects.rds"))$derivation_dat
pearson_screen <- read.csv(file.path(.integration_leisure,
  "01_Spearman_all_799_metabolites_full_nested.csv"), check.names = FALSE)
for (i in seq_len(nrow(pearson_screen))) {
  ct <- cor.test(pearson_dat$nonoccupational_PA_MET_h_day,
                 pearson_dat[[pearson_screen$COMPOUND.Name[i]]], method = "pearson")
  pearson_screen$Estimate[i] <- unname(ct$estimate)
  pearson_screen$P.value[i] <- ct$p.value
}
pearson_screen$p_fdr <- p.adjust(pearson_screen$P.value, "BH")
write.csv(pearson_screen, file.path(pearson_dir, "01_Pearson_all_799_metabolites_full_nested.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

for(i in seq_len(nrow(jobs))) {
  d <- file.path(root,jobs$rel[i]); m <- jobs$method[i]
  draw_one(
    file.path(d,paste0("01_",m,"_all_799_metabolites_full_nested.csv")),
    file.path(d,paste0("03_",m,"_nominal_P_Manhattan.pdf")),
    file.path(d,paste0("03_",m,"_nominal_P_Manhattan.png")),
    file.path(d,paste0("03a_",m,"_nominal_P_Manhattan_plot_data.csv")),
    file.path(d,paste0("03b_",m,"_nominal_P_Manhattan_pathway_summary.csv")),
    jobs$seed[i],jobs$criterion[i])

  # Export FDR files for total and occupational PA.
  if (jobs$criterion[i] == "FDR") {
    file.copy(file.path(d,paste0("03_",m,"_nominal_P_Manhattan.pdf")),
              file.path(d,paste0("03_",m,"_FDR_Manhattan.pdf")),overwrite=TRUE)
    file.copy(file.path(d,paste0("03_",m,"_nominal_P_Manhattan.png")),
              file.path(d,paste0("03_",m,"_FDR_Manhattan.png")),overwrite=TRUE)
    file.copy(file.path(d,paste0("03a_",m,"_nominal_P_Manhattan_plot_data.csv")),
              file.path(d,paste0("03a_",m,"_FDR_Manhattan_plot_data.csv")),overwrite=TRUE)
    file.copy(file.path(d,paste0("03b_",m,"_nominal_P_Manhattan_pathway_summary.csv")),
              file.path(d,paste0("03b_",m,"_FDR_Manhattan_pathway_summary.csv")),overwrite=TRUE)
  }
}

})
