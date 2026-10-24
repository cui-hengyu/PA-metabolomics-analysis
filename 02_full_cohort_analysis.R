# Full-cohort analysis.

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
.integration_input_root <- dirname(.integration_code_root)
.integration_data <- file.path(.integration_input_root, "data")
.integration_results <- file.path(.integration_code_root, "results")
.integration_analysis <- file.path(.integration_results, "analysis")
.integration_nested <- file.path(.integration_analysis, "nested")
.integration_total <- file.path(.integration_nested, "total_PA")
.integration_leisure <- file.path(.integration_nested, "nonoccupational_PA")
dir.create(.integration_results, recursive = TRUE, showWarnings = FALSE)
if (.Platform$OS.type == "windows" && getRversion() >= "4.3.0" &&
    getRversion() < "4.4.0" && dir.exists(file.path(.integration_input_root, "_library"))) {
  .libPaths(c(file.path(.integration_input_root, "_library"), .libPaths()))
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

# 1. Baseline characteristics.
.integration_run_stage("14_baseline_table_19519.R",
  file.path(.integration_analysis, "13_full_cohort_PA_RCS_min_reference"), {
# Rebuild from raw Excel; exclude CRC diagnosed within one year.
# Retain unclassified smoking responses as Unknown.

Sys.setlocale("LC_CTYPE", "Chinese_China.utf8")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(openxlsx))

root <- .integration_analysis
out_dir <- file.path(root, "13_full_cohort_PA_RCS_min_reference")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_out <- function(x, filename) {
  write.csv(x, file.path(out_dir, filename), row.names = FALSE,
            fileEncoding = "UTF-8", na = "")
}

# Cohort preparation.
source(file.path(.integration_code_root, "01_data_preprocessing.R"),
       local = TRUE, encoding = "UTF-8")
stopifnot(exists("data1"))

num0 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
}

data1$total_PA_MET_h_day <- suppressWarnings(as.numeric(data1$MET小时))
data1$occupational_PA_MET_h_day <-
  (8 * num0(data1$E2) * num0(data1$E33) +
     4 * num0(data1$E5) * num0(data1$E66)) / 60 / 7
data1$transport_PA_MET_h_day <-
  4 * num0(data1$E8) * num0(data1$E99) / 60 / 7
data1$leisure_PA_MET_h_day <-
  (8 * num0(data1$E11) * num0(data1$E122) +
     4 * num0(data1$E14) * num0(data1$E155)) / 60 / 7
data1$nonoccupational_PA_MET_h_day <-
  data1$leisure_PA_MET_h_day + data1$transport_PA_MET_h_day

smoking_levels <- c(
  "Never smoke", "Former smoke, <20 pack-years",
  "Former smoke, >=20 pack-years", "Current smoke, <20 pack-years",
  "Current smoke, >=20 pack-years", "Unknown"
)
data1$smoking5 <- with(data1, ifelse(
  B1 == 2, "Never smoke",
  ifelse(B3 == 3 & packyear1 < 20, "Former smoke, <20 pack-years",
  ifelse(B3 == 3 & packyear1 >= 20, "Former smoke, >=20 pack-years",
  ifelse(B3 %in% c(1, 2) & packyear1 < 20, "Current smoke, <20 pack-years",
  ifelse(B3 %in% c(1, 2) & packyear1 >= 20,
         "Current smoke, >=20 pack-years", "Unknown"))))
))
data1$smoking5[is.na(data1$smoking5)] <- "Unknown"
data1$smoking5 <- factor(data1$smoking5, levels = smoking_levels)

analysis_vars <- c(
  "ID", "colo_cancerfb", "endtime5", "A2", "AGE", "BMI", "education",
  "smoking5", "drink", "DMTIME", "H1H2", "redmeat", "fruit", "veg",
  "total_PA_MET_h_day", "occupational_PA_MET_h_day",
  "nonoccupational_PA_MET_h_day"
)
dat <- droplevels(data1[complete.cases(data1[, analysis_vars, drop = FALSE]), ,
                        drop = FALSE])
lag1_dat <- droplevels(dat[!(dat$colo_cancerfb == 1L & dat$endtime5 <= 1), ,
                           drop = FALSE])
stopifnot(nrow(lag1_dat) == 19519L)

lag1_dat$sex <- factor(
  ifelse(lag1_dat$A2 == 1, "Male",
         ifelse(lag1_dat$A2 == 2, "Female", NA_character_)),
  levels = c("Male", "Female")
)
lag1_dat$education_label <- factor(
  ifelse(lag1_dat$education == 1, "Illiterate",
  ifelse(lag1_dat$education == 2, "Primary school",
  ifelse(lag1_dat$education == 3, "Middle school",
  ifelse(lag1_dat$education == 4, "High school and above", NA_character_)))),
  levels = c("Illiterate", "Primary school", "Middle school",
             "High school and above")
)
lag1_dat$drink_label <- factor(
  ifelse(lag1_dat$drink == 1, "Ever drink",
         ifelse(lag1_dat$drink == 2, "Never drink", NA_character_)),
  levels = c("Ever drink", "Never drink")
)
lag1_dat$medication_label <- factor(
  ifelse(lag1_dat$H1H2 == 1, "Yes",
         ifelse(lag1_dat$H1H2 == 2, "No", NA_character_)),
  levels = c("Yes", "No")
)
lag1_dat$crc_label <- factor(
  ifelse(lag1_dat$colo_cancerfb == 1, "Yes", "No"),
  levels = c("Yes", "No")
)

male_n <- sum(lag1_dat$sex == "Male", na.rm = TRUE)
female_n <- sum(lag1_dat$sex == "Female", na.rm = TRUE)
stopifnot(male_n + female_n == nrow(lag1_dat))

# Summaries and tests.
fmt_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}
fmt_mean_sd <- function(x) sprintf("%.2f (%.2f)", mean(x), sd(x))
fmt_cat <- function(v, g, level, sex_value) {
  denom <- sum(g == sex_value, na.rm = TRUE)
  n <- sum(v == level & g == sex_value, na.rm = TRUE)
  sprintf("%d (%.1f)", n, 100 * n / denom)
}
wilcox_p <- function(x, g) {
  suppressWarnings(wilcox.test(x ~ g, exact = FALSE)$p.value)
}
chisq_p <- function(v, g) {
  tab <- table(v, g, useNA = "ifany")
  if (nrow(tab) < 2L || ncol(tab) < 2L) return(NA_real_)
  suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
}

rows <- list()
add_row <- function(label, male = "", female = "", p = NA_real_) {
  rows[[length(rows) + 1L]] <<- data.frame(
    Characteristics = label,
    Male = male,
    Female = female,
    `P value` = fmt_p(p),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
add_section <- function(label, p) add_row(label, p = p)

add_row("Participants, n", as.character(male_n), as.character(female_n))
add_row("Age, years", fmt_mean_sd(lag1_dat$AGE[lag1_dat$sex == "Male"]),
        fmt_mean_sd(lag1_dat$AGE[lag1_dat$sex == "Female"]),
        wilcox_p(lag1_dat$AGE, lag1_dat$sex))
add_row("BMI, kg/m2", fmt_mean_sd(lag1_dat$BMI[lag1_dat$sex == "Male"]),
        fmt_mean_sd(lag1_dat$BMI[lag1_dat$sex == "Female"]),
        wilcox_p(lag1_dat$BMI, lag1_dat$sex))

edu_levels <- levels(lag1_dat$education_label)
add_section("Education, n (%)", chisq_p(lag1_dat$education_label, lag1_dat$sex))
for (z in edu_levels) {
  add_row(paste0("   ", z),
          fmt_cat(lag1_dat$education_label, lag1_dat$sex, z, "Male"),
          fmt_cat(lag1_dat$education_label, lag1_dat$sex, z, "Female"))
}

add_section("Cigarette smoking status, n (%)",
            chisq_p(lag1_dat$smoking5, lag1_dat$sex))
for (z in smoking_levels) {
  display_z <- sub("^Former smoke", "Former smoker", z)
  display_z <- sub("^Current smoke", "Current smoker", display_z)
  display_z <- sub("^Never smoke$", "Never smoker", display_z)
  add_row(paste0("   ", display_z),
          fmt_cat(lag1_dat$smoking5, lag1_dat$sex, z, "Male"),
          fmt_cat(lag1_dat$smoking5, lag1_dat$sex, z, "Female"))
}

add_section("Alcohol drinking status, n (%)",
            chisq_p(lag1_dat$drink_label, lag1_dat$sex))
for (z in levels(lag1_dat$drink_label)) {
  add_row(paste0("   ", z),
          fmt_cat(lag1_dat$drink_label, lag1_dat$sex, z, "Male"),
          fmt_cat(lag1_dat$drink_label, lag1_dat$sex, z, "Female"))
}

add_row("Diabetes duration, years",
        fmt_mean_sd(lag1_dat$DMTIME[lag1_dat$sex == "Male"]),
        fmt_mean_sd(lag1_dat$DMTIME[lag1_dat$sex == "Female"]),
        wilcox_p(lag1_dat$DMTIME, lag1_dat$sex))

add_section("Oral antidiabetic medication, n (%)",
            chisq_p(lag1_dat$medication_label, lag1_dat$sex))
for (z in levels(lag1_dat$medication_label)) {
  add_row(paste0("   ", z),
          fmt_cat(lag1_dat$medication_label, lag1_dat$sex, z, "Male"),
          fmt_cat(lag1_dat$medication_label, lag1_dat$sex, z, "Female"))
}

continuous <- list(
  list(label = "Red meat consumption, g/day", x = lag1_dat$redmeat),
  list(label = "Fruit consumption, g/day", x = lag1_dat$fruit),
  list(label = "Vegetable consumption, g/day", x = lag1_dat$veg),
  list(label = "Total physical activity, MET-h/d", x = lag1_dat$total_PA_MET_h_day),
  list(label = "Occupational physical activity, MET-h/d", x = lag1_dat$occupational_PA_MET_h_day),
  list(label = "Non-occupational physical activity, MET-h/d", x = lag1_dat$nonoccupational_PA_MET_h_day)
)
for (item in continuous) {
  add_row(item$label,
          fmt_mean_sd(item$x[lag1_dat$sex == "Male"]),
          fmt_mean_sd(item$x[lag1_dat$sex == "Female"]),
          wilcox_p(item$x, lag1_dat$sex))
}

add_section("Incident CRC, n (%)", chisq_p(lag1_dat$crc_label, lag1_dat$sex))
for (z in levels(lag1_dat$crc_label)) {
  add_row(paste0("   ", z),
          fmt_cat(lag1_dat$crc_label, lag1_dat$sex, z, "Male"),
          fmt_cat(lag1_dat$crc_label, lag1_dat$sex, z, "Female"))
}

table_out <- do.call(rbind, rows)
write_out(table_out, "JBPCD队列基线表_剔除1年_19519.csv")

# Missingness and metadata.
missing_n <- c(
  Education = sum(lag1_dat$education == 999, na.rm = TRUE),
  `Alcohol drinking status` = sum(lag1_dat$drink == 999, na.rm = TRUE),
  `Oral antidiabetic medication use` = sum(lag1_dat$H1H2 == 999, na.rm = TRUE)
)
missing_out <- data.frame(
  Variable = names(missing_n), Missing_n = as.integer(missing_n),
  Definition = c("education == 999", "drink == 999", "H1H2 == 999"),
  check.names = FALSE
)
write_out(missing_out, "JBPCD队列基线表_剔除1年_19519_缺失情况.csv")

metadata <- data.frame(
  Item = c("Table title", "Participants", "Male", "Female", "CRC cases",
           "Excluded CRC cases diagnosed within first year", "Continuous-variable test",
           "Categorical-variable test", "Smoking missing/invalid handling", "Footnote"),
  Value = c(
    "Baseline characteristics of patients with T2D in the JBPCD cohort after excluding CRC cases diagnosed within the first year",
    nrow(lag1_dat), male_n, female_n, sum(lag1_dat$colo_cancerfb == 1L),
    sum(dat$colo_cancerfb == 1L & dat$endtime5 <= 1),
    "Wilcoxon rank-sum test", "Chi-square test",
    "Retained as Unknown in the five-category smoking variable",
    paste0("Continuous variables were described as means (standard deviations), and categorical variables were described as percentages. P values were calculated using Wilcoxon rank-sum tests for continuous variables and chi-square tests for categorical variables. Percentages may not sum to 100% due to missing data (education, n = ", missing_n[["Education"]], "; alcohol drinking status, n = ", missing_n[["Alcohol drinking status"]], "; oral antidiabetic medication use, n = ", missing_n[["Oral antidiabetic medication use"]], ").")
  ),
  stringsAsFactors = FALSE
)
write_out(metadata, "JBPCD队列基线表_剔除1年_19519_元数据.csv")

message("Baseline table generated for N = ", nrow(lag1_dat),
        " (male = ", male_n, ", female = ", female_n,
        "; CRC cases = ", sum(lag1_dat$colo_cancerfb == 1L), ").")
})

# 2. PA associations: Cox models and restricted cubic splines.
.integration_run_stage("code_full_cohort_PA_RCS_min_reference.R",
  file.path(.integration_analysis, "13_full_cohort_PA_RCS_min_reference"), {
# Use the minimum observed PA as the spline reference.

Sys.setlocale("LC_CTYPE", "Chinese_China.utf8")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(rms)
  library(gridExtra)
})

root <- .integration_analysis
out_dir <- file.path(root, "13_full_cohort_PA_RCS_min_reference")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_out <- function(x, filename, where = out_dir) {
  write.csv(x, file.path(where, filename), row.names = FALSE, fileEncoding = "UTF-8")
}
fmt_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}
extract_cox <- function(fit, term) {
  tab <- summary(fit)$coefficients
  if (!term %in% rownames(tab)) return(c(HR = NA, Lower = NA, Upper = NA, P = NA))
  b <- tab[term, "coef"]
  se <- tab[term, "se(coef)"]
  c(HR = exp(b), Lower = exp(b - 1.96 * se),
    Upper = exp(b + 1.96 * se), P = tab[term, "Pr(>|z|)"])
}

# Cohort preparation.
source(file.path(.integration_code_root, "01_data_preprocessing.R"),
       local = TRUE, encoding = "UTF-8")
stopifnot(exists("data1"))

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
data1$smoking5[is.na(data1$smoking5)] <- "Unknown"
data1$smoking5 <- factor(data1$smoking5, levels = smoking_levels)

num0 <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  x
}
data1$total_PA_MET_h_day <- as.numeric(data1$MET小时)
data1$occupational_PA_MET_h_day <-
  (8 * num0(data1$E2) * num0(data1$E33) +
     4 * num0(data1$E5) * num0(data1$E66)) / 60 / 7
data1$transport_PA_MET_h_day <- 4 * num0(data1$E8) * num0(data1$E99) / 60 / 7
data1$leisure_PA_MET_h_day <-
  (8 * num0(data1$E11) * num0(data1$E122) +
     4 * num0(data1$E14) * num0(data1$E155)) / 60 / 7
data1$nonoccupational_PA_MET_h_day <-
  data1$leisure_PA_MET_h_day + data1$transport_PA_MET_h_day

data1$colo_cancerfb <- as.integer(as.character(data1$colo_cancerfb))
data1$endtime5 <- as.numeric(data1$endtime5)
data1$A2 <- factor(data1$A2)
data1$education <- factor(data1$education)
data1$drink <- factor(data1$drink)
data1$H1H2 <- factor(data1$H1H2)

exposures <- c(
  total = "total_PA_MET_h_day",
  occupational = "occupational_PA_MET_h_day",
  nonoccupational = "nonoccupational_PA_MET_h_day"
)
labels <- c(
  total = "Total physical activity",
  occupational = "Occupational physical activity",
  nonoccupational = "Non-occupational physical activity (leisure-time + transport)"
)
plot_labels <- c(
  total = "Total physical activity",
  occupational = "Occupational physical activity",
  nonoccupational = "Non-occupational physical activity"
)
covariates <- c(
  "AGE", "A2", "BMI", "education", "smoking5", "drink", "DMTIME",
  "H1H2", "redmeat", "fruit", "veg"
)
analysis_vars <- c("ID", "colo_cancerfb", "endtime5", covariates,
                   unname(exposures))
dat <- droplevels(data1[complete.cases(data1[, analysis_vars, drop = FALSE]), , drop = FALSE])
lag1_dat <- droplevels(dat[!(dat$colo_cancerfb == 1L & dat$endtime5 <= 1), , drop = FALSE])
if (any(lag1_dat$endtime5 <= 0) || anyNA(lag1_dat[, c("colo_cancerfb", "endtime5", covariates, unname(exposures))])) {
  stop("Missing or non-positive values remain in lag-1 full-cohort analysis variables.")
}

write_out(data.frame(
  Original_cohort_N = nrow(data1), Baseline_analysis_N = nrow(dat),
  Baseline_CRC_cases = sum(dat$colo_cancerfb),
  Lag1_analysis_N = nrow(lag1_dat), Lag1_CRC_cases = sum(lag1_dat$colo_cancerfb),
  Excluded_CRC_cases_within_1_year = sum(dat$colo_cancerfb == 1L & dat$endtime5 <= 1)
), "00_sample_counts.csv")

# Fully adjusted Cox spline models.
adjustment <- paste(covariates, collapse = " + ")
rcs_stats <- list()
plots <- list()
tertile_guides <- list()
rcs_datadist_vars <- unique(c("endtime5", "colo_cancerfb", covariates, unname(exposures)))
dd_lag1 <- rms::datadist(lag1_dat[, rcs_datadist_vars, drop = FALSE])
assign("dd_lag1", dd_lag1, envir = .GlobalEnv)
old_datadist <- options(datadist = "dd_lag1")
on.exit(options(old_datadist), add = TRUE)

for (nm in names(exposures)) {
  x <- exposures[[nm]]
  reference_value <- min(lag1_dat[[x]], na.rm = TRUE)
  # Use exactly the same type-7 cutpoints and right-closed groups as Table 2.
  # Non-occupational PA has cut1 = 0: T1 consists of zero values only.
  tertile_cuts <- unname(quantile(lag1_dat[[x]], c(1 / 3, 2 / 3),
                                  type = 7, na.rm = TRUE))
  stopifnot(length(tertile_cuts) == 2L, all(is.finite(tertile_cuts)),
            tertile_cuts[1] < tertile_cuts[2])
  tertile_guides[[nm]] <- data.frame(
    Variable = labels[[nm]], Unit = "MET-h/day",
    T1_T2_boundary = tertile_cuts[1], T2_T3_boundary = tertile_cuts[2],
    T1_N = sum(lag1_dat[[x]] <= tertile_cuts[1]),
    T2_N = sum(lag1_dat[[x]] > tertile_cuts[1] &
                 lag1_dat[[x]] <= tertile_cuts[2]),
    T3_N = sum(lag1_dat[[x]] > tertile_cuts[2]),
    Boundary_rule = "T1: x <= cut1; T2: cut1 < x <= cut2; T3: x > cut2"
  )
  dd_lag1$limits["Adjust to", x] <- reference_value
  assign("dd_lag1", dd_lag1, envir = .GlobalEnv)

  z <- paste0(nm, "_z")
  lag1_dat[[z]] <- as.numeric(scale(lag1_dat[[x]]))
  fit_linear <- coxph(as.formula(paste0("Surv(endtime5, colo_cancerfb) ~ ", z,
                                       " + ", adjustment)),
                      data = lag1_dat, ties = "efron")
  p_linear <- extract_cox(fit_linear, z)["P"]

  rcs_term <- paste0("rcs(", x, ", 3)")
  fit_rcs <- rms::cph(as.formula(paste0(
    "Surv(endtime5, colo_cancerfb) ~ ", rcs_term, " + ", adjustment)),
    data = lag1_dat, x = TRUE, y = TRUE, surv = TRUE)
  knots <- as.numeric(fit_rcs$Design$parms[[x]])
  an <- anova(fit_rcs)
  p_col <- grep("^P$", colnames(an), value = TRUE)[1]
  rn <- trimws(rownames(an))
  overall_idx <- which(rn == x)[1]
  nonlinear_idx <- which(rn == "Nonlinear")[1]
  p_overall <- if (is.finite(overall_idx)) as.numeric(an[overall_idx, p_col]) else NA_real_
  p_nonlinear <- if (is.finite(nonlinear_idx)) as.numeric(an[nonlinear_idx, p_col]) else NA_real_

  # Include the observed minimum in the prediction grid.
  q99 <- as.numeric(quantile(lag1_dat[[x]], 0.99, na.rm = TRUE, type = 7))
  grid <- sort(unique(c(reference_value, seq(reference_value, q99, length.out = 200))))
  pred <- as.data.frame(do.call(rms::Predict,
                                c(list(fit_rcs), setNames(list(grid), x),
                                  list(ref.zero = TRUE, fun = exp))))
  pred_out <- data.frame(Exposure_MET_h_day = pred[[x]],
                         HR = pred$yhat, Lower = pred$lower, Upper = pred$upper)
  stub <- c(total = "总体PA", occupational = "职业PA", nonoccupational = "非职业PA")[nm]
  write_out(pred_out, paste0("RCS_剔除1年_", stub, "_预测值.csv"))

  iref <- which.min(abs(pred_out$Exposure_MET_h_day - reference_value))
  if (abs(pred_out$HR[iref] - 1) > 1e-8) {
    stop("Minimum-reference check failed for ", labels[[nm]],
         ": HR at reference = ", pred_out$HR[iref])
  }
  rcs_stats[[nm]] <- data.frame(
    Variable = labels[[nm]], Knots = paste(sprintf("%.4f", knots), collapse = ", "),
    P_linear = unname(p_linear), P_overall = p_overall,
    P_nonlinear = p_nonlinear, Reference = reference_value,
    N = nrow(lag1_dat), CRC_cases = sum(lag1_dat$colo_cancerfb), Status = "Fitted"
  )

  plots[[nm]] <- ggplot(pred_out, aes(Exposure_MET_h_day, HR)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#5B9BD5", alpha = 0.25) +
    geom_vline(xintercept = tertile_cuts, linetype = "dotted",
               colour = "grey45", linewidth = 0.45) +
    # Display rounded PA values; keep the lines at the exact cutpoints.
    annotate("text", x = tertile_cuts, y = Inf,
             label = sprintf("%.2f", tertile_cuts),
             hjust = -0.12, vjust = 1.3, size = 3.5,
             colour = "grey25", family = "Arial") +
    geom_line(colour = "#2C6E9B", linewidth = 1) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
    annotate("text", x = Inf, y = Inf,
             label = paste0("P for linear: ", fmt_p(p_linear)),
             hjust = 1.05, vjust = 1.6, size = 3.5, family = "Arial") +
    annotate("text", x = Inf, y = Inf,
             label = paste0("P for non-linear: ", fmt_p(p_nonlinear)),
             hjust = 1.05, vjust = 3.2, size = 3.5, family = "Arial") +
    labs(title = plot_labels[[nm]],
         subtitle = "CRC risk after excluding cases diagnosed within 1 year",
         x = paste0(plot_labels[[nm]], " (MET-h/day)"), y = "Hazard ratio (95% CI)") +
    theme_classic(base_family = "Arial", base_size = 11) +
    theme(plot.title = element_text(size = 14, family = "Arial"),
          plot.subtitle = element_text(size = 10, family = "Arial"),
          axis.title = element_text(size = 11, family = "Arial"),
          axis.text = element_text(family = "Arial"),
          plot.margin = margin(5.5, 8, 5.5, 5.5))
  ggsave(file.path(out_dir, paste0("RCS_剔除1年_", stub, ".pdf")),
         plots[[nm]], width = 7, height = 5.5, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0("RCS_剔除1年_", stub, ".png")),
         plots[[nm]], width = 7, height = 5.5, dpi = 600)
}

rcs_statistics <- do.call(rbind, rcs_stats)
rownames(rcs_statistics) <- NULL
write_out(rcs_statistics, "12_剔除1年内发病_PA_RCS统计量.csv")
write_out(data.frame(
  Variable = rcs_statistics$Variable,
  Reference = rcs_statistics$Reference,
  Reference_definition = "Minimum observed exposure in the lag-1 full cohort"
), "RCS_reference_values.csv")

combined <- arrangeGrob(plots$total + labs(tag = "(A)"),
                        plots$occupational + labs(tag = "(B)"),
                        plots$nonoccupational + labs(tag = "(C)"), ncol = 3)
ggsave(file.path(root, "RCS_剔除1年_PA及构成与CRC.pdf"), combined,
       width = 12.7, height = 3.33, device = cairo_pdf)
ggsave(file.path(root, "RCS_剔除1年_PA及构成与CRC.png"), combined,
       width = 12.7, height = 3.33, dpi = 600)
ggsave(file.path(out_dir, "RCS_剔除1年_PA及构成与CRC_min_reference.pdf"), combined,
       width = 12.7, height = 3.33, device = cairo_pdf)
ggsave(file.path(out_dir, "RCS_剔除1年_PA及构成与CRC_min_reference.png"), combined,
       width = 12.7, height = 3.33, dpi = 600)

# Manuscript Figure S2: preserve the existing three-panel order, axis wording,
# curve, confidence band, reference and P values; add only the tertile guides.
# The paper's "leisure-time" panel uses non-occupational PA (leisure + transport),
# consistently with the current manuscript and Table 2 definition.
paper_axis_labels <- c(total = "Total physical activity",
                       occupational = "Occupational physical activity",
                       nonoccupational = "Leisure-time physical activity")
paper_panels <- lapply(seq_along(exposures), function(i) {
  nm <- names(exposures)[i]
  plots[[nm]] +
    labs(title = NULL, subtitle = NULL, tag = NULL,
         x = paste0(paper_axis_labels[[nm]], " (MET-h/day)")) +
    theme(plot.margin = margin(8, 8, 8, 8))
})
tagged_panels <- lapply(seq_along(paper_panels), function(i) {
  arrangeGrob(paper_panels[[i]],
    top = grid::textGrob(paste0("(", LETTERS[i], ")"),
      x = grid::unit(0.02, "npc"), hjust = 0,
      gp = grid::gpar(fontfamily = "Arial", fontsize = 14)),
    padding = grid::unit(8, "pt"))
})
paper_s2 <- arrangeGrob(grobs = tagged_panels, ncol = 3)
ggsave(file.path(.integration_code_root, "Fig S2_tertile_lines_600dpi.png"),
       paper_s2, width = 10000 / 600, height = 2615 / 600, units = "in",
       dpi = 600, bg = "white")
ggsave(file.path(.integration_code_root, "Fig S2_tertile_lines.pdf"),
       paper_s2, width = 10000 / 600, height = 2615 / 600, units = "in",
       device = cairo_pdf, bg = "white")
write_out(do.call(rbind, tertile_guides), "Fig S2_tertile_cutpoints.csv",
          where = .integration_code_root)

message("Full-cohort PA RCS completed. All references are the minimum observed exposure.")

# Cox associations by tertile and per SD.
# Use the spline covariates for adjustment; retain Unknown smoking.
make_tertiles <- function(x) {
  cuts <- unname(quantile(x, c(1 / 3, 2 / 3), type = 7, na.rm = TRUE))
  if (length(cuts) != 2L || !all(is.finite(cuts)) || cuts[1] >= cuts[2]) {
    stop("Two distinct tertile cutpoints could not be defined.")
  }
  group <- cut(x, breaks = c(-Inf, cuts, Inf), labels = c("T1", "T2", "T3"),
               include.lowest = TRUE, right = TRUE)
  ranges <- c(
    sprintf("[%.2f, %.2f]", min(x, na.rm = TRUE), cuts[1]),
    sprintf("(%.2f, %.2f]", cuts[1], cuts[2]),
    sprintf("(%.2f, %.2f]", cuts[2], max(x, na.rm = TRUE))
  )
  list(group = group, cuts = cuts, ranges = ranges)
}

cox_exposure_groups <- lapply(exposures, function(v) make_tertiles(lag1_dat[[v]]))
for (nm in names(exposures)) {
  lag1_dat[[paste0(nm, "_tertile")]] <- cox_exposure_groups[[nm]]$group
  lag1_dat[[paste0(nm, "_z")]] <- as.numeric(scale(lag1_dat[[exposures[[nm]]]]))
}

fit_cox_row <- function(nm, adjusted) {
  exposure <- exposures[[nm]]
  zvar <- paste0(nm, "_z")
  group_var <- paste0(nm, "_tertile")
  extra <- if (adjusted) adjustment else character()
  f_sd <- as.formula(paste("Surv(endtime5, colo_cancerfb) ~",
                           paste(c(zvar, extra), collapse = " + ")))
  f_t <- as.formula(paste("Surv(endtime5, colo_cancerfb) ~",
                          paste(c(group_var, extra), collapse = " + ")))
  m_sd <- survival::coxph(f_sd, data = lag1_dat, ties = "efron")
  m_t <- survival::coxph(f_t, data = lag1_dat, ties = "efron")
  e_sd <- extract_cox(m_sd, zvar)
  e_t2 <- extract_cox(m_t, paste0(group_var, "T2"))
  e_t3 <- extract_cox(m_t, paste0(group_var, "T3"))
  data.frame(
    Variable = labels[[nm]],
    Model = if (adjusted) "Fully adjusted" else "Unadjusted",
    T1 = "1.00 (ref.)",
    T2 = sprintf("%.2f (%.2f-%.2f)", e_t2["HR"], e_t2["Lower"], e_t2["Upper"]),
    T3 = sprintf("%.2f (%.2f-%.2f)", e_t3["HR"], e_t3["Lower"], e_t3["Upper"]),
    `Per 1-SD` = sprintf("%.2f (%.2f-%.2f)", e_sd["HR"], e_sd["Lower"], e_sd["Upper"]),
    `P for per 1-SD` = unname(e_sd["P"]),
    SD_MET_h_day = sd(lag1_dat[[exposure]]),
    N = nrow(lag1_dat), CRC_cases = sum(lag1_dat$colo_cancerfb == 1L),
    check.names = FALSE
  )
}

cox_results <- do.call(rbind, lapply(names(exposures), function(nm) {
  rbind(fit_cox_row(nm, FALSE), fit_cox_row(nm, TRUE))
}))
rownames(cox_results) <- NULL
write_out(cox_results, "14_PA_CRC_Cox_tertiles_and_perSD.csv")

cox_cutpoints <- do.call(rbind, lapply(names(exposures), function(nm) {
  z <- cox_exposure_groups[[nm]]
  data.frame(
    Variable = labels[[nm]], Unit = "MET-h/day",
    T1 = z$ranges[1], T2 = z$ranges[2], T3 = z$ranges[3],
    Cut1 = z$cuts[1], Cut2 = z$cuts[2], SD = sd(lag1_dat[[exposures[[nm]]]]),
    Minimum = min(lag1_dat[[exposures[[nm]]]], na.rm = TRUE),
    Maximum = max(lag1_dat[[exposures[[nm]]]], na.rm = TRUE),
    N = nrow(lag1_dat), CRC_cases = sum(lag1_dat$colo_cancerfb == 1L)
  )
}))
rownames(cox_cutpoints) <- NULL
write_out(cox_cutpoints, "15_PA_Cox_tertile_cutpoints_and_SD.csv")

write_out(data.frame(
  Item = c("Analysis population", "Participants", "CRC cases",
           "Outcome", "Unadjusted model", "Fully adjusted model",
           "Covariates"),
  Value = c(
    "Full JBPCD cohort after excluding CRC cases diagnosed within the first year",
    nrow(lag1_dat), sum(lag1_dat$colo_cancerfb == 1L),
    "Surv(endtime5, colo_cancerfb)", "PA exposure only",
    "PA exposure + AGE + A2 + BMI + education + smoking5 + drink + DMTIME + H1H2 + redmeat + fruit + veg",
    paste(covariates, collapse = ", ")
  ), stringsAsFactors = FALSE
), "16_PA_CRC_Cox_model_metadata.csv")

message("Full-cohort PA RCS and Cox analyses completed for N = ", nrow(lag1_dat),
        " and ", sum(lag1_dat$colo_cancerfb == 1L), " CRC cases.")
})

# 3. Observed follow-up: median, Q1, and Q3 in the final analytic cohort.
.integration_run_stage("calc_median_followup.R",
  file.path(.integration_results, "follow_up"), {
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(openxlsx))

source(file.path(.integration_code_root, "01_data_preprocessing.R"),
       local = TRUE, encoding = "UTF-8")

num0 <- function(x) { x <- suppressWarnings(as.numeric(x)); x[is.na(x)] <- 0; x }
data1$total_PA_MET_h_day <- suppressWarnings(as.numeric(data1$MET小时))
data1$occupational_PA_MET_h_day <-
  (8 * num0(data1$E2) * num0(data1$E33) +
     4 * num0(data1$E5) * num0(data1$E66)) / 60 / 7
data1$transport_PA_MET_h_day <- 4 * num0(data1$E8) * num0(data1$E99) / 60 / 7
data1$leisure_PA_MET_h_day <-
  (8 * num0(data1$E11) * num0(data1$E122) +
     4 * num0(data1$E14) * num0(data1$E155)) / 60 / 7
data1$nonoccupational_PA_MET_h_day <-
  data1$leisure_PA_MET_h_day + data1$transport_PA_MET_h_day

smoking_levels <- c(
  "Never smoke", "Former smoke, <20 pack-years",
  "Former smoke, >=20 pack-years", "Current smoke, <20 pack-years",
  "Current smoke, >=20 pack-years", "Unknown")
data1$smoking5 <- with(data1, ifelse(
  B1 == 2, "Never smoke",
  ifelse(B3 == 3 & packyear1 < 20, "Former smoke, <20 pack-years",
  ifelse(B3 == 3 & packyear1 >= 20, "Former smoke, >=20 pack-years",
  ifelse(B3 %in% c(1, 2) & packyear1 < 20, "Current smoke, <20 pack-years",
  ifelse(B3 %in% c(1, 2) & packyear1 >= 20,
         "Current smoke, >=20 pack-years", "Unknown"))))))
data1$smoking5[is.na(data1$smoking5)] <- "Unknown"
data1$smoking5 <- factor(data1$smoking5, levels = smoking_levels)

analysis_vars <- c(
  "ID", "colo_cancerfb", "endtime5", "A2", "AGE", "BMI", "education",
  "smoking5", "drink", "DMTIME", "H1H2", "redmeat", "fruit", "veg",
  "total_PA_MET_h_day", "occupational_PA_MET_h_day",
  "nonoccupational_PA_MET_h_day")
dat <- droplevels(data1[complete.cases(data1[, analysis_vars, drop = FALSE]), , drop = FALSE])
lag1_dat <- droplevels(dat[!(dat$colo_cancerfb == 1L & dat$endtime5 <= 1), , drop = FALSE])

cat("Complete baseline cohort:", nrow(dat), "\n")
cat("Final lag-1 cohort:", nrow(lag1_dat), "\n")
cat("CRC cases in final cohort:", sum(lag1_dat$colo_cancerfb == 1L), "\n")
cat("Excluded CRC cases diagnosed within 1 year:",
    sum(dat$colo_cancerfb == 1L & dat$endtime5 <= 1), "\n")
cat("Follow-up duration unit: years (endtime5)\n")
stopifnot(nrow(lag1_dat) == 19519L,
          sum(lag1_dat$colo_cancerfb == 1L) == 149L,
          !anyDuplicated(lag1_dat$ID),
          all(is.finite(lag1_dat$endtime5)),
          all(lag1_dat$endtime5 >= 0))

# Preserve the original definition: observed duration from baseline to CRC
# diagnosis or the non-CRC endpoint in 01_data_preprocessing.R. This is not a
# reverse Kaplan-Meier estimate. Excluding early cases does not reset time zero.
followup_years <- lag1_dat$endtime5
followup_quantiles <- quantile(followup_years, probs = c(0.25, 0.50, 0.75),
                              na.rm = TRUE, type = 7, names = FALSE)
followup_summary <- data.frame(
  analysis_population = "Full JBPCD cohort; CRC diagnoses within first year excluded",
  participants = nrow(lag1_dat),
  crc_cases = sum(lag1_dat$colo_cancerfb == 1L),
  complete_cohort_before_lag_exclusion = nrow(dat),
  crc_cases_excluded_within_first_year =
    sum(dat$colo_cancerfb == 1L & dat$endtime5 <= 1),
  q1_years = followup_quantiles[1L],
  median_years = followup_quantiles[2L],
  q3_years = followup_quantiles[3L],
  iqr_width_years = followup_quantiles[3L] - followup_quantiles[1L],
  mean_years = mean(followup_years),
  minimum_years = min(followup_years),
  maximum_years = max(followup_years),
  total_person_years = sum(followup_years),
  method = "Observed endtime5; quantile type 7; baseline time origin retained",
  stringsAsFactors = FALSE
)
followup_sentence <- sprintf(
  paste0("During a median follow-up of %.2f years (IQR, %.2f-%.2f years), ",
         "%d incident CRC cases were identified."),
  followup_summary$median_years, followup_summary$q1_years,
  followup_summary$q3_years, followup_summary$crc_cases
)
write.csv(followup_summary,
          file.path(.integration_results, "follow_up", "followup_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
writeLines(c(
  "Population: final full-cohort analytic sample, not the nested case-control sample.",
  "Method: observed follow-up duration (endtime5); quantile type 7.",
  "Time origin: baseline; early CRC cases excluded without subtracting one year.",
  "Q1 and Q3 are the 25th and 75th percentiles. IQR width is Q3 minus Q1.",
  capture.output(print(followup_summary, digits = 12, row.names = FALSE)),
  "", "Response sentence (rounded to two decimal places):", followup_sentence
), file.path(.integration_results, "follow_up", "followup_summary.txt"),
useBytes = TRUE)
cat("Q1 follow-up:", format(followup_quantiles[1L], digits = 12), "years\n")
cat("Median follow-up:", format(followup_quantiles[2L], digits = 12), "years\n")
cat("Q3 follow-up:", format(followup_quantiles[3L], digits = 12), "years\n")
cat("IQR width:", format(followup_summary$iqr_width_years, digits = 12), "years\n")
cat("Mean follow-up:", mean(followup_years), "years\n")
cat(followup_sentence, "\n")
print(summary(followup_years))
})
