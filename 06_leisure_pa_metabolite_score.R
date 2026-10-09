# Leisure-PA 16-metabolite score.
# Select the seed by same-sample adjusted CRC OR; estimates are exploratory.
# Run the selected-seed analysis once at the search endpoint.

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

.integration_require_output("05_leisure_pa_selection.R",
  file.path(.integration_leisure, "forced_18_opposite_direction_56met_LASSO",
            "19g_VIF_18_forced_metabolites_plus_covariates.csv"))

# 1. Search seeds 1-300.
.integration_run_stage("code_seed_search_1_300.R",
  file.path(.integration_leisure, "forced_15_seed1_300_min_fullmodel_OR"), {
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
})

# Use package-local inputs and results.
out_dir <- getwd()
base_dir <- dirname(out_dir)
quadrant_dir <- file.path(base_dir, "four_quadrant_plot_56_metabolites")
previous_vif_dir <- file.path(base_dir, "forced_18_opposite_direction_56met_LASSO")

quadrant <- read.csv(file.path(quadrant_dir, "01_quadrant_plot_data.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
screen <- read.csv(file.path(base_dir, "02_Spearman_nominal_P_lt_0.05_full_nested.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
obj <- readRDS(file.path(.integration_nested, "05_complete_analysis_objects.rds"))
derivation <- obj$derivation_dat
lag1 <- obj$lag1_dat
previous_vif <- read.csv(file.path(previous_vif_dir,
                          "19g_VIF_18_forced_metabolites_plus_covariates.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
previous_vif <- subset(previous_vif,
                       Predictor_set == "18 forced metabolites + covariates" & VIF >= 5)
excluded <- previous_vif$Metabolite
opposite <- (quadrant$Spearman_rho > 0 & quadrant$log_OR_CRC < 0) |
            (quadrant$Spearman_rho < 0 & quadrant$log_OR_CRC > 0)
forced <- quadrant$Metabolite[opposite & !(quadrant$Metabolite %in% excluded)]
candidates <- setdiff(screen$COMPOUND.Name, excluded)
stopifnot(nrow(quadrant) == 56L, nrow(screen) == 56L,
          nrow(derivation) == 288L, nrow(lag1) == 256L,
          length(unique(derivation$Match)) == 144L,
          length(unique(lag1$Match)) == 128L,
          length(excluded) == 3L, length(forced) == 15L,
          length(candidates) == 53L, all(forced %in% candidates),
          all(candidates %in% names(derivation)), all(candidates %in% names(lag1)))

x_met <- as.matrix(derivation[, candidates, drop = FALSE])
x_cov <- model.matrix(~ AGE + BMI + A2 + smoking5, data = derivation)[, -1, drop = FALSE]
pen_met <- ifelse(candidates %in% forced, 0, 1)
x <- cbind(x_met, x_cov)
penalty <- c(pen_met, rep(0, ncol(x_cov)))
stopifnot(!anyNA(x), !anyNA(derivation$nonoccupational_PA_MET_h_day))
covars <- c("AGE", "BMI", "education", "smoking5", "drink", "DMTIME",
            "H1H2", "redmeat", "fruit", "veg")
full_formula <- as.formula(paste("Group ~ signature_z +",
                                 paste(covars, collapse = " + "),
                                 "+ strata(Match)"))

fit_one_seed <- function(seed) {
  tryCatch({
    set.seed(seed)
    cv <- cv.glmnet(x, derivation$nonoccupational_PA_MET_h_day,
                    family = "gaussian", alpha = 1, nfolds = 10,
                    penalty.factor = penalty)
    fit <- glmnet(x, derivation$nonoccupational_PA_MET_h_day,
                  family = "gaussian", alpha = 1, lambda = cv$lambda.min,
                  penalty.factor = penalty)
    cf <- as.matrix(coef(fit))
    selected <- intersect(rownames(cf)[abs(cf[, 1]) > 1e-10], candidates)
    if (!all(forced %in% selected)) {
      return(data.frame(Seed = seed, Status = "forced_metabolite_zero",
                        Selected_count = length(selected),
                        Forced_retained_count = sum(forced %in% selected),
                        Lambda_min = cv$lambda.min, OR = NA_real_,
                        CI_lower = NA_real_, CI_upper = NA_real_, P = NA_real_,
                        Pearson_r = NA_real_, Pearson_r2 = NA_real_,
                        Selected_metabolites = paste(selected, collapse = " | "),
                        stringsAsFactors = FALSE))
    }
    weights <- as.numeric(cf[selected, 1])
    names(weights) <- selected
    train_score <- as.numeric(as.matrix(derivation[, selected, drop = FALSE]) %*% weights)
    test_score <- as.numeric(as.matrix(lag1[, selected, drop = FALSE]) %*% weights)
    train_sd <- sd(train_score)
    if (!is.finite(train_sd) || train_sd <= 0) stop("Invalid derivation score SD")
    d <- lag1
    d$signature_z <- (test_score - mean(train_score)) / train_sd
    model <- suppressWarnings(clogit(full_formula, data = d, method = "exact",
                                     na.action = na.fail))
    z <- summary(model)$coefficients["signature_z", ]
    beta <- unname(z["coef"])
    se <- unname(z["se(coef)"])
    r <- unname(cor(train_score, derivation$nonoccupational_PA_MET_h_day,
                    method = "pearson"))
    data.frame(Seed = seed, Status = "ok",
               Selected_count = length(selected),
               Forced_retained_count = sum(forced %in% selected),
               Lambda_min = cv$lambda.min, OR = exp(beta),
               CI_lower = exp(beta - 1.96 * se),
               CI_upper = exp(beta + 1.96 * se),
               P = unname(z["Pr(>|z|)"]), Pearson_r = r,
               Pearson_r2 = r^2,
               Selected_metabolites = paste(selected, collapse = " | "),
               stringsAsFactors = FALSE)
  }, error = function(e) {
    data.frame(Seed = seed, Status = paste0("error: ", conditionMessage(e)),
               Selected_count = NA_integer_, Forced_retained_count = NA_integer_,
               Lambda_min = NA_real_, OR = NA_real_, CI_lower = NA_real_,
               CI_upper = NA_real_, P = NA_real_, Pearson_r = NA_real_,
               Pearson_r2 = NA_real_, Selected_metabolites = NA_character_,
               stringsAsFactors = FALSE)
  })
}

scan <- do.call(rbind, lapply(1:300, function(seed) {
  result <- fit_one_seed(seed)
  if (seed %% 25L == 0L) cat(sprintf("Completed %d / 300 seeds\n", seed))
  result
}))
write.csv(scan, file.path(out_dir, "00_seed_scan_all_1_to_300.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
eligible <- subset(scan, Status == "ok" & is.finite(OR))
if (nrow(eligible) == 0L) stop("No valid seed produced a full-model estimate.")
eligible <- eligible[order(eligible$OR, eligible$Seed), ]
winner <- eligible[1, , drop = FALSE]
winner$Selection_rule <- "Minimum full-model per-SD OR among seeds 1-300"
winner$Valid_seed_count <- nrow(eligible)
winner$Failed_seed_count <- nrow(scan) - nrow(eligible)
write.csv(winner, file.path(out_dir, "00_selected_seed_minimum_fullmodel_OR.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Selected seed %d: OR %.6f, P %.6g, selected metabolites %d\n",
            winner$Seed, winner$OR, winner$P, winner$Selected_count))
Sys.setenv(LEISURE_PA_LASSO_SEED = as.character(winner$Seed))
# 1.1 Complete analysis using the selected seed.
# Evaluate in the original global environment.
eval(quote({
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(Hmisc)
  library(ggplot2)
})
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))

# Use this run's full nested analysis object.
out_dir <- getwd()
base_dir <- dirname(out_dir)
quadrant_dir <- file.path(base_dir, "four_quadrant_plot_56_metabolites")
previous_vif_dir <- file.path(base_dir, "forced_18_opposite_direction_56met_LASSO")

quadrant <- read.csv(file.path(quadrant_dir, "01_quadrant_plot_data.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
screen <- read.csv(file.path(base_dir, "02_Spearman_nominal_P_lt_0.05_full_nested.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
obj <- readRDS(file.path(.integration_nested, "05_complete_analysis_objects.rds"))
derivation <- obj$derivation_dat
lag1 <- obj$lag1_dat
stopifnot(nrow(quadrant) == 56L, nrow(screen) == 56L,
          nrow(derivation) == 288L, nrow(lag1) == 256L,
          length(unique(derivation$Match)) == 144L,
          length(unique(lag1$Match)) == 128L)

# CRC directions: adjusted lag-1 conditional logistic models; outcome-guided selection.
opposite <- (quadrant$Spearman_rho > 0 & quadrant$log_OR_CRC < 0) |
            (quadrant$Spearman_rho < 0 & quadrant$log_OR_CRC > 0)
previous_vif <- read.csv(file.path(previous_vif_dir,
                          "19g_VIF_18_forced_metabolites_plus_covariates.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
previous_vif <- subset(previous_vif,
                       Predictor_set == "18 forced metabolites + covariates" & VIF >= 5)
excluded <- previous_vif$Metabolite
forced <- quadrant$Metabolite[opposite & !(quadrant$Metabolite %in% excluded)]
candidates <- setdiff(screen$COMPOUND.Name, excluded)
stopifnot(length(excluded) == 3L, length(forced) == 15L,
          !anyDuplicated(forced), !anyDuplicated(candidates),
          setequal(candidates, setdiff(quadrant$Metabolite, excluded)),
          all(forced %in% candidates), all(excluded %in% screen$COMPOUND.Name),
          all(screen$P.value < 0.05),
          all(candidates %in% names(derivation)), all(candidates %in% names(lag1)))
excluded_table <- merge(
  previous_vif[, c("Metabolite", "VIF")],
  quadrant[, c("Metabolite", "HMDB_ID", "CRC_OR_per_1SD", "CRC_P")],
  by = "Metabolite", all.x = TRUE, sort = FALSE
)
excluded_table$Abs_log_OR <- abs(log(excluded_table$CRC_OR_per_1SD))
excluded_table$Exclusion_reason <- "VIF >= 5 in the 18-forced-metabolite model"
write.csv(excluded_table, file.path(out_dir, "00_excluded_VIF_ge5_metabolites.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
forced_table <- quadrant[opposite & !(quadrant$Metabolite %in% excluded),
                         c("Metabolite", "HMDB_ID", "Super_pathway",
                           "Spearman_rho", "Spearman_P",
                           "CRC_OR_per_1SD", "CRC_P", "log_OR_CRC")]
forced_table$Quadrant <- ifelse(forced_table$Spearman_rho > 0,
                                "PA positive, CRC inverse",
                                "PA inverse, CRC positive")
write.csv(forced_table, file.path(out_dir, "01_forced_15_metabolites.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Candidates: 53 metabolites in N=288; force 15 with zero penalty.
# Unpenalized covariates adjust the model and are excluded from the score.
x_met <- as.matrix(derivation[, candidates, drop = FALSE])
x_cov <- model.matrix(~ AGE + BMI + A2 + smoking5, data = derivation)[, -1, drop = FALSE]
pen_met <- ifelse(candidates %in% forced, 0, 1)
x <- cbind(x_met, x_cov)
penalty <- c(pen_met, rep(0, ncol(x_cov)))
stopifnot(!anyNA(x), !anyNA(derivation$nonoccupational_PA_MET_h_day))
analysis_seed <- as.integer(Sys.getenv("LEISURE_PA_LASSO_SEED", "222"))
stopifnot(is.finite(analysis_seed), analysis_seed >= 1L)
set.seed(analysis_seed)
cv <- cv.glmnet(x, derivation$nonoccupational_PA_MET_h_day,
                family = "gaussian", alpha = 1, nfolds = 10,
                penalty.factor = penalty)
fit <- glmnet(x, derivation$nonoccupational_PA_MET_h_day,
              family = "gaussian", alpha = 1, lambda = cv$lambda.min,
              penalty.factor = penalty)
cf <- as.matrix(coef(fit))
selected <- intersect(rownames(cf)[abs(cf[, 1]) > 1e-10], candidates)
weights <- setNames(as.numeric(cf[selected, 1]), selected)
if (!all(forced %in% selected)) {
  stop("At least one of the 15 required metabolites has a zero coefficient; ",
       "the signature would not satisfy the forced-inclusion requirement.")
}
coefs <- data.frame(
  Metabolite = selected,
  HMDB_ID = quadrant$HMDB_ID[match(selected, quadrant$Metabolite)],
  Super_pathway = quadrant$Super_pathway[match(selected, quadrant$Metabolite)],
  Sub_pathway = quadrant$Sub_pathway[match(selected, quadrant$Metabolite)],
  LASSO_coefficient = as.numeric(weights),
  Penalty_factor = pen_met[match(selected, candidates)],
  Forced_opposite_direction = selected %in% forced
)
write.csv(coefs, file.path(out_dir, "02_signature_metabolites_and_coefficients.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
direction_check <- merge(
  coefs,
  quadrant[, c("Metabolite", "Spearman_rho", "Spearman_P")],
  by = "Metabolite", all.x = TRUE, sort = FALSE
)
direction_check$Spearman_direction <- ifelse(direction_check$Spearman_rho > 0,
                                              "positive", "negative")
direction_check$LASSO_direction <- ifelse(direction_check$LASSO_coefficient > 0,
                                           "positive", "negative")
direction_check$Direction_concordant <- sign(direction_check$Spearman_rho) ==
  sign(direction_check$LASSO_coefficient)
direction_check <- direction_check[
  order(!direction_check$Forced_opposite_direction,
        !direction_check$Direction_concordant,
        direction_check$Metabolite), ]
write.csv(direction_check,
          file.path(out_dir, "02b_LASSO_vs_Spearman_direction_check.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
# Forced-metabolite correlations, before and after covariate residualization.
x_forced <- x_met[, forced, drop = FALSE]
raw_cor <- cor(x_forced, use = "pairwise.complete.obs", method = "pearson")
cov_design_for_resid <- cbind(intercept = 1, x_cov)
resid_met_for_cor <- apply(x_forced, 2L, function(v) {
  as.numeric(lm.fit(cov_design_for_resid, v)$residuals)
})
adjusted_cor <- cor(resid_met_for_cor, use = "pairwise.complete.obs",
                    method = "pearson")
write.csv(raw_cor,
          file.path(out_dir, "19a_Pearson_correlation_matrix_retained_15.csv"),
          fileEncoding = "UTF-8")
write.csv(adjusted_cor,
          file.path(out_dir, "19b_covariate_residualized_correlation_matrix_retained_15.csv"),
          fileEncoding = "UTF-8")

pair_index <- t(combn(forced, 2L))
pair_cor <- data.frame(
  Metabolite_1 = pair_index[, 1],
  Metabolite_2 = pair_index[, 2],
  Pearson_r_raw = NA_real_,
  Pearson_r_covariate_residualized = NA_real_,
  N = nrow(x_met),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(pair_index))) {
  a <- pair_index[i, 1]
  b <- pair_index[i, 2]
  pair_cor$Pearson_r_raw[i] <- raw_cor[a, b]
  pair_cor$Pearson_r_covariate_residualized[i] <- adjusted_cor[a, b]
}
pair_cor$Abs_r_raw <- abs(pair_cor$Pearson_r_raw)
pair_cor$Abs_r_covariate_residualized <-
  abs(pair_cor$Pearson_r_covariate_residualized)
pair_cor$Metabolite_1_selected <- pair_cor$Metabolite_1 %in% selected
pair_cor$Metabolite_2_selected <- pair_cor$Metabolite_2 %in% selected
pair_cor$Metabolite_1_forced <- pair_cor$Metabolite_1 %in% forced
pair_cor$Metabolite_2_forced <- pair_cor$Metabolite_2 %in% forced
pair_cor$Contains_direction_reversal <- pair_cor$Metabolite_1 %in%
  direction_check$Metabolite[!direction_check$Direction_concordant] |
  pair_cor$Metabolite_2 %in%
  direction_check$Metabolite[!direction_check$Direction_concordant]
pair_cor <- pair_cor[order(-pair_cor$Abs_r_covariate_residualized), ]
write.csv(pair_cor,
          file.path(out_dir, "19c_pairwise_correlations_retained_15_ranked.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(subset(pair_cor, Abs_r_covariate_residualized >= 0.70),
          file.path(out_dir, "19c_high_correlations_retained_15_abs_r_ge_0.70.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

draw_corr_heatmap <- function(mat, file_stub, subtitle, width = 18, height = 17,
                               axis_size = 6) {
  ord <- hclust(as.dist(1 - abs(mat)), method = "average")$order
  labs <- colnames(mat)[ord]
  d <- expand.grid(Metabolite_x = labs, Metabolite_y = rev(labs),
                   stringsAsFactors = FALSE)
  d$Correlation <- as.vector(mat[rev(labs), labs, drop = FALSE])
  d$Metabolite_x <- factor(d$Metabolite_x, levels = labs)
  d$Metabolite_y <- factor(d$Metabolite_y, levels = rev(labs))
  g <- ggplot(d, aes(x = Metabolite_x, y = Metabolite_y, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    coord_fixed() +
    labs(title = "Pairwise metabolite correlation matrix", subtitle = subtitle,
         x = NULL, y = NULL) +
    theme_minimal(base_family = "Arial", base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = axis_size),
          axis.text.y = element_text(size = axis_size),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9),
          legend.position = "right")
  ggsave(file.path(out_dir, paste0(file_stub, ".png")), g,
         width = width, height = height, dpi = 320, bg = "white")
  ggsave(file.path(out_dir, paste0(file_stub, ".pdf")), g,
         width = width, height = height, device = cairo_pdf, bg = "white")
}
draw_corr_heatmap(raw_cor, "19d_correlation_heatmap_raw_retained_15",
                   "The 15 retained metabolites; full nested derivation (N = 288)",
                   width = 15, height = 15, axis_size = 7)
draw_corr_heatmap(adjusted_cor,
                   "19e_correlation_heatmap_covariate_residualized_retained_15",
                   "The 15 retained metabolites; adjusted for age, BMI, sex, and smoking",
                   width = 15, height = 15, axis_size = 7)

calc_metabolite_vif <- function(target, predictor_set, model_name) {
  other_metabolites <- setdiff(predictor_set, target)
  design <- cbind(intercept = 1,
                  x_met[, other_metabolites, drop = FALSE], x_cov)
  y <- x_met[, target]
  fit <- lm.fit(design, y)
  tss <- sum((y - mean(y))^2)
  r_squared <- 1 - sum(fit$residuals^2) / tss
  vif <- if (r_squared >= 1 - 1e-12) Inf else 1 / (1 - r_squared)
  data.frame(Metabolite = target, Predictor_set = model_name,
             R_squared_from_other_predictors = r_squared, VIF = vif,
             stringsAsFactors = FALSE)
}
vif_15 <- do.call(rbind, lapply(forced, calc_metabolite_vif,
                               predictor_set = forced,
                               model_name = "15 retained forced metabolites + covariates"))
vif_53 <- do.call(rbind, lapply(forced, calc_metabolite_vif,
                               predictor_set = candidates,
                               model_name = "53 remaining candidates + covariates"))
vif_table <- rbind(vif_15, vif_53)
vif_table$VIF_category <- ifelse(vif_table$VIF >= 10, "At least 10",
                                 ifelse(vif_table$VIF >= 5, "5 to <10", "Below 5"))
vif_table <- merge(vif_table,
                   direction_check[, c("Metabolite", "Spearman_rho",
                                        "LASSO_coefficient",
                                        "Direction_concordant")],
                   by = "Metabolite", all.x = TRUE, sort = FALSE)
vif_table <- vif_table[order(vif_table$Predictor_set, -vif_table$VIF), ]
write.csv(vif_table, file.path(out_dir, "19f_VIF_retained_15_metabolites.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(subset(vif_table, Predictor_set == "15 retained forced metabolites + covariates"),
          file.path(out_dir, "19g_VIF_retained_15_plus_covariates.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(data.frame(Feature = colnames(x_cov), Penalty_factor = 0,
                     Coefficient = as.numeric(cf[colnames(x_cov), 1])),
          file.path(out_dir, "03_nonpenalized_covariates.csv"), row.names = FALSE)
write.csv(data.frame(Metabolite = candidates,
                     Spearman_rho = screen$Estimate[match(candidates, screen$COMPOUND.Name)],
                     Spearman_P = screen$P.value[match(candidates, screen$COMPOUND.Name)],
                     Forced_opposite_direction = candidates %in% forced,
                     Penalty_factor = pen_met,
                     Selected_in_signature = candidates %in% selected),
          file.path(out_dir, "04_all_53_remaining_candidate_penalties.csv"), row.names = FALSE)
write.csv(data.frame(Seed = analysis_seed, Derivation_N = nrow(derivation),
                     Derivation_pairs = length(unique(derivation$Match)),
                     CRC_N = nrow(lag1), CRC_pairs = length(unique(lag1$Match)),
                     Initial_candidate_count = nrow(screen),
                     Candidate_count = length(candidates),
                     Excluded_VIF_ge5_count = length(excluded),
                     Forced_count = length(forced), Selected_count = length(selected),
                     Forced_retained_count = sum(forced %in% selected),
                     Lambda_min = cv$lambda.min, Lambda_1se = cv$lambda.1se,
                     Cross_validation_folds = 10L),
          file.path(out_dir, "05_analysis_specification.csv"), row.names = FALSE)

score <- function(d) as.numeric(as.matrix(d[, selected, drop = FALSE]) %*% weights)
train_score <- score(derivation)
valid_score <- score(lag1)
train_mean <- mean(train_score)
train_sd <- sd(train_score)
stopifnot(is.finite(train_sd), train_sd > 0)
write.csv(data.frame(ID = derivation$ID, Group = derivation$Group,
                     Match = derivation$Match, Signature_score = train_score,
                     Leisure_PA_MET_h_day = derivation$nonoccupational_PA_MET_h_day),
          file.path(out_dir, "06_signature_scores_full_nested.csv"), row.names = FALSE)
write.csv(data.frame(ID = lag1$ID, Group = lag1$Group,
                     Match = lag1$Match, Signature_score = valid_score,
                     Leisure_PA_MET_h_day = lag1$nonoccupational_PA_MET_h_day),
          file.path(out_dir, "07_signature_scores_lag1.csv"), row.names = FALSE)
pearson <- function(s, pa, label) {
  test <- cor.test(s, pa, method = "pearson")
  r <- unname(test$estimate)
  data.frame(Sample = label, N = length(s), Pearson_r = r,
             CI_lower = test$conf.int[1], CI_upper = test$conf.int[2],
             P = test$p.value, Pearson_r_squared = r^2,
             Variance_explained_percent = 100 * r^2)
}
cor_table <- rbind(
  pearson(train_score, derivation$nonoccupational_PA_MET_h_day,
          "Full nested derivation"),
  pearson(valid_score, lag1$nonoccupational_PA_MET_h_day,
          "Lag-1 CRC analysis")
)
write.csv(cor_table, file.path(out_dir, "08_signature_leisure_PA_Pearson_r_r2.csv"),
          row.names = FALSE)

fit_clogit <- function(formula, data) suppressWarnings(
  clogit(formula, data = data, method = "exact", na.action = na.fail))
extract_or <- function(fit, term) {
  z <- summary(fit)$coefficients
  beta <- unname(z[term, "coef"])
  se <- unname(z[term, "se(coef)"])
  c(OR = exp(beta), Lower = exp(beta - 1.96 * se),
    Upper = exp(beta + 1.96 * se), P = unname(z[term, "Pr(>|z|)"]))
}
fmt <- function(z) sprintf("%.2f (%.2f-%.2f)", z["OR"], z["Lower"], z["Upper"])
covars <- c("AGE", "BMI", "education", "smoking5", "drink", "DMTIME",
            "H1H2", "redmeat", "fruit", "veg")
cov_text <- paste(covars, collapse = " + ")
d <- lag1
d$signature_score <- valid_score
d$signature_z <- (valid_score - train_mean) / train_sd
cuts <- as.numeric(quantile(train_score, c(1/3, 2/3), type = 2, names = FALSE))
stopifnot(cuts[1] < cuts[2])
d$signature_tertile <- cut(valid_score, breaks = c(-Inf, cuts, Inf),
                           labels = c("T1", "T2", "T3"), include.lowest = TRUE)
ranges <- c(sprintf("[%.2f, %.2f]", min(train_score), cuts[1]),
            sprintf("(%.2f, %.2f]", cuts[1], cuts[2]),
            sprintf("(%.2f, %.2f]", cuts[2], max(train_score)))
write.csv(data.frame(Tertile = c("T1", "T2", "T3"),
                     Range = ranges, Lower = c(min(train_score), cuts[1], cuts[2]),
                     Upper = c(cuts[1], cuts[2], max(train_score)),
                     Signature_SD_full_nested = train_sd,
                     Reference_sample = "Full nested derivation, N=288"),
          file.path(out_dir, "09_signature_tertile_ranges_and_SD.csv"),
          row.names = FALSE)

models <- list(
  Unadjusted_continuous = fit_clogit(Group ~ signature_z + strata(Match), d),
  Full_continuous = fit_clogit(as.formula(paste(
    "Group ~ signature_z +", cov_text, "+ strata(Match)")), d),
  Unadjusted_tertile = fit_clogit(Group ~ signature_tertile + strata(Match), d),
  Full_tertile = fit_clogit(as.formula(paste(
    "Group ~ signature_tertile +", cov_text, "+ strata(Match)")), d)
)
signature_rows <- lapply(c("Unadjusted", "Full"), function(model) {
  cont <- extract_or(models[[paste0(model, "_continuous")]], "signature_z")
  tfit <- models[[paste0(model, "_tertile")]]
  data.frame(Model = model, T1_range = ranges[1], T1 = "1 (ref)",
             T2_range = ranges[2], T2 = fmt(extract_or(tfit, "signature_tertileT2")),
             T3_range = ranges[3], T3 = fmt(extract_or(tfit, "signature_tertileT3")),
             Per_1SD_OR_95CI = fmt(cont), OR = cont["OR"],
             CI_lower = cont["Lower"], CI_upper = cont["Upper"],
             P = cont["P"], Signature_SD = train_sd, N = nrow(d),
             Matched_pairs = length(unique(d$Match)))
})
signature_table <- do.call(rbind, signature_rows)
write.csv(signature_table, file.path(out_dir, "10_signature_CRC_conditional_logistic.csv"),
          row.names = FALSE)

# Fit each metabolite separately per lag-1 sample SD.
single <- do.call(rbind, lapply(selected, function(met) {
  one <- d
  one$met_z <- as.numeric(scale(one[[met]]))
  fit <- fit_clogit(as.formula(paste("Group ~ met_z +", cov_text,
                                       "+ strata(Match)")), one)
  est <- extract_or(fit, "met_z")
  data.frame(Metabolite = met,
             HMDB_ID = quadrant$HMDB_ID[match(met, quadrant$Metabolite)],
             Forced_opposite_direction = met %in% forced,
             Spearman_rho = quadrant$Spearman_rho[match(met, quadrant$Metabolite)],
             CRC_OR_per_1SD = est["OR"], CI_lower = est["Lower"],
             CI_upper = est["Upper"], OR_95CI = fmt(est), P = est["P"])
}))
single$FDR_P_within_selected <- p.adjust(single$P, method = "BH")
single$FDR_P_within_original_56 <- quadrant$CRC_FDR_P_within_56[
  match(single$Metabolite, quadrant$Metabolite)]
write.csv(single, file.path(out_dir, "11_signature_metabolites_CRC_full_model.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Three-knot spline; reference = minimum score.
rd <- d[, c("Group", "Match", "signature_score", covars), drop = FALSE]
basis <- Hmisc::rcspline.eval(rd$signature_score, nk = 3, inclx = TRUE)
knots <- attr(basis, "knots")
colnames(basis) <- paste0("rcs", seq_len(ncol(basis)))
rd <- cbind(rd, as.data.frame(basis))
spline_terms <- colnames(basis)
spline_fit <- fit_clogit(as.formula(paste("Group ~",
  paste(spline_terms, collapse = " + "), "+", cov_text,
  "+ strata(Match)")), rd)
linear_fit <- fit_clogit(as.formula(paste(
  "Group ~ signature_score +", cov_text, "+ strata(Match)")), rd)
null_fit <- fit_clogit(as.formula(paste(
  "Group ~", cov_text, "+ strata(Match)")), rd)
lrt <- function(small, large) {
  statistic <- 2 * (as.numeric(logLik(large)) - as.numeric(logLik(small)))
  df <- attr(logLik(large), "df") - attr(logLik(small), "df")
  pchisq(statistic, df = df, lower.tail = FALSE)
}
p_linear <- signature_table$P[signature_table$Model == "Full"]
p_overall <- lrt(null_fit, spline_fit)
p_nonlinear <- lrt(linear_fit, spline_fit)
reference <- min(rd$signature_score)
xpred <- seq(reference, max(rd$signature_score), length.out = 300)
pb <- Hmisc::rcspline.eval(xpred, knots = knots, inclx = TRUE)
rb <- Hmisc::rcspline.eval(reference, knots = knots, inclx = TRUE)
contrast <- sweep(pb, 2, as.numeric(rb), "-")
beta <- coef(spline_fit)[spline_terms]
variance <- vcov(spline_fit)[spline_terms, spline_terms, drop = FALSE]
lp <- as.numeric(contrast %*% beta)
se <- sqrt(rowSums((contrast %*% variance) * contrast))
pred <- data.frame(Signature_score = xpred, OR = exp(lp),
                   Lower = exp(lp - 1.96 * se), Upper = exp(lp + 1.96 * se))
write.csv(pred, file.path(out_dir, "12_RCS_predictions.csv"), row.names = FALSE)
write.csv(data.frame(Knot = seq_along(knots), Value = as.numeric(knots)),
          file.path(out_dir, "13_RCS_knots.csv"), row.names = FALSE)
rcs_stats <- data.frame(Test = c("Linear per 1-SD", "Overall spline", "Nonlinear"),
                        P = c(p_linear, p_overall, p_nonlinear),
                        Reference_score = reference, Knots = 3L,
                        N = nrow(rd), Matched_pairs = length(unique(rd$Match)),
                        Model = "Full covariates; exact conditional logistic")
write.csv(rcs_stats, file.path(out_dir, "14_RCS_statistics.csv"), row.names = FALSE)

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
p <- ggplot(pred, aes(Signature_score, OR)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#5B9BD5", alpha = 0.25) +
  geom_line(color = "#2C78B8", linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.35) +
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 50)) +
  annotate("text", x = quantile(xpred, 0.25), y = Inf,
           label = paste0("P for linear: ", fmt_p(p_linear)),
           vjust = 1.2, size = 3.4, family = "Arial") +
  annotate("text", x = quantile(xpred, 0.25), y = Inf,
           label = paste0("P for non-linear: ", fmt_p(p_nonlinear)),
           vjust = 3.0, size = 3.4, family = "Arial") +
  theme_classic(base_family = "Arial", base_size = 11) +
  labs(x = "Leisure-time physical activity-related metabolic signature score",
       y = "Odds ratio (log scale)")
ggsave(file.path(out_dir, "15_signature_CRC_RCS_3knots.png"), p,
       width = 6, height = 5.5, dpi = 600, bg = "white")
ggsave(file.path(out_dir, "15_signature_CRC_RCS_3knots.pdf"), p,
       width = 6, height = 5.5, device = cairo_pdf, bg = "white")
saveRDS(list(cv = cv, fit = fit, coefficients = coefs, correlation = cor_table,
             signature_CRC = signature_table, constituents_CRC = single,
             RCS_model = spline_fit, RCS_statistics = rcs_stats),
        file.path(out_dir, "16_analysis_objects.rds"))
cat("Forced:", length(forced), "Selected:", length(selected), "\n")
print(cor_table)
print(signature_table)
print(rcs_stats)
}), envir = .GlobalEnv)

# Record seed selection and exploratory interpretation.
write.csv(data.frame(Seed_search = "1-300",
                     Selection_criterion = "Minimum fully adjusted per-SD OR for signature-CRC association",
                     Forced_metabolites = length(forced),
                     Remaining_candidate_metabolites = length(candidates),
                     Derivation_N = nrow(derivation),
                     CRC_analysis_N = nrow(lag1),
                     Valid_seed_count = nrow(eligible),
                     Selected_seed = winner$Seed,
                     Selected_OR = winner$OR,
                     Selection_bias_note = "Seed selected using CRC outcome in this same lag-1 sample; reported association and P value are exploratory and optimistically biased."),
          file.path(out_dir, "00_seed_search_provenance.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
})

# 2. LASSO cross-validation and coefficient paths.
.integration_run_stage("code_LASSO_regression_metabolite_selection_AB.R",
  file.path(.integration_leisure, "forced_15_seed1_300_min_fullmodel_OR"), {
# Plot the saved LASSO model.
# Panel A: 10-fold CV means +/- 1 SE; mark lambda.min and lambda.1se.
# Panel B: saved coefficient paths versus L1 norm.
# Top-axis counts include covariates; the score includes metabolites only.
# Export 600-dpi PNG, vector PDF, and model audit CSV.

suppressPackageStartupMessages(library(glmnet))
analysis_dir <- file.path(.integration_leisure, "forced_15_seed1_300_min_fullmodel_OR")

obj <- readRDS(file.path(analysis_dir, '16_analysis_objects.rds'))
spec <- read.csv(file.path(analysis_dir, '05_analysis_specification.csv'), check.names = FALSE)
coefficients <- read.csv(file.path(analysis_dir, '02_signature_metabolites_and_coefficients.csv'),
                         check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = 'UTF-8')
penalties <- read.csv(file.path(analysis_dir, '04_all_53_remaining_candidate_penalties.csv'),
                     check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = 'UTF-8')
stopifnot(inherits(obj$cv, 'cv.glmnet'), inherits(obj$fit, 'glmnet'), nrow(spec) == 1L)
cv <- obj$cv
path_fit <- cv$glmnet.fit
stopifnot(inherits(path_fit, 'glmnet'), length(path_fit$lambda) > 1L,
          isTRUE(all.equal(cv$lambda.min, spec$Lambda_min[1], tolerance = 1e-10)))

# Check selected coefficients against the saved coefficient table.
final_cf <- as.matrix(coef(obj$fit))[, 1]
stopifnot(all(coefficients$Metabolite %in% names(final_cf)))
coefficient_difference <- max(abs(final_cf[coefficients$Metabolite] - coefficients$LASSO_coefficient))
stopifnot(coefficient_difference < 1e-9, nrow(coefficients) == spec$Selected_count[1])

draw_AB <- function() {
  par(mfrow = c(1, 2), family = 'Arial', mar = c(4.6, 4.7, 4.5, 1.0),
      mgp = c(2.8, 0.8, 0), cex = 1, cex.axis = 0.95, cex.lab = 1,
      las = 0, tcl = -0.3, bty = 'l')
  plot(cv, xlab = expression(log(lambda)), ylab = 'Mean-squared error')
  mtext('(A)', side = 3, line = 3.1, adj = 0, cex = 1.2, font = 1)
  # Use the CV object for the full path; obj$fit contains only lambda.min.
  plot(path_fit, xvar = 'norm', label = FALSE, xlab = 'L1 norm', ylab = 'Coefficients')
  mtext('(B)', side = 3, line = 3.1, adj = 0, cex = 1.2, font = 1)
}

stem <- file.path(analysis_dir, '17_LASSO_regression_metabolite_selection_AB')
png(paste0(stem, '.png'), width = 12, height = 5.5, units = 'in', res = 600,
    pointsize = 11, type = 'cairo', family = 'Arial', bg = 'white')
draw_AB()
invisible(dev.off())
cairo_pdf(paste0(stem, '.pdf'), width = 12, height = 5.5, pointsize = 11,
          family = 'Arial', bg = 'white')
draw_AB()
invisible(dev.off())

audit <- data.frame(
  Seed = spec$Seed[1], Derivation_N = spec$Derivation_N[1],
  Initial_candidates = spec$Initial_candidate_count[1],
  Candidate_count = spec$Candidate_count[1], Forced_count = spec$Forced_count[1],
  Selected_metabolites = nrow(coefficients), CV_folds = spec$Cross_validation_folds[1],
  Lambda_min = cv$lambda.min, Lambda_1se = cv$lambda.1se,
  Path_lambda_count = length(path_fit$lambda),
  Saved_fit_lambda_count = length(obj$fit$lambda),
  Max_coefficient_difference_from_existing_table = coefficient_difference,
  Model_refitted = FALSE, Font = 'Arial', Width_inches = 12, Height_inches = 5.5,
  PNG_dpi = 600, PDF_is_vector = TRUE,
  Plot_source = '16_analysis_objects.rds: cv and cv$glmnet.fit')
write.csv(audit, file.path(analysis_dir, '18_LASSO_AB_plot_model_audit.csv'), row.names = FALSE)
cat('Existing saved model verified: seed=', spec$Seed[1],
    '; selected metabolites=', nrow(coefficients), '; path lambdas=', length(path_fit$lambda), '\n', sep = '')
cat(paste0(stem, c('.png', '.pdf')), sep = '\n')
})
