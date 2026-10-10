# Leisure-PA metabolic signature.

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

# 1. LASSO and signature analysis.
.integration_run_stage("code_seed257_complete_analysis.R",
  file.path(.integration_leisure, "forced_15_seed257_complete_analysis"), {
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(Hmisc)
  library(ggplot2)
})
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))

# Prepare model data.
out_dir <- getwd()
source(file.path(.integration_code_root, "01_data_preprocessing.R"),
       local = TRUE, encoding = "UTF-8")
smoking_levels <- c("Never", "Former, <20 pack-years", "Former, >=20 pack-years",
                    "Current, <20 pack-years", "Current, >=20 pack-years", "Unknown")
data1$smoking5 <- with(data1, ifelse(
  B1 == 2, "Never",
  ifelse(B3 == 3 & packyear1 < 20, "Former, <20 pack-years",
  ifelse(B3 == 3 & packyear1 >= 20, "Former, >=20 pack-years",
  ifelse(B3 %in% c(1, 2) & packyear1 < 20, "Current, <20 pack-years",
  ifelse(B3 %in% c(1, 2) & packyear1 >= 20, "Current, >=20 pack-years", "Unknown"))))))
data1$smoking5[is.na(data1$smoking5)] <- "Unknown"
data1$smoking5 <- factor(data1$smoking5, levels = smoking_levels)
z0 <- function(x) ifelse(is.na(x), 0, x)
data1$leisure_only_PA_MET_h_day <-
  (8 * z0(data1$E11) * data1$E122 + 4 * z0(data1$E14) * data1$E155) / 60 / 7
data1$transport_PA_MET_h_day <- 4 * z0(data1$E8) * data1$E99 / 60 / 7
data1$nonoccupational_PA_MET_h_day <-
  data1$leisure_only_PA_MET_h_day + data1$transport_PA_MET_h_day
dat <- metarowdata_0.2_imputed1
idx <- match(dat$ID, data1$ID)
stopifnot(!anyNA(idx), !anyDuplicated(data1$ID))
update_vars <- c("AGE", "A2", "BMI", "education", "smoking5", "drink", "DMTIME",
                 "H1H2", "redmeat", "fruit", "veg", "nonoccupational_PA_MET_h_day")
for (v in update_vars) dat[[v]] <- data1[[v]][idx]
metabolite_original <- colnames(scaled_metarowdata)
mi <- match(metabolite_original, names(dat))
stopifnot(!anyNA(mi))
names(dat)[mi] <- gsub("\\.", " ", metabolite_original)
dat$Group <- as.integer(as.character(dat$Group))
dat$Match <- factor(dat$Match)
dat$education <- factor(dat$education)
dat$smoking5 <- droplevels(factor(dat$smoking5, levels = smoking_levels))
dat$drink <- factor(dat$drink)
dat$H1H2 <- factor(dat$H1H2)
dat$A2 <- factor(dat$A2)
derivation <- droplevels(dat)
early_case <- derivation$Group == 1 & as.numeric(derivation$endtime5) <= 1
removed_matches <- unique(as.character(derivation$Match[early_case]))
lag1 <- droplevels(derivation[!(as.character(derivation$Match) %in% removed_matches), ])
check_pairs <- function(d, expected_pairs) {
  pt <- table(d$Match, d$Group)
  stopifnot(nrow(pt) == expected_pairs, ncol(pt) == 2L,
            all(pt == 1L), nrow(d) == 2L * expected_pairs)
}
check_pairs(derivation, 144L)
check_pairs(lag1, 128L)
annotation <- data.frame(
  Metabolite = gsub("\\.", " ", metarowdata2$COMPOUND.Name),
  HMDB_ID = metarowdata2$HMDB.ID,
  Super_pathway = metarowdata2$SUPER.META.PATHWAY,
  Sub_pathway = metarowdata2$SUB.META.PATHWAY,
  stringsAsFactors = FALSE
)
annotation$HMDB_ID[is.na(annotation$HMDB_ID)] <- ""

# Model inputs.
candidates <- c(
  "4-hydroxychlorothalonil",
  "2-hydroxy-3-methylvalerate",
  "phenylpyruvate",
  "acisoga",
  "C-glycosyltryptophan",
  "kynurenate",
  "dopamine 3-O-sulfate",
  "gamma-CEHC glucuronide*",
  "tetrahydrocortisol glucuronide",
  "tetrahydrocortisone glucuronide (5)",
  "docosadioate (C22-DC)",
  "eicosanedioate (C20-DC)",
  "hexadecenedioate (C16:1-DC)*",
  "octadecenedioate (C18:1-DC)*",
  "3-hydroxyoleate*",
  "1-linolenoylglycerol (18:3)",
  "linolenate [alpha or gamma; (18:3n3 or 6)]",
  "3beta,7alpha-dihydroxy-5-cholestenoate",
  "5,6-dihydrouridine",
  "sedoheptulose",
  "2R,3R-dihydroxybutyrate",
  "orotidine",
  "hydroxyasparagine**",
  "gamma-carboxyglutamate",
  "N-methylpipecolate",
  "5-hydroxylysine",
  "hydroxy-N6,N6,N6-trimethyllysine*",
  "N-acetylputrescine",
  "tryptophan betaine",
  "argininate*",
  "homocitrulline",
  "N,N,N-trimethyl-alanylproline betaine (TMAP)",
  "2-methylmalonylcarnitine (C4-DC)",
  "N2,N2-dimethylguanosine",
  "isoleucylhydroxyproline*",
  "gamma-glutamyl-alpha-lysine",
  "gamma-glutamyl-epsilon-lysine",
  "N,N-dimethyl-pro-pro",
  "linolenoylcarnitine (C18:3)*",
  "octadecenedioylcarnitine (C18:1-DC)*",
  "1-linolenoyl-GPC (18:3)*",
  "1-linoleoyl-GPE (18:2)*",
  "1-oleoyl-GPE (18:1)",
  "1-oleoyl-2-linoleoyl-GPC (18:1/18:2)*",
  "1-palmitoleoyl-2-linolenoyl-GPC (16:1/18:3)*",
  "1-palmitoyl-2-oleoyl-GPC (16:0/18:1)",
  "1-stearoyl-2-oleoyl-GPC (18:0/18:1)",
  "1-oleoyl-2-arachidonoyl-GPE (18:1/20:4)*",
  "1-oleoyl-2-linoleoyl-GPE (18:1/18:2)*",
  "sphingomyelin (d18:1/21:0, d17:1/22:0, d16:1/23:0)*",
  "sphingomyelin (d18:1/22:1, d18:2/22:0, d16:1/24:1)*",
  "sphingomyelin (d18:1/25:0, d19:0/24:1, d20:1/23:0, d19:1/24:0)*",
  "sphingomyelin (d18:2/24:2)*"
)
pen_met <- c(1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1)
unpenalized <- candidates[pen_met == 0]
stopifnot(length(candidates) == length(pen_met), !anyDuplicated(candidates),
          all(candidates %in% names(derivation)), all(candidates %in% names(lag1)),
          all(candidates %in% annotation$Metabolite))
write.csv(annotation[match(unpenalized, annotation$Metabolite), ],
          file.path(out_dir, "01_unpenalized_metabolites.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Unpenalized covariates adjust the model and are excluded from the score.
x_met <- as.matrix(derivation[, candidates, drop = FALSE])
x_cov <- model.matrix(~ AGE + BMI + A2 + smoking5, data = derivation)[, -1, drop = FALSE]
x <- cbind(x_met, x_cov)
penalty <- c(pen_met, rep(0, ncol(x_cov)))
stopifnot(!anyNA(x), !anyNA(derivation$nonoccupational_PA_MET_h_day))
analysis_seed <- 257L
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
stopifnot(all(unpenalized %in% selected))
coefs <- data.frame(
  Metabolite = selected,
  HMDB_ID = annotation$HMDB_ID[match(selected, annotation$Metabolite)],
  Super_pathway = annotation$Super_pathway[match(selected, annotation$Metabolite)],
  Sub_pathway = annotation$Sub_pathway[match(selected, annotation$Metabolite)],
  LASSO_coefficient = as.numeric(weights),
  Penalty_factor = pen_met[match(selected, candidates)],
  Unpenalized = selected %in% unpenalized
)
write.csv(coefs, file.path(out_dir, "02_signature_metabolites_and_coefficients.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(data.frame(Feature = colnames(x_cov), Penalty_factor = 0,
                     Coefficient = as.numeric(cf[colnames(x_cov), 1])),
          file.path(out_dir, "03_nonpenalized_covariates.csv"), row.names = FALSE)
write.csv(data.frame(Metabolite = candidates, Penalty_factor = pen_met,
                     Selected_in_signature = candidates %in% selected),
          file.path(out_dir, "04_candidate_penalties.csv"), row.names = FALSE)
write.csv(data.frame(Seed = analysis_seed, Derivation_N = nrow(derivation),
                     Derivation_pairs = length(unique(derivation$Match)),
                     CRC_N = nrow(lag1), CRC_pairs = length(unique(lag1$Match)),
                     Candidate_count = length(candidates),
                     Unpenalized_count = length(unpenalized), Selected_count = length(selected),
                     Unpenalized_retained_count = sum(unpenalized %in% selected),
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
             signature_CRC = signature_table,
             RCS_model = spline_fit, RCS_statistics = rcs_stats),
        file.path(out_dir, "16_analysis_objects.rds"))
cat("Unpenalized:", length(unpenalized), "Selected:", length(selected), "\n")
print(cor_table)
print(signature_table)
print(rcs_stats)


})

# 2. LASSO cross-validation and coefficient paths.
.integration_run_stage("code_LASSO_regression_metabolite_selection_AB.R",
  file.path(.integration_leisure, "forced_15_seed257_complete_analysis"), {
# Plot cross-validation and coefficient paths.

suppressPackageStartupMessages(library(glmnet))
analysis_dir <- file.path(.integration_leisure, "forced_15_seed257_complete_analysis")

obj <- readRDS(file.path(analysis_dir, '16_analysis_objects.rds'))
spec <- read.csv(file.path(analysis_dir, '05_analysis_specification.csv'), check.names = FALSE)
coefficients <- read.csv(file.path(analysis_dir, '02_signature_metabolites_and_coefficients.csv'),
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
  Candidate_count = spec$Candidate_count[1], Unpenalized_count = spec$Unpenalized_count[1],
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
