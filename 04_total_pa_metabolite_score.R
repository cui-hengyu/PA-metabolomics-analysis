# Total-PA 18-metabolite score.

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

.integration_require_output("03_nested_sample_screening.R",
  file.path(.integration_nested, "05_complete_analysis_objects.rds"))

# 1. CRC ring plot for 181 screened metabolites.
.integration_run_stage("code_ring_plot_181_screened_metabolites_CRC.R",
  file.path(.integration_total, "ring_plot_181_screened_metabolites"), {
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(gtable)
  library(grid)
})
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))

write_csv <- function(x, f) write.csv(x, f, row.names = FALSE, fileEncoding = "UTF-8")

# Screening: N=288, Spearman BH-FDR <0.05; CRC: lag-1 N=256 (128 pairs).
obj <- readRDS(file.path(.integration_nested, "05_complete_analysis_objects.rds"))
d <- obj$lag1_dat
screen <- read.csv(
  file.path(.integration_nested, "1_total_PA_signature", "02_Spearman_FDR_lt_0.05_full_nested.csv"),
  check.names = FALSE, fileEncoding = "UTF-8"
)
stopifnot(nrow(screen) == 181L, nrow(d) == 256L,
          length(unique(d$Match)) == 128L,
          all(screen$COMPOUND.Name %in% names(d)))

crc_covars <- c("AGE", "BMI", "education", "smoking5", "drink", "DMTIME",
                "H1H2", "redmeat", "fruit", "veg")
crc_results <- do.call(rbind, lapply(seq_len(nrow(screen)), function(i) {
  m <- screen$COMPOUND.Name[i]
  f <- as.formula(paste0("Group ~ scale(`", m, "`) + ",
                         paste(crc_covars, collapse = " + "),
                         " + strata(Match)"))
  fit <- suppressWarnings(clogit(f, data = d, method = "exact", na.action = na.fail))
  z <- summary(fit)$coefficients[1, ]
  data.frame(
    Metabolite = m,
    HMDB_ID = screen$HMDB.ID[i],
    SUPER.META.PATHWAY = screen$SUPER.META.PATHWAY[i],
    SUB.META.PATHWAY = screen$SUB.META.PATHWAY[i],
    OR = exp(z["coef"]),
    CI_lower = exp(z["coef"] - 1.96 * z["se(coef)"]),
    CI_upper = exp(z["coef"] + 1.96 * z["se(coef)"]),
    P = z["Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}))
crc_results$FDR_P <- p.adjust(crc_results$P, method = "BH")
crc_results$SUPER.META.PATHWAY[
  is.na(crc_results$SUPER.META.PATHWAY) | crc_results$SUPER.META.PATHWAY == ""
] <- "Partially Characterized Molecules"
write_csv(crc_results, "01_181_screened_metabolites_CRC_full_model.csv")
write_csv(crc_results[crc_results$P < 0.05, , drop = FALSE],
          "02_181_screened_metabolites_CRC_nominal_P_lt_0.05.csv")

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
association_levels <- c("Positive association (OR > 1)",
                        "Negative association (OR < 1)",
                        "Non-significant (P >= 0.05)")
association_colors <- c(
  "Positive association (OR > 1)" = "#D75B4E",
  "Negative association (OR < 1)" = "#3F7FB3",
  "Non-significant (P >= 0.05)" = "#C8CDD2"
)

ring <- crc_results
ring$Class <- factor(ring$SUPER.META.PATHWAY, levels = class_order)
if (anyNA(ring$Class)) stop("Unrecognized metabolite class in ring-plot data.")
ring <- ring[order(ring$Class, ring$P, ring$OR), , drop = FALSE]
ring$class_number <- as.integer(ring$Class)
gap_slots <- 2L
ring$index <- seq_len(nrow(ring)) + (ring$class_number - 1L) * gap_slots
total_slots <- max(ring$index) + gap_slots
ring$logP <- -log10(pmax(ring$P, .Machine$double.xmin))
ring$Association <- ifelse(
  ring$P < .05 & ring$OR > 1, "Positive association (OR > 1)",
  ifelse(ring$P < .05 & ring$OR < 1, "Negative association (OR < 1)",
         "Non-significant (P >= 0.05)")
)
ring$Association <- factor(ring$Association, levels = association_levels)

class_ring_inner <- 2.46; class_ring_outer <- 2.72; bar_base <- 2.82
ring$bar_top <- bar_base + ring$logP
p05_height <- bar_base - log10(.05)
sector_split <- split(ring$index, ring$Class, drop = TRUE)
sectors <- data.frame(
  Class = factor(names(sector_split), levels = class_order),
  xmin = vapply(sector_split, min, numeric(1)) - .44,
  xmax = vapply(sector_split, max, numeric(1)) + .44,
  mid = vapply(sector_split, function(x) mean(range(x)), numeric(1)),
  row.names = NULL
)

sig <- ring[ring$P < .05, , drop = FALSE]
if (nrow(sig)) {
  sig <- sig[order(sig$index), , drop = FALSE]
  sig$label_cluster <- cumsum(c(Inf, diff(sig$index)) > 3)
  sig$label_radius <- sig$bar_top + .16
  for (cluster_id in unique(sig$label_cluster)) {
    ii <- which(sig$label_cluster == cluster_id)
    if (length(ii) >= 2L) {
      oi <- ii[order(sig$label_radius[ii], sig$index[ii])]
      spacing <- if (length(ii) == 2L) 1.20 else .55
      for (j in 2:length(oi)) {
        sig$label_radius[oi[j]] <- max(sig$label_radius[oi[j]],
                                       sig$label_radius[oi[j - 1L]] + spacing)
      }
    }
  }
  raw_angle <- 90 - 360 * (sig$index - .5) / total_slots
  sig$label_angle <- ifelse(raw_angle < -90, raw_angle + 180, raw_angle)
  sig$label_hjust <- ifelse(raw_angle < -90, 1, 0)
}

max_logp <- max(c(ring$logP, -log10(.05)), na.rm = TRUE)
outer_limit <- max(c(ring$bar_top, if (nrow(sig)) sig$label_radius else NA_real_),
                   na.rm = TRUE) + .85
grid_values <- seq_len(max(1L, ceiling(max_logp)))
grid_values <- grid_values[bar_base + grid_values <= outer_limit]
grid_data <- data.frame(y = bar_base + grid_values)
grid_labels <- data.frame(x = .5 + total_slots / 4, y = bar_base + grid_values,
                          label = as.character(grid_values))

p <- ggplot() +
  geom_hline(data = grid_data, aes(yintercept = y), linewidth = .28,
             color = "#E3E6E8") +
  geom_rect(data = sectors,
            aes(xmin = xmin, xmax = xmax, ymin = class_ring_inner,
                ymax = class_ring_outer, fill = Class),
            color = "white", linewidth = .45, inherit.aes = FALSE) +
  geom_hline(yintercept = bar_base, linewidth = .30, color = "#B8BDC2") +
  geom_segment(data = ring,
               aes(x = index, xend = index, y = bar_base, yend = bar_top,
                   color = Association),
               linewidth = 1.55, lineend = "butt") +
  geom_hline(yintercept = p05_height, linetype = 2, linewidth = .55,
             color = "#D75B4E") +
  geom_label(data = grid_labels, aes(x = x, y = y, label = label),
             family = "Arial", size = 3, color = "#686D72", fill = "white",
             linewidth = 0, label.padding = unit(.06, "lines"),
             label.r = unit(0, "lines"), inherit.aes = FALSE) +
  {if (nrow(sig)) geom_text(data = sig,
     aes(index, label_radius, label = Metabolite, angle = label_angle,
         hjust = label_hjust), size = 2.55, family = "Arial", color = "#202020",
     inherit.aes = FALSE)} +
  coord_polar(start = 0, clip = "off") +
  scale_x_continuous(limits = c(.5, total_slots + .5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, outer_limit), expand = c(0, 0)) +
  scale_color_manual(values = association_colors, drop = FALSE, guide = "none") +
  scale_fill_manual(values = class_colors, labels = class_labels, drop = TRUE) +
  guides(fill = guide_legend(order = 1, ncol = 2, byrow = TRUE,
                             override.aes = list(color = NA))) +
  labs(fill = "Metabolite class",
       caption = "Bars show -log10(P); dashed circle indicates P = 0.05.") +
  theme_void(base_family = "Arial") +
  theme(
    legend.position = c(.50, .50), legend.justification = c(.50, .50),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 11.5, face = "bold", hjust = .5,
                                family = "Arial"),
    legend.text = element_text(size = 9.5, lineheight = .95, family = "Arial"),
    legend.key.height = unit(.40, "cm"), legend.key.width = unit(.47, "cm"),
    plot.caption = element_text(size = 8, hjust = .5, color = "#4B4B4B",
                                family = "Arial", margin = margin(t = 4)),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(28, 34, 16, 34)
  )

direction_key <- grobTree(
  textGrob("Effect direction", x = .79, y = .94, just = "left",
           gp = gpar(fontfamily = "Arial", fontsize = 14, fontface = "bold")),
  segmentsGrob(x0 = .79, x1 = .83, y0 = .90, y1 = .90,
               gp = gpar(col = association_colors["Negative association (OR < 1)"], lwd = 5)),
  textGrob("OR < 1", x = .85, y = .90, just = "left",
           gp = gpar(fontfamily = "Arial", fontsize = 12)),
  segmentsGrob(x0 = .79, x1 = .83, y0 = .865, y1 = .865,
               gp = gpar(col = association_colors["Positive association (OR > 1)"], lwd = 5)),
  textGrob("OR > 1", x = .85, y = .865, just = "left",
           gp = gpar(fontfamily = "Arial", fontsize = 12))
)
pg <- ggplotGrob(p)
pg <- gtable_add_grob(pg, direction_key, t = 1, l = 1,
                      b = nrow(pg), r = ncol(pg), clip = "off")

ggsave("181_total_PA_screened_metabolites_CRC_ring.png", pg,
       width = 14, height = 13.2, dpi = 600, limitsize = FALSE, bg = "white")
ggsave("181_total_PA_screened_metabolites_CRC_ring.pdf", pg,
       width = 14, height = 13.2, limitsize = FALSE, device = cairo_pdf,
       bg = "white")
write_csv(ring, "03_181_metabolites_ring_plot_data.csv")
write_csv(sectors, "04_181_metabolites_ring_sector_positions.csv")
write_csv(data.frame(
  N_screened_metabolites = nrow(ring), N_nominal_P_lt_0.05 = sum(ring$P < .05),
  N_FDR_P_lt_0.05 = sum(ring$FDR_P < .05),
  N_OR_lt_1_among_nominal = sum(ring$P < .05 & ring$OR < 1),
  N_OR_gt_1_among_nominal = sum(ring$P < .05 & ring$OR > 1),
  Sample = "Lag-1 nested case-control sample: N=256; 128 matched pairs",
  Model = paste("Exact conditional logistic regression, adjusted for",
                paste(crc_covars, collapse = "; "))
), "00_ring_plot_analysis_summary.csv")

cat("Ring plot completed.\n")
print(read.csv("00_ring_plot_analysis_summary.csv", check.names = FALSE))
})

# 2. score and complete analysis.
.integration_run_stage("code_seed20_GPA_only_complete_analysis.R",
  .integration_total, {
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(ggplot2)
  library(Hmisc)
})

write_csv <- function(x, file) write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8")
fmt_est <- function(x) sprintf("%.2f (%.2f-%.2f)", x["OR"], x["Lower"], x["Upper"])
fmt_p <- function(x) ifelse(x < .001, "<0.001", sprintf("%.3f", x))
fit_clogit <- function(f, d) suppressWarnings(clogit(f, data = d, method = "exact", na.action = na.fail))
extract_or <- function(fit, term) {
  z <- summary(fit)$coefficients
  b <- z[term, "coef"]; se <- z[term, "se(coef)"]
  c(OR = exp(b), Lower = exp(b - 1.96 * se), Upper = exp(b + 1.96 * se),
    P = z[term, "Pr(>|z|)"])
}

# Score derivation: N=288; CRC models: lag-1 N=256 (128 pairs).
obj <- readRDS(file.path(.integration_nested, "05_complete_analysis_objects.rds"))
derivation_dat <- obj$derivation_dat
lag1_dat <- obj$lag1_dat
# Initialize from this run's screening and annotation.
screen <- read.csv(file.path(.integration_nested, "1_total_PA_signature",
                            "02_Spearman_FDR_lt_0.05_full_nested.csv"),
                   check.names = FALSE, fileEncoding = "UTF-8")

pa_var <- "total_PA_MET_h_day"
lasso_covars <- c("AGE", "BMI", "A2", "smoking5")
crc_covars <- c("AGE", "BMI", "education", "smoking5", "drink", "DMTIME",
                "H1H2", "redmeat", "fruit", "veg")
candidates <- screen$COMPOUND.Name
gpa <- candidates[grepl("^1-linoleoyl-GPA \\(18:2\\)\\*?$", candidates)]
stopifnot(length(candidates) == 181L, length(gpa) == 1L,
          nrow(derivation_dat) == 288L, nrow(lag1_dat) == 256L,
          length(unique(derivation_dat$Match)) == 144L,
          length(unique(lag1_dat$Match)) == 128L)

# Force GPA with zero penalty; unpenalized covariates are excluded from the score.
met_x <- as.matrix(derivation_dat[, candidates, drop = FALSE])
unpen_x <- model.matrix(~ AGE + BMI + A2 + smoking5, data = derivation_dat)[, -1, drop = FALSE]
met_penalty <- ifelse(colnames(met_x) == gpa, 0, 1)
X <- cbind(met_x, unpen_x)
penalty <- c(met_penalty, rep(0, ncol(unpen_x)))

set.seed(20)
cvfit <- cv.glmnet(X, derivation_dat[[pa_var]], alpha = 1, family = "gaussian",
                   penalty.factor = penalty, nfolds = 10)
fit <- glmnet(X, derivation_dat[[pa_var]], alpha = 1, family = "gaussian",
              penalty.factor = penalty, lambda = cvfit$lambda.min)
cf <- as.matrix(coef(fit))
selected <- intersect(rownames(cf)[cf[, 1] != 0], candidates)
weights <- as.numeric(cf[selected, 1]); names(weights) <- selected
stopifnot(length(selected) == 18L, gpa %in% selected)

annot <- screen[match(selected, screen$COMPOUND.Name),
                c("COMPOUND.Name", "SUPER.META.PATHWAY", "SUB.META.PATHWAY",
                  "HMDB.ID", "KEGG.ID")]
coef_tab <- data.frame(
  Metabolite = selected,
  Super_pathway = annot$SUPER.META.PATHWAY,
  Sub_pathway = annot$SUB.META.PATHWAY,
  HMDB_ID = annot$HMDB.ID,
  KEGG_ID = annot$KEGG.ID,
  LASSO_coefficient = weights,
  Penalty_factor = met_penalty[match(selected, candidates)],
  Forced_inclusion = selected == gpa,
  check.names = FALSE
)
write_csv(coef_tab, "01_signature_metabolite_coefficients.csv")
write_csv(data.frame(Feature = colnames(unpen_x), Penalty_factor = 0),
          "02_LASSO_nonpenalized_covariates.csv")
write_csv(data.frame(Seed = 20, N_candidates = length(candidates),
                     N_selected = length(selected), Forced_metabolite = gpa,
                     Lambda_min = cvfit$lambda.min, Lambda_1se = cvfit$lambda.1se),
          "03_LASSO_specification.csv")

train_score <- as.numeric(as.matrix(derivation_dat[, selected, drop = FALSE]) %*% weights)
lag1_score <- as.numeric(as.matrix(lag1_dat[, selected, drop = FALSE]) %*% weights)
train_mean <- mean(train_score); train_sd <- sd(train_score)
sig_cuts <- as.numeric(quantile(train_score, c(1/3, 2/3), type = 2, names = FALSE))
sig_ranges <- c(sprintf("[%.2f, %.2f]", min(train_score), sig_cuts[1]),
                sprintf("(%.2f, %.2f]", sig_cuts[1], sig_cuts[2]),
                sprintf("(%.2f, %.2f]", sig_cuts[2], max(train_score)))

cor_row <- function(score, pa, sample) {
  z <- cor.test(score, pa, method = "pearson")
  data.frame(Sample = sample, N = length(score), Pearson_r = unname(z$estimate),
             Pearson_95CI_lower = z$conf.int[1], Pearson_95CI_upper = z$conf.int[2],
             Pearson_P = z$p.value, R_squared = unname(z$estimate)^2,
             Variance_explained_percent = 100 * unname(z$estimate)^2)
}
corr <- rbind(cor_row(train_score, derivation_dat[[pa_var]], "Full nested derivation"),
              cor_row(lag1_score, lag1_dat[[pa_var]], "Lag-1 CRC analysis"))
write_csv(corr, "04_signature_total_PA_Pearson_r_and_r2.csv")

d <- lag1_dat
d$signature_score <- lag1_score
d$signature_z <- (lag1_score - train_mean) / train_sd
d$signature_tertile <- cut(lag1_score, breaks = c(-Inf, sig_cuts, Inf),
                           labels = c("T1", "T2", "T3"), include.lowest = TRUE)
d$pa_z <- as.numeric(scale(d[[pa_var]]))
pa_sd <- sd(d[[pa_var]])
pa_cuts <- as.numeric(quantile(d[[pa_var]], c(1/3, 2/3), type = 2, names = FALSE))
d$pa_tertile <- cut(d[[pa_var]], breaks = c(-Inf, pa_cuts, Inf),
                    labels = c("T1", "T2", "T3"), include.lowest = TRUE)
pa_ranges <- c(sprintf("[%.2f, %.2f]", min(d[[pa_var]]), pa_cuts[1]),
               sprintf("(%.2f, %.2f]", pa_cuts[1], pa_cuts[2]),
               sprintf("(%.2f, %.2f]", pa_cuts[2], max(d[[pa_var]])))

sig_sd0 <- fit_clogit(Group ~ signature_z + strata(Match), d)
sig_t0 <- fit_clogit(Group ~ signature_tertile + strata(Match), d)
sig_sd1 <- fit_clogit(as.formula(paste("Group ~ signature_z +",
                                      paste(crc_covars, collapse = " + "), "+ strata(Match)")), d)
sig_t1 <- fit_clogit(as.formula(paste("Group ~ signature_tertile +",
                                     paste(crc_covars, collapse = " + "), "+ strata(Match)")), d)
sig_row <- function(model, msd, mt) {
  es <- extract_or(msd, "signature_z")
  e2 <- extract_or(mt, "signature_tertileT2")
  e3 <- extract_or(mt, "signature_tertileT3")
  data.frame(Variable = "Total PA-related metabolic signature", Model = model,
             T1_range = sig_ranges[1], T1 = "1.00 (ref.)", T2_range = sig_ranges[2],
             T2 = fmt_est(e2), T3_range = sig_ranges[3], T3 = fmt_est(e3),
             `Per 1-SD` = fmt_est(es), P = es["P"], Signature_SD = train_sd,
             Range_and_SD_reference = "Full nested derivation sample (N=288)",
             check.names = FALSE)
}
sig_crc <- rbind(sig_row("Unadjusted", sig_sd0, sig_t0),
                 sig_row("Full model", sig_sd1, sig_t1))
write_csv(sig_crc, "05_signature_CRC_tertiles_and_perSD.csv")
write_csv(data.frame(Signature_SD = train_sd, Tertile_cut1 = sig_cuts[1],
                     Tertile_cut2 = sig_cuts[2], T1_range = sig_ranges[1],
                     T2_range = sig_ranges[2], T3_range = sig_ranges[3],
                     Reference_sample = "Full nested derivation sample (N=288)"),
          "06_signature_SD_and_tertile_ranges.csv")

# Fully adjusted CRC associations per metabolite SD.
individual_crc <- do.call(rbind, lapply(selected, function(m) {
  d$met_z <- as.numeric(scale(d[[m]]))
  f <- fit_clogit(as.formula(paste("Group ~ met_z +",
                                  paste(crc_covars, collapse = " + "),
                                  "+ strata(Match)")), d)
  e <- extract_or(f, "met_z")
  data.frame(Metabolite = m, HMDB_ID = annot$HMDB.ID[match(m, annot$COMPOUND.Name)],
             OR = e["OR"], CI_lower = e["Lower"], CI_upper = e["Upper"],
             `Per 1-SD OR (95% CI)` = fmt_est(e), P = e["P"], check.names = FALSE)
}))
individual_crc$FDR_P <- p.adjust(individual_crc$P, method = "BH")
write_csv(individual_crc, "07_signature_metabolites_CRC_full_model.csv")

# Mutually adjusted PA and score models: tertiles and per SD.
mutual_sd <- fit_clogit(as.formula(paste("Group ~ pa_z + signature_z +",
                                        paste(crc_covars, collapse = " + "),
                                        "+ strata(Match)")), d)
mutual_pa_t <- fit_clogit(as.formula(paste("Group ~ pa_tertile + signature_z +",
                                          paste(crc_covars, collapse = " + "),
                                          "+ strata(Match)")), d)
mutual_sig_t <- fit_clogit(as.formula(paste("Group ~ signature_tertile + pa_z +",
                                           paste(crc_covars, collapse = " + "),
                                           "+ strata(Match)")), d)
mutual_row <- function(variable, term, mt, prefix, ranges, sd_value, sd_ref) {
  es <- extract_or(mutual_sd, term); e2 <- extract_or(mt, paste0(prefix, "T2"));
  e3 <- extract_or(mt, paste0(prefix, "T3"))
  data.frame(Variable = variable, Model = "Mutually adjusted full model",
             T1_range = ranges[1], T1 = "1.00 (ref.)", T2_range = ranges[2],
             T2 = fmt_est(e2), T3_range = ranges[3], T3 = fmt_est(e3),
             `Per 1-SD` = fmt_est(es), P = es["P"], SD = sd_value,
             Range_and_SD_reference = sd_ref, check.names = FALSE)
}
mutual <- rbind(
  mutual_row("Total physical activity", "pa_z", mutual_pa_t, "pa_tertile",
             pa_ranges, pa_sd, "Lag-1 CRC analysis sample (N=256)"),
  mutual_row("Metabolic signature", "signature_z", mutual_sig_t,
             "signature_tertile", sig_ranges, train_sd,
             "Full nested derivation sample (N=288)")
)
write_csv(mutual, "08_PA_and_signature_mutually_adjusted_full_model.csv")

# Change-in-estimate mediation; use Pearson correlation for PA versus score.
pa_total_fit <- fit_clogit(as.formula(paste("Group ~ pa_z +",
                                           paste(crc_covars, collapse = " + "),
                                           "+ strata(Match)")), d)
et <- extract_or(pa_total_fit, "pa_z")
ed <- extract_or(mutual_sd, "pa_z")
em <- extract_or(mutual_sd, "signature_z")
ct <- cor.test(d[[pa_var]], d$signature_score, method = "pearson")
# Adjusted a-path: include sex explicitly in the linear model.
lm_covars <- c("A2", crc_covars)
pa_sig_lm <- lm(as.formula(paste("signature_z ~ pa_z +",
                                  paste(lm_covars, collapse = " + "))),
                data = d, na.action = na.fail)
lm_beta <- unname(coef(pa_sig_lm)["pa_z"])
lm_ci <- unname(confint(pa_sig_lm, "pa_z", level = 0.95)[1, ])
lm_p <- unname(summary(pa_sig_lm)$coefficients["pa_z", "Pr(>|t|)"])
med_prop <- 100 * (coef(pa_total_fit)["pa_z"] - coef(mutual_sd)["pa_z"]) /
  coef(pa_total_fit)["pa_z"]
mediation <- data.frame(
  Exposure = "Total physical activity",
  Mediator = "Total PA-related metabolic signature",
  N = nrow(d), Matched_pairs = length(unique(d$Match)),
  Exposure_mediator_Pearson_r = unname(ct$estimate),
  Pearson_95CI_lower = ct$conf.int[1], Pearson_95CI_upper = ct$conf.int[2],
  Pearson_P = ct$p.value,
  Exposure_mediator_LM_beta_sigSD_per_1SD_PA = lm_beta,
  LM_beta_95CI_lower = lm_ci[1], LM_beta_95CI_upper = lm_ci[2],
  LM_beta_P = lm_p,
  LM_beta_unit = "Derivation-sample signature SD units per 1 lag-1-sample SD total PA",
  Signature_SD_derivation = train_sd,
  PA_SD_lag1_MET_h_day = pa_sd,
  LM_adjustment = paste(lm_covars, collapse = "; "),
  Total_effect_OR = et["OR"], Total_effect_CI_lower = et["Lower"],
  Total_effect_CI_upper = et["Upper"], Total_effect_P = et["P"],
  Direct_effect_OR = ed["OR"], Direct_effect_CI_lower = ed["Lower"],
  Direct_effect_CI_upper = ed["Upper"], Direct_effect_P = ed["P"],
  Signature_OR_adjusted_for_PA = em["OR"],
  Signature_CI_lower_adjusted_for_PA = em["Lower"],
  Signature_CI_upper_adjusted_for_PA = em["Upper"],
  Signature_P_adjusted_for_PA = em["P"],
  Change_in_estimate_percent = unname(med_prop),
  Method_note = "Change-in-estimate estimate; not a formal causal mediation estimand",
  check.names = FALSE
)
write_csv(mediation, "09_mediation_change_in_estimate_full_model.csv")

# Three-knot spline; reference = minimum score.
rd <- d[, c("Group", "Match", "signature_score", crc_covars), drop = FALSE]
basis <- Hmisc::rcspline.eval(rd$signature_score, nk = 3, inclx = TRUE)
knots <- attr(basis, "knots"); colnames(basis) <- paste0("rcs", seq_len(ncol(basis)))
rd <- cbind(rd, as.data.frame(basis)); spline_terms <- colnames(basis)
f_spline <- fit_clogit(as.formula(paste("Group ~", paste(spline_terms, collapse = " + "),
                                        "+", paste(crc_covars, collapse = " + "),
                                        "+ strata(Match)")), rd)
f_linear_raw <- fit_clogit(as.formula(paste("Group ~ signature_score +",
                                            paste(crc_covars, collapse = " + "),
                                            "+ strata(Match)")), rd)
f_null <- fit_clogit(as.formula(paste("Group ~", paste(crc_covars, collapse = " + "),
                                      "+ strata(Match)")), rd)
lrt_p <- function(small, large) {
  stat <- 2 * (as.numeric(logLik(large)) - as.numeric(logLik(small)))
  df <- attr(logLik(large), "df") - attr(logLik(small), "df")
  pchisq(stat, df = df, lower.tail = FALSE)
}
xpred <- seq(min(rd$signature_score), max(rd$signature_score), length.out = 300)
pb <- Hmisc::rcspline.eval(xpred, knots = knots, inclx = TRUE)
rb <- Hmisc::rcspline.eval(min(rd$signature_score), knots = knots, inclx = TRUE)
xd <- sweep(pb, 2, as.numeric(rb), "-")
b <- coef(f_spline)[spline_terms]
v <- vcov(f_spline)[spline_terms, spline_terms, drop = FALSE]
lp <- as.numeric(xd %*% b); se <- sqrt(rowSums((xd %*% v) * xd))
pred <- data.frame(Signature_score = xpred, OR = exp(lp),
                   Lower = exp(lp - 1.96 * se), Upper = exp(lp + 1.96 * se))
linear_p <- extract_or(sig_sd1, "signature_z")["P"]
nonlinear_p <- lrt_p(f_linear_raw, f_spline)
overall_p <- lrt_p(f_null, f_spline)
p <- ggplot(pred, aes(Signature_score, OR)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#5B9BD5", alpha = .25) +
  geom_line(color = "#2C78B8", linewidth = .9) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = .35) +
  theme_classic(base_size = 11) +
  labs(x = "Total physical activity-related metabolic signature score",
       y = "Odds ratio") +
  annotate("text", x = quantile(xpred, .25), y = Inf,
           label = paste0("P for linear: ", fmt_p(linear_p)), vjust = 1.2, size = 3.4) +
  annotate("text", x = quantile(xpred, .25), y = Inf,
           label = paste0("P for non-linear: ", fmt_p(nonlinear_p)), vjust = 3, size = 3.4)
ggsave("10_signature_CRC_RCS_3knots.png", p, width = 6, height = 5.5, dpi = 600, bg = "white")
ggsave("10_signature_CRC_RCS_3knots.pdf", p, width = 6, height = 5.5, device = cairo_pdf)
write_csv(pred, "10a_signature_CRC_RCS_predictions.csv")
write_csv(data.frame(Test = c("Linear association (per 1-SD)",
                              "Overall spline association", "Non-linear association"),
                     P_value = c(linear_p, overall_p, nonlinear_p),
                     Reference = "Minimum signature score", Knots = 3,
                     N = nrow(d), Matched_pairs = length(unique(d$Match))),
          "10b_signature_CRC_RCS_statistics.csv")
write_csv(data.frame(Knot = seq_along(knots), Knot_value = knots),
          "10c_signature_CRC_RCS_knots.csv")

pdf("11_LASSO_cross_validation.pdf", width = 7, height = 6); plot(cvfit); dev.off()
png("11_LASSO_cross_validation.png", width = 2100, height = 1800, res = 300); plot(cvfit); dev.off()
pdf("12_LASSO_coefficient_paths.pdf", width = 7, height = 6); plot(cvfit$glmnet.fit); dev.off()
png("12_LASSO_coefficient_paths.png", width = 2100, height = 1800, res = 300); plot(cvfit$glmnet.fit); dev.off()

write_csv(data.frame(Sample = "Full nested derivation", ID = derivation_dat$ID,
                     Group = derivation_dat$Group, Match = derivation_dat$Match,
                     Signature_score = train_score), "13_signature_scores_full_nested.csv")
write_csv(data.frame(Sample = "Lag-1 CRC analysis", ID = lag1_dat$ID,
                     Group = lag1_dat$Group, Match = lag1_dat$Match,
                     Signature_score = lag1_score), "14_signature_scores_lag1.csv")
saveRDS(list(seed = 20, candidates = candidates, forced_metabolite = gpa,
             selected = selected, weights = weights, cvfit = cvfit, fit = fit,
             train_score = train_score, lag1_score = lag1_score, correlation = corr,
             signature_crc = sig_crc, individual_crc = individual_crc,
             mutual = mutual, mediation = mediation, rcs_fit = f_spline,
             rcs_predictions = pred), "15_complete_analysis_objects.rds")

writeLines(c(
  "Signature derivation: complete nested case-control sample, N=288 (144 matched pairs).",
  "Candidate metabolites: 181 identified by full-nested Spearman BH-FDR P<0.05 screening.",
  "LASSO random seed: 20; only 1-linoleoyl-GPA (18:2) was forced (penalty factor=0).",
  "AGE, BMI, A2, and smoking5 dummy variables were non-penalized covariates and were excluded from the score.",
  "CRC analyses: lag-1 sample, N=256 (128 matched pairs), using exact conditional logistic regression.",
  paste0("Full-model covariates: ", paste(crc_covars, collapse = "; "), "."),
  "RCS: three knots; minimum observed signature score used as reference.",
  "Mediation table uses Pearson correlation and a change-in-estimate proportion; it is not a formal causal mediation estimate."
), "README_analysis_specification.txt", useBytes = TRUE)

cat("Completed seed=20 GPA-only signature analysis.\n")
print(corr)
print(sig_crc)
print(mutual)
print(mediation)
})

# 3. Adjusted linear-model mediation.
.integration_run_stage("code_update_mediation_adjusted_LM.R",
  .integration_total, {
options(stringsAsFactors = FALSE)
source_data <- readRDS(file.path(.integration_nested, "05_complete_analysis_objects.rds"))
analysis <- readRDS("15_complete_analysis_objects.rds")
d <- source_data$lag1_dat
d$signature_score <- analysis$lag1_score
train_sd <- sd(analysis$train_score)
d$signature_z <- (d$signature_score - mean(analysis$train_score)) / train_sd
pa_var <- "total_PA_MET_h_day"
pa_sd <- sd(d[[pa_var]])
d$pa_z <- as.numeric(scale(d[[pa_var]]))
stopifnot(nrow(d) == 256L, length(unique(d$Match)) == 128L,
          length(d$signature_score) == 256L,
          !anyNA(d[, c("signature_z", "pa_z")]), is.finite(train_sd), train_sd > 0,
          is.finite(pa_sd), pa_sd > 0)

crc_covars <- c("AGE", "BMI", "education", "smoking5", "drink",
                "DMTIME", "H1H2", "redmeat", "fruit", "veg")
lm_covars <- c("A2", crc_covars)
lm_formula <- as.formula(paste("signature_z ~ pa_z +",
                               paste(lm_covars, collapse = " + ")))
lm_fit <- lm(lm_formula, data = d, na.action = na.fail)
beta <- unname(coef(lm_fit)["pa_z"])
ci <- unname(confint(lm_fit, "pa_z", level = 0.95)[1, ])
p <- unname(summary(lm_fit)$coefficients["pa_z", "Pr(>|t|)"])
stopifnot(nobs(lm_fit) == 256L, is.finite(beta), all(is.finite(ci)),
          is.finite(p))

file <- "09_mediation_change_in_estimate_full_model.csv"
out <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(out) == 1L, out$N[1] == 256L)
remove_cols <- c("Exposure_mediator_LM_beta_per_1SD_PA",
                 "LM_beta_per_1MET_h_day", "LM_beta_standardized_signature",
                 "LM_beta_standardized_95CI_lower",
                 "LM_beta_standardized_95CI_upper", "PA_SD_MET_h_day",
                 "Signature_SD_derivation", "Exposure_mediator_LM_beta_per_1MET_h_day",
                 "Exposure_mediator_LM_beta_sigSD_per_1MET_h_day",
                 "Exposure_mediator_LM_beta_sigSD_per_1SD_PA", "PA_SD_lag1_MET_h_day")
out <- out[, !(names(out) %in% remove_cols), drop = FALSE]
out$Exposure_mediator_LM_beta_sigSD_per_1SD_PA <- beta
out$LM_beta_95CI_lower <- ci[1]
out$LM_beta_95CI_upper <- ci[2]
out$LM_beta_P <- p
out$LM_beta_unit <- "Derivation-sample signature SD units per 1 lag-1-sample SD total PA"
out$Signature_SD_derivation <- train_sd
out$PA_SD_lag1_MET_h_day <- pa_sd
out$LM_adjustment <- paste(lm_covars, collapse = "; ")
write.csv(out, file, row.names = FALSE, fileEncoding = "UTF-8")
cat("Adjusted lm beta per 1 SD total PA:", beta,
    "95% CI:", ci[1], ci[2], "P:", p, "N:", nobs(lm_fit), "\n")
})

# 4. Separate and mutually adjusted models.
.integration_run_stage("table_S5_mutually_adjusted.R",
  .integration_total, {
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(survival))

if (!file.exists(file.path(.integration_nested, "05_complete_analysis_objects.rds"))) {
  stop("Run this script from its analysis directory.")
}

# Score derivation: N=288; CRC analysis: lag-1 N=256 (128 pairs).
source_obj <- readRDS(file.path(.integration_nested, "05_complete_analysis_objects.rds"))
d <- source_obj$lag1_dat
score_file <- read.csv("14_signature_scores_lag1.csv", check.names = FALSE)
sd_file <- read.csv("06_signature_SD_and_tertile_ranges.csv", check.names = FALSE)

key <- function(z) paste(z$ID, z$Match, z$Group, sep = "|")
score_index <- match(key(d), key(score_file))
stopifnot(nrow(d) == 256L, nrow(score_file) == 256L,
          length(unique(d$Match)) == 128L,
          all(table(d$Match) == 2L),
          all(tapply(d$Group, d$Match, sum) == 1L),
          !anyNA(score_index), !anyDuplicated(key(score_file)),
          nrow(sd_file) == 1L)

pa_sd <- sd(d$total_PA_MET_h_day)
signature_sd <- sd_file$Signature_SD[[1]]
d$pa_z <- (d$total_PA_MET_h_day - mean(d$total_PA_MET_h_day)) / pa_sd
d$signature_z <- score_file$Signature_score[score_index] / signature_sd

covariates <- c("AGE", "BMI", "education", "smoking5", "drink",
                "DMTIME", "H1H2", "redmeat", "fruit", "veg")
model_vars <- c("Group", "Match", "pa_z", "signature_z", covariates)
stopifnot(all(model_vars %in% names(d)), all(complete.cases(d[, model_vars])))

model_formula <- function(exposures) {
  as.formula(paste("Group ~", paste(c(exposures, covariates), collapse = " + "),
                   "+ strata(Match)"))
}
fit_model <- function(exposures) {
  suppressWarnings(clogit(model_formula(exposures), data = d,
                         method = "exact", na.action = na.fail))
}

# Model 1: separate, fully adjusted exposures.
pa_model1 <- fit_model("pa_z")
signature_model1 <- fit_model("signature_z")
# Model 2: both exposures, fully adjusted.
mutual_model2 <- fit_model(c("pa_z", "signature_z"))

estimate <- function(fit, term, model, exposure) {
  z <- summary(fit)$coefficients
  beta <- unname(z[term, "coef"])
  se <- unname(z[term, "se(coef)"])
  data.frame(Model = model, Exposure = exposure, N = nrow(d),
             Matched_pairs = length(unique(d$Match)),
             Beta = beta, SE = se, OR = exp(beta),
             CI_lower = exp(beta - 1.96 * se),
             CI_upper = exp(beta + 1.96 * se),
             P = unname(z[term, "Pr(>|z|)"]),
             Exposure_SD = if (term == "pa_z") pa_sd else signature_sd,
             SD_reference = if (term == "pa_z")
               "Lag-1 analysis sample (N=256)" else
               "Full nested signature derivation sample (N=288)")
}
long <- rbind(
  estimate(pa_model1, "pa_z", "Model 1", "Physical activity"),
  estimate(signature_model1, "signature_z", "Model 1", "Metabolic signature"),
  estimate(mutual_model2, "pa_z", "Model 2", "Physical activity"),
  estimate(mutual_model2, "signature_z", "Model 2", "Metabolic signature")
)

fmt_or <- function(x) sprintf("%.2f (%.2f-%.2f)", x$OR, x$CI_lower, x$CI_upper)
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
wide <- data.frame(
  Model = c("Model 1", "Model 2"),
  `Physical activity: OR (95% CI)` = c(fmt_or(long[1, ]), fmt_or(long[3, ])),
  `Physical activity: P` = fmt_p(long$P[c(1, 3)]),
  `Metabolic signature: OR (95% CI)` = c(fmt_or(long[2, ]), fmt_or(long[4, ])),
  `Metabolic signature: P` = fmt_p(long$P[c(2, 4)]),
  check.names = FALSE
)

write.csv(wide, "TableS5_PA_signature_CRC_model1_model2.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(long, "TableS5_PA_signature_CRC_model1_model2_numeric.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

caption <- paste(
  "Table S5. Associations of physical activity and the metabolic signature with",
  "colorectal cancer risk in the nested case-control study of JBPCD"
)
footnote <- paste0(
  "ORs are per 1-SD increment. All models used exact conditional logistic regression",
  " stratified by the 128 matched case-control pairs (N=256), after excluding CRC",
  " cases diagnosed within one year and their matched controls. Model 1 included",
  " age at baseline, BMI, education, cigarette smoking status, alcohol drinking",
  " status, diabetes duration, oral antidiabetic medication use, red meat",
  " consumption, fruit consumption, and vegetable consumption. Model 2 further",
  " adjusted for the metabolic signature in the physical activity column and for",
  " physical activity in the metabolic signature column. Physical activity",
  " denotes total physical activity. Controls were matched by sex and age",
  " (within 2 years); age was also included as a covariate. One SD was ",
  sprintf("%.2f", pa_sd), " MET-h/day for physical",
  " activity in the lag-1 sample and ", sprintf("%.2f", signature_sd),
  " score units for the metabolic signature in the full nested derivation sample."
)
md <- c(
  caption, "",
  "| Model | Physical activity, OR (95% CI) | P | Metabolic signature, OR (95% CI) | P |",
  "|:--|:--|:--|:--|:--|",
  sprintf("| %s | %s | %s | %s | %s |",
          wide$Model, wide[[2]], wide[[3]], wide[[4]], wide[[5]]),
  "", footnote,
  "", "Abbreviations: BMI, body mass index; CI, confidence interval; CRC, colorectal cancer;",
  "JBPCD, Jiangsu Biobank for the Prevention and Control of Diabetes;",
  "OR, odds ratio; SD, standard deviation."
)
writeLines(md, "TableS5_PA_signature_CRC_model1_model2.md", useBytes = TRUE)
cat(paste(md, collapse = "\n"), "\n")
})

# 5. MSEA for three pathways.
.integration_run_stage("code_MSEA_three_pathways_total_PA_CRC.R",
  file.path(.integration_total, "MSEA_three_pathways"), {
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(clusterProfiler)
  library(enrichplot)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(BiocParallel)
})
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))
BiocParallel::register(BiocParallel::SerialParam())

set.seed(20260922)
out_dir <- "."
# Use local SMPDB copies for Windows path compatibility.
smpdb_dir <- file.path(.integration_data, "pathways")
write_csv <- function(x, f) write.csv(x, file.path(out_dir, f), row.names = FALSE,
                                      fileEncoding = "UTF-8")

targets <- data.frame(
  Panel = c("A", "B", "C"),
  SMPDB_ID = c("SMP0000456", "SMP0000035", "SMP0000071"),
  User_ID = c("SMP00456", "SMP00035", "SMP00071"),
  Pathway = c("Fatty Acid Biosynthesis", "Bile Acid Biosynthesis",
              "Ketone Body Metabolism"),
  stringsAsFactors = FALSE
)

# Use saved derivation and lag-1 samples.
obj <- readRDS(file.path(.integration_nested, "05_complete_analysis_objects.rds"))
full_dat <- obj$derivation_dat
lag1_dat <- obj$lag1_dat
stopifnot(nrow(full_dat) == 288L, length(unique(full_dat$Match)) == 144L,
          nrow(lag1_dat) == 256L, length(unique(lag1_dat$Match)) == 128L)

# Metabolite annotation.
annotation_file <- file.path(.integration_nested, "1_total_PA_signature",
                             "01_Spearman_all_799_metabolites_full_nested.csv")
ann <- read.csv(annotation_file, check.names = FALSE, fileEncoding = "UTF-8")
stopifnot(nrow(ann) == 799L)
metabolites <- ann$COMPOUND.Name
stopifnot(all(metabolites %in% names(full_dat)), all(metabolites %in% names(lag1_dat)))

# PA ranking: N=288; CRC ranking: fully adjusted lag-1 conditional logistic model.
pa_rank_all <- data.frame(
  COMPOUND.Name = metabolites,
  HMDB.ID = ann$HMDB.ID,
  Spearman_rho = vapply(metabolites, function(m) {
    suppressWarnings(cor(full_dat$total_PA_MET_h_day, full_dat[[m]],
                         method = "spearman", use = "complete.obs"))
  }, numeric(1)),
  stringsAsFactors = FALSE
)

crc_covars <- c("AGE", "BMI", "education", "smoking5", "drink", "DMTIME",
                "H1H2", "redmeat", "fruit", "veg")
crc_rank_all <- data.frame(
  COMPOUND.Name = metabolites, HMDB.ID = ann$HMDB.ID,
  logOR_per_SD = NA_real_, OR_per_SD = NA_real_, CI_lower = NA_real_,
  CI_upper = NA_real_, P = NA_real_, stringsAsFactors = FALSE
)
for (i in seq_along(metabolites)) {
  m <- metabolites[i]
  f <- as.formula(paste0("Group ~ scale(`", m, "`) + ",
                         paste(crc_covars, collapse = " + "), " + strata(Match)"))
  fit <- suppressWarnings(clogit(f, data = lag1_dat, method = "exact",
                                 na.action = na.fail))
  z <- summary(fit)$coefficients[1, ]
  crc_rank_all$logOR_per_SD[i] <- z["coef"]
  crc_rank_all$OR_per_SD[i] <- exp(z["coef"])
  crc_rank_all$CI_lower[i] <- exp(z["coef"] - 1.96 * z["se(coef)"])
  crc_rank_all$CI_upper[i] <- exp(z["coef"] + 1.96 * z["se(coef)"])
  crc_rank_all$P[i] <- z["Pr(>|z|)"]
}
crc_rank_all$FDR_P <- p.adjust(crc_rank_all$P, method = "BH")
write_csv(pa_rank_all, "01_all_metabolites_total_PA_Spearman_full_nested.csv")
write_csv(crc_rank_all, "02_all_metabolites_CRC_full_model_lag1.csv")

# Retain the reference GSEA duplicate-HMDB rules.
make_unique_rank <- function(x, stat) {
  z <- x[!is.na(x$HMDB.ID) & nzchar(x$HMDB.ID) & is.finite(x[[stat]]),
         c("COMPOUND.Name", "HMDB.ID", stat), drop = FALSE]
  z <- z[z$COMPOUND.Name != "N-delta-acetylornithine", , drop = FALSE]
  z <- z[!z$HMDB.ID %in% c("HMDB0240429", "HMDB0010383"), , drop = FALSE]
  if (anyDuplicated(z$HMDB.ID)) stop("Duplicate HMDB IDs remain in ranked list.")
  names(z)[3] <- "Rank_statistic"
  z[order(z$Rank_statistic, decreasing = TRUE), ]
}
pa_rank <- make_unique_rank(pa_rank_all, "Spearman_rho")
crc_rank <- make_unique_rank(crc_rank_all, "logOR_per_SD")
write_csv(pa_rank, "03_total_PA_ranked_unique_HMDB.csv")
write_csv(crc_rank, "04_CRC_ranked_unique_HMDB.csv")

pathway_members <- do.call(rbind, lapply(seq_len(nrow(targets)), function(i) {
  f <- file.path(smpdb_dir, paste0(targets$SMPDB_ID[i], "_metabolites.csv"))
  if (!file.exists(f)) stop("SMPDB file not found: ", f)
  x <- as.data.frame(fread(f, encoding = "UTF-8", fill = TRUE, quote = "\"",
                           showProgress = FALSE), check.names = FALSE)
  req <- c("SMPDB ID", "Pathway Name", "Metabolite Name", "HMDB ID")
  stopifnot(all(req %in% names(x)))
  x <- x[!is.na(x[["HMDB ID"]]) & nzchar(x[["HMDB ID"]]), req, drop = FALSE]
  x <- x[!duplicated(x[["HMDB ID"]]), , drop = FALSE]
  names(x) <- c("SMPDB_ID", "Pathway", "SMPDB_metabolite", "HMDB_ID")
  x$SMPDB_ID <- targets$SMPDB_ID[i]
  x$Pathway <- targets$Pathway[i]
  x
}))
term2gene <- unique(pathway_members[, c("SMPDB_ID", "HMDB_ID")])
names(term2gene) <- c("term", "gene")
term2name <- targets[, c("SMPDB_ID", "Pathway")]
names(term2name) <- c("term", "name")

measured <- unique(pa_rank[, c("COMPOUND.Name", "HMDB.ID")])
names(measured) <- c("Measured_metabolite", "HMDB_ID")
overlap_members <- merge(pathway_members, measured, by = "HMDB_ID")
overlap_members$Order <- match(overlap_members$SMPDB_ID, targets$SMPDB_ID)
overlap_members <- overlap_members[order(overlap_members$Order, overlap_members$HMDB_ID),
                                   setdiff(names(overlap_members), "Order")]
write_csv(overlap_members, "05_pathway_members_measured_in_JBPCD.csv")
overlap <- do.call(rbind, lapply(seq_len(nrow(targets)), function(i) {
  id <- targets$SMPDB_ID[i]
  data.frame(Panel = targets$Panel[i], SMPDB_ID = id, User_ID = targets$User_ID[i],
             Pathway = targets$Pathway[i],
             SMPDB_members = length(unique(pathway_members$HMDB_ID[pathway_members$SMPDB_ID == id])),
             Measured_members = length(unique(overlap_members$HMDB_ID[overlap_members$SMPDB_ID == id])))
}))
write_csv(overlap, "06_pathway_overlap_summary.csv")

run_msea <- function(rank_df) {
  gl <- rank_df$Rank_statistic; names(gl) <- rank_df$HMDB.ID
  gl <- sort(gl, decreasing = TRUE)
  set.seed(20260922)
  GSEA(geneList = gl, TERM2GENE = term2gene, TERM2NAME = term2name,
       minGSSize = 2, maxGSSize = 500, pvalueCutoff = 1,
       pAdjustMethod = "none", eps = 0, verbose = FALSE, seed = TRUE, by = "fgsea")
}
pa_gsea <- run_msea(pa_rank)
crc_gsea <- run_msea(crc_rank)

tidy <- function(g, label) {
  x <- as.data.frame(g)
  x <- x[match(targets$SMPDB_ID, x$ID), , drop = FALSE]
  data.frame(Analysis = label, Panel = targets$Panel, SMPDB_ID = targets$SMPDB_ID,
             User_ID = targets$User_ID, Pathway = targets$Pathway,
             Set_size = x$setSize, Enrichment_score = x$enrichmentScore,
             NES = x$NES, P = x$pvalue, P_adjusted = x$p.adjust,
             Rank_at_max = x$rank, Leading_edge = x$leading_edge,
             Core_enrichment = x$core_enrichment, stringsAsFactors = FALSE)
}
pa_res <- tidy(pa_gsea, "Total PA (Spearman rho; full nested N=288)")
crc_res <- tidy(crc_gsea, "CRC (fully adjusted log OR; lag-1 N=256)")
write_csv(pa_res, "07_MSEA_total_PA_three_pathways.csv")
write_csv(crc_res, "08_MSEA_CRC_three_pathways.csv")
write_csv(rbind(pa_res, crc_res), "09_MSEA_three_pathways_combined_results.csv")

make_plot <- function(gsea_obj, res, i, analysis_stub) {
  id <- targets$SMPDB_ID[i]; nm <- targets$Pathway[i]; panel <- targets$Panel[i]
  nes <- res$NES[i]; pv <- res$P[i]
  pexpr <- if (pv < .001) "italic(P) < 0.001" else sprintf("italic(P) == %.3f", pv)
  ann_expr <- sprintf("atop('NES = %.2f', %s)", nes, pexpr)
  parts <- gseaplot2(gsea_obj, geneSetID = id, title = nm,
                     base_size = 12, pvalue_table = FALSE)
  parts[[1]] <- parts[[1]] +
    labs(title = paste0("(", panel, ") ", nm)) +
    annotate("text", x = Inf, y = Inf, label = ann_expr, parse = TRUE,
             hjust = 1.08, vjust = 1.10, size = 3.4, family = "Arial",
             lineheight = 1.05) +
    theme(text = element_text(family = "Arial"),
          plot.title = element_text(family = "Arial", face = "plain", size = 11,
                                    hjust = 0),
          plot.margin = margin(6, 12, 0, 6))
  parts[[2]] <- parts[[2]] + theme(text = element_text(family = "Arial"))
  parts[[3]] <- parts[[3]] + theme(text = element_text(family = "Arial"))
  p <- wrap_plots(parts, ncol = 1, heights = c(1.45, .45, 1))
  stub <- paste0(analysis_stub, "_", targets$User_ID[i], "_",
                 gsub("[^A-Za-z0-9]+", "_", nm))
  ggsave(paste0(stub, ".png"), p, width = 8.5, height = 7.2,
         dpi = 600, bg = "white")
  ggsave(paste0(stub, ".pdf"), p, width = 8.5, height = 7.2,
         device = cairo_pdf)
  p
}

pa_plots <- lapply(seq_len(nrow(targets)), function(i) make_plot(pa_gsea, pa_res, i, "Total_PA"))
crc_plots <- lapply(seq_len(nrow(targets)), function(i) make_plot(crc_gsea, crc_res, i, "CRC"))
pa_combined <- wrap_plots(pa_plots, nrow = 1)
crc_combined <- wrap_plots(crc_plots, nrow = 1)
ggsave("MSEA_total_PA_three_pathways_combined.png", pa_combined,
       width = 18, height = 6.4, dpi = 600, bg = "white")
ggsave("MSEA_total_PA_three_pathways_combined.pdf", pa_combined,
       width = 18, height = 6.4, device = cairo_pdf)
ggsave("MSEA_CRC_three_pathways_combined.png", crc_combined,
       width = 18, height = 6.4, dpi = 600, bg = "white")
ggsave("MSEA_CRC_three_pathways_combined.pdf", crc_combined,
       width = 18, height = 6.4, device = cairo_pdf)

write_csv(data.frame(
  Item = c("Total PA ranking sample", "CRC ranking sample", "Total PA statistic",
           "CRC statistic", "CRC covariates", "Matching", "MSEA minGSSize",
           "MSEA maxGSSize", "P-value type", "Font"),
  Value = c("Full nested N=288 (144 matched pairs)",
            "Lag-1 N=256 (128 matched pairs)", "Spearman rho",
            "Fully adjusted conditional-logistic log OR per 1 SD",
            paste(crc_covars, collapse = "; "), "strata(Match)", "2", "500",
            "Nominal P; targeted exploratory MSEA", "Arial")
), "00_MSEA_analysis_metadata.csv")

saveRDS(list(targets = targets, pathway_members = pathway_members,
             overlap_members = overlap_members, pa_rank = pa_rank, crc_rank = crc_rank,
             pa_gsea = pa_gsea, crc_gsea = crc_gsea, pa_results = pa_res,
             crc_results = crc_res), "10_MSEA_analysis_objects.rds")
cat("MSEA completed.\n"); print(pa_res[, c("Pathway", "Set_size", "NES", "P")]);
print(crc_res[, c("Pathway", "Set_size", "NES", "P")])
})
