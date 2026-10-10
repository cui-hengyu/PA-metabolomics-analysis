# Smoking GVIF checks.
# Rebuild inputs from raw Excel; retain the legacy smoking model.

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

# 1. Smoking covariate GVIF.
.integration_run_stage("calculate_original_smoking_gvif.R",
  file.path(.integration_results, "gvif"), {
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(survival)
  library(car)
})

source(file.path(.integration_code_root, "01_data_preprocessing.R"),
       local = TRUE, encoding = "UTF-8")
stopifnot(exists("data1"), exists("metarowdata_0.2_imputed1"))

dat <- metarowdata_0.2_imputed1
idx <- match(dat$ID, data1$ID)
stopifnot(!anyNA(idx), !anyDuplicated(data1$ID))

# Smoking status and pack-years.
data1$smoke_original <- NA_real_
data1$smoke_original[!is.na(data1$B3)] <- 1
data1$smoke_original[data1$B1 == 2] <- 2
data1$smoke_original[is.na(data1$smoke_original)] <- 999

for (v in c("AGE", "BMI", "education", "drink", "DMTIME", "H1H2",
            "redmeat", "fruit", "veg", "packyear1")) {
  dat[[v]] <- data1[[v]][idx]
}
dat$smoke_original <- factor(data1$smoke_original[idx])
dat$education <- factor(dat$education)
dat$drink <- factor(dat$drink)
dat$H1H2 <- factor(dat$H1H2)

required <- c("Group", "Match", "MET", "AGE", "BMI", "education",
              "smoke_original", "packyear1", "drink", "DMTIME", "H1H2",
              "redmeat", "fruit", "veg")
cat("N before complete-case filtering:", nrow(dat), "\n")
cat("Missing counts:\n")
print(colSums(is.na(dat[, required])))
dat <- droplevels(dat[complete.cases(dat[, required]), ])
cat("N used:", nrow(dat), "; matched sets:", length(unique(dat$Match)), "\n")
cat("Smoking status by pack-year summary:\n")
print(aggregate(packyear1 ~ smoke_original, dat, function(x)
  c(N = length(x), Mean = mean(x), SD = sd(x), Min = min(x), Max = max(x))))

fit <- clogit(
  Group ~ scale(MET) + AGE + BMI + education + smoke_original + packyear1 +
    drink + DMTIME + H1H2 + redmeat + fruit + veg + strata(Match),
  data = dat, method = "exact", x = TRUE, model = TRUE
)

v <- car::vif(fit)
cat("\ncar::vif(clogit model):\n")
print(v)

if (is.matrix(v)) {
  gvif <- data.frame(
    Variable = rownames(v), GVIF = v[, "GVIF"], Df = v[, "Df"],
    GVIF_adjusted = v[, "GVIF^(1/(2*Df))"], row.names = NULL,
    check.names = FALSE
  )
} else {
  gvif <- data.frame(
    Variable = names(v), GVIF = as.numeric(v), Df = 1,
    GVIF_adjusted = sqrt(as.numeric(v)), row.names = NULL,
    check.names = FALSE
  )
}
gvif$VIF_equivalent <- gvif$GVIF_adjusted^2
gvif$Flag_adjusted_ge_sqrt5 <- gvif$GVIF_adjusted >= sqrt(5)
gvif$Flag_adjusted_ge_5 <- gvif$GVIF_adjusted >= 5

cat("\nStandardized GVIF table:\n")
print(gvif, digits = 6, row.names = FALSE)
cat("\nSmoking rows:\n")
print(gvif[grepl("smoke|packyear", gvif$Variable, ignore.case = TRUE), ],
      digits = 6, row.names = FALSE)

write.csv(gvif,
          file.path(.integration_results, "gvif", "original_smoking_GVIF_results.csv"),
          row.names = FALSE)
})
