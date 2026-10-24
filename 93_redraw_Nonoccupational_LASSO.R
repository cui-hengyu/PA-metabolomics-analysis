# Draw only the non-occupational LASSO AB figure from the saved model.
# No model fitting, re-screening, or seed changes.
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(args) == 1L)
run_root <- dirname(normalizePath(sub("^--file=", "", args), winslash = "/"))
.libPaths(c(file.path(dirname(run_root), "_library"), .libPaths()))
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))
exprs <- parse(file.path(run_root, "05_leisure_pa_metabolite_score.R"), encoding = "UTF-8")
is_stage <- vapply(exprs, function(x)
  is.call(x) && identical(x[[1L]], as.name(".integration_run_stage")), logical(1))
for (i in seq_len(which(is_stage)[1L]-1L)) eval(exprs[[i]], envir = .GlobalEnv)
target <- which(vapply(exprs, function(x) is.call(x) &&
  identical(x[[1L]], as.name(".integration_run_stage")) &&
  identical(x[[2L]], "code_LASSO_regression_metabolite_selection_AB.R"), logical(1)))
stopifnot(length(target) == 1L)
eval(exprs[[target]], envir = .GlobalEnv)
stem <- file.path(run_root, "results/analysis/nested/nonoccupational_PA/forced_15_seed257_complete_analysis/17_LASSO_regression_metabolite_selection_AB")
stopifnot(file.copy(paste0(stem, ".pdf"), file.path(run_root, "figures/Nonoccupational_PA_LASSO_AB.pdf"), overwrite = TRUE),
          file.copy(paste0(stem, ".png"), file.path(run_root, "figures/Nonoccupational_PA_LASSO_AB_600dpi.png"), overwrite = TRUE))
cat("Added panel A top frame line; final PDF and 600-dpi PNG synchronized.\n")
