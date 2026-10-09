# Data preprocessing.

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

# Reuse only preprocessing generated from the bundled inputs and current code.
.prep_file <- file.path(.integration_results, "preprocessing", "prepared_inputs.rds")
.prep_inputs <- file.path(.integration_data, c("cohort.xlsx", "metabolomics.xlsx"))
if (!all(file.exists(.prep_inputs))) stop("Missing bundled raw Excel inputs.")
.prep_signature <- unname(tools::md5sum(c(.prep_inputs,
  file.path(.integration_code_root, "01_data_preprocessing.R"))))
.prep_cache <- if (file.exists(.prep_file)) readRDS(.prep_file) else NULL
if (!is.null(.prep_cache) && identical(.prep_cache$signature, .prep_signature)) {
  list2env(.prep_cache$objects, envir = environment())
  cohort_file <- .prep_inputs[1L]
  metabolomics_file <- .prep_inputs[2L]
  message("[Preparation] Using current-input preprocessing cache.")
} else {
# 1. Rebuild cohort and metabolomics data from raw Excel files.
.integration_run_stage("00_从原始Excel准备数据.R",
  file.path(.integration_results, "preprocessing"), {
# Apply sample exclusions and screening rules in downstream analyses.
suppressPackageStartupMessages({library(openxlsx);library(lubridate)})
local({
  preparation_dir <- file.path(.integration_results, "preprocessing")
  cohort_file <- file.path(.integration_data, "cohort.xlsx")
  metabolomics_file <- file.path(.integration_data, "metabolomics.xlsx")
  stopifnot(file.exists(cohort_file), file.exists(metabolomics_file))
  message("Reading raw cohort Excel: ", cohort_file)
  message("Reading raw metabolomics Excel: ", metabolomics_file)
  # 1.1 Cohort cleaning.
  # Retain original variable definitions and fixed imputation values.
  # MET hours: original total PA; MET: values outside mean +/- 6 SD set to missing.
  if (!exists("cohort_file", inherits = FALSE)) {
    cohort_file <- file.path(.integration_data, "cohort.xlsx")
  }
  data <- openxlsx::read.xlsx(cohort_file)
  data1 <- subset(data,cancer_基线==0)  # Exclude baseline cancer.
  # Cancer: ICD-10 C00-C97, excluding C44.
  data1$total_cancerfb <- ifelse(data1$ICD10_cancer == "D18" | data1$ICD10_cancer == "D18.0" | data1$ICD10_cancer == "D32.0" | data1$ICD10_cancer == "D32.9" | 
                                   data1$ICD10_cancer == "D43.2" | data1$ICD10_cancer == "D46.9" | data1$ICD10_cancer == "D61.9" | data1$ICD10_cancer == "D64.9" | 
                                   data1$ICD10_cancer == "C44.1" | data1$ICD10_cancer == "C44.3" | data1$ICD10_cancer == "C44.4" | data1$ICD10_cancer == "C44.5" |
                                   data1$ICD10_cancer == "C44.7" | data1$ICD10_cancer == "C44.9" | is.na(data1$ICD10_cancer),0,1)
  # Colorectal cancer: C18-C20.
  data1$colo_cancerfb <- ifelse(data1$ICD10_cancer == "C18.0" | data1$ICD10_cancer == "C18.1" | data1$ICD10_cancer == "C18.2" | data1$ICD10_cancer == "C18.4" | 
                                  data1$ICD10_cancer == "C18.6" | data1$ICD10_cancer == "C18.7" | data1$ICD10_cancer == "C18.9" | data1$ICD10_cancer=="C19" | 
                                  data1$ICD10_cancer == "C20",1,0)
  data1$colo_cancerfb[is.na(data1$colo_cancerfb)] <- 0
  # Liver cancer: C22.
  data1$liver_cancerfb <- ifelse(data1$ICD10_cancer == "C22" | data1$ICD10_cancer == "C22.0" | data1$ICD10_cancer == "C22.1" | data1$ICD10_cancer == "C22.9",1,0)
  data1$liver_cancerfb[is.na(data1$liver_cancerfb)] <- 0
  # Lung cancer: C34.
  data1$lung_cancerfb <- ifelse(data1$ICD10_cancer == "C34" | data1$ICD10_cancer == "C34.0" | data1$ICD10_cancer == "C34.1" | data1$ICD10_cancer == "C34.2" | 
                                  data1$ICD10_cancer == "C34.3" | data1$ICD10_cancer == "C34.8" | data1$ICD10_cancer == "C34.9",1,0)
  data1$lung_cancerfb[is.na(data1$lung_cancerfb)] <- 0
  # Stomach cancer: C16.
  data1$gastric_cancerfb <- ifelse(data1$ICD10_cancer == "C16" | data1$ICD10_cancer == "C16.0" | data1$ICD10_cancer == "C16.1" | data1$ICD10_cancer == "C16.2"|
                                     data1$ICD10_cancer == "C16.3" | data1$ICD10_cancer == "C16.5" | data1$ICD10_cancer == "C16.8" | data1$ICD10_cancer == "C16.9",1,0)
  data1$gastric_cancerfb[is.na(data1$gastric_cancerfb)] <- 0
  # Esophageal cancer: C15.
  data1$eso_cancerfb <-  ifelse(data1$ICD10_cancer == "C15" | data1$ICD10_cancer == "C15.1" | data1$ICD10_cancer == "C15.3" | data1$ICD10_cancer == "C15.4" | 
                                  data1$ICD10_cancer == "C15.5" | data1$ICD10_cancer == "C15.8" | data1$ICD10_cancer == "C15.9",1,0) 
  data1$eso_cancerfb[is.na(data1$eso_cancerfb)] <- 0
  # Breast cancer: C50.
  data1$breast_cancerfb <- ifelse(data1$ICD10_cancer == "C50" | data1$ICD10_cancer == "C50.0" | data1$ICD10_cancer == "C50.2" | data1$ICD10_cancer == "C50.3" |
                                    data1$ICD10_cancer == "C50.4" | data1$ICD10_cancer == "C50.9",1,0)
  data1$breast_cancerfb[is.na(data1$breast_cancerfb)] <- 0
  # Pancreatic cancer: C25.
  data1$pancreas_cancerfb <- ifelse(data1$ICD10_cancer == "C25.0" | data1$ICD10_cancer == "C25.1" | data1$ICD10_cancer == "C25.2" | data1$ICD10_cancer == "C25.4" | data1$ICD10_cancer == "C25.9",1,0)
  data1$pancreas_cancerfb[is.na(data1$pancreas_cancerfb)] <- 0

  # Follow-up times.
  data1$cancer_发病时间 <- as.Date(data1$cancer_发病时间, origin = "1899-12-30")
  data1$DEA1 <- as.Date(data1$DEA1, origin = "1899-12-30")
  data1$surveydate2013 <- as.Date(data1$surveydate2013, origin = "1899-12-30")

  data1$cancer_发病时间 <- as.character(data1$cancer_发病时间)
  data1$DEA1 <- as.character(data1$DEA1)

  data1$total_cancerfb_time <- ifelse(data1$total_cancerfb == 1,data1$cancer_发病时间,data1$DEA1)
  data1$total_cancerfb_time <- as.Date(data1$total_cancerfb_time)
  data1$endtime1 <- time_length(interval(data1$surveydate2013,data1$total_cancerfb_time),'year')
  data1$lung_cancerfb_time <- ifelse(data1$lung_cancerfb == 1,data1$cancer_发病时间,data1$DEA1)
  data1$lung_cancerfb_time <- as.Date(data1$lung_cancerfb_time)
  data1$endtime2 <- time_length(interval(data1$surveydate2013,data1$lung_cancerfb_time),'year')
  data1$gastric_cancerfb_time <- ifelse(data1$gastric_cancerfb == 1,data1$cancer_发病时间,data1$DEA1)
  data1$gastric_cancerfb_time <- as.Date(data1$gastric_cancerfb_time)
  data1$endtime3 <- time_length(interval(data1$surveydate2013,data1$gastric_cancerfb_time),'year')
  data1$eso_cancerfb_time <- ifelse(data1$eso_cancerfb == 1,data1$cancer_发病时间,data1$DEA1)
  data1$eso_cancerfb_time <- as.Date(data1$eso_cancerfb_time)
  data1$endtime4 <- time_length(interval(data1$surveydate2013,data1$eso_cancerfb_time),'year')
  data1$colo_cancerfb_time <- ifelse(data1$colo_cancerfb == 1,data1$cancer_发病时间,data1$DEA1)
  data1$colo_cancerfb_time <- as.Date(data1$colo_cancerfb_time)
  data1$endtime5 <- time_length(interval(data1$surveydate2013,data1$colo_cancerfb_time),'year')
  data1$liver_cancerfb_time <- ifelse(data1$liver_cancerfb == 1,data1$cancer_发病时间,data1$DEA1)
  data1$liver_cancerfb_time <- as.Date(data1$liver_cancerfb_time)
  data1$endtime6 <- time_length(interval(data1$surveydate2013,data1$liver_cancerfb_time),'year')
  data1$breast_cancerfb_time <- ifelse(data1$breast_cancerfb == 1,data1$cancer_发病时间,data1$DEA1)
  data1$breast_cancerfb_time <- as.Date(data1$breast_cancerfb_time)
  data1$endtime7 <- time_length(interval(data1$surveydate2013,data1$breast_cancerfb_time),'year')
  data1$pancreas_cancerfb_time <- ifelse(data1$pancreas_cancerfb == 1,data1$cancer_发病时间,data1$DEA1)
  data1$pancreas_cancerfb_time <- as.Date(data1$pancreas_cancerfb_time)
  data1$endtime8 <- time_length(interval(data1$surveydate2013,data1$pancreas_cancerfb_time),'year')

  # Legacy MET variable.
  sd_MET <- sd(data1$MET小时)
  data1$MET <- ifelse(data1$MET小时 >= mean(data1$MET小时)-6*sd_MET&data1$MET小时 <= mean(data1$MET小时)+6*sd_MET,data1$MET小时,NA)

  # Covariates.
  # BMI: impute missing values with sex-specific medians.
  data1$BMI[is.na(data1$BMI) & (data1$A2==1)] <- "25.15590"
  data1$BMI[is.na(data1$BMI) & (data1$A2==2)] <- "25.11721"
  data1$BMI <- as.numeric(data1$BMI)

  # Education: none, primary, middle, high school or above.
  data1$education[data1$A4 == 1] <- 1
  data1$education[data1$A4 == 2 | data1$A4 == 3] <- 2
  data1$education[data1$A4 == 4] <- 3
  data1$education[data1$A4 == 5 | data1$A4 == 6 | data1$A4 == 7 | data1$A4 == 8] <- 4
  data1$education[is.na(data1$education)] <- 999

  # Smoking status.
  data1$smoke[!is.na(data1$B3)] <- 1
  data1$smoke[data1$B1 == 2] <- 2
  data1$smoke[is.na(data1$smoke)] <- 999

  # Pack-years: packs/day x years smoked; impute with sex-specific medians.
  data1$smoke_pack <- data1$xyl/20
  data1$smoke_year <- data1$AGE - data1$B2
  data1$packyear <- data1$smoke_pack * data1$smoke_year
  data1$packyear1 <- ifelse(data1$smoke == 2,0,data1$packyear)
  data1$packyear1[is.na(data1$packyear1) & (data1$A2==1)] <- "27.75"
  data1$packyear1[is.na(data1$packyear1) & (data1$A2==2)] <- "18.55"
  data1$packyear1 <- as.numeric(data1$packyear1)
  # Alcohol use.
  data1$drink[data1$C1 == 1|data1$C1 == 2] <- 1
  data1$drink[data1$C1 == 3] <- 2
  data1$drink[is.na(data1$drink)] <- 999

  # Diabetes medication.
  data1$H1H2[is.na(data1$H1H2)] <- 999

  # Red meat intake (g/day).

  data1$redmeat <- NA
  data1$redmeat<- ifelse(is.na(data1$redmeat)==T,data1$F4B*data1$F4F*50,data1$redmeat)
  data1$redmeat<- ifelse(is.na(data1$redmeat)==T,data1$F4C*data1$F4F*50/7,data1$redmeat)
  data1$redmeat<- ifelse(is.na(data1$redmeat)==T,data1$F4D*data1$F4F*50/30,data1$redmeat)
  data1$redmeat<- ifelse(is.na(data1$redmeat)==T,data1$F4E*data1$F4F*50/365,data1$redmeat)
  data1$redmeat[is.na(data1$redmeat) & (data1$A2==1)&(data1$F4A)==1] <- "35.71429"
  data1$redmeat[is.na(data1$redmeat) & (data1$A2==2)&(data1$F4A)==1] <- '28.57143'
  data1$redmeat <- ifelse(data1$F4A==2,'0',data1$redmeat)
  data1$redmeat <- as.numeric(data1$redmeat)
  data1$redmeat[is.na(data1$redmeat)&(data1$A2==1)] <- '35.71429'
  data1$redmeat[is.na(data1$redmeat)&(data1$A2==2)] <- '25.00000'

  # Fruit intake (g/day).
  data1$fruit <- NA

  data1$fruit<- ifelse(is.na(data1$fruit)==T,data1$F9B*data1$F9F*50,data1$fruit)
  data1$fruit<- ifelse(is.na(data1$fruit)==T,data1$F9C*data1$F9F*50/7,data1$fruit)
  data1$fruit<- ifelse(is.na(data1$fruit)==T,data1$F9D*data1$F9F*50/30,data1$fruit)
  data1$fruit<- ifelse(is.na(data1$fruit)==T,data1$F9E*data1$F9F*50/365,data1$fruit)
  data1$fruit[is.na(data1$fruit) & (data1$A2==1)&(data1$F9A)==1] <- "28.57143"
  data1$fruit[is.na(data1$fruit) & (data1$A2==2)&(data1$F9A)==1] <- '21.42857'
  data1$fruit <- ifelse(data1$F9A==2,'0',data1$fruit)
  data1$fruit <- as.numeric(data1$fruit)
  data1$fruit[is.na(data1$fruit)&(data1$A2==1)] <- '7.142857'
  data1$fruit[is.na(data1$fruit)&(data1$A2==2)] <- '7.142857'

  # Vegetable intake (g/day).
  data1$veg <- NA

  data1$veg<- ifelse(is.na(data1$veg)==T,data1$F10B*data1$F10F*50,data1$veg)
  data1$veg<- ifelse(is.na(data1$veg)==T,data1$F10C*data1$F10F*50/7,data1$veg)
  data1$veg<- ifelse(is.na(data1$veg)==T,data1$F10D*data1$F10F*50/30,data1$veg)
  data1$veg<- ifelse(is.na(data1$veg)==T,data1$F10E*data1$F10F*50/365,data1$veg)
  data1$veg[is.na(data1$veg) & (data1$A2==1)&(data1$F10A)==1] <- "250"
  data1$veg[is.na(data1$veg) & (data1$A2==2)&(data1$F10A)==1] <- '250'
  data1$veg <- ifelse(data1$F10A==2,'0',data1$veg)
  data1$veg <- as.numeric(data1$veg)
  data1$veg[is.na(data1$veg)&(data1$A2==1)] <- '250'
  data1$veg[is.na(data1$veg)&(data1$A2==2)] <- '250'

  data1$MET小时 <- as.numeric(data1$MET小时)
  data1$MET <- as.numeric(data1$MET)
  data1$AGE <- as.numeric(data1$AGE)
  data1$A2 <- as.character(data1$A2)
  data1$BMI <- as.numeric(data1$BMI)
  data1$education <- as.character(data1$education)
  data1$smoke <- as.character(data1$smoke)
  data1$packyear1 <- as.numeric(data1$packyear1)
  data1$drink <- as.character(data1$drink)
  data1$DMTIME <- as.numeric(data1$DMTIME)
  data1$H1H2 <- as.character(data1$H1H2)
  data1$redmeat <- as.numeric(data1$redmeat)
  data1$fruit <- as.numeric(data1$fruit)
  data1$veg <- as.numeric(data1$veg)
  # 1.2 Metabolomics preprocessing and existing matches.
  # Preprocess all 288 participants; exclude matched pairs downstream.
  if (!exists("metabolomics_file", inherits = FALSE)) {
    metabolomics_file <- file.path(.integration_data, "metabolomics.xlsx")
  }
  metarowdata <- openxlsx::read.xlsx(metabolomics_file, sheet="1.凯莱谱-检测原始峰值")
  metarowdata1 <- metarowdata[,c("COMPOUND.Name","SUPER.META.PATHWAY")]
  metarowdata2 <- metarowdata[,c("COMPOUND.Name","SUPER.META.PATHWAY","SUB.META.PATHWAY","HMDB.ID","KEGG.ID")]
  metadata <- openxlsx::read.xlsx(metabolomics_file,sheet="4.合并人群数据库（用于分析）")
  metarowdata_0.2 <- openxlsx::read.xlsx(metabolomics_file,sheet="2.数据导入表")
  stopifnot(!anyDuplicated(data1$ID),!anyDuplicated(metadata$ID))

  # Retain metabolites with missingness <20%.
  metabolite_names <- setdiff(names(metarowdata_0.2),c("Sample","label"))
  missing_ratio <- colMeans(is.na(metarowdata_0.2[,metabolite_names,drop=FALSE]))
  selected_columns <- names(missing_ratio)[missing_ratio < .20]
  metarowdata_0.2_clean <- metarowdata_0.2[,selected_columns,drop=FALSE]

  # Impute with half the minimum observed peak value, then apply natural log.
  metarowdata_0.2_imputed <- metarowdata_0.2_clean
  for (v in selected_columns) {
    minimum <- min(metarowdata_0.2_imputed[[v]],na.rm=TRUE)
    stopifnot(is.finite(minimum),minimum > 0)
    metarowdata_0.2_imputed[[v]][is.na(metarowdata_0.2_imputed[[v]])] <- minimum/2
  }
  metarowdata_0.2_imputed <- log(metarowdata_0.2_imputed)

  # Median-center and divide by IQR; retain column order.
  median_values <- apply(metarowdata_0.2_imputed,2,median)
  iqr_values <- apply(metarowdata_0.2_imputed,2,IQR)
  stopifnot(all(iqr_values > 0))
  scaled_metarowdata <- sweep(sweep(metarowdata_0.2_imputed,2,median_values,"-"),2,iqr_values,"/")

  # Join existing Excel match and group keys.
  # Label sign: positive = case, negative = control; absolute value = Match.
  source_match <- abs(metarowdata_0.2$label)
  source_group <- as.integer(metarowdata_0.2$label > 0)
  source_key <- paste(source_match,source_group,sep="_")
  target_key <- paste(metadata$Match,metadata$Group,sep="_")
  stopifnot(!anyDuplicated(source_key),!anyDuplicated(target_key))
  met_order <- match(target_key,source_key)
  cohort_order <- match(metadata$ID,data1$ID)
  stopifnot(!anyNA(met_order),!anyNA(cohort_order))
  sample_number <- as.integer(sub("sample-","",metarowdata_0.2$Sample))
  stopifnot(all(sample_number[met_order] == seq_len(nrow(metadata))))
  scaled_metarowdata1 <- data.frame(Match=metadata$Match,Sample=sample_number[met_order],
    label=metarowdata_0.2$label[met_order],scaled_metarowdata[met_order,,drop=FALSE],check.names=FALSE)
  covariate_names <- c("AGE","A2","BMI","education","smoke","packyear1","drink","DMTIME",
                      "H1H2","redmeat","fruit","veg","colo_cancerfb","endtime5","MET小时","MET")
  data2 <- data1[,c("ID",covariate_names)]
  metarowdata_0.2_imputed1 <- data.frame(metadata[,c("ID","Group")],
    data1[cohort_order,covariate_names,drop=FALSE],scaled_metarowdata1,check.names=FALSE)
  rownames(metarowdata_0.2_imputed1) <- NULL
  stopifnot(nrow(metarowdata_0.2_imputed1)==288L,ncol(scaled_metarowdata)==799L,
            all(table(metadata$Match,metadata$Group)==1L))
  exports <- c("data","data1","data2","metarowdata","metarowdata1","metarowdata2","metadata",
    "metarowdata_0.2","metarowdata_0.2_clean","metarowdata_0.2_imputed","metarowdata_0.2_imputed1",
    "scaled_metarowdata","scaled_metarowdata1","median_values","iqr_values","selected_columns",
    "cohort_file","metabolomics_file")
  list2env(mget(exports),envir=parent.env(environment()))
})
})

  .prep_objects <- c("data", "data1", "data2", "metarowdata", "metarowdata1", "metarowdata2", "metadata", "metarowdata_0.2", "metarowdata_0.2_clean", "metarowdata_0.2_imputed", "metarowdata_0.2_imputed1", "scaled_metarowdata", "scaled_metarowdata1", "median_values", "iqr_values", "selected_columns")
  saveRDS(list(signature = .prep_signature, objects = mget(.prep_objects)), .prep_file)
}
