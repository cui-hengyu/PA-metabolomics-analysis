# Nested sample, baseline table, and metabolite associations.

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

# 1. Matched samples and metabolite associations.
.integration_run_stage("nested_samples_and_spearman.R", .integration_nested, {
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(tableone); library(ggplot2)})
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))
output_root <- .integration_nested
write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
# Prepare matched samples.
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

# Signature derivation sample.
derivation_dat <- droplevels(dat)
early_case <- derivation_dat$Group == 1 & as.numeric(derivation_dat$endtime5) <= 1
removed_matches <- unique(as.character(derivation_dat$Match[early_case]))
# CRC analysis sample.
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
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
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
  stopifnot(nrow(out) == 794L)
  write_csv(out, file.path(out_dir,"01_Spearman_all_794_metabolites_full_nested.csv"))
  selected_name <- if (criterion == "FDR") "02_Spearman_FDR_lt_0.05_full_nested.csv" else "02_Spearman_nominal_P_lt_0.05_full_nested.csv"
  write_csv(selected, file.path(out_dir,selected_name))
  list(all=out, selected=selected)
}


total_screen <- screen_spearman("total_PA_MET_h_day",
  file.path(output_root, "1_total_PA_signature"), "FDR")
leisure_screen <- screen_spearman("nonoccupational_PA_MET_h_day",
  .integration_leisure, "nominal_P")
stopifnot(nrow(total_screen$all) == 794L, nrow(leisure_screen$all) == 794L,
  nrow(total_screen$selected) > 0L, nrow(leisure_screen$selected) > 0L)
write_csv(data.frame(Exposure = c("Total PA", "Non-occupational PA (original leisure label)"),
  N_tested = 794L, Criterion = c("BH-FDR <0.05", "Nominal P <0.05"),
  N_selected = c(nrow(total_screen$selected), nrow(leisure_screen$selected))),
  file.path(output_root, "02_screening_counts_794.csv"))
saveRDS(list(derivation_dat = derivation_dat, lag1_dat = lag1_dat, flow = flow),
        file.path(output_root, "05_complete_analysis_objects.rds"))

# 2. Nested baseline characteristics.
dirs <- c(baseline_pa = file.path(.integration_analysis, "01_baseline_and_PA"))
dir.create(dirs["baseline_pa"], recursive = TRUE, showWarnings = FALSE)
pa_exposures <- c(total = "total_PA_MET_h_day",
                  occupational = "occupational_PA_MET_h_day",
                  nonoccupational = "nonoccupational_PA_MET_h_day")
# Baseline table.
baseline_vars <- c(
  "A2", "AGE", "BMI", "education", "smoking5", "drink", "DMTIME", "H1H2",
  "redmeat", "fruit", "veg", unname(pa_exposures)
)
baseline_factor_vars <- c("A2", "education", "smoking5", "drink", "H1H2")
# Continuous variables: mean (SD).
baseline_continuous_vars <- c(
  "AGE", "BMI", "DMTIME", "redmeat", "fruit", "veg",
  "total_PA_MET_h_day", "occupational_PA_MET_h_day",
  "nonoccupational_PA_MET_h_day"
)
stopifnot(setequal(baseline_continuous_vars,
                   setdiff(baseline_vars, baseline_factor_vars)))
# Smoking categories.
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
# Wilcoxon and chi-square tests.
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
# Include empty smoking categories.
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


})

# 3. Manhattan plots with inset tables.
.integration_run_stage("manhattan_plots.R",
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
  # One font size for every header, class label, and count. A taller header
  # accommodates the three-line final heading without shrinking its text.
  table_text_size <- 2.85
  table_header_rows <- 1.65
  table_row_h <- (table_top-table_bottom)/(nrow(table_body)+table_header_rows)
  table_header_h <- table_header_rows*table_row_h
  table_col_edges <- c(table_left,700,746,table_right)
  table_hlines <- data.frame(x=table_left,xend=table_right,
    y=c(table_top, table_top-table_header_h-(0:nrow(table_body))*table_row_h),
    yend=c(table_top, table_top-table_header_h-(0:nrow(table_body))*table_row_h))
  table_vlines <- data.frame(x=table_col_edges,xend=table_col_edges,
                             y=table_bottom,yend=table_top)
  table_body$y <- table_top-table_header_h-table_row_h*(seq_len(nrow(table_body))-0.5)

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
    annotate("text",x=table_left+10,y=table_top-table_header_h/2,label="Super pathway",
             hjust=0,family="Arial",size=table_text_size,color="black")+
    annotate("text",x=mean(table_col_edges[2:3]),y=table_top-table_header_h/2,
             label="No. of\nmetabolites",lineheight=.86,family="Arial",size=table_text_size,color="black")+
    annotate("text",x=mean(table_col_edges[3:4]),y=table_top-table_header_h/2,
             label="No. of\nsignificant\nmetabolites",lineheight=.86,family="Arial",size=table_text_size,color="black")+
    geom_point(data=table_body[table_body$SUPER.META.PATHWAY!="Total",],
               aes(x=table_left+6,y=y,color=SUPER.META.PATHWAY),inherit.aes=FALSE,
               size=3.2,show.legend=FALSE)+
    geom_text(data=table_body,aes(x=table_left+11,y=y,label=display_name),inherit.aes=FALSE,
              hjust=0,family="Arial",size=table_text_size,color="black")+
    geom_text(data=table_body,aes(x=mean(table_col_edges[2:3]),y=y,label=count_all),
              inherit.aes=FALSE,family="Arial",size=table_text_size,color="black")+
    geom_text(data=table_body,aes(x=mean(table_col_edges[3:4]),y=y,label=count_significant),
              inherit.aes=FALSE,family="Arial",size=table_text_size,color="black")+
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
  directory = c(file.path(.integration_nested, "1_total_PA_signature"), .integration_leisure),
  stem = c("03_Spearman_FDR_Manhattan", "03_Spearman_nominal_P_Manhattan"),
  seed = c(20260921L, 20260923L),
  criterion = c("FDR", "nominal_P"), stringsAsFactors = FALSE
)
for (i in seq_len(nrow(jobs))) {
  directory <- jobs$directory[i]
  draw_one(
    file.path(directory, "01_Spearman_all_794_metabolites_full_nested.csv"),
    file.path(directory, paste0(jobs$stem[i], ".pdf")),
    file.path(directory, paste0(jobs$stem[i], ".png")),
    file.path(directory, sub("^03_", "03a_", paste0(jobs$stem[i], "_plot_data.csv"))),
    file.path(directory, sub("^03_", "03b_", paste0(jobs$stem[i], "_pathway_summary.csv"))),
    jobs$seed[i], jobs$criterion[i]
  )
}
})
