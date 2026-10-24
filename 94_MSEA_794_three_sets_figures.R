# Updated manuscript Figures S6 and S7, following the supplied Supplementary
# Materials layout. Compute from the current 794-feature saved samples, without
# altering the original 799-feature scripts, data or results.
options(stringsAsFactors = FALSE)
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(arg)==1L)
run_root <- dirname(normalizePath(sub("^--file=","",arg),winslash="/"))
input_root <- dirname(run_root)
.libPaths(c(file.path(input_root,"_library"),.libPaths()))
suppressPackageStartupMessages({
  library(survival); library(clusterProfiler); library(enrichplot)
  library(data.table); library(ggplot2); library(patchwork); library(BiocParallel)
})
BiocParallel::register(BiocParallel::SerialParam())
if(.Platform$OS.type=="windows") windowsFonts(Arial=windowsFont("Arial"))
if(!"Arial"%in%names(pdfFonts())) pdfFonts(Arial=pdfFonts("Helvetica")[[1L]])
options(device=function(...) grDevices::pdf(file=NULL,...))
nested <- file.path(run_root,"results/analysis/nested")
out <- file.path(nested,"total_PA/MSEA_794_three_sets")
fig <- file.path(run_root,"figures")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(fig,showWarnings=FALSE)
write_csv <- function(d,name) write.csv(d,file.path(out,name),row.names=FALSE,fileEncoding="UTF-8")
if("--plots-only"%in%commandArgs(TRUE)) {
  cache <- readRDS(file.path(out,"07_MSEA_794_full_precision_objects.rds"))
  list2env(cache,envir=.GlobalEnv)
} else {
obj <- readRDS(file.path(nested,"05_complete_analysis_objects.rds"))
full_dat <- obj$derivation_dat
lag1_dat <- obj$lag1_dat
ann <- read.csv(file.path(nested,"1_total_PA_signature/01_Spearman_all_794_metabolites_full_nested.csv"),
                check.names=FALSE,fileEncoding="UTF-8")
metabolites <- ann$COMPOUND.Name
stopifnot(nrow(ann)==794L,!anyDuplicated(metabolites),nrow(full_dat)==288L,nrow(lag1_dat)==256L,
          length(unique(full_dat$Match))==144L,length(unique(lag1_dat$Match))==128L,
          all(metabolites%in%names(full_dat)),all(metabolites%in%names(lag1_dat)))
removed <- c("3-amino-2-piperidone","4-hydroxyphenylacetylglutamine","N,N-dimethyl-pro-pro",
             "O-acetylhomoserine","alpha-tocopherol")
stopifnot(!any(metabolites%in%removed))

# Keep full-precision double statistics in memory/RDS through GSEA. Do not
# reconstruct ranks from rounded CSVs, which can change tail P calculations.
pa_all <- data.frame(COMPOUND.Name=metabolites,HMDB.ID=ann$HMDB.ID,
  Spearman_rho=vapply(metabolites,function(m)suppressWarnings(cor(full_dat$total_PA_MET_h_day,
    full_dat[[m]],method="spearman",use="complete.obs")),numeric(1)))
covars <- c("AGE","BMI","education","smoking5","drink","DMTIME","H1H2","redmeat","fruit","veg")
crc_all <- data.frame(COMPOUND.Name=metabolites,HMDB.ID=ann$HMDB.ID,
                      logOR_per_SD=NA_real_,OR_per_SD=NA_real_,CI_lower=NA_real_,CI_upper=NA_real_,P=NA_real_)
for(i in seq_along(metabolites)) {
  f <- as.formula(paste0("Group ~ scale(`",metabolites[i],"`) + ",paste(covars,collapse=" + ")," + strata(Match)"))
  model <- suppressWarnings(clogit(f,data=lag1_dat,method="exact",na.action=na.fail))
  z <- summary(model)$coefficients[1,]
  crc_all$logOR_per_SD[i] <- z["coef"]
  crc_all$OR_per_SD[i] <- exp(z["coef"])
  crc_all$CI_lower[i] <- exp(z["coef"]-1.96*z["se(coef)"])
  crc_all$CI_upper[i] <- exp(z["coef"]+1.96*z["se(coef)"])
  crc_all$P[i] <- z["Pr(>|z|)"]
  if(i%%100L==0L) {cat("CRC rankings:",i,"/",length(metabolites),"\n");flush.console()}
}
crc_all$FDR_P <- p.adjust(crc_all$P,method="BH")
make_unique_rank <- function(x,stat) {
  z <- x[!is.na(x$HMDB.ID)&nzchar(x$HMDB.ID)&is.finite(x[[stat]]),
         c("COMPOUND.Name","HMDB.ID",stat),drop=FALSE]
  z <- z[z$COMPOUND.Name!="N-delta-acetylornithine",,drop=FALSE]
  z <- z[!z$HMDB.ID%in%c("HMDB0240429","HMDB0010383"),,drop=FALSE]
  stopifnot(!anyDuplicated(z$HMDB.ID))
  names(z)[3] <- "Rank_statistic"
  z[order(z$Rank_statistic,decreasing=TRUE),,drop=FALSE]
}
pa_rank <- make_unique_rank(pa_all,"Spearman_rho")
crc_rank <- make_unique_rank(crc_all,"logOR_per_SD")
stopifnot(nrow(pa_rank)==638L,nrow(crc_rank)==638L)
targets <- data.frame(Panel=c("A","B","C"),SMPDB_ID=c("SMP0000456","SMP0000035","SMP0000071"),
  Pathway=c("Fatty Acid Biosynthesis","Bile Acid Biosynthesis","Ketone Body Metabolism"))
members <- do.call(rbind,lapply(seq_len(nrow(targets)),function(i) {
  p <- file.path(input_root,"data/pathways",paste0(targets$SMPDB_ID[i],"_metabolites.csv"))
  x <- as.data.frame(fread(p,encoding="UTF-8",fill=TRUE,quote="\"",showProgress=FALSE),check.names=FALSE)
  fields <- c("SMPDB ID","Pathway Name","Metabolite Name","HMDB ID")
  stopifnot(all(fields%in%names(x)))
  x <- x[!is.na(x[["HMDB ID"]])&nzchar(x[["HMDB ID"]]),fields,drop=FALSE]
  x <- x[!duplicated(x[["HMDB ID"]]),,drop=FALSE]
  names(x) <- c("SMPDB_ID","Pathway","SMPDB_metabolite","HMDB_ID")
  x$SMPDB_ID <- targets$SMPDB_ID[i];x$Pathway <- targets$Pathway[i]
  x
}))
term2gene <- unique(members[,c("SMPDB_ID","HMDB_ID")]);names(term2gene)<-c("term","gene")
term2name <- targets[,c("SMPDB_ID","Pathway")];names(term2name)<-c("term","name")
run_msea <- function(ranks) {
  gl <- ranks$Rank_statistic;names(gl)<-ranks$HMDB.ID
  gl <- sort(gl,decreasing=TRUE)
  set.seed(20260922)
  GSEA(geneList=gl,TERM2GENE=term2gene,TERM2NAME=term2name,minGSSize=2,maxGSSize=500,
    pvalueCutoff=1,pAdjustMethod="none",eps=0,verbose=FALSE,seed=TRUE,by="fgsea")
}
pa_gsea <- run_msea(pa_rank)
crc_gsea <- run_msea(crc_rank)
tidy <- function(g,label) {
  x<-as.data.frame(g);x<-x[match(targets$SMPDB_ID,x$ID),,drop=FALSE]
  stopifnot(identical(as.integer(x$setSize),c(8L,12L,2L)))
  data.frame(Analysis=label,Panel=targets$Panel,SMPDB_ID=targets$SMPDB_ID,Metabolite_set=targets$Pathway,
    Set_size=x$setSize,ES=x$enrichmentScore,NES=x$NES,P=x$pvalue,
    P_type="Nominal; pAdjustMethod=none",Rank_at_max=x$rank,Core_enrichment=x$core_enrichment)
}
pa_results <- tidy(pa_gsea,"Total PA: Spearman rho; full nested N=288")
crc_results <- tidy(crc_gsea,"CRC: fully adjusted log OR per SD; lag-1 N=256")
write_csv(pa_all,"01_all_794_total_PA_Spearman.csv")
write_csv(crc_all,"02_all_794_CRC_full_model.csv")
write_csv(pa_rank,"03_total_PA_638_unique_HMDB_ranks.csv")
write_csv(crc_rank,"04_CRC_638_unique_HMDB_ranks.csv")
write_csv(members,"05_SMPDB_set_members.csv")
write_csv(rbind(pa_results,crc_results),"06_MSEA_794_three_sets_results.csv")
saveRDS(list(targets=targets,pathway_members=members,pa_rank=pa_rank,crc_rank=crc_rank,
  pa_gsea=pa_gsea,crc_gsea=crc_gsea,pa_results=pa_results,crc_results=crc_results),
  file.path(out,"07_MSEA_794_full_precision_objects.rds"))
}

make_panel <- function(g,res,i) {
  pval <- res$P[i]
  # Literal text avoids Symbol-font mathematical '=' glyphs and preserves
  # trailing zeros. Only the P is italic, using Arial Italic.
  ux <- grid::unit(.90,"npc")
  uy <- grid::unit(1,"npc")-grid::unit(5.5,"mm")
  p_suffix <- grid::textGrob(if(pval<.001)" < 0.001" else sprintf(" = %.3f",pval),
    x=ux,y=uy,just=c("right","top"),gp=grid::gpar(fontfamily="Arial",fontsize=10))
  annotation <- grid::grobTree(
    grid::textGrob(sprintf("NES = %.2f",res$NES[i]),
      x=ux,y=grid::unit(1,"npc")-grid::unit(1,"mm"),just=c("right","top"),
      gp=grid::gpar(fontfamily="Arial",fontsize=10)),
    grid::textGrob("P",x=ux-grid::grobWidth(p_suffix),y=uy,just=c("right","top"),
      gp=grid::gpar(fontfamily="Arial",fontface="italic",fontsize=10)),p_suffix)
  pieces <- gseaplot2(g,geneSetID=targets$SMPDB_ID[i],title=targets$Pathway[i],
    color="green",base_size=12,pvalue_table=FALSE)
  pieces[[1]] <- pieces[[1]] +
    labs(title=targets$Pathway[i],tag=paste0("(",targets$Panel[i],")")) +
    annotation_custom(annotation,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf) +
    theme(text=element_text(family="Arial",colour="black"),
      plot.title=element_text(family="Arial",face="plain",size=16,hjust=.5),
      plot.title.position="plot",
      plot.tag=element_text(family="Arial",face="plain",size=16,hjust=0,vjust=1),
      plot.tag.position=c(0,1),plot.margin=margin(7,9,0,9))
  pieces[[2]] <- pieces[[2]] + theme(text=element_text(family="Arial",colour="black"))
  pieces[[3]] <- pieces[[3]] + theme(text=element_text(family="Arial",colour="black"))
  wrap_plots(pieces,ncol=1,heights=c(1.45,.45,1))
}
for(domain in c("PA","CRC")) {
  g <- if(domain=="PA")pa_gsea else crc_gsea
  res <- if(domain=="PA")pa_results else crc_results
  panels <- lapply(seq_len(nrow(targets)),function(i)make_panel(g,res,i))
  plot <- wrap_plots(panels,nrow=1)
  stem <- if(domain=="PA")"Figure_S6_MSEA_total_PA" else "Figure_S7_MSEA_CRC"
  ggsave(file.path(fig,paste0(stem,".pdf")),plot,width=18,height=6.4,units="in",device=cairo_pdf,bg="white")
  ggsave(file.path(fig,paste0(stem,"_600dpi.png")),plot,width=18,height=6.4,units="in",dpi=600,bg="white")
}
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
cat("Updated MSEA results (nominal P):\n");print(rbind(pa_results,crc_results),row.names=FALSE)
cat("Two manuscript figures exported to",fig,"\n")
