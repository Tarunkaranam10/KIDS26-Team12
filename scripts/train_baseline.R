# Run from repository root: Rscript scripts/train_baseline.R beta.tsv master.tsv results/run01
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop("Usage: Rscript scripts/train_baseline.R beta.tsv master.tsv output_dir")
source("R/model.R")
if (!requireNamespace("data.table",quietly=TRUE)) stop("Install dependencies with scripts/setup.R")
beta <- data.table::fread(args[1],data.table=FALSE,check.names=FALSE)
if (anyDuplicated(beta[[1]])) stop("Duplicate probe IDs")
x <- t(as.matrix(beta[,-1,drop=FALSE])); storage.mode(x)<-"double";colnames(x)<-beta[[1]]
if(any(is.finite(x)&(x<0|x>1))) stop("Expected beta values in [0,1]")
meta <- read.delim(args[2],check.names=FALSE,stringsAsFactors=FALSE)
if(anyDuplicated(meta$sample_id)||anyDuplicated(meta$patient_id)) stop("Duplicate patient/specimen")
if(any(!rownames(x)%in%meta$sample_id)) stop("Missing sample metadata")
meta <- meta[match(rownames(x),meta$sample_id),]
if(any(!is.finite(meta$HRDsum))||any(meta$HRDsum<0)) stop("Invalid target")
if(!all(abs(meta$HRDsum-meta$HRD_LOH-meta$LST-meta$TAI)<1e-6)) stop("Label sum mismatch")
if (!all(meta$quality_annotation=="published_450K_no_exclusion")) stop("Unadjudicated QC")
# Never silently run a smoke-probe fixture as scientific evidence.
prov <- sub("\\.[^.]+$",".provenance.json",args[1])
if(!file.exists(prov)) stop("Missing matrix provenance JSON")
if(!requireNamespace("jsonlite",quietly=TRUE)) stop("Install jsonlite")
pr <- jsonlite::fromJSON(prov)
if(isTRUE(pr$engineering_only)) stop("Engineering-only matrix. Build a real technical probe allowlist first.")
if(is.null(pr$probe_allowlist_sha256)) stop("Missing technical allowlist hash")
out <- args[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
cns <- meta$cancer_type %in% c("GBM","LGG")
if(length(unique(meta$cancer_type[!cns]))<3) stop("Need >=3 development cancers; six-patient smoke subset is insufficient")
if(sum(!cns)<50) stop("Need >=50 development patients for this protocol; do not claim transfer from a tiny fixture")
all_predictions <- list();all_metrics<-list()
for (type in sort(unique(meta$cancer_type[!cns]))) {
 tr <- !cns & meta$cancer_type!=type;te <- !cns & meta$cancer_type==type
 b <- fit_en(x[tr,,drop=FALSE],meta$HRDsum[tr],meta$patient_id[tr],meta$cancer_type[tr])
 p <- predict_en(b,x[te,,drop=FALSE]);p$actual<-meta$HRDsum[te];p$cancer_type<-type;p$patient_id<-meta$patient_id[te]
 p$null_prediction<-b$training_mean;p$split<-"development_LOCO"
 all_predictions[[type]]<-p
 mm<-metrics(p$actual,p$predicted_reference_HRDsum);mm$cancer_type<-type;mm$train_n<-sum(tr);mm$selected_features<-length(b$preprocess$features);mm$alpha<-b$alpha;mm$lambda<-b$lambda;mm$null_MAE<-mean(abs(p$actual-b$training_mean));mm$abstention_rate<-mean(!p$reportable)
 all_metrics[[type]]<-mm
 saveRDS(b,file.path(out,paste0("loco_",type,".rds")))
 write.table(b$inner_folds,file.path(out,paste0("inner_folds_",type,".tsv")),sep="\t",row.names=FALSE,quote=FALSE)
}
write.table(do.call(rbind,all_predictions),file.path(out,"loco_predictions.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
mt<-do.call(rbind,all_metrics);write.table(mt,file.path(out,"loco_metrics.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines(paste("macro_MAE",mean(mt$MAE)),file.path(out,"macro_metrics.txt"))
# Deterministic calibration reservation within non-CNS cancers; no refit on calibration.
set.seed(260910);dev<-which(!cns);cal<-unlist(lapply(split(dev,meta$cancer_type[dev]),function(ii) sample(ii,max(1L,floor(length(ii)*0.2)))))
tr<-setdiff(dev,cal)
b<-fit_en(x[tr,,drop=FALSE],meta$HRDsum[tr],meta$patient_id[tr],meta$cancer_type[tr])
cp<-predict_en(b,x[cal,,drop=FALSE]);b$interval_q<-conformal_q(meta$HRDsum[cal],cp$predicted_reference_HRDsum)
b$calibration_ids<-meta$patient_id[cal];b$training_ids<-meta$patient_id[tr];b$probe_provenance<-pr
b$interval_note<-"Split-conformal residual interval; no pediatric/domain-shift coverage guarantee"
saveRDS(b,file.path(out,"frozen_nonCNS.rds"))
cc<-as.matrix(coef(b$model,s=b$lambda))
write.table(data.frame(feature=rownames(cc),coefficient=as.numeric(cc[,1]),stringsAsFactors=FALSE),file.path(out,"frozen_coefficients.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(data.frame(sample_id=meta$sample_id,patient_id=meta$patient_id,cancer_type=meta$cancer_type,role=ifelse(cns,"locked_CNS",ifelse(seq_len(nrow(meta))%in%cal,"calibration","training"))),file.path(out,"final_partitions.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
# Deliberately do not evaluate CNS here: architecture selection must precede locked testing.
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
cat("Development LOCO complete; frozen_nonCNS.rds saved. Review development and lock before external prediction.\n")
