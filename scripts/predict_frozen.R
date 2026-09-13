# Rscript scripts/predict_frozen.R model.rds beta.tsv output.tsv
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=3)stop("Expected model, beta matrix, output")
source("R/model.R");if(!requireNamespace("data.table",quietly=TRUE))stop("Install data.table")
b<-readRDS(args[1]);d<-data.table::fread(args[2],data.table=FALSE,check.names=FALSE)
if(anyDuplicated(d[[1]]))stop("Duplicate probes")
x<-t(as.matrix(d[,-1,drop=FALSE]));storage.mode(x)<-"double";colnames(x)<-d[[1]]
if(any(is.finite(x)&(x<0|x>1)))stop("Invalid beta")
p<-predict_en(b,x);q<-b$interval_q
if(is.null(q))q<-Inf
p$lower<-p$predicted_reference_HRDsum-q;p$upper<-p$predicted_reference_HRDsum+q
p$interval_status<-if(is.finite(q))"empirical_no_domain_shift_guarantee" else "unavailable_insufficient_calibration"
p$warning<-ifelse(p$qc_fail,"excessive_missing_features",ifelse(p$ood,"outside_training_distance_reference","research_prediction_unvalidated_domain"))
p$estimate_for_display<-ifelse(p$reportable,p$predicted_reference_HRDsum,NA_real_)
p$lower[!p$reportable|!is.finite(q)]<-NA_real_;p$upper[!p$reportable|!is.finite(q)]<-NA_real_
p$HRD_high_probability<-NA_real_
p$model_md5<-unname(tools::md5sum(args[1]))
p$provenance<-"frozen_research_model;domain_validity_requires_review"
dir.create(dirname(args[3]),recursive=TRUE,showWarnings=FALSE)
write.table(p,args[3],sep="\t",row.names=FALSE,quote=FALSE)
