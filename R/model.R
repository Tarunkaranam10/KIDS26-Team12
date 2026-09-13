# Fold-safe elastic-net helpers. All learned preprocessing belongs to a fit object.
fit_preprocess <- function(x, max_features=5000L, max_missing=0.05) {
  stopifnot(is.matrix(x), !is.null(colnames(x)), !anyDuplicated(colnames(x)))
  keep <- colMeans(!is.finite(x)) <= max_missing
  if (!any(keep)) stop("No features pass train-only missingness")
  z <- x[,keep,drop=FALSE]
  med <- apply(z,2,median,na.rm=TRUE)
  for (j in seq_len(ncol(z))) z[!is.finite(z[,j]),j] <- med[j]
  v <- apply(z,2,var); ii <- which(is.finite(v) & v>0)
  ii <- ii[order(-v[ii],names(v)[ii])]
  ii <- head(ii,max_features)
  if (length(ii)<2) stop("Fewer than two nonconstant features")
  z <- z[,ii,drop=FALSE]; mu <- colMeans(z); s <- apply(z,2,sd)
  list(features=colnames(z), median=med[colnames(z)],center=mu,scale=s,
       max_missing=max_missing)
}
apply_preprocess <- function(x,pp, reject_missing=TRUE) {
  if (anyDuplicated(colnames(x))) stop("Duplicate input features")
  z <- matrix(NA_real_,nrow(x),length(pp$features),dimnames=list(rownames(x),pp$features))
  common <- intersect(colnames(x),pp$features);z[,common] <- x[,common,drop=FALSE]
  miss <- rowMeans(!is.finite(z))
  if (reject_missing && any(miss>pp$max_missing)) stop("Too many missing model features")
  for (j in seq_len(ncol(z))) z[!is.finite(z[,j]),j] <- pp$median[j]
  z <- sweep(sweep(z,2,pp$center,"-"),2,pp$scale,"/")
  attr(z,"missing_fraction") <- miss;z
}
inner_folds <- function(patient,cancer,seed=260910L) {
  if (anyDuplicated(patient)) stop("Use one specimen per patient for this MVP")
  if (length(unique(cancer))>=3) return(match(cancer,sort(unique(cancer))))
  if (length(patient)<9) stop("Need >=9 training patients for patient-fold fallback")
  set.seed(seed); as.integer(sample(rep(1:3,length.out=length(patient))))
}
fit_en <- function(x,y,patient,cancer,max_features=5000L,seed=260910L) {
  if (!requireNamespace("glmnet",quietly=TRUE)) stop("Install glmnet via scripts/setup.R")
  stopifnot(length(y)==nrow(x),all(is.finite(y)),length(patient)==length(y))
  if (sd(y)==0) stop("Constant training target")
  folds <- inner_folds(patient,cancer,seed)
  grid <- expand.grid(alpha=c(0.1,0.5,1),lambda=c(0.01,0.1,1,10,100))
  losses <- matrix(NA_real_,nrow(grid),length(unique(folds)))
  for (f in sort(unique(folds))) {
    tr <- folds!=f;va <- !tr
    pp <- fit_preprocess(x[tr,,drop=FALSE],max_features)
    ztr <- apply_preprocess(x[tr,,drop=FALSE],pp,FALSE)
    zva <- apply_preprocess(x[va,,drop=FALSE],pp,FALSE)
    for (a in unique(grid$alpha)) {
      mod <- glmnet::glmnet(ztr,y[tr],alpha=a,lambda=sort(unique(grid$lambda),decreasing=TRUE),standardize=FALSE)
      ids <- which(grid$alpha==a)
      for (g in ids) {
        pred <- as.numeric(predict(mod,zva,s=grid$lambda[g]))
        # Equal weight to each inner-validation cancer, independent of size.
        losses[g,match(f,sort(unique(folds)))] <- mean(tapply(abs(pred-y[va]),cancer[va],mean))
      }
    }
  }
  grid$inner_macro_mae <- rowMeans(losses)
  best <- which.min(grid$inner_macro_mae); pp <- fit_preprocess(x,max_features)
  z <- apply_preprocess(x,pp,FALSE)
  mod <- glmnet::glmnet(z,y,alpha=grid$alpha[best],lambda=sort(unique(grid$lambda),decreasing=TRUE),standardize=FALSE)
  # Heuristic score only: not a calibrated domain-shift test.
  dist <- sqrt(rowMeans(z^2));ood_cut <- as.numeric(quantile(dist,0.99,names=FALSE))
  list(preprocess=pp,model=mod,alpha=grid$alpha[best],lambda=grid$lambda[best],
       tuning=grid,inner_folds=data.frame(patient_id=patient,cancer_type=cancer,fold=folds),
       ood_cut=ood_cut,seed=seed,training_n=length(y),training_mean=mean(y),target="reference_HRDsum")
}
predict_en <- function(bundle,x) {
  z <- apply_preprocess(x,bundle$preprocess,FALSE)
  raw <- as.numeric(predict(bundle$model,z,s=bundle$lambda))
  score <- sqrt(rowMeans(z^2));missing <- attr(z,"missing_fraction")
  fail <- missing>bundle$preprocess$max_missing
  ood <- score>bundle$ood_cut
  data.frame(sample_id=rownames(x),predicted_reference_HRDsum=raw,
             missing_fraction=missing,ood_score=score,ood=ood,qc_fail=fail,
             reportable=!(fail|ood),stringsAsFactors=FALSE)
}
conformal_q <- function(y,pred,coverage=0.95) {
  good <- is.finite(y)&is.finite(pred);res <- sort(abs(y[good]-pred[good]));n <- length(res)
  k <- ceiling((n+1)*coverage)
  if (!n || k>n) return(Inf)
  res[k]
}
metrics <- function(y,p) {
  ok <- is.finite(y)&is.finite(p);y<-y[ok];p<-p[ok];n<-length(y)
  if (!n) return(data.frame(n=0,MAE=NA,RMSE=NA,R2=NA,Pearson=NA,Spearman=NA,calibration_intercept=NA,calibration_slope=NA,bias=NA))
  variable <- n>2 && sd(y)>0 && sd(p)>0
  cal <- if (variable) coef(lm(y~p)) else c(NA,NA)
  data.frame(n=n,MAE=mean(abs(y-p)),RMSE=sqrt(mean((y-p)^2)),
             R2=if(n>1 && sum((y-mean(y))^2)>0) 1-sum((y-p)^2)/sum((y-mean(y))^2) else NA,
             Pearson=if(variable) cor(y,p) else NA,Spearman=if(variable) cor(y,p,method="spearman") else NA,
             calibration_intercept=unname(cal[1]),calibration_slope=unname(cal[2]),bias=mean(p-y))
}
