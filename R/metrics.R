## Metrics usable in cross-validation

## Utility to calculate area under a performance curve
auc <- function(perf){
  x <- perf@x.values[[1]]
  y <- perf@y.values[[1]]
  auc <- 0.0
  for(j in 2:length(x)){
    i <- j-1
    w <- x[j] - x[i]
    if(w > 0 && !is.nan(y[i])){
      h <- y[j] - y[i]
      auc <- auc + y[i]*w + (0.5 * h * w)
    }
  }
  auc
}

## Compute aucs for both ROC and PR
rocpr.aucs <- function(p,t,...){
  library(ROCR)
  pred <- prediction(p$values,t)
  roc <- performance(pred,"tpr","fpr")
  pr  <- performance(pred,"prec","rec")
  c(roc=auc(roc),pr=auc(pr))
}


## Computes the rank of the single positive class
ranks <- function(preds, class, pids){
  prots <- unique(pids)
  rnks <- sapply(prots,function(i){
    idxs <- which(pids == i)
    cls <- class[idxs]
    prd <- preds[idxs]
    nat <- which(cls == max(cls))
    ord <- order(prd, decreasing=TRUE)
    which(ord == nat)
  })
  rnks
}
## Computes the z-score: separation between decoys and native
zscores <- function(preds, class, pids){
  prots <- unique(pids)
  zs <- sapply(prots,function(i){
    idxs <- which(pids == i)
    cls <- class[idxs]
    prd <- preds[idxs]
    nat <- which(cls == max(cls))
    e.decoys <- prd[-nat]
    e.native <- prd[nat]
    z <- (mean(e.decoys) - e.native) / sd(e.decoys)
  })
  zs
}

## Produces a function that can be used in cross validation to
## evaluate various statistics assocaited with decoys
decoy.metric <- function(pids){
   function(preds, class, test,...){
    rp <- rocpr.aucs(preds,class)
    rnks <- ranks(preds$values,class,pids[test])
    c(rp,
      rank=mean(rnks),
      tops=sum(rnks==1)/length(rnks),
      zscore=mean(zscores(preds$values,class,pids[test])))
  }
}
