svm <-
function (x, ...)
    UseMethod ("svm")

predict.svmmodel <- function (model, newdata, probability=FALSE, ...){
    stopifnot(class(model) == "svmmodel",
              !missing(newdata),
              class(newdata) == "svmproblem",
              model$nCols == newdata$nCols,
              model$tot.nSV > 0)

    ret <- list(pred = double(newdata$nRows),
                dec = double(newdata$nRows * model$nclasses * (model$nclasses - 1) / 2),
                prob = double(newdata$nRows * model$nclasses))
    
    .External("svmpredict",
              train = newdata$x,
              nRows = newdata$nRows,
              probability = probability,
              model = model,
              ret = ret)

    ## PACKAGE = "e1071"
    ## )
    ## class(ret) <- "svm.prediction"
    ret
}

print.svmproblem <-  function(x,...) {
    cat("\nSVM Problem Data\n")
    cat("Rows:",x$nRows,"\n")
    cat("nCols:", x$nCols,"\n")
    cat("y: ",if(is.null(x$y)){"Missing"}else{"Present"},"\n\n")
}
                 
setMethod("[", "svmproblem", #signature(x = "svmproblem", i = "index", j = "missing", drop = "logical"),
          function(x, i, j, ..., drop){
              x$x <- x$x[i]
              if(!is.null(x$y)){
                  x$y <- y[i]
              }
              x
          }
)

"[.svmproblem" <- function(x, i, j, ..., drop){
    x$x <- x$x[i]
    if(!is.null(x$y)){
        x$y <- x$y[i]
    }
    x$nRows <- length(i)
    x
}




print.svmmodel <-
function (x, ...)
{
    cat("\nCall:", deparse(x$call, 0.8 * getOption("width")), "\n", sep="\n")
    cat("Parameters:\n")
    cat("   SVM-Type: ", c("C-classification",
                           "nu-classification",
                           "one-classification",
                           "eps-regression",
                           "nu-regression")[x$type+1], "\n")
    cat(" SVM-Kernel: ", c("linear",
                           "polynomial",
                           "radial",
                           "sigmoid")[x$kernel+1], "\n")
    if (x$type==0 || x$type==3 || x$type==4)
        cat("       cost: ", x$cost, "\n")
    if (x$kernel==1)
        cat("     degree: ", x$degree, "\n")
    cat("      gamma: ", x$gamma, "\n")
    if (x$kernel==1 || x$kernel==3)
        cat("     coef.0: ", x$coef0, "\n")
    if (x$type==1 || x$type==2 || x$type==4)
        cat("         nu: ", x$nu, "\n")
    if (x$type==3) {
        cat("    epsilon: ", x$epsilon, "\n\n")
        if (x$compprob)
            cat("Sigma: ", x$sigma, "\n\n")
    }

    cat("\nNumber of Support Vectors: ", x$tot.nSV)
    cat("\n\n")

}

summary.svm <-
function(object, ...)
    structure(object, class="summary.svm")

print.summary.svm <-
function (x, ...)
{
    print.svm(x)
    if (x$type<2) {
        cat(" (", x$nSV, ")\n\n")
        cat("\nNumber of Classes: ", x$nclasses, "\n\n")
        cat("Levels:", if(is.numeric(x$levels)) "(as integer)", "\n", x$levels)
    }
    cat("\n\n")
    if (x$type==2) cat("\nNumber of Classes: 1\n\n\n")

    if ("MSE" %in% names(x)) {
        cat(length (x$MSE), "-fold cross-validation on training data:\n\n", sep="")
        cat("Total Mean Squared Error:", x$tot.MSE, "\n")
        cat("Squared Correlation Coefficient:", x$scorrcoef, "\n")
        cat("Mean Squared Errors:\n", x$MSE, "\n\n")
    }
    if ("accuracies" %in% names(x)) {
        cat(length (x$accuracies), "-fold cross-validation on training data:\n\n", sep="")
        cat("Total Accuracy:", x$tot.accuracy, "\n")
        cat("Single Accuracies:\n", x$accuracies, "\n\n")
    }
    cat("\n\n")
}

scale.data.frame <-
function(x, center = TRUE, scale = TRUE)
{
    i <- sapply(x, is.numeric)
    if (ncol(x[, i, drop = FALSE])) {
        x[, i] <- tmp <- scale.default(x[, i, drop = FALSE], na.omit(center), na.omit(scale))
        if(center || !is.logical(center))
            attr(x, "scaled:center")[i] <- attr(tmp, "scaled:center")
        if(scale || !is.logical(scale))
            attr(x, "scaled:scale")[i]  <- attr(tmp, "scaled:scale")
    }
    x
}

plot.svm <-
function(x, data, formula = NULL, fill = TRUE,
         grid = 50, slice = list(), symbolPalette = palette(),
         svSymbol = "x", dataSymbol = "o", ...)
{
    if (x$type < 3) {
        if (is.null(formula) && ncol(data) == 3) {
            formula <- formula(delete.response(terms(x)))
            formula[2:3] <- formula[[2]][2:3]
        }
        if (is.null(formula))
            stop("missing formula.")
        if (fill) {
            sub <- model.frame(formula, data)
            xr <- seq(min(sub[, 2]), max(sub[, 2]), length = grid)
            yr <- seq(min(sub[, 1]), max(sub[, 1]), length = grid)
            l <- length(slice)
            if (l < ncol(data) - 3) {
                slnames <- names(slice)
                slice <- c(slice, rep(list(0), ncol(data) - 3 -
                                      l))
                names <- labels(delete.response(terms(x)))
                names(slice) <- c(slnames, names[!names %in%
                                                 c(colnames(sub), slnames)])
            }
            for (i in names(which(sapply(data, is.factor))))
                if (!is.factor(slice[[i]])) {
                    levs <- levels(data[[i]])
                    lev <- if (is.character(slice[[i]])) slice[[i]] else levs[1]
                    fac <- factor(lev, levels = levs)
                    if (is.na(fac))
                        stop(paste("Level", dQuote(lev), "could not be found in factor", sQuote(i)))
                    slice[[i]] <- fac
                }

            lis <- c(list(yr), list(xr), slice)
            names(lis)[1:2] <- colnames(sub)
            new <- expand.grid(lis)[, labels(terms(x))]
            preds <- predict(x, new)
            filled.contour(xr, yr,
                           matrix(as.numeric(preds),
                                  nrow = length(xr), byrow = TRUE),
                           plot.axes = {
                               axis(1)
                               axis(2)
                               colind <- as.numeric(model.response(model.frame(x, data)))
                               dat1 <- data[-x$index,]
                               dat2 <- data[x$index,]
                               coltmp1 <- symbolPalette[colind[-x$index]]
                               coltmp2 <- symbolPalette[colind[x$index]]
                               points(formula, data = dat1, pch = dataSymbol, col = coltmp1)
                               points(formula, data = dat2, pch = svSymbol, col = coltmp2)
                           },
                           levels = 1:(length(levels(preds)) + 1),
                           key.axes = axis(4, 1:(length(levels(preds))) + 0.5,
                           labels = levels(preds),
                           las = 3),
                           plot.title = title(main = "SVM classification plot",
                           xlab = names(lis)[2], ylab = names(lis)[1]),
                           ...)
        }
        else {
            plot(formula, data = data, type = "n", ...)
            colind <- as.numeric(model.response(model.frame(x,
                                                            data)))
            dat1 <- data[-x$index,]
            dat2 <- data[x$index,]
            coltmp1 <- symbolPalette[colind[-x$index]]
            coltmp2 <- symbolPalette[colind[x$index]]
            points(formula, data = dat1, pch = dataSymbol, col = coltmp1)
            points(formula, data = dat2, pch = svSymbol, col = coltmp2)
            invisible()
        }
    }
}

write.svm <-
function (object, svm.file="Rdata.svm", scale.file = "Rdata.scale",
          yscale.file = "Rdata.yscale")
{

    ret <- .C ("svmwrite",
               ## model
               as.double  (if (object$sparse) object$SV@ra else t(object$SV)),
               as.integer (nrow(object$SV)), as.integer(ncol(object$SV)),
               as.integer (if (object$sparse) object$SV@ia else 0),
               as.integer (if (object$sparse) object$SV@ja else 0),
               as.double  (as.vector(object$coefs)),
               as.double  (object$rho),
               as.double  (object$probA),
               as.double  (object$probB),
               as.integer (object$nclasses),
               as.integer (object$tot.nSV),
               as.integer (object$labels),
               as.integer (object$nSV),
               as.integer (object$sparse),

               ## parameter
               as.integer (object$type),
               as.integer (object$kernel),
               as.integer (object$degree),
               as.double  (object$gamma),
               as.double  (object$coef0),

               ## filename
               as.character(svm.file),

               PACKAGE = "e1071"
               )$ret

    write.table(data.frame(center = object$x.scale$"scaled:center",
                           scale  = object$x.scale$"scaled:scale"),
                file=scale.file, col.names=FALSE, row.names=FALSE)

    if (!is.null(object$y.scale))
        write.table(data.frame(center = object$y.scale$"scaled:center",
                               scale  = object$y.scale$"scaled:scale"),
                    file=yscale.file, col.names=FALSE, row.names=FALSE)
}

read.svmproblem <-
  function(filename){
    stopifnot(is.character(filename),
              file.exists(filename))
    prob <- .Call("svm_read_problem",
                  filename
                                        #PACKAGE="libsvmR"
                  )
    ## Should have a list of 4 things
    ## length, targetvalues, feature pointers, node space
    names(prob) <- c("nRows","nCols","y","x","x.space")
    class(prob) <- "svmproblem"
    prob
  }



set.svmdebug <- function(d){
  .C("svm_set_debug",as.integer(d))
}

svm.default <-
function (p,
          type        = NULL,
          kernel      = "radial",
          degree      = 3L,
          gamma       = 1 / p$nCols,
          coef0       = 0,
          cost        = 1,
          nu          = 0.5,
          class.weights = NULL,
          cachesize   = 40,
          tolerance   = 0.001,
          epsilon     = 0.1,
          shrinking   = TRUE,
          probability = FALSE,
          fitted      = TRUE,
          seed        = 1L,
          ...)
{
  stopifnot(class(p) == "svmproblem")
  y <- p$y

  ## NULL parameters?
  if(is.null(degree)) stop(sQuote("degree"), " must not be NULL!")
  if(is.null(gamma)) stop(sQuote("gamma"), " must not be NULL!")
  if(is.null(coef0)) stop(sQuote("coef0"), " must not be NULL!")
  if(is.null(cost)) stop(sQuote("cost"), " must not be NULL!")
  if(is.null(nu)) stop(sQuote("nu"), " must not be NULL!")
  if(is.null(epsilon)) stop(sQuote("epsilon"), " must not be NULL!")
  if(is.null(tolerance)) stop(sQuote("tolerance"), " must not be NULL!")

  ## determine model type
  if (is.null(type)) type <-
    if (is.null(y)) "one-classification"
    else if (is.factor(y)) "C-classification"
    else "eps-regression"

  type <- pmatch(type, c("C-classification",
                         "nu-classification",
                         "one-classification",
                         "eps-regression",
                         "nu-regression"), 99) - 1

  if (type > 10) stop("wrong type specification!")

  kernel.name <- kernel
  kernel <- pmatch(kernel, c("linear",
                             "polynomial",
                             "radial",
                             "sigmoid"), 99) - 1

  if (kernel > 10) stop("wrong kernel specification!")
  type <- as.integer(type)
  kernel <- as.integer(kernel)  

  ## further parameter checks
  nr <- p$nRows

  if (!is.vector(y) && !is.factor (y) && type != 2)
    stop("y must be a vector or a factor.")
  if (type != 2 && length(y) != nr)
    stop("x and y don't match.")

  if (cachesize < 0.1)
    cachesize <- 0.1
  
  if (type > 2 && !is.numeric(y))
    stop("Need numeric dependent variable for regression.")

  lev <- NULL
  weightlabels <- NULL

  ## in case of classification: transform factors into integers
  if (type == 2) # one class classification --> set dummy
    y <- rep(1, nr)
  else
    if (is.factor(y)) {
        lev <- levels(y)
        y <- as.integer(y)
        if (!is.null(class.weights)) {
            if (is.null(names(class.weights)))
              stop("Weights have to be specified along with their according level names !")
            weightlabels <- match (names(class.weights), lev)
            if (any(is.na(weightlabels)))
              stop("At least one level name is missing or misspelled.")
        }
    } else {
        if (type < 3) {
            if(any(as.integer(y) != y))
              stop("dependent variable has to be of factor or integer type for classification mode.")
            y <- as.factor(y)
            lev <- levels(y)
            y <- as.integer(y)
        } else lev <- unique(y)
    }

  nclass <- 2
  if (type < 2) nclass <- length(lev)

  if (type > 1 && length(class.weights) > 0) {
      class.weights <- NULL
      warning(sQuote("class.weights"), " are set to NULL for regression mode. For classification, use a _factor_ for ", sQuote("y"),
              ", or specify the correct ", sQuote("type"), " argument.")
  }

  ret <- list(## Allocate results ahead of time and modify in C call
              nclasses = integer  (1),
              tot.nSV  = integer  (1),  # for all classes
              SV       = NULL,          # External pointer for SVs
              labels   = integer  (nclass),
              nSV      = integer  (nclass),
              rho      = double   (nclass * (nclass - 1) / 2),
              coefs    = double   (nr * (nclass - 1)),
              sigma    = if (probability){double (nclass * (nclass - 1) / 2)} else NULL,
              probA    = if (probability){double (nclass * (nclass - 1) / 2)} else NULL,
              probB    = if (probability){double (nclass * (nclass - 1) / 2)} else NULL,
              x.space  = p$x.space)

  .External("svmtrain",
            nRows        = nr,
            x            = p$x,
            y            = as.numeric(y),
            type         = as.integer(type),
            kernel       = kernel,
            degree       = degree,
            gamma        = gamma,
            coef0        = coef0,
            cost         = cost,
            nu           = nu,
            weightlabels = weightlabels,
            weights      = class.weights,
            nWeights     = length(class.weights),
            cachesize    = cachesize,
            tolerance    = tolerance,
            epsilon      = epsilon,
            shrinking    = shrinking,
            probability  = probability,
            seed         = seed,
            ret          = ret)
  ## PACKAGE = "e1071")

  ret <- append(ret,
                list(
                     call        = match.call(),
                     type        = type,
                     kernel.name = kernel.name,
                     kernel      = kernel,
                     cost        = cost,
                     degree      = degree,
                     gamma       = gamma,
                     coef0       = coef0,
                     nu          = nu,
                     epsilon     = epsilon,
                     levels      = lev,
                     nCols       = p$nCols))

  class (ret) <- "svmmodel"
  ## if (fitted) {
  ##       ret$fitted <- na.action(predict(ret, xhold,
  ##                                       decision.values = TRUE))
  ##       ret$decision.values <- attr(ret$fitted, "decision.values")
  ##       attr(ret$fitted, "decision.values") <- NULL
  ##       if (type > 1) ret$residuals <- y - ret$fitted
  ##   }
  ret
}
