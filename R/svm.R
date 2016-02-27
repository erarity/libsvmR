## Types of problem
typenames <-  c("C-classification",
                "nu-classification",
                "one-classification",
                "eps-regression",
                "nu-regression")
## Kernels
kernelnames <- c("linear",
                 "polynomial",
                 "radial",
                 "sigmoid")

predict.svmmodel <- function (model, newdata,
                              probability=FALSE, multicore=FALSE, ...){
    stopifnot(class(model) == "svmmodel",
              !missing(newdata),
              class(newdata) == "svmproblem",
              model$nCols == newdata$nCols,
              model$tot.nSV > 0)

    ret <- list(class = double(newdata$nRows),
                values = double(newdata$nRows * model$nclasses * (model$nclasses - 1) / 2),
                prob = double(newdata$nRows * model$nclasses))
    
    .External("svmpredict",
              train = newdata$x,
              nRows = newdata$nRows,
              probability = as.integer(probability),
              model = model,
              ret = ret,
              PACKAGE = "libsvmR")
    ret
}

print.svmproblem <-  function(x,...) {
    cat("\nSVM Problem Data\n")
    cat("nRows:", x$nRows,"\n")
    cat("nCols:", x$nCols,"\n")
    cat("y: ",if(is.null(x$y)){"Missing"}else{"Present"},"\n\n")
}
                 
"[.svmproblem" <- function(x, i, j, ..., drop){
    x$x <- x$x[i]
    stopifnot(!sapply(x$x,is.null))
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
    cat("   SVM-Type: ", typenames[x$type+1], "\n")
    cat(" SVM-Kernel: ", kernelnames[x$kernel+1], "\n")
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
        if (x$probability)
            cat("Sigma: ", x$sigma, "\n\n")
    }

    cat("\nNumber of Support Vectors: ", x$tot.nSV)
    cat("\n\n")

}

svm.debuglvl <- function(d){
  invisible(.C("svm_set_debug",as.integer(d)))
}


## Compute the hyperplane, only works for a linear model
compute.hyperplane <- function(mod){
  Reduce(
         function(cur,i){
           cur + mod$coefs[[i]] * node2vec(mod$SV[[i]],mod$nCols)
         },
         1:mod$tot.nSV,
         numeric(mod$nCols))
}

## Train an svm
svm <- function (p,
	  skips	      = NULL,	#MODTAG
          y           = NULL,
          type        = NULL,
          kernel      = "linear",
          degree      = 3,
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
          max.classes = 2,              # This or fewer classes automatically becomes classification
          ...)
{
  stopifnot(class(p) == "svmproblem")
  if(is.null(y)){
    stopifnot(!is.null(p$y))
    y <- p$y
  }

  #MODTAG
  # Is the list captured properly?
  #cat("Captured skip list (R) as: ",skips)

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
    else if (length(unique(y)) <= max.classes) "C-classification"
    else "eps-regression"

  type <- pmatch(type, c(typenames,99))-1
  ## ("C-classification",
  ##                        "nu-classification",
  ##                        "one-classification",
  ##                        "eps-regression",
  ##                        "nu-regression"), 99) - 1

  if (type > 10) stop("wrong type specification!")

  kernel.name <- kernel
  kernel <- pmatch(kernel, c(kernelnames,99))-1
  ## c("linear",
  ##                            "polynomial",
  ##                            "radial",
  ##                            "sigmoid"), 99) - 1

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
            y <- as.ordered(y)
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
              SVidx    = NULL,          # Indices to support vectors
              labels   = integer  (nclass),
              nSV      = integer  (nclass),
              rho      = double   (nclass * (nclass - 1) / 2),
              coefs    = double   (nr * (nclass - 1)),
              sigma    = if (probability){double (nclass * (nclass - 1) / 2)} else NULL,
              probA    = if (probability){double (nclass * (nclass - 1) / 2)} else NULL,
              probB    = if (probability){double (nclass * (nclass - 1) / 2)} else NULL,
              x.space  = p$x.space)

  .External("svmtrain",
	    skips	 = as.integer(skips), 	#MODTAG
	    numskips	 = length(skips),	#MODTAG
            nRows        = nr,
            x            = p$x,
            y            = as.numeric(y),
            type         = as.integer(type),
            kernel       = kernel,
            degree       = as.integer(degree),
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
            ret          = ret,
            PACKAGE = "libsvmR")
  ret <- append(ret,
                list(
                     call        = match.call(),
                     type        = type,
                     type.name   = typenames[type+1],
                     kernel.name = kernel.name,
                     kernel      = kernel,
                     cost        = cost,
                     degree      = as.integer(degree),
                     gamma       = gamma,
                     coef0       = coef0,
                     nu          = nu,
                     epsilon     = epsilon,
                     levels      = lev,
                     probability  = probability,
                     hyperplane  = NULL,
                     nCols       = p$nCols))

  if(kernel.name=="linear"){
    ret$hyperplane <- compute.hyperplane(ret)
  }

  class(ret) <- "svmmodel"
  ret
}
            

             
read.svmproblem <- function(filename, max.classes=2){
  stopifnot(is.character(filename),
            file.exists(filename))
  filename <- path.expand(filename)
  prob <- .Call("svm_read_problem",
                filename,
                PACKAGE="libsvmR"
                )
  ## Should have a list of 4 things
  ## length, targetvalues, feature pointers, node space
  names(prob) <- c("nRows","nCols","y","x","x.space")
  if (length(unique(prob$y)) <= max.classes) {    #2-class problem
    prob$y <- as.ordered(prob$y)
  }

  class(prob) <- "svmproblem"
  prob
}

save.svmmodel <- function(filename,model){
  stopifnot(class(model) == "svmmodel",  
            is.character(filename))
  filename <- path.expand(filename)
  invisible(.Call("svm_save_model",
                  filename,
                  model,
                  PACKAGE = "libsvmR"))
}

load.svmmodel <- function(filename){
  stopifnot(is.character(filename),
            file.exists(filename))
  ret <- list( ## Allocate results ahead of time and modify in C call
              type     = integer (1),
              kernel   = integer (1),
              nclasses = integer (1),
              degree = integer (1),
              gamma = numeric (1),
              coef0 = numeric (1),
              tot.nSV  = integer (1),  # for all classes
              SV       = NULL,          # External pointer for SVs
              labels   = NULL,
              nSV      = NULL,
              rho      = NULL,
              coefs    = NULL,
              probA    = NULL,
              probB    = NULL,
              x.space  = NULL)
  .Call("svm_load_model_R",
        filename,
        ret,
        PACKAGE = "libsvmR")
  
  ret <- append(ret,
                list(call = match.call(),
                     type.name   = typenames[ret$type+1],
                     kernel.name = kernelnames[ret$kernel+1]
                     ))
  class(ret) <- "svmmodel"
  ret
}

crossvalidate <- function(learn,        #Function to do learning
                          data,         #Features
                          targets,      #target values/classes
                          metric,

                          folds = 10,   #Number of folds or indices
                          dopredict = predict,      #Function to do prediction
                          seed = 1L,
                          cv.iter = lapply,
                          ...){
  
  set.seed(as.integer(seed))
  folds <- as.integer(folds)
  nfold <- max(unique(folds))
  foldids <-
    if(length(folds) == 1){               #nway-folds, construct ids
      ord <- sample(1:length(targets))    #Random ordering
      rep(1:folds,length.out=length(targets))[ord]
    }
    else {
      folds
    }
  
  cv <- cv.iter(1:nfold,function(i){
    cat("Fold ",i,"\n")
    train <- which(foldids!=i)
    test <-  which(foldids==i)
    model <- learn(data[train,], targets[train], ...)
    pred <- dopredict(model,data[test,])
    res <- metric(pred, targets[test], testids=test,...)
    cat("\n")
    res
  })

  t(simplify2array(as.matrix(cv)))
}
               

## Make a set of parameters for an SVM
make.svmparams <- function(...){
  p <- list(...)
  expand.grid(p,stringsAsFactors=FALSE)
}  
                           



tune.svm <- function(x, params, metric, tune.iter = lapply, ...){
  tuned <-
    tune.iter(1:nrow(params), function(p){
      par <- lapply(params[p,],unlist)
      c("--------------------------\n")
      
      cv <- do.call(crossvalidate,
                    c(list(learn=svm, data=x, targets=x$y, metric=metric),
                      par,
                      list(...)))
      means <- apply(cv,2,mean)
      names(means) <- paste("mean.",names(means),sep="")
      sds <- apply(cv,2,sd)
      names(sds) <- paste("sd.",names(sds),sep="")
      res <- c(par,means,sds)
      invisible(print(format(c("Results: ",res))))
      res
    }
  )
  data.frame(t(simplify2array(tuned)))
}

## Convert a node to a dense vector.  Dangerous as they types can't be
## checked
node2vec <- function(node,ncols){
  stopifnot(mode(node) == "externalptr")
  vec <- numeric(ncols)
  .Call("svm_node_as_vec",
        node,
        vec)
  vec
}
