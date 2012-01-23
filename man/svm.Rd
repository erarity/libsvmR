\name{svm}
\alias{svm}
\alias{summary.svm}
\alias{print.svm}
\alias{svmmodel}
\title{Support Vector Machines}
\description{
\code{svm} is used to train a support vector machine. It can be used to carry
out general regression and classification (of nu and epsilon-type), as
well as density-estimation.
}
\usage{
\method{svm}{default}(x, y = NULL, type = NULL, kernel =
"radial", degree = 3, gamma = 1 / ncol(x),
coef0 = 0, cost = 1, nu = 0.5,
class.weights = NULL, cachesize = 40, tolerance = 0.001, epsilon = 0.1,
shrinking = TRUE, probability = FALSE, seed = 1L)
}
\arguments{
  \item{x}{an svmproblem created with \code{\link{read.svmproblem}}}
  \item{y}{a response vector with one label for each row/component of
    \code{x}. Can be either a factor (for classification tasks)
    or a numeric vector (for regression).  If not provided, x$y is
    used.}
  \item{type}{\code{svm} can be used as a classification
    machine, as a regression machine, or for novelty detection.
    Depending of whether \code{y} is
    a factor or not, the default setting for \code{type} is \code{C-classification} or \code{eps-regression}, respectively, but may be overwritten by setting an explicit value.\cr
    Valid options are:
    \itemize{
      \item \code{C-classification}
      \item \code{nu-classification}
      \item \code{one-classification} (for novelty detection)
      \item \code{eps-regression}
      \item \code{nu-regression}
    }
  }
  \item{kernel}{the kernel used in training and predicting. You
    might consider changing some of the following parameters, depending
    on the kernel type.\cr
    \describe{
      \item{linear:}{\eqn{u'v}{u'*v}}
      \item{polynomial:}{\eqn{(\gamma u'v + coef0)^{degree}}{(gamma*u'*v + coef0)^degree}}
      \item{radial basis:}{\eqn{e^(-\gamma |u-v|^2)}{exp(-gamma*|u-v|^2)}}
      \item{sigmoid:}{\eqn{tanh(\gamma u'v + coef0)}{tanh(gamma*u'*v + coef0)}}
      }
    }
  \item{degree}{parameter needed for kernel of type \code{polynomial} (default: 3)}
  \item{gamma}{parameter needed for all kernels except \code{linear}
    (default: 1/(data dimension))}
  \item{coef0}{parameter needed for kernels of type \code{polynomial}
    and \code{sigmoid} (default: 0)}
  \item{cost}{cost of constraints violation (default: 1)---it is the
    \sQuote{C}-constant of the regularization term in the Lagrange formulation.}
  \item{nu}{parameter needed for \code{nu-classification},
    \code{nu-regression}, and \code{one-classification}}
  \item{class.weights}{a named vector of weights for the different
    classes, used for asymmetric class sizes. Not all factor levels have
    to be supplied (default weight: 1). All components have to be named.}
  \item{cachesize}{cache memory in MB (default 40)}
  \item{tolerance}{tolerance of termination criterion (default: 0.001)}
  \item{epsilon}{epsilon in the insensitive-loss function (default: 0.1)}
  \item{shrinking}{option whether to use the shrinking-heuristics
    (default: \code{TRUE})}
  \item{probability}{logical indicating whether the model should
    allow for probability predictions.}
  \item{seed}{integer seed for libsvm (used for cross-validation and
    probability prediction models).}
}

\value{
  An object of class \code{"svmmodel"} containing the fitted model, including:
  \item{SV}{Pointers to the support vectors}
  \item{coefs}{The corresponding coefficients times the training labels.}
  \item{rho}{The negative intercept.}
  \item{sigma}{In case of a probabilistic regression model, the scale
    parameter of the hypothesized (zero-mean) laplace distribution estimated by
    maximum likelihood.}
  \item{probA, probB}{numeric vectors of length k(k-1)/2, k number of
    classes, containing the parameters of the logistic distributions fitted to
    the decision values of the binary classifiers (1 / (1 + exp(a x + b))).}
}
\details{
  For multiclass-classification with k levels, k>2, \code{libsvm} uses the
  \sQuote{one-against-one}-approach, in which k(k-1)/2 binary classifiers are
  trained; the appropriate class is found by a voting scheme.
  
  If the predictor variables include factors, the formula interface must be used to get a
  correct model matrix.

  \code{plot.svm} allows a simple graphical
  visualization of classification models.

  The probability model for classification fits a logistic distribution
  using maximum likelihood to the decision values of all binary
  classifiers, and computes the a-posteriori class probabilities for the
  multi-class problem using quadratic optimization. The probabilistic
  regression model assumes (zero-mean) laplace-distributed errors for the
  predictions, and estimates the scale parameter using maximum likelihood.
}
\note{
Data are scaled internally, usually yielding better results.

Parameters of SVM-models usually \emph{must} be tuned to yield sensible results!
}
\references{
  \itemize{
    \item
      Chang, Chih-Chung and Lin, Chih-Jen:\cr
      \emph{LIBSVM: a library for Support Vector Machines}\cr
      \url{http://www.csie.ntu.edu.tw/~cjlin/libsvm}

     \item
      Exact formulations of models, algorithms, etc. can be found in the
      document:\cr
      Chang, Chih-Chung and Lin, Chih-Jen:\cr
      \emph{LIBSVM: a library for Support Vector Machines}\cr
      \url{http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.ps.gz}
    
      \item
      More implementation details and speed benchmarks can be found on:
      Rong-En Fan and Pai-Hsune Chen and Chih-Jen Lin:\cr
      \emph{Working Set Selection Using the Second Order Information for Training SVM}\cr
      \url{http://www.csie.ntu.edu.tw/~cjlin/papers/quadworkset.pdf}
    
  }
}
\author{
  Chris Kauffman, David Meyer (based on C/C++-code by Chih-Chung Chang and Chih-Jen Lin)\cr
  \email{kauffman@cs.umn.edu}
}
\seealso{
  \code{\link{predict.svmmodel}}
  \code{\link{read.svmproblem}}
  \code{\link{save.svmmodel}}
  \code{\link{load.svmmodel}}
}
\keyword{neural}
\keyword{nonlinear}
\keyword{classif}









