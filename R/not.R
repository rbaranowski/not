#' @title Narrowest-Over-Threshold Change-Point Detection
#' @description Implements the Narrowest-Over-Threshold approach for general multiple change-point 
#' detection in one-dimensional data following 'deterministic signal + noise' model. Scenarios that are currently implemented are: piecewise-constant signal, piecewise-constant signal with a heavy tailed noise, piecewise-linear signal, piecewise-quadratic signal, piecewise-constant signal and with piecewise-constant standard deviation of the noise. The main routines of the package are \code{\link{not}} and \code{\link{features}}.	
#' @references R. Baranowski, Y. Chen, and P. Fryzlewicz (2016). Narrowest-Over-Threshold Change-Point Detection.  (\url{http://personal.lse.ac.uk/baranows/not/not.pdf})
#' @docType package
#' @useDynLib not
#' @name not-package
#' @importFrom  graphics plot lines title
#' @importFrom stats cov lm.fit logLik predict residuals sd mad quantile runif ts.plot var
#' @importFrom splines bs


NULL

#' @title Narrowest-Over-Threshold Change-Point Detection
#' @description Identifies potential locations of the change-points in the data following 'deterministic signal + noise' model (see details below) in a number of different scenarios.
#' The object returned by this routine can be further passed to the \code{\link{features}} function,  which finds the final estimate of the change-points based on a chosen stopping criterion. 
#' It can be also passed to \code{\link{plot}}, \code{\link{predict}} and \code{\link{residuals}} routines.
#' @details The data points provided in \code{x} are assumed to follow
#'  \deqn{Y_{t} = f_{t}+\sigma_{t}\varepsilon_{t},}{Y_t= f_t + sigma_t varepsilon_t,}
#'  for \eqn{t=1,\ldots,n}{t=1,...,n}, where \eqn{n}{n} is the number of observations in \code{x}, the signal \eqn{f_{t}}{f_t} and the standard deviation \eqn{\sigma_{t}}{sigma_{t}} 
#'  are non-stochastic with structural breaks at unknown locations in time \eqn{t}{t}. Currently, thefollowing scenarios for \eqn{f_{t}}{f_t} and \eqn{\sigma_{t}}{sigma_t} are implemented:
#'  \itemize{
#'    \item Piecewise-constant signal with a Gaussian noise and constant standard deviation. 
#'    
#'    Use \code{contrast="pcwsConstMean"} here.
#'    \item Piecewise-constant mean with a heavy-tailed noise and constant standard deviation. 
#'    
#'    Use \code{contrast="pcwsConstMeanHT"} here.
#'    \item Piecewise-linear continuous signal with Gaussian noise and constant standard deviation. 
#'    
#'    Use \code{contrast="pcwsLinContMean"} here.
#'    \item Piecewise-linear signal with Gaussian noise and constant standard deviation.
#'    
#'      Use \code{contrast="pcwsLinMean"} here.
#'    \item Piecewise-quadratic signal with Gaussian noise and constant standard deviation. 
#'    
#'    Use \code{contrast="pcwsQuadMean"} here.
#'    \item Piecewise-constant signal and piecewise-constant standard deviation of the Gaussian noise. 
#'    
#'    Use \code{contrast="pcwsConstMeanVar"} here.
#'  }
#' @param x A numeric vector with data points.
#' @param M A number of intervals drawn in the procedure. 
#' @param method  Choice of "not" (recommended) and "max". If \code{method="not"}, the Narrowest-Over-Threshold intervals are used in the algorithm. 
#' If \code{method="max"}, the intervals corresponding to the largest contrast function are used. For an explanation, see the references.
#' @param contrast A type of the contrast function used in the NOT algorithm. 
#' Choice of \code{"pcwsConstMean"}, \code{"pcwsConstMeanHT"}, \code{"pcwsLinContMean"}, \code{"pcwsLinMean"}, \code{"pcwsQuadMean"}, \code{"pcwsConstMeanVar"}. 
#' For the explanation, see details below. 
#' @param rand.intervals A logical variable. If \code{rand.intervals=TRUE} intervals used in the procedure are drawn uniformly using the \code{\link{random.intervals}} routine. 
#' If \code{rand.intervals=FALSE}, the intervals need to be passed using the \code{intervals} argument.
#' @param parallel A logical variable. If TRUE some of computations are run in parallel using OpenMP framework. Currently this option is not supported on Windows.
#' @param augmented A logical variable. if TRUE, the entire data are considered when the NOT segmentation tree is constructed (see the solution path algorithm in the references).
#' @param intervals  A 2-column matrix with the intervals considered in the algorithm, with start- and end- points of the intervals in, respectively, the first and  the second column. 
#' The intervals are used only if \code{rand.intervals=FALSE}. 
#' @param ... Not in use.
#' @rdname not
#' @export
#' @examples 
#' # **** Piecewisce-constant mean with Gaussian noise.
#' # *** signal
#' pcws.const.sig <- c(rep(0, 100), rep(1,100))
#' # *** data vector
#' x <- pcws.const.sig + rnorm(100)
#' # *** identify potential locations of the change-points
#' w <- not(x, contrast = "pcwsConstMean") 
#' # *** some examples of how the w object can be used
#' plot(w)
#' plot(residuals(w))
#' plot(predict(w))
#' # *** this is how to extract the change-points
#' fo <- features(w)
#' fo$cpt
#' 
#' # **** Piecewisce-constant mean with a heavy-tailed noise.
#' # *** data vector, signal the same as in the previous example, but heavy tails
#' x <- pcws.const.sig + rt(100, 3) 
#' # *** identify potential locations of the change-points, 
#' # using a contrast taylored to heavy-tailed data
#' w <- not(x, contrast = "pcwsConstMeanHT") 
#' plot(w)
#' 
#' # **** Piecewisce-constant mean and piecewise-constant variance
#' # *** signal's standard deviation
#' pcws.const.sd <- c(rep(2, 50), rep(1,150))
#' # *** data vector with pcws-const mean and variance
#' x <- pcws.const.sig + pcws.const.sd * rnorm(100)
#' # *** identify potential locations of the change-points in this model
#' w <- not(x, contrast = "pcwsConstMeanVar") 
#' # *** extracting locations of the change-points
#' fo <- features(w)
#' fo$cpt
#' 
#' # **** Piecewisce-linear coninuous mean
#' # *** signal with a change in slope
#' pcws.lin.cont.sig <- cumsum(c(rep(-1/50, 100), rep(1/50,100)))
#' # *** data vector 
#' x <- pcws.lin.cont.sig +  rnorm(100)
#' # *** identify potential locations of the change-points in the slope coefficient
#' w <- not(x, contrast = "pcwsLinContMean") 
#' # *** ploting the results
#' plot(w)
#' # *** location(s) of the change-points
#' fo <- features(w)
#' fo$cpt
#' 
#' # **** Piecewisce-linear mean with jumps
#' # *** signal with a change in slope and jumpe
#' pcws.lin.sig <- pcws.lin.cont.sig + pcws.const.sig
#' # *** data vector 
#' x <- pcws.lin.sig +  rnorm(100)
#' # *** identify potential locations of the change-points in the slope coefficient and the intercept
#' w <- not(x, contrast = "pcwsLinMean") 
#' # *** ploting the results
#' plot(w)
#' # *** location(s) of the change-points
#' fo <- features(w)
#' fo$cpt
#' 
#' # **** Piecewisce-quadratic mean with jumps
#' # *** Piecewise-quadratic signal
#' pcws.quad.sig <- 2*c((1:50)^2 /1000, rep(2, 100), 1:50 / 50 )
#' # *** data vector 
#' x <- pcws.quad.sig +  rnorm(100)
#' # *** identify potential locations of the change-points in the slope coefficient and the intercept
#' w <- not(x, contrast = "pcwsQuadMean") 
#' # *** ploting the results
#' plot(w)
#' # *** location(s) of the change-points
#' fo <- features(w)
#' fo$cpt

not <- function(x, ...)  UseMethod("not")

#' @method not default
#' @export 
#' @rdname not
#' @return An object of class "not", which contains the following fields:
#' \item{x}{The input vector.}
#' \item{n}{The length of \code{x}.}
#' \item{contrast}{A scenario for the change-points.}
#' \item{contrasts}{A 5-column matrix with the values of the contrast function, where 's' and 'e' denote start-
#' end points of the intervals in which change-points candidates 'arg.max' have been found; 'length' shows the length of the intervals drawn,
#' column 'max.contrast' contains corresponding value of the contrast statistic.}
#' \item{solution.path}{A list with the solution path of the NOT algorithm (see the references) containing three fields of the same length: \code{cpt} - a list with consecutive solutions, i.e. s the sets of change-point candidates, 
#' \code{th} - a vector of thresholds corresponding to the solutions, \code{n.cpt} - a vector with the number of change-points for each solution.}
#' @references R. Baranowski, Y. Chen, and P. Fryzlewicz (2016). Narrowest-Over-Threshold Change-Point Detection.  (\url{http://personal.lse.ac.uk/baranows/not/not.pdf})

not.default <- function(x, 
                        M=10000,  
                        method = c("not", "max"),
                        contrast = c("pcwsConstMean", 
                                          "pcwsConstMeanHT", 
                                          "pcwsLinContMean",
                                          "pcwsLinMean",
                                          "pcwsQuadMean",
                                          "pcwsConstMeanVar"),
                        rand.intervals = TRUE, 
                        parallel = FALSE,
                        augmented = FALSE,
                        intervals,
                        ...){
	
  results <- list()	
  
  #veryfing the input parameters - x
	results$x <- as.numeric(x)
	storage.mode(results$x) <- "double"
	results$n <- length(results$x)
	if(results$n < 2) stop("Data vector x should contain at least two elements.")
	if(any(is.na(results$x))) stop("Data vector x vector cannot contain NA's")
	if(var(x) <= sqrt(.Machine$double.eps)) stop("Data vector x is essentially a constant vector, change-point detection is not needed")
	
	#veryfing the input parameters - M
	M <- as.integer(M)
	if(any(is.na(M))) stop("M cannot be NA")
	if(length(M)> 1)  stop("M should be a single integer.")
	if(M<0)  stop("M should be an integer > 0.")
  
	#veryfing the input parameters - method
	method  <- match.arg(method, c("not", "max"))
	
	#veryfing the input parameters - contrast
	results$contrast <- match.arg(contrast, 
	                                   c("pcwsConstMean", 
	                                     "pcwsConstMeanHT", 
	                                     "pcwsLinContMean",
	                                     "pcwsLinMean",
	                                     "pcwsQuadMean",
	                                     "pcwsConstMeanVar"))
	
	#veryfing the input parameters - rand.intervals
	rand.intervals <- as.logical(rand.intervals[1])
	
	#veryfing the input parameters - parallel
	parallel <- as.integer(parallel[1])
  if(parallel != 0 ) parallel <- as.integer(1)
	
	if(.Platform$OS.type != "unix") 
	  if(parallel == 1) {
	    warning("Running computations in parallel on non-unix systems is currently not supported. Using a single core.")
	    parallel <- as.integer(0) 
	  }
	
	#veryfing the input parameters - augmented
	augmented <- as.logical(augmented[1])
	
	#drawing the intervals over which the contrast function will be computed
	if(rand.intervals){
	  intervals <-  matrix(random.intervals(results$n,M, ...), ncol=2)
	} else if(missing(intervals)){
	  stop("intervals need to be specified when rand.intervals==FALSE")
	}else{
	  intervals <- matrix(as.integer(intervals), ncol=2)
	  M <- nrow(intervals)
	  
	  if(any(intervals[,1] > intervals[,2])) stop("intervals[,2] must be > than intervals[,1]")
	  if(any(intervals < 1 || intervals > length(x))) stop("endpoints in intervals must be between 1 and length(x)")
	  
	}
    
	
	#making sure that not_r_wrapper gets the right type
	storage.mode(intervals) <- "integer"
	
	
	#choosing the way the contrast functions are aggregated
	if(method  == "not") method <- as.integer(0)
	else if(method  == "max") method <- as.integer(1)
	
	#should we additionally compute the contrast function while constructing the segmenation tree?
	if(augmented == FALSE) augmented <- as.integer(0)
	else augmented <- as.integer(1)
	
	#picking the correct type of the contrast function
	if(results$contrast == "pcwsConstMean") contrast <- as.integer(0)
	else if (results$contrast == "pcwsLinContMean") contrast <- as.integer(1)
	else if (results$contrast == "pcwsLinMean") contrast <- as.integer(2)
	else if (results$contrast == "pcwsQuadMean") contrast <- as.integer(3)
	else if (results$contrast == "pcwsConstMeanVar") contrast <- as.integer(4)
	else if (results$contrast == "pcwsConstMeanHT") contrast <- as.integer(5)
	
	#calling the C code
	tmp <-  .Call("not_r_wrapper", results$x, intervals, method, contrast, parallel, augmented)
	
	#combining the results with the info on how the method has been run
	results <- c(results, tmp)
  class(results) <- "not"
  
	return(results)
	
}





