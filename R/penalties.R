#' @rdname sic.penalty
#' @title Schwarz Information Criterion penalty
#' @description The function evaluates the penalty term for Schwarz Information Criterion. 
#' If \code{alpha} is greater than 1,  the strengthen SIC proposed proposed in Fryzlewicz (2014) is calculated. This routine is typically not called directly by the user; 
#' its name can be passed as an argument to \code{\link{features}}.
#' @param n The number of observations.
#' @param n.param The number of parameters in the model for which the penalty is evaluated.
#' @param alpha A scalar greater or equal than one. 
#' @param ... Not in use.
#' @export sic.penalty
#' @return the penalty term \eqn{\code{n.param}\times(\log(n))^{\code{alpha}}}{n.param * (log(n))^(alpha)}.
#' @references 
#' R. Baranowski, Y. Chen, and P. Fryzlewicz (2016). Narrowest-Over-Threshold Change-Point Detection.  (\url{http://personal.lse.ac.uk/baranows/not/not.pdf})
#' 
#' P. Fryzlewicz (2014). Wild Binary Segmentation for multiple change-point detection. Annals of Statistics. (\url{http://stats.lse.ac.uk/fryzlewicz/wbs/wbs.pdf}) 
#' @examples
#' #*** a simple example how to use the AIC penalty
#' x <- rnorm(300) + c(rep(1,50),rep(0,250))
#' w <- not(x)
#' w.cpt <- features(w, penalty="sic")
#' w.cpt$cpt[[1]]

sic.penalty <- function(n, n.param, alpha=1.00, ...){
	alpha <- as.numeric(alpha)
	
	pen <- log(n)^alpha
	
	return(n.param*pen)
}


#' @title Akaike Information Criterion penalty
#' @description The function evaluates the penalty term for Akaike Information Criterion. 
#' This routine is typically not called directly by the user; its name can be passed as an argument to \code{\link{features}}.
#' @param n The number of observations.
#' @param n.param The number of parameters in the model for which the penalty is evaluated.
#' @param ... Not in use.
#' @return The penalty term \eqn{2 \times \code{n.param}}{2*n.param}.
#' @references 
#' R. Baranowski, Y. Chen, and P. Fryzlewicz (2016). Narrowest-Over-Threshold Change-Point Detection.  (\url{http://personal.lse.ac.uk/baranows/not/not.pdf})
#' @examples
#' #*** a simple example how to use the AIC penalty
#' x <- rnorm(300) + c(rep(1,50),rep(0,250))
#' w <- not(x)
#' w.cpt <- features(w, penalty="aic")
#' w.cpt$cpt[[1]]

aic.penalty <- function(n, n.param,  ...){
  return(2*n.param)
}