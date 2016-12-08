#' @title Extract likelihood from a 'not' object
#' @description Calculates the Gaussian log-likelihood for the signal estimated using \code{\link{predict.not}} with the change-points at \code{cpt}. The type of the signal depends on 
#' on the value of \code{contrast} that has been passed to \code{\link{not}} (see \code{\link{predict.not}}).
#' @param object An object of class 'not', returned by \code{\link{not}}.
#' @param cpt An integer vector with locations of the change-points.
#' If missing, the \code{\link{features}} is called internally to extract the change-points from \code{object}.
#' @param  ... Further parameters that can be passed to \code{\link{predict.not}} and \code{\link{features}}.
#' @method logLik not
#' @export
#' @rdname loglik.not
#' @examples
#' #' # **** Piecewisce-constant mean with Gaussian noise.
#' x <- c(rep(0, 100), rep(1,100)) + rnorm(100)
#' # *** identify potential locations of the change-points
#' w <- not(x, contrast = "pcwsConstMean")
#' # *** log-likelihood for the model with the change-point estimated  via 'not'
#' logLik(w)
#' # *** log-likelihood for the model with the change-point at 100
#' logLik(w, cpt=100)

logLik.not <- function(object, cpt, ...){
  
  if(missing(cpt)) cpt <- features(object, ...)$cpt
  if(!is.null(cpt)) if(any(is.na(cpt))) cpt <- cpt[!is.na(cpt)]
  
  
  n.cpt <- length(cpt)
  

  if(object$contrast == "pcwsConstMean" ||
     object$contrast == "pcwsConstMeanHT" ||
     object$contrast == "pcwsLinMean" ||
     object$contrast == "pcwsLinContMean" ||
     object$contrast == "pcwsQuadMean"){
    
    res <- - length(object$x)/2 * log(mean(residuals(object, cpt=cpt)^2))
    
  }else if (object$contrast == "pcwsConstMeanVar") {
    
    fit <- predict(object, cpt=cpt)
    int.len <- diff(sort(unique(c(cpt,0,length(object$x)))))
    res <- -sum((int.len) * (log(fit[c(cpt, length(object$x)),2])))
    
  }
  
    
  if(object$contrast == "pcwsConstMean" || object$contrast == "pcwsConstMeanHT") attr(res, "df") <- 2*n.cpt + 2
  else if (object$contrast == "pcwsLinMean") attr(res, "df") <- 3*n.cpt + 3
  else if (object$contrast == "pcwsLinContMean") attr(res, "df") <- 2*n.cpt + 3
  else if (object$contrast == "pcwsQuadMean") attr(res, "df") <- 4*n.cpt + 4
  else if (object$contrast == "pcwsConstMeanVar") attr(res, "df") <- 3*n.cpt + 2
  
  
  class(res) <- "logLik"
  
  return(res)

}