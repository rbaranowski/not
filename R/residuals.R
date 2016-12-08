#' @title Extract residuals from a 'not' object
#' @description Returns a difference between \code{x} in \code{object} and the estimated signal with change-points at \code{cpt}.
#' Type of the signal depends on the value of \code{contrast} that has been passed to \code{\link{not}} in order to construct \code{object} (see details of \code{\link{predict.not}}).
#' @param object An object of class 'not', returned by \code{\link{not}}.
#' @param cpt An integer vector with locations of the change-points.
#' If missing, the \code{\link{features}} is called internally to extract the change-points from \code{object}.
#' @param type Choice of "raw" and "standardised".
#' @param  ... Further parameters that can be passed to \code{\link{predict.not}} and \code{\link{features}}.
#' @method residuals not
#' @rdname residuals.not
#' @examples
#' pcws.const.sig <- c(rep(0, 100), rep(1,100))
#' x <- pcws.const.sig + rnorm(100)
#' w <- not(x, contrast = "pcwsConstMean")
#' # *** plot residuals obtained via fitting piecewise-constant function with estimated change-points
#' plot(residuals(w))
#' # *** plot residuals with obtained via fitting piecewise-constant function with true change-point
#' plot(residuals(w, cpt=100))
#'# *** plot standardised residuals
#' plot(residuals(w, type="standardised"))
#' @export
#' @return If \code{type="raw"}, the difference between the data and the estimated signal. If \code{type="standardised"}, the difference between the data and the estimated signal, divided by the estimated standard deviation.

residuals.not <-
  function(object,
           cpt,
           type = c("raw", "standardised"),
           ...) {
    type <- match.arg(type, c("raw", "standardised"))
    
    if (missing(cpt))
      cpt <- features(object, ...)$cpt
    if (!is.null(cpt))
      if (any(is.na(cpt)))
        cpt <- cpt[!is.na(cpt)]
      
      
      if (object$contrast == "pcwsConstMean" ||
          object$contrast == "pcwsConstMeanHT" ||
          object$contrast == "pcwsLinMean" ||
          object$contrast == "pcwsLinContMean" ||
          object$contrast == "pcwsQuadMean") {
        res <- object$x -  predict(object, cpt = cpt)
        
        if (type == "raw")
          return(res)
        else
          return(res / sd(res))
        
      } else if (object$contrast == "pcwsConstMeanVar") {
        fit <- predict(object, cpt = cpt)
        int.len <- diff(sort(unique(c(
          cpt, 0, length(object$x)
        ))))
        
        if (type == "raw")
          return((object$x -  fit[, 1]))
        else
          return((object$x -  fit[, 1]) / fit[, 2])
        
      }
      
  }