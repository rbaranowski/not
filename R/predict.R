#' @title Estimate signal for a 'not' object.
#' @description Estimates signal in \code{object$x} with change-points at \code{cpt}. The type of the signal depends on 
#' on the value of \code{contrast} that has been passed to \code{\link{not}} (see details below).
#' @details The data points provided in \code{object$x} are assumed to follow
#'  \deqn{Y_{t} = f_{t}+\sigma_{t}\varepsilon_{t},}{Y_t= f_t + sigma_t varepsilon_t,}
#'  for \eqn{t=1,\ldots,n}{t=1,...,n}, where \eqn{n}{n} is the number of observations in \code{object$x}, the signal \eqn{f_{t}}{f_t} and the standard deviation \eqn{\sigma_{t}}{sigma_{t}} 
#'  are non-stochastic with change-points at locations given in \code{cpt} and \eqn{\varepsilon_{t}}{varepsilon_t} is a white-noise. Denote by \eqn{\tau_{1}, \ldots, \tau_{q}}{tau_1, ..., tau_q} 
#'  the elements in \code{cpt} and set \eqn{\tau_{0}=0}{tau_0=0} and \eqn{\tau_{q+1}=T}{tau_q+1=T}. Depending on the value of \code{contrast} that has been passed to \code{\link{not}} to construct \code{object},  the returned value is calculated as follows.
#'  \itemize{
#'    \item For \code{contrast="pcwsConstantMean"} and \code{contrast="pcwsConstantMeanHT"}, in each  segment  \eqn{(\tau_{j}+1, \tau_{j+1})}{(tau_j +1, tau_(j+1))},
#'    \eqn{f_{t}}{f_t} for \eqn{t\in(\tau_{j}+1, \tau_{j+1})}{t in (tau_j +1, tau_(j+1))} is approximated by the mean of \eqn{Y_{t}}{Y_t} calculated over \eqn{t\in(\tau_{j}+1, \tau_{j+1})}{t in (tau_j +1, tau_(j+1))}. 
#'    \item For \code{contrast="pcwsLinContMean"}, \eqn{f_{t}}{f_{t}} is approximated by the linear spline fit with knots at \eqn{\tau_{1}, \ldots, \tau_{q}}{tau_1, ..., tau_q} minimising the l2 distance between the fit and the data.
#'    \item For \code{contrast="pcwsLinMean"} in each  segment  \eqn{(\tau_{j}+1, \tau_{j+1})}{(tau_j +1, tau_(j+1))}, the signal
#'    \eqn{f_{t}}{f_t} for \eqn{t\in(\tau_{j}+1, \tau_{j+1})}{t in (tau_j +1, tau_(j+1))} is approximated by the line \eqn{\alpha_{j} + \beta_{j} t}{alpha_j + beta_j j}, where the regression coefficients are 
#'    found using the least squares method.
#'    \item For \code{contrast="pcwsQuad"}, the signal
#'    \eqn{f_{t}}{f_t} for \eqn{t\in(\tau_{j}+1, \tau_{j+1})}{t in (tau_j +1, tau_(j+1))} is approximated by the curve \eqn{\alpha_{j} + \beta_{j} t + \gamma_{j} t^2}{alpha_j + beta_j j + gamma_j^2}, where the regression coefficients are 
#'    found using the least squares method.
#'    \item For \code{contrast="pcwsConstMeanVar"},  in each  segment  \eqn{(\tau_{j}+1, \tau_{j+1})}{(tau_j +1, tau_(j+1))}, 
#'    \eqn{f_{t}}{f_t} and \eqn{\sigma_{t}}{sigma_t} for \eqn{t\in(\tau_{j}+1, \tau_{j+1})}{t in (tau_j +1, tau_(j+1))} are approximated by, respectively, the mean and the standard deviation of \eqn{Y_{t}}{Y_t}, both calculated over \eqn{t\in(\tau_{j}+1, \tau_{j+1})}{t in (tau_j +1, tau_(j+1))}.     
#'  }
#' @param object An object of class 'not', returned by \code{\link{not}}.
#' @param cpt An integer vector with locations of the change-points.
#' If missing, the \code{\link{features}} is called internally to extract the change-points from \code{object}.
#' @param  ... Further parameters that can be passed to \code{\link{predict.not}} and \code{\link{features}}.
#' @method predict not
#' @export
#' @rdname predict.not
#' @seealso \code{\link{not}}
#' @examples
#' # **** Piecewisce-constant mean with Gaussian noise.
#' x <- c(rep(0, 100), rep(1,100)) + rnorm(100)
#' # *** identify potential locations of the change-points
#' w <- not(x, contrast = "pcwsConstMean")
#' # *** when 'cpt' is omitted, 'features' function is used internally 
#' # to choose change-points locations
#' signal.est <- predict(w)
#' # *** estimate the signal specifying the location of the change-point
#' signal.est.known.cpt <- predict(w, cpt=100)
#' # *** pass arguments of the 'features' function through 'predict'.
#' signal.est.aic <- predict(w, penalty.type="aic")
#'
#' # **** Piecewisce-constant mean and variance with Gaussian noise.
#' x <- c(rep(0, 100), rep(1,100)) + c(rep(2, 100), rep(1,100)) * rnorm(100)
#' # *** identify potential locations of the change-points
#' w <- not(x, contrast = "pcwsConstMeanVar")
#' # *** here signal is two-dimensional
#' signal.est <- predict(w)
#' @return A vector wit the estimated signal or a two-column matrix with the estimated estimated signal and standard deviation if \code{contrast="pcwsConstMeanVar"} was used to construct \code{object}.

predict.not <- function(object, cpt, ...) {
  if (missing(cpt))
    cpt <- features(object, ...)$cpt
  if (!is.null(cpt))
    if (any(is.na(cpt)))
      cpt <- cpt[!is.na(cpt)]
    
    
    cpt <- sort(unique(c(cpt, 0, length(object$x))))
    
    
    if (object$contrast == "pcwsConstMean" ||
        object$contrast == "pcwsConstMeanHT") {
      fit <- rep(0, length(object$x))
      
      for (i in 1:(length(cpt) - 1)) {
        fit[(cpt[i] + 1):cpt[i + 1]] <- mean(object$x[(cpt[i] + 1):cpt[i + 1]])
        
      }
      
    } else if (object$contrast == "pcwsLinMean") {
      fit <- rep(0, length(object$x))
      
      for (i in 1:(length(cpt) - 1)) {
        y <- object$x[(cpt[i] + 1):cpt[i + 1]]
        
        if (length(y) == 1)
          fit[(cpt[i] + 1):cpt[i + 1]] <- y
        else{
          n <- length(y)
          x <- 1:n
          beta <- cov(x, y) / (1 / 12 * (-1 + n ^ 2))
          alpha <- mean(y) - (n + 1) * beta / 2
          
          
          fit[(cpt[i] + 1):cpt[i + 1]] <- alpha + beta * x
          
          
        }
        
      }
      
      
    } else if (object$contrast == "pcwsQuadMean") {
      fit <- rep(0, length(object$x))
      
      for (i in 1:(length(cpt) - 1)) {
        y <- object$x[(cpt[i] + 1):cpt[i + 1]]
        
        if (length(y) < 2)
          fit[(cpt[i] + 1):cpt[i + 1]] <- mean(y)
        else{
          n <- length(y)
          x <- cbind(rep(1, n), 1:n, (1:n) ^ 2)
          
          
          
          fit[(cpt[i] + 1):cpt[i + 1]] <- lm.fit(x, y)$fitted.values
          
          
        }
        
      }
    } else if (object$contrast == "pcwsLinContMean") {
      fit <- rep(0, length(object$x))
      cpt <- setdiff(cpt, c(0, length(object$x)))
      X <-
        bs(
          1:length(object$x),
          knots = cpt,
          degree = 1,
          intercept = TRUE
        )
      
      fit <- lm.fit(X, object$x)$fitted.values
      
    } else if (object$contrast == "pcwsConstMeanVar") {
      fit <- matrix(0, nrow = length(object$x), ncol = 2)
      colnames(fit) <- c("mean", "volatility")
      
      for (i in 1:(length(cpt) - 1)) {
        y <- object$x[(cpt[i] + 1):cpt[i + 1]]
        
        fit[(cpt[i] + 1):cpt[i + 1], 1] <- mean(y)
        fit[(cpt[i] + 1):cpt[i + 1], 2] <- sqrt(mean((y - mean(y)) ^ 2))
        
        
      }
      
    }
    
    
    return(fit)
    
}