#' @title Plot a 'not' object
#' @description Plots the input vector used to generate 'not' object \code{x} with the signal fitted with \code{\link{predict.not}}.
#' @method plot not
#' @export 
#' @param x An object of class 'not', returned by \code{\link{not}}.
#' @param ... Further parameters which may be passed to \code{\link{predict.not}} and \code{\link{features}}.  
#' @seealso \code{\link{predict.not}} \code{\link{not}}  \code{\link{features}}
#' @examples 
#' # **** Piecewisce-constant mean with Gaussian noise.
#' x <- c(rep(0, 100), rep(1,100)) + rnorm(100)
#' # *** identify potential locations of the change-points
#' w <- not(x, contrast = "pcwsConstMean")
#' # *** when 'cpt' is omitted, 'features' function is used internally 
#' # to choose change-points locations
#' plot(w)
#' # *** estimate and plot the signal specifying the location of the change-point
#' plot(w, cpt=100)

plot.not <- function(x,...){
	
  
  if(x$contrast == "pcwsConstMean" ||
     x$contrast == "pcwsConstMeanHT" ||
     x$contrast == "pcwsLinMean" ||
     x$contrast == "pcwsLinContMean" ||
     x$contrast == "pcwsQuadMean"){
    
    plot(x$x,ylab="x",type="l")
    lines(x=predict(x,...), type="l",col="red")  
    title("Data and the fitted signal")
    
  }else if(x$contrast == "pcwsConstMeanVar"){
    
    fit <- predict(x, ...)
    
    plot(x$x,ylab="x",type="l")
    lines(x=fit[,1], type="l",col="red")  
    title("Data and the fitted mean function")

    readline(prompt="Press [enter] to continue\n")
    
    plot(abs(x$x-fit[,1]),ylab="x",type="l")
    lines(x=fit[,2], type="l",col="red")  
    title("Centered data and the fitted volatility function")
    
  }
 
}

