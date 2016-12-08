#' Generate random intervals
#'
#' The function generates \code{M} intervals of the length smaller or equal than \code{max.length}, whose endpoints are 
#' are drawn uniformly without replacements from \code{1},\code{2},..., \code{n}. This routine can be
#' used inside \code{\link{not}} function and is typically not called directly by the user.
#' @param n a number of endpoints to choose from
#' @param M a number of intervals to generate
#' @param max.length an integer specifying maximum interval length
#' @param min.length an integer specifying minimum interval length
#' @param ... not in use
#' @return a \code{M} by 2 matrix with start (first column) and end (second column) points of an interval in each row
#' @examples
#' #*** draw 100 intervals with the endpoints in 1,...,100
#' intervals <- random.intervals(50, 100)
#' @export random.intervals
#' @seealso \code{\link{not}}


random.intervals <-	function(n, M, min.length=1, max.length=n, ...) {

	n <- as.integer(n)
	M <- as.integer(M)
	
	if( max.length <= 3 || max.length > n) stop("max.length must satisfy 3 < max.lenght <= n")
	
	intervals <- matrix(0,nrow=M,ncol=2)
	
	intervals[,1] <- sample(n-min.length, M, replace=TRUE)
	intervals[,2] <- sample(n-min.length, M, replace=TRUE)
	
	intervals <- apply(intervals, 1, function(x){
	  a <- min(x)
	  b <- max(x)
	  x[1] <- a
	  x[2] <- b+min.length
	  x
	})

	t(intervals)
	
}
