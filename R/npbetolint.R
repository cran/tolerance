npbetol.int <- function(x, Beta = 0.95, side = 1, upper = NULL, lower = NULL){
  n <- length(x)
  x <- sort(x)
  n.e <- min(ceiling(Beta*(n+1)),n)
  n.e2 <- max(floor((n-n.e)/2),1)
  if(side==1){
    if(is.null(upper)) upper <- x[n.e]
    if(is.null(lower)) lower <- x[max(n-n.e+1,1)]
  } else{
    if(is.null(upper)) upper <- x[min(n.e+n.e2,n)]
    if(is.null(lower)) lower <- x[n.e2]
  }
  temp <- data.frame(cbind(Beta, lower, upper))
  colnames(temp) <- c("beta", paste(side,"-sided.lower",sep=""), 
                      paste(side,"-sided.upper",sep=""))
  temp
}
