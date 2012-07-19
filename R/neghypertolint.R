#This gives tolerance limits for the number of future draws you can expect to make to get m successes.

neghypertol.int <- function (x, n, N, m = NULL, alpha = 0.05, P = 0.99, side = 1, 
	method = c("LS", "CC")) 
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (side == 2) {
        alpha <- alpha/2
        P <- (P + 1)/2
    }
	rate <- x/N
	if (rate < 0.05) warning("Sampling rate < 0.05.  Results may not be accurate!", call. = FALSE)
    method <- match.arg(method)
    if (length(x) > 1) 
        x <- sum(x)
    nu.hat <- n/x
    k <- qnorm(1 - alpha)
    if(is.null(m)) m <- n
	fpc <- sqrt((N-x)/(N-1))
    se.nu.hat <- sqrt((nu.hat^2*(1 - nu.hat))/x)*fpc
    if(method=="LS"){
        lower.p <- max(0.0000001, nu.hat - k * se.nu.hat,na.rm=TRUE)
        upper.p <- min(nu.hat + k * se.nu.hat, 1,na.rm=TRUE)
    }
    if (method == "CC") {
        lower.p <- nu.hat - k * se.nu.hat - 1/(2*x)
        upper.p <- nu.hat + k * se.nu.hat + 1/(2*x)
    }
    lower.p <- max(0,lower.p)
    upper.p <- min(upper.p,1)
	Mu <- min(ceiling(N * upper.p), N)
	Ml <- max(0, floor(N * lower.p))
    lower <- qnhyper(1-P, m = Mu, n = N, k = m)
    upper <- qnhyper(P, m = Ml, n = N, k = m)    
	if (side == 2) {
        alpha <- 2 * alpha
        P <- (2 * P) - 1
    }
    temp <- data.frame(cbind(alpha, P, rate, nu.hat, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "rate", "p.hat", "2-sided.lower", 
            "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "rate", "p.hat", "1-sided.lower", 
            "1-sided.upper")
    }
    temp
}



