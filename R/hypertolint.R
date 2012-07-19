hypertol.int <- function (x, n, N, m = NULL, alpha = 0.05, P = 0.99, side = 1, 
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
	rate <- n/N
	if (rate < 0.05) warning("Sampling rate < 0.05.  Results may not be accurate!", call. = FALSE)
    method <- match.arg(method)
    if (length(x) > 1) 
        x <- sum(x)
    p.hat <- x/n
    k <- qnorm(1 - alpha)
    if(is.null(m)) m <- n
	fpc <- sqrt((N-n)/(N-1))
    if (method == "LS") {
        lower.p <- p.hat - k * sqrt(p.hat * (1 - p.hat)/n) * fpc
        upper.p <- p.hat + k * sqrt(p.hat * (1 - p.hat)/n) * fpc
    }
    if (method == "CC") {
        lower.p <- p.hat - k * (sqrt(p.hat * (1 - p.hat)/n) * fpc) - 1/(2*n)
        upper.p <- p.hat + k * (sqrt(p.hat * (1 - p.hat)/n) * fpc) + 1/(2*n)
    }
    lower.p <- max(0,lower.p)
    upper.p <- min(upper.p,1)
	Mu <- min(ceiling(N * upper.p), N)
	Ml <- max(0, floor(N * lower.p))
    lower <- qhyper(1-P, m = Ml, n = N - Mu, k = m)
    upper <- qhyper(P, m = Mu, n = N - Ml, k = m)    
	if (side == 2) {
        alpha <- 2 * alpha
        P <- (2 * P) - 1
    }
    temp <- data.frame(cbind(alpha, P, rate, p.hat, lower, upper))
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



