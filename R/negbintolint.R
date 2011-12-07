negbintol.int <- function (x, n = NULL, N, m = NULL, alpha = 0.05, P = 0.99, side = 1) 
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (side == 2) {
        alpha <- alpha/2
        P <- (P + 1)/2
    }
	if(is.null(n)) n <- length(x)
	if(n==1) warning("MLE is calculated assuming n=1.  Specify a value for n.", call. = FALSE)
	if(length(x)>1) x <- sum(x)
	if(is.null(m)) m <- N
	  mu.hat <- x/n
	  nu.hat <- N/(N+mu.hat)
	  mu.hat <- m*(1/nu.hat-1)

	  se.mu.hat <- sqrt((mu.hat+(mu.hat^2)/m)/n)

        lower.mu <- max(0,mu.hat-qnorm(1-alpha)*se.mu.hat)
        upper.mu <- mu.hat+qnorm(1-alpha)*se.mu.hat


    f1 <- function(J, N, mu.hat, P) pnbinom((J - 1), size = N, mu = mu.hat, lower.tail = FALSE) - P
    f2 <- function(J, N, mu.hat, P) pnbinom(J, size = N, mu = mu.hat) - P
    lower <- suppressWarnings(try(floor(uniroot(f1, interval = c(0, 1e+101), N = m, mu.hat = lower.mu, P = P)$root), 
		silent = TRUE))
    upper <- suppressWarnings(try(floor(uniroot(f2, interval = c(0, 1e+101), N = m, mu.hat = upper.mu, P = P)$root), 
		silent = TRUE))
    if (class(lower) == "try-error") {
        lower <- 0
    }
    else {
        J1.temp <- lower + c(-15:15)
        J1 <- cbind(J1.temp, suppressWarnings(f1(J1.temp, m, lower.mu, P)))
        J1.ind <- (J1[, 2] >= 0 & J1[, 1] >= 0)
        J1 <- matrix(J1[J1.ind, ], ncol = 2)
        lower <- J1[which.min(J1[, 2]), 1]
    }
    if (class(upper) == "try-error") {
        upper <- 0
    }
    else {
        J2.temp <- upper + c(-15:15)
        J2 <- cbind(J2.temp, suppressWarnings(f2(J2.temp, m, upper.mu, P)))
        J2.ind <- (J2[, 2] >= 0 & J2[, 1] >= 0)
        J2 <- matrix(J2[J2.ind, ], ncol = 2)
        upper <- J2[which.min(J2[, 2]), 1]
    }
    if (side == 2) {
        alpha <- 2 * alpha
        P <- (2 * P) - 1
    }
    temp <- data.frame(cbind(alpha, P, nu.hat, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "nu.hat", "2-sided.lower", 
            "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "nu.hat", "1-sided.lower", 
            "1-sided.upper")
    }
    temp
}
