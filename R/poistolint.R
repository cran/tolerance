poistol.int <- function (x, n, m, alpha = 0.05, P = 0.99, side = 1, method = c("TAB", 
    "LS", "SC")) 
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    method <- match.arg(method)
    if(length(x) > 1) x <- sum(x)
    if (side == 2) {
        alpha <- alpha/2
        P <- (P + 1)/2
    }
    if (method == "TAB") {
        lower.lambda <- 0.5 * qchisq(alpha, df = (2 * x))/n
        upper.lambda <- 0.5 * qchisq(1 - alpha, df = (2 * x + 
            2))/n
    }
    if (method == "LS") {
        lower.lambda <- (x/n) - (qnorm(1 - alpha) * sqrt(x))/n
        upper.lambda <- (x/n) + (qnorm(1 - alpha) * sqrt(x))/n
    }
    if (method == "SC") {
	  k <- qnorm(1 - alpha)
        lower.lambda <- (x/n) + (k^2/(2*n)) - (k/sqrt(n))*sqrt((x/n)+(k^2/(4*n)))
        upper.lambda <- (x/n) + (k^2/(2*n)) + (k/sqrt(n))*sqrt((x/n)+(k^2/(4*n)))
    }
    f1 <- function(J, m, P, lambda1) ppois((J - 1), lambda = (m * 
        lambda1), lower.tail = FALSE) - P
    f2 <- function(J, m, P, lambda1) ppois(J, lambda = (m * lambda1)) - 
        P
    lower <- try(floor(uniroot(f1, interval = c(0, 1e+101), m = m, 
        P = P, lambda1 = lower.lambda)$root), silent = TRUE)
    upper <- try(floor(uniroot(f2, interval = c(0, 1e+101), m = m, 
        P = P, lambda1 = upper.lambda)$root), silent = TRUE)
    if (class(lower) == "try-error") {
        lower <- 0
    }
    else {
        J1.temp <- lower + c(-15:15)
        J1 <- cbind(J1.temp, f1(J1.temp, m, P, lower.lambda))
        J1.ind <- (J1[, 2] >= 0 & J1[, 1] >= 0)
        J1 <- matrix(J1[J1.ind, ], ncol = 2)
        lower <- J1[which.min(J1[, 2]), 1]
    }
    if (class(upper) == "try-error") {
        upper <- 0
    }
    else {
        J2.temp <- upper + c(-15:15)
        J2 <- cbind(J2.temp, f2(J2.temp, m, P, upper.lambda))
        J2.ind <- (J2[, 2] >= 0 & J2[, 1] >= 0)
        J2 <- matrix(J2[J2.ind, ], ncol = 2)
        upper <- J2[which.min(J2[, 2]), 1]
    }
    if (side == 2) {
        alpha <- 2 * alpha
        P <- (2 * P) - 1
    }
    temp <- data.frame(cbind(alpha, P, x/n, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "lambda.hat", "2-sided.lower", 
            "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "lambda.hat", "1-sided.lower", 
            "1-sided.upper")
    }
    temp
}