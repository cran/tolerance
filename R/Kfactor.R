K.factor <- function (n, f = NULL, alpha = 0.05, P = 0.99, side = 1, method = c("HE", 
    "WBE", "ELL")) 
{
    if (is.null(f)) 
        f <- n - 1
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    method <- match.arg(method)
    if (side == 1) {
        z.p <- qnorm(P)
        ncp <- sqrt(n) * z.p
        t.a <- suppressWarnings(qt(1 - alpha, df = f, ncp = ncp))
        K <- t.a/sqrt(n)
    }
    else {
        chi.a <- qchisq(alpha, f)
        if (method == "HE") {
            z.p <- qnorm((1 + P)/2)
            z.a <- qnorm((2 - alpha)/2)
            df.cut <- n^2 * (1 + 1/z.a^2)
            V <- 1 + z.a^2/n + ((3 - z.p^2) * z.a^4)/(6 * n^2)
            K.1 <- suppressWarnings(z.p * sqrt(V * (1 + (n * V/(2 * 
                f)) * (1 + 1/z.a^2))))
            G <- (f - 2 - chi.a)/(2 * (n + 1)^2)
            K.2 <- suppressWarnings(z.p * sqrt(((f * (1 + 1/n))/(chi.a)) * 
                (1 + G)))
            if (f > df.cut) {
                K <- K.1
            }
            else {
                K <- K.2
                if (is.na(K)) 
                  K <- 0
            }
        }
        else if (method == "WBE") {
            r <- 0.5
            delta <- 1
            while (abs(delta) > 1e-08) {
                P.new <- pnorm(1/sqrt(n) + r) - pnorm(1/sqrt(n) - 
                  r)
                delta <- P.new - P
                diff <- dnorm(1/sqrt(n) + r) + dnorm(1/sqrt(n) - 
                  r)
                r <- r - delta/diff
            }
            K <- r * sqrt(f/chi.a)
        }
        else if (method == "ELL") {
            if (f < (n^2)) 
                warning("The Ellison method should only be used for f appreciably larger than n^2.", 
                  call. = FALSE)
            r <- 0.5
            delta <- 1
            z.p <- qnorm((1 + P)/2)
            while (abs(delta) > 1e-08) {
                P.new <- pnorm(z.p/sqrt(n) + r) - pnorm(z.p/sqrt(n) - 
                  r)
                delta <- P.new - P
                diff <- dnorm(z.p/sqrt(n) + r) + dnorm(z.p/sqrt(n) - 
                  r)
                r <- r - delta/diff
            }
            K <- r * sqrt(f/chi.a)
        }
    }
    K
}
