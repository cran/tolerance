mvtol.region <- function (x, alpha = 0.05, P = 0.99, B = 1000) 
{
    n <- nrow(x)
    p <- ncol(x)
    q.squared <- matrix(rchisq(p * B, df = 1), ncol = p)/n
    L <- t(sapply(1:B, function(i) eigen(rwishart(n - 1, p))$values))
    c1 <- apply((1 + q.squared)/L, 1, sum)
    c2 <- apply((1 + 2 * q.squared)/L^2, 1, sum)
    c3 <- apply((1 + 3 * q.squared)/L^3, 1, sum)
    a <- (c2^3)/(c3^2)
    T <- matrix(sapply(1:length(P), function(i) (n - 1) * (sqrt(c2/a) * 
        (qchisq(P[i], a) - a) + c1)), ncol = length(P))
    tol <- matrix(sapply(1:length(alpha), function(i) apply(T, 
        2, quantile, 1 - alpha[i])), ncol = length(alpha))
    colnames(tol) <- alpha
    rownames(tol) <- P
    tol
}
