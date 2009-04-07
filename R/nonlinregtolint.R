nlregtol.int <- function (formula, xy.data = data.frame(), x.new = NULL, side = 1, 
    alpha = 0.05, P = 0.99, maxiter = 50, ...) 
{
    n <- nrow(xy.data)
    form <- as.formula(formula)
    test.sig <- "try-error"
    while (test.sig == "try-error") {
        out <- try(suppressWarnings(nls(formula = form, data = xy.data, 
            control = list(maxiter = maxiter, warnOnly = TRUE), 
            ...)), silent = TRUE)
        test.sig <- class(try(summary(out)$sigma, silent = TRUE))
        maxiter <- ceiling(maxiter/2)
        if (maxiter <= 1) 
            test.sig <- "quit"
    }
    sigma <- summary(out)$sigma
    beta.hat <- coef(out)
    beta.names <- names(beta.hat)
    temp <- data.frame(matrix(beta.hat, ncol = length(beta.hat)))
    colnames(temp) <- beta.names
    attach(temp)
    pars <- length(beta.hat)
    (fx <- deriv(form, beta.names))
    P.mat <- attr(eval(fx), "gradient")
    PTP <- t(P.mat) %*% P.mat
    PTP2 <- try(solve(PTP), silent = TRUE)
    test.PTP <- class(PTP2)
    if (test.PTP == "try-error") {
        PTP0 <- PTP
        while (test.PTP == "try-error") {
            PTP3 <- PTP0 + diag(rep(min(diag(PTP))/1000, length(diag(PTP))))
            PTP.new <- try(solve(PTP3), silent = TRUE)
            test.PTP <- class(PTP > new)
            PTP0 <- PTP3
        }
        PTP <- PTP.new
    }
    else PTP <- PTP2
    if (is.null(x.new) == FALSE) {
        x.temp <- cbind(NA, x.new)
        colnames(x.temp) <- colnames(xy.data)
        xy.data <- rbind(xy.data, x.temp)
        P.mat <- attr(eval(fx, xy.data), "gradient")
    }
    y.hat <- predict(out, newdata = xy.data)
    n.star <- rep(NULL, nrow(xy.data))
    for (i in 1:nrow(xy.data)) {
        n.star[i] <- c(as.numeric((t(P.mat[i, ]) %*% PTP %*% t(t(P.mat[i, 
            ])))))
    }
    n.star <- n.star^(-1)
    detach(temp)
    if (side == 1) {
        z.p <- qnorm(P)
        delta <- sqrt(n.star) * z.p
        t.delta <- suppressWarnings(qt(1 - alpha, df = n - pars, 
            ncp = delta))
        t.delta[is.na(t.delta)] <- Inf
        K <- t.delta/sqrt(n.star)
        K[is.na(K)] <- Inf
        upper <- y.hat + sigma * K
        lower <- y.hat - sigma * K
        temp <- data.frame(cbind(alpha, P, y.hat, xy.data[, 1], 
            lower, upper))
        colnames(temp) <- c("alpha", "P", "y.hat", "y", "1-sided.lower", 
            "1-sided.upper")
    }
    else {
        z.p <- qnorm((1 + P)/2)
        z.a <- qnorm((2 - alpha)/2)
        df.cut <- n.star^2 * (1 + 1/z.a^2)
        K <- NULL
        df <- n - pars
        for (i in 1:length(n.star)) {
            V <- 1 + z.a^2/n.star[i] + ((3 - z.p^2) * z.a^4)/(6 * 
                n.star[i]^2)
            K.1 <- suppressWarnings(z.p * sqrt(V * (1 + (n.star[i] * 
                V/(2 * df)) * (1 + 1/z.a^2))))
            chi.a <- qchisq(alpha, df = df)
            K.2 <- suppressWarnings(z.p * sqrt(((df * (1 + 1/n.star[i]))/(chi.a)) * 
                (1 + (df - 2 - chi.a)/(2 * (n.star[i] + 1)^2))))
            if (df > df.cut[i]) {
                K[i] <- K.1
            }
            else {
                K[i] <- K.2
                if (is.na(K[i])) 
                  K[i] <- 0
            }
        }
        upper <- y.hat + sigma * K
        lower <- y.hat - sigma * K
        temp <- data.frame(cbind(alpha, P, y.hat, xy.data[, 1], 
            lower, upper))
        colnames(temp) <- c("alpha", "P", "y.hat", "y", "2-sided.lower", 
            "2-sided.upper")
    }
    index <- which(names(temp) == "y")
    temp <- data.matrix(temp[order(temp[, index]), ], rownames.force = FALSE)
    temp
}
