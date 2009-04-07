regtol.int <- function (reg, new.x = NULL, side = 1, alpha = 0.05, P = 0.99) 
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (class(reg) != "lm") {
        stop(paste("Input must be of class 'lm'.", "\n"))
    }
    if (is.vector(new.x)) 
        new.x <- matrix(new.x, ncol = 1)
    n <- length(reg$res)
    pars <- length(reg$coef)
    x <- data.matrix(reg$model[, -1], rownames.force = FALSE)
    new.length <- 0
    if (is.null(new.x) == FALSE) {
        new.length <- nrow(new.x)
        x <- rbind(x, new.x)
    }
    x <- data.frame(x)
    y <- c(reg$model[, 1], rep(NA, new.length))
    names(x) <- names(reg$coef[-1])
    est <- predict(reg, newdata = x, se.fit = TRUE)
    y.hat <- est$fit
    se.y <- est$se.fit
    a.out <- anova(reg)
    MSE <- a.out$"Mean Sq"[length(a.out$"Mean Sq")]
    df <- a.out$Df[length(a.out$Df)]
    n.star <- MSE/se.y^2
    if (side == 1) {
        z.p <- qnorm(P)
        delta <- sqrt(n.star) * z.p
        t.delta <- suppressWarnings(qt(1 - alpha, df = n - pars, 
            ncp = delta))
        K <- t.delta/sqrt(n.star)
        upper <- y.hat + sqrt(MSE) * K
        lower <- y.hat - sqrt(MSE) * K
        temp <- data.frame(cbind(alpha, P, y, y.hat, lower, upper))
        colnames(temp) <- c("alpha", "P", "y", "y.hat", "1-sided.lower", 
            "1-sided.upper")
    }
    else {
        z.p <- qnorm((1 + P)/2)
        z.a <- qnorm((2 - alpha)/2)
        df.cut <- n.star^2 * (1 + 1/z.a^2)
        K <- NULL
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
        upper <- y.hat + sqrt(MSE) * K
        lower <- y.hat - sqrt(MSE) * K
        temp <- data.frame(cbind(alpha, P, y, y.hat, lower, upper))
        colnames(temp) <- c("alpha", "P", "y", "y.hat", "2-sided.lower", 
            "2-sided.upper")
    }
    index <- which(names(temp) == "y.hat")
    temp <- data.matrix(temp[order(temp[, index]), ], rownames.force = FALSE)
    temp
}

