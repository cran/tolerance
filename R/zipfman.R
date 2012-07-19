dzipfman <- function (x, s, b = NULL, N = NULL, log = FALSE) 
{
    if((is.null(b)) & (N < Inf)){
    if (s <= 0) {
        stop(paste("Invalid value for s!", "\n"))
    }
    x <- ceiling(x)
    out <- rep(0, length(x))
    out[(x > 0) & (x <= N)] <- x[(x > 0) & (x <= N)]^(-s)/(sum((1:N)^(-s)))
    } else if(!is.null(b)){
    if (s <= 0 | b < 0) {
        stop(paste("Invalid value for s and/or b!", "\n"))
    }
    if (N == Inf) {
        stop(paste("N must be finite!", "\n"))
    }
    x <- ceiling(x)
    out <- rep(0, length(x))
    out[(x > 0) & (x <= N)] <- (x[(x > 0) & (x <= N)] + b)^(-s)/(sum((c(1:N) + 
        b)^(-s)))
    } else{
    if (s <= 1) {
        stop(paste("Invalid value for s!", "\n"))
    }
    x <- ceiling(x)
    out <- rep(0, length(x))
    out[x > 0] <- x[x > 0]^(-s)/zeta.fun(s)
    }
    if (log) 
        out <- log(out)
    out
}


rzipfman <- function (n, s, b = NULL, N = NULL) 
{
    if((is.null(b)) & (N < Inf)){
    if (s <= 0) {
        stop(paste("Invalid value for s!", "\n"))
    }
    out <- qzipfman(p = runif(n), s = s, N = N)
    } else if(!is.null(b)){
    if (s <= 0 | b < 0) {
        stop(paste("Invalid value for s and/or b!", "\n"))
    }
    if (N == Inf) {
        stop(paste("N must be finite!", "\n"))
    }
    out <- qzipfman(runif(n), s = s, b = b, N = N)
    } else{
    if (s <= 1) {
        stop(paste("Invalid value for s!", "\n"))
    }
    out <- qzipfman(runif(n), s = s, N = Inf)
    }
    out
}


qzipfman <- function (p, s, b = NULL, N = NULL, lower.tail = TRUE, log.p = FALSE) 
{
    if (log.p) 
        p <- exp(p)
    if (lower.tail == FALSE) 
        p <- 1 - p
    if((is.null(b)) & (N < Inf)){
    if (s <= 0) {
        stop(paste("Invalid value for s!", "\n"))
    }
    all.p <- NULL
    for (i in 1:length(p)) {
        if (pzipfman(q = 1, s = s, N = N) >= p[i]) {
            all.p <- c(all.p, 1)
        }
        else {
            temp <- try(ceiling(uniroot(function(n) (pzipfman(q = n, 
                s = s, N = N) - p[i]), c(1, N))$root), silent = TRUE)
            if (class(temp) == "try-error") {
                temp <- NaN
            }
            else if (temp > 1) {
                temp.n_1 <- pzipfman(q = temp - 1, s = s, N = N)
                if (temp.n_1 > p[i]) 
                  temp <- temp - 1
            }
            all.p <- c(all.p, temp)
        }
    }
    } else if(!is.null(b)){
    if (s <= 0 | b < 0) {
        stop(paste("Invalid value for s and/or b!", "\n"))
    }
    if (N == Inf) {
        stop(paste("N must be finite!", "\n"))
    }
    all.p <- NULL
    for (i in 1:length(p)) {
        if (pzipfman(q = 1, s = s, b = b, N = N) >= p[i]) {
            all.p <- c(all.p, 1)
        }
        else {
            temp <- try(ceiling(uniroot(function(n) (pzipfman(q = n, 
                s = s, b = b, N = N) - p[i]), c(1, N))$root), 
                silent = TRUE)
            if (class(temp) == "try-error") {
                temp <- NaN
            }
            else if (temp > 2) {
                temp.n_1 <- pzipfman(q = temp - 1, s = s, b = b, 
                  N = N)
                if (temp.n_1 > p[i]) 
                  temp <- temp - 1
            }
            all.p <- c(all.p, temp)
        }
    }
    } else{
    if (s <= 1) {
        stop(paste("Invalid value for s!", "\n"))
    }
    all.p <- rep(NA, length(p))
    temp <- (1/zeta.fun(s)) * cumsum((1:1e+06)^(-s))
    temp.1 <- (1/zeta.fun(s))
    temp.max <- max(temp)
    TEMP1 <- which(p <= temp.1)
    if (sum(TEMP1) != 0) 
        all.p[TEMP1] <- 1
    TEMP2 <- which(p > temp.max)
    if (sum(TEMP2) != 0) 
        all.p[TEMP2] <- 1e+06
    TEMP3 <- which(p > temp.1 & p <= temp.max)
    if (sum(TEMP3) != 0) 
        all.p[TEMP3] <- sapply(1:length(TEMP3), function(i) min(which(temp >= 
            p[TEMP3[i]])))
    }
    all.p
}

pzipfman <- function (q, s, b = NULL, N = NULL, lower.tail = TRUE, log.p = FALSE) 
{
    q <- ceiling(q)
    temp <- rep(0, length(q))
    if((is.null(b)) & (N < Inf)){
    if (s <= 0) {
        stop(paste("Invalid value for s!", "\n"))
    }
    for (i in 1:length(q)) {
        if (q[i] <= 0) {
            temp[i] <- 0
        }
        else if (q[i] > N) {
            temp[i] <- 1
        }
        else temp[i] <- sum(dzipfman(x = 1:q[i], s = s, N = N))
    }
    } else if(!is.null(b)){
    if (s <= 0 | b < 0) {
        stop(paste("Invalid value for s and/or b!", "\n"))
    }
    if (N == Inf) {
        stop(paste("N must be finite!", "\n"))
    }
    for (i in 1:length(q)) {
        if (q[i] <= 0) {
            temp[i] <- 0
        }
        else if (q[i] > N) {
            temp[i] <- 1
        }
        else temp[i] <- sum(dzipfman(x = 1:q[i], s = s, b = b, N = N))
    }
    } else{
    if (s <= 1) {
        stop(paste("Invalid value for s!", "\n"))
    }
    for (i in 1:length(q)) {
        if (q[i] <= 0) {
            temp[i] <- 0
        }
        else if (q[i] == Inf) {
            temp[i] <- 1
        }
        else temp[i] <- sum(dzipfman(x = 1:q[i], s = s, N = Inf))
    }
    }
    if (lower.tail == FALSE) 
        temp <- 1 - temp
    if (log.p) 
        temp <- log(temp)
    temp
}

