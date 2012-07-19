dnhyper <- function(x, m, n, k, log = FALSE){
    if (k > m) {
        stop(paste("k cannot be larger than m!", 
            "\n"))
    }
	p <- choose(x - 1, k - 1) * choose(n - x, m - k) / choose(n, m)
    if (log) 
        p <- log(p)
    p
}

pnhyper <- function(q, m, n, k, lower.tail = TRUE, log.p = FALSE){
    if (k > m) {
        stop(paste("k cannot be larger than m!", 
            "\n"))
    }
    q <= ceiling(q)
    temp <- rep(0, length(q))
    for (i in 1:length(q)) {
        if (q[i] < k) {
            temp[i] <- 0
        }
        else if (q[i] >= (n - m + k)) {
            temp[i] <- 1
        }
        else temp[i] <- sum(dnhyper(k:q[i], m = m, n = n, k = k))
    }
    if (lower.tail == FALSE) 
        temp <- 1 - temp
    if (log.p) 
        temp <- log(temp)
    temp
}

qnhyper <- function(p, m, n, k, lower.tail = TRUE, log.p = FALSE){
    if (k > m) {
        stop(paste("k cannot be larger than m!", 
            "\n"))
    }
    if (log.p) p <- exp(p)
    all.p <- NULL
	temp <- (k - 1):(n - m + k + 1) 
	if(lower.tail){
		temp.out <- cbind(temp, pnhyper(temp, m, n, k))
	} else temp.out <- cbind(temp, pnhyper(temp, m, n, k, lower.tail = FALSE))
    for (i in 1:length(p)) {
	if(lower.tail){
		all.p <- c(all.p, min(temp.out[which(temp.out[,2]>=p[i]),1]))
	} else all.p <- c(all.p, min(temp.out[which(temp.out[,2]<p[i]),1]))
   	}
    all.p
}

rnhyper <- function(nn, m, n, k){
    if (k > m) {
        stop(paste("k cannot be larger than m!", 
            "\n"))
    }
	qnhyper(runif(nn), m, n, k)
}
