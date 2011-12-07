dzeta <- function(x, s, log = FALSE){
	if(s <= 1){
		stop(paste("Invalid value for s!", "\n"))
	}
    x <= ceiling(x)
    out <- rep(0, length(x))
    out[x>0] <- x[x>0]^(-s)/zeta.fun(s)
    if(log) out <- log(out)
	out
}

pzeta <- function(q, s, lower.tail = TRUE, log.p = FALSE){
	if(s <= 1){
		stop(paste("Invalid value for s!", "\n"))
	}
    q <= ceiling(q)
	temp <- rep(0,length(q))
    for(i in 1:length(q)){ 
        if(q[i] <= 0){
        	temp[i] <- 0
        } else if(q[i] == Inf){
							temp[i] <- 1
        				} else temp[i] <- sum(dzeta(1:q[i], s = s))
    }
    if (lower.tail == FALSE) temp <- 1 - temp
    if (log.p) temp <- log(temp)
	temp
}

qzeta <- function(p, s, lower.tail = TRUE, log.p = FALSE){
	if(s <= 1){
		stop(paste("Invalid value for s!", "\n"))
	}
    if (log.p) p <- exp(p)
    if (lower.tail == FALSE) p <- 1 - p
	all.p <- rep(NA,length(p))
	temp <- (1/zeta.fun(s))*cumsum((1:1e6)^(-s))
	temp.1 <- (1/zeta.fun(s))
	temp.max <- max(temp)
	TEMP1 <- which(p<=temp.1)
	if(sum(TEMP1)!=0) all.p[TEMP1] <- 1
	TEMP2 <- which(p>temp.max)
	if(sum(TEMP2)!=0) all.p[TEMP2] <- 1e6
	TEMP3 <- which(p>temp.1&p<=temp.max)
	if(sum(TEMP3)!=0) all.p[TEMP3] <- sapply(1:length(TEMP3),function(i) min(which(temp>=p[TEMP3[i]])))
all.p
}

rzeta <- function(n, s){
	if(s <= 1){
		stop(paste("Invalid value for s!", "\n"))
	}
	qzeta(runif(n), s = s)
}





