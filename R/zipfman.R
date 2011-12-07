dzipfman <- function(x, s, b, N, log = FALSE){
	if(s <= 0|b < 0){
		stop(paste("Invalid value for s and/or b!", "\n"))
	}
    x <= ceiling(x)
    out <- rep(0, length(x))
    out[(x>0)&(x<=N)] <- (x[(x>0)&(x<=N)]+b)^(-s)/(sum((c(1:N)+b)^(-s)))
    if(log) out <- log(out)
	out
}

pzipfman <- function(q, s, b, N, lower.tail = TRUE, log.p = FALSE){
	if(s <= 0|b < 0){
		stop(paste("Invalid value for s and/or b!", "\n"))
	}
    q <= ceiling(q)
	temp <- rep(0,length(q))
    for(i in 1:length(q)){ 
        if(q[i] <= 0){ 
            temp[i] <- 0 
        } else if(q[i] > N){
            temp[i] <- 1
        } else temp[i] <- sum(dzipfman(1:q[i], s = s, b = b, N = N))
    }
    if (lower.tail == FALSE) temp <- 1 - temp
    if (log.p) temp <- log(temp)
	temp
}

qzipfman <- function(p, s, b, N, lower.tail = TRUE, log.p = FALSE){
	if(s <= 0|b < 0){
		stop(paste("Invalid value for s and/or b!", "\n"))
	}
    if (log.p) p <- exp(p)
    if (lower.tail == FALSE) p <- 1 - p
	all.p <- NULL
	for(i in 1:length(p)){
		if(pzipfman(1, s = s, b = b, N = N) >= p[i]){
			all.p <- c(all.p, 1)
		} else{
			temp <- try(ceiling(uniroot(function(n) 
(pzipfman(n, s = s, b = b, N = N)-p[i]), c(1, N))$root), silent = TRUE)
				if(class(temp)=="try-error"){
					temp <- NaN
				} else if(temp > 2){
						temp.n_1 <- pzipfman(temp-1, s = s, b = b, N = N)
						if(temp.n_1 > p[i]) temp <- temp-1
					}
			all.p <- c(all.p, temp)			
		}
	}
	all.p
}


rzipfman <- function(n, s, b, N){
	if(s <= 0|b < 0){
		stop(paste("Invalid value for s and/or b!", "\n"))
	}
	qzipfman(runif(n), s = s, b = b, N = N)
}




