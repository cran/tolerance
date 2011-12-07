dzipf <- function(x, s, N, log = FALSE){
	if(s <= 0){
		stop(paste("Invalid value for s!", "\n"))
	}
    x <= ceiling(x)
    out <- rep(0, length(x))
    out[(x>0)&(x<=N)] <- x[(x>0)&(x<=N)]^(-s)/(sum((1:N)^(-s)))
    if(log) out <- log(out)
	out
}

pzipf <- function(q, s, N, lower.tail = TRUE, log.p = FALSE){
	if(s <= 0){
		stop(paste("Invalid value for s!", "\n"))
	}
    q <= ceiling(q)
	temp <- rep(0,length(q))
    for(i in 1:length(q)){ 
        if(q[i] <= 0){ 
            temp[i] <- 0 
        } else if(q[i] > N){
            temp[i] <- 1
        } else temp[i] <- sum(dzipf(1:q[i], s = s, N = N))
    }
    if (lower.tail == FALSE) temp <- 1 - temp
    if (log.p) temp <- log(temp)
	temp
}

qzipf <- function(p, s, N, lower.tail = TRUE, log.p = FALSE){
	if(s <= 0){
		stop(paste("Invalid value for s!", "\n"))
	}
    if (log.p) p <- exp(p)
    if (lower.tail == FALSE) p <- 1 - p
	all.p <- NULL
	for(i in 1:length(p)){
		if(pzipf(1,s,N) >= p[i]){
			all.p <- c(all.p, 1)
		} else{
			temp <- try(ceiling(uniroot(function(n) 
(pzipf(n, s, N)-p[i]), c(1, N))$root), silent = TRUE)
				if(class(temp)=="try-error"){
					temp <- NaN
				} else if(temp > 1){
						temp.n_1 <- pzipf(temp-1, s = s, N = N)
						if(temp.n_1 > p[i]) temp <- temp-1
					}
			all.p <- c(all.p, temp)			
		}
	}
	all.p
}

rzipf <- function(n, s, N){
	if(s <= 0){
		stop(paste("Invalid value for s!", "\n"))
	}
	qzipf(p = runif(n), s = s, N = N)
}
