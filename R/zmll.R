zm.ll <- function (x, N=NULL, s = 1, b = 1, dist = c("Zipf", "Zipf-Man", "Zeta"), exact = TRUE, ...) {
	dist <- match.arg(dist)
	if(class(x)!="table"){ 
		x <- table(x)
	}
	if(sum(is.na(suppressWarnings(as.numeric(names(x)))))>0) names(x)=1:length(x)
	x.labs <- as.numeric(names(x))
	N.temp <- max(x.labs)
	if(dist == "Zeta"){
	N <- N.temp
	}
	if(is.null(N)){
		N <- N.temp
	}
	if(N<N.temp) stop(paste("N cannot be smaller than the maximum observed value in x!",
		"\n"))
	N.seq <- 1:N
	if(length(N.seq)!=length(x.labs)){
	x.0 <- N.seq[as.logical(1-(N.seq%in%x.labs))]
	temp.table <- as.table(rep(0,length(x.0)))
	names(temp.table) <- as.character(x.0)
	x <- c(x,temp.table)
	x <- x[order(as.numeric(names(x)))]
	}
	if(exact==FALSE) x <- x[order(x,decreasing=TRUE)]; names(x)=1:N 
	if(dist=="Zipf"){
		ll.zipf <- function(s) sum(x*(s*log(N.seq)+log(sum(1/(N.seq)^s))))
		fit <- suppressWarnings(mle(ll.zipf,start=list(s=s),lower=0,...))
	}
	if(dist=="Zipf-Man"){
		ll.zima <- function(s,b) sum(x*(s*log(N.seq+b)+log(sum(1/(N.seq+b)^s))))
		fit <- suppressWarnings(mle(ll.zima,start=list(s=s,b=b),lower=c(0,0),...))
	}
	if(dist=="Zeta"){
		ll.zeta <- function(s) sum(x*(s*log(N.seq)+log(zeta.fun(s))))
		fit <- suppressWarnings(mle(ll.zeta,start=list(s=s),lower=1+1e-14,...))
	}
	fit
}