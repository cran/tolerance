zipftol.int <- function (x, N = NULL, alpha = 0.05, P = 0.99, side = 1, 
s = 1, b = 1, dist = c("Zipf", "Zipf-Man", "Zeta"), exact = TRUE, ...) {
    dist = match.arg(dist)
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (side == 2) {
        alpha <- alpha/2
        P <- (P + 1)/2
    }

fit <- zm.ll(x = x, N = N, s = s, b = b, dist = dist, exact = exact, ...)

if(class(x)!="table"){ 
x=table(x)
}
if(sum(is.na(suppressWarnings(as.numeric(names(x)))))>0) names(x)=1:length(x)
x.labs=as.numeric(names(x))
N.temp=max(x.labs)
if(is.null(N)){
N=N.temp
}

if(dist=="Zipf"){
s.hat <- as.numeric(coef(fit))
s.se <- sqrt(as.numeric(vcov(fit)))
CI <- s.hat + c(-1,1)*qnorm(1-alpha)*s.se
lower.s <- max(CI[1],0)
upper.s <- CI[2]

    f1 <- function(J, N, P, s.hat) 1 - pzipfman(q = J, N = N, s = s.hat) - P
    f2 <- function(J, N, P, s.hat) pzipfman(q = J, N = N, s = s.hat) - P

lower <- try(floor(uniroot(f1, interval = c(1, N), N = N, 
        P = P, s.hat = upper.s)$root), silent = TRUE)
upper <- try(floor(uniroot(f2, interval = c(1, N), N = N, 
        P = P, s.hat = lower.s)$root), silent = TRUE)
} 

if(dist=="Zipf-Man"){
s.hat <- as.numeric(coef(fit)[1])
s.se <- sqrt(as.numeric(vcov(fit)[1,1]))
b.hat <- as.numeric(coef(fit)[2])
b.se <- sqrt(as.numeric(vcov(fit)[2,2]))
    if (b.hat==0) warning("MLE for b is 0! Consider fitting a Zipf distribution.")
s.CI <- s.hat + c(-1,1)*qnorm(1-alpha)*s.se #These intervals might need to be divided by 2 for joint confidence.
b.CI <- b.hat + c(-1,1)*qnorm(1-alpha)*b.se
lower.s <- max(s.CI[1],0)
upper.s <- s.CI[2]
lower.b <- max(b.CI[1],0)
upper.b <- b.CI[2]
    f3 <- function(J, N, P, s.hat, b.hat) 1 - pzipfman(q = J, N = N, s = s.hat, b = b.hat) - P
    f4 <- function(J, N, P, s.hat, b.hat) pzipfman(q = J, N = N, s = s.hat, b = b.hat) - P

lower <- try(floor(uniroot(f3, interval = c(1, N), N = N, 
        P = P, s.hat = upper.s, b.hat = lower.b)$root), silent = TRUE)
upper <- try(floor(uniroot(f4, interval = c(1, N), N = N, 
        P = P, s.hat = lower.s, b.hat = upper.b)$root), silent = TRUE)
}

if(dist=="Zeta"){
N <- N.temp
s.hat <- as.numeric(coef(fit))
s.se <- sqrt(as.numeric(vcov(fit)))
CI <- s.hat + c(-1,1)*qnorm(1-alpha)*s.se
lower.s <- max(CI[1],0)
upper.s <- CI[2]

    f5 <- function(J, P, s.hat) 1 - pzipfman(q = J, s = s.hat, N = Inf) - P
    f6 <- function(J, P, s.hat) pzipfman(q = J, s = s.hat, N = Inf) - P
lower <- try(floor(uniroot(f5, interval = c(1, 1000*N), 
        P = P, s.hat = upper.s)$root), silent = TRUE)
upper <- try(floor(uniroot(f6, interval = c(1, 1000*N),  
        P = P, s.hat = lower.s)$root), silent = TRUE)
}

    if (class(lower) == "try-error") {
        lower <- 1
    }
    if (class(upper) == "try-error") {
        upper <- N
    }
    if (side == 2) {
        alpha <- 2 * alpha
        P <- (2 * P) - 1
    }
    
    if(dist=="Zipf-Man"){
    temp <- data.frame(cbind(alpha, P, s.hat, b.hat, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "s.hat", "b.hat", "2-sided.lower", 
            "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "s.hat", "b.hat", "1-sided.lower", 
            "1-sided.upper")
    }
    } else{
    temp <- data.frame(cbind(alpha, P, s.hat, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "s.hat", "2-sided.lower", 
            "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "s.hat", "1-sided.lower", 
            "1-sided.upper")
    }
    }
    temp
}