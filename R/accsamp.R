acc.samp=function(n,N,alpha=0.05,P=0.99,AQL){
RQL=1-P
if(RQL-AQL<4e-8){stop(paste("RQL must be greater than AQL!"))}
D=floor((1-P)*N)
m.h=D
n.h=N-D
ff=function(k,c,m,n) (alpha)-phyper(c,m,n,k)
c=try(floor(uniroot(ff,interval=c(0,D),k=n,m=m.h,n=n.h)$root),silent=TRUE)
if(class(c)=="try-error") c=0
if(phyper(c,m=m.h,n=n.h,k=n)>alpha) c=max(c-1,0)
prob.rej.good=1-phyper(c,floor(AQL*N),N-floor(AQL*N),n)
prob.rej.bad=1-phyper(c,m=m.h,n=n.h,k=n)
temp=c(round(c,0),round(N,0),round(RQL,4),1-round(alpha,4),round(AQL,4),round(n,0),round(prob.rej.good,4),round((1-prob.rej.bad),4))
temp=as.matrix(temp)
rownames(temp)=c("acceptance.limit","lot.size","RQL","confidence","AQL","sample.size","prod.risk","cons.risk")
colnames(temp)=""
print(temp)
if(temp[8,]>alpha) cat("Warning: Desired confidence level not attained!","\n")
}


