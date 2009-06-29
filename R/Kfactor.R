K.factor=function(n,alpha=0.05,P=0.99,side=1,method=c("HE","WBE")){
if(side!=1 && side!=2){
	stop(paste("Must specify a one-sided or two-sided procedure!","\n"))
	}
method=match.arg(method)

if(side==1){
z.p=qnorm(P)
ncp=sqrt(n)*z.p
t.a=suppressWarnings(qt(1-alpha,df=n-1,ncp=ncp))
K=t.a/sqrt(n)
} else{
chi.a=qchisq(alpha,n-1)
if(method=="HE"){
u=qnorm((1+P)/2)*sqrt(1+1/n)
v=sqrt((n-1)/(chi.a))
w=sqrt(1+(n-3-chi.a)/(2*(n+1)^2))
K=u*v*w
}
else if(method=="WBE"){
r=.5
delta=1
while(abs(delta)>1e-8){
P.new=pnorm(1/sqrt(n)+r)-pnorm(1/sqrt(n)-r)
delta=P.new-P
diff=dnorm(1/sqrt(n)+r)+dnorm(1/sqrt(n)-r)
r=r-delta/diff
}
K=r*sqrt((n-1)/chi.a)
}
}
K
}
