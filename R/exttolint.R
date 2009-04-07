exttol.int<-function(x,alpha=0.05,P=0.99,dist=c("Weibull","Gumbel"),NR.delta=1e-8){
n=length(x)
dist=match.arg(dist)
if(dist=="Weibull") x=log(x)
delta=sqrt((mean(x^2)-mean(x)^2)*6/pi^2)
xi=mean(x)+digamma(1)*delta
theta.old=c(xi,delta)
#NR to find MLE's;MOM estimates used as starting values
diff=1
while(sum(diff>NR.delta)>0){
f=sum(x*exp(x/delta))
f.1=-sum(x^2*exp(x/delta))/(delta^2)
g=sum(exp(x/delta))
g.1=-f/(delta^2)
d=delta+mean(x)-(f/g)
d.1=1-(g*f.1-f*g.1)/(g^2)
delta.new=delta-d/d.1
xi.new=-delta.new*log(n/sum(exp(x/delta.new)))
delta.old=delta
xi.old=xi
delta=delta.new
xi=xi.new
if(is.na(xi) | is.na(delta) | delta<0){
xi=theta.old[1]
delta=theta.old[2]
diff=NR.delta/5
} else diff=c(abs(delta.new-delta.old),abs(xi.new-xi.old))
}
theta=c(delta,xi)
if(dist=="Weibull"){
x=exp(x)
theta.new=c(1/theta[1],exp(theta[2]))
} else theta.new=theta
lambda=function(P) log(-log(P))
k.t=function(x1,x2) suppressWarnings(qt(1-x1,df=(n-1),ncp=(-sqrt(n)*lambda(x2))))
lower=xi-delta*k.t(alpha,P)/sqrt(n-1)
upper=xi-delta*k.t(1-alpha,1-P)/sqrt(n-1)
if(dist=="Weibull"){
lower=exp(lower)
upper=exp(upper)
}
temp=data.frame(cbind(alpha,P,theta.new[1],theta.new[2],lower,upper))
colnames(temp)=c("alpha","P","shape.1","shape.2","1-sided.lower","1-sided.upper")
temp
}

