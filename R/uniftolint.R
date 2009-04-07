uniftol.int=function(x,alpha=0.05,P=0.99){
n=length(x)
x.n=max(x)
lower=x.n*(1-P)/(1-alpha)^(1/n)
upper=x.n*P/(alpha)^(1/n)
temp=data.frame(cbind(alpha,P,lower,upper))
colnames(temp)=c("alpha","P","1-sided.lower","1-sided.upper")
temp}