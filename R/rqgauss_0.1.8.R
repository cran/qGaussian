rqgauss <-function(n,q=0,mu=0,sig=1,meth="Box-Muller"){
nsam <- n
qv <- q
if(qv >= 3 || qv == 1) stop("q value must be < 3 or != 1")
if(meth != "Quantile") if(meth != "Chaotic") 
if(meth != "Box-Muller") stop("invalid method")
qPDF <- array(1:nsam)
qPDF[] <-0
if (meth == "Box-Muller"){ 
logq <-
function(x,qva)
 {
 mde<-.Machine$double.exponent-2
 if(abs(qva-1)<10^-mde){
 a<-log(x)
 }else{
 b<-as.double(x/x^qva/(1-qva))
 a<-as.double(b-1/(1-qva))
 }
 return(a)
 }
 qgen<-as.double((1+qv)/(3-qv))
 for(i in 1:nsam){
 u1<-as.double(runif(1,0,1))
 u2<-as.double(runif(1,0,1))
 qPDF[i]<-mu+sig*as.double(sqrt(abs(2*logq(u1,qgen)))*sin(2*pi*u2))
 }
}else{
if (meth == "Quantile") qPDF<-cqgauss(runif(nsam),qv,mu,sig) 
}
if (meth == "Chaotic") qPDF<-Chaotic(nsam,qv,mu,sig)
return(qPDF)
}
