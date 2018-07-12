qbymc<-function(x){
#library(robustbase)
N<-length(x)
yyy<-cut((x), c(min(x)-1,median(x), max(x)),label<-c(0,1))
mmm<-data.frame(x,yyy)
vmcl<-abs(by(as.double(mmm[,1]),factor(mmm[,2]),mc)[2])
vmcr<-abs(by(as.double(mmm[,1]),factor(mmm[,2]),mc)[1])
vmct<-(vmcl+vmcr)/2
if (vmct > 0.348) j<-c(0.1797145,.38767097,-.00837164)
else j<-c(0.2017775,.28213917,.083140836)
Z<-(-j[2]+sqrt(j[2]^2-4*j[3]*(j[1]-atanh(vmct))))/(2*j[3])
dqdmq<-cosh(j[1]+j[2]*Z+j[3]*Z^2)^2/((j[2]+2*j[3]*Z)*(1+Z)^2)
qv<-(3*Z+1)/(Z+1)
if (qv < 5/3) dm<-exp(.5)/sqrt(N)
else  dm<-exp(.5)/sqrt(N)*.5^(qv-5/3)
dq<-dm*dqdmq*sqrt(2)
attr(qv,"names") <- "Estimate"
attr(dq,"names") <- "Std. Error"
return(c(qv,dq))
}
