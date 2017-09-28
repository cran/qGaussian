
#include <Rcpp.h>
//#include <math.h>
//#include <iostream.h>

using namespace Rcpp;
double qnormal_x,qnormal_y,qnormal_z;
double expq(double q, double w){
if(q==1.0){
return(exp(w));
}
else{
return (expl(log(1.0+(1.0-q)*w)/(1.0-q)));
}
}
double lnq(double q, double w){
if(q==1.0){
return(log(w));
}
else{
return ((exp(log(w)*(1.0-q))-1.0)/(1.0-q));
}
}
void setseed_qnormal(double v0, double z0){
qnormal_x = sqrt(1-v0*v0);
qnormal_y = v0;
qnormal_z = z0;
}
double Q8(double w, double v){
return(8*w*v*(((16.0*w*w-24.0)*w*w+10.0)*w*w-1.0));
}
double P8(double w){
return((((128.0*w*w-256.0)*w*w+160.0)*w*w-32.0)*w*w+1.0);
}
double f(double z){
return(1.0-fabs(1.0-1.99999*z));
}
void qnormal(double q){
double qq;
qnormal_y = Q8(qnormal_x,qnormal_y);
qnormal_x = P8(qnormal_x);
qq = (q+1.0)/(3.0-q);
qnormal_z = f(expq(qq,-qnormal_z*qnormal_z*0.5));
qnormal_z = sqrt(-2.0*lnq(qq,qnormal_z));
}  

// [[Rcpp::export]]
NumericVector Chaotic(int n,double q,double v0, double z0){
  double eta;
  double xi;
  NumericVector ux(n);
setseed_qnormal(v0,z0);
for(int i=0;i<n;i++){
qnormal(q);
xi = qnormal_x*qnormal_z;
eta = qnormal_y*qnormal_z;
ux[i]= xi;
}
return ux;
}   
