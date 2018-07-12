#' @title The q-gaussian Distribution
#' @name pqgauss
#'
#' @description Density, distribution function, quantile function and 
#'   random generation for the q-gaussian distribution with mu equal to 
#'   q-mean and sig equal to q-variance. 
#'
#' @param x vector of quantiles. 
#' @param p vector of probabilities.
#' @param q entropic index.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param mu a value for q-mean.
#' @param sig a value for q-variance. 
pqgauss <- function(x,q=0,mu=0,sig=1, lower.tail = TRUE){
qv <- q
Z<-(qv-1)/(3-qv)
if(qv >= 3 || qv == 1) {
stop("q value must be < 3 or != 1)")
}else{
if (qv < 1){
supp<-sig/sqrt(-Z)
if (max(x-mu)> supp || min(x-mu) < -supp){
stop("|x-mu| > Support")
}}}
#library(zipfR) # library de Rbeta
X<-1+Z/(sig^2)*(x-mu)^2
if (qv < 1){
xl<-X
a <- 1-1/Z 
}else{
a <- 1/Z
xl<-1/X
}
P <- array(seq_len(length(xl)))
for(i in seq_len(length(xl))){
if (lower.tail == "TRUE"){
if(x[i]<mu){
P[i] <-Rbeta(xl[i],a/2,1/2)/2
}else{
P[i] <-1-Rbeta(xl[i],a/2,1/2)/2
}
}else{
if (lower.tail == "FALSE"){
if(x[i]>mu){
P[i] <-Rbeta(xl[i],a/2,1/2)/2
}else{
P[i] <-1-Rbeta(xl[i],a/2,1/2)/2
}
}
}

}
return(P)
}