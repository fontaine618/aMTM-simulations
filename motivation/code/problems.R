#############################################
# Definition of the target density

# dimension
d = 2
# number of mixture components
M = 2
# weight of the components
w.t = c(0.3, 0.7)
# matrix of means
mu.t = matrix(c(20, 0, 0, 8), M, d, T)
# array of covariances
Sigma.t <- array(c(9, 0, 0 ,1 ,1 , 0, 0 ,9), dim = c(d,d,M))


# precomputation to speed up the evaluations
S <- sapply(seq(M), function(i)chol(Sigma.t[,,i]), simplify = 'array')
Sinv <- sapply(seq(M), function(i)solve(Sigma.t[,,i]), simplify = 'array')
denum <- sapply(seq(M), function(i)sqrt(det(Sigma.t[,,i]* 2*3.14159)) )
B <- 10^max(round(log(denum, 10)))
denum <- denum/B
#target function
parms <- list(M=M,mu=t(mu.t),Sinv=Sinv,w=w.t,denum=denum)
p <- function(x, pa){
   x <- as.matrix(x)
   apply(x,1,function(x){
      log(sum(sapply(seq(pa$M), function(i){
         pa$w[i] * exp(-0.5*(x-pa$mu[,i]) %*% pa$Sinv[,,i] %*% (x-pa$mu[,i]))/pa$denum[i]
      })))})
}

# check that it actually works
# p(t(parms$mu),parms)

# precompile for faster evaluation
logp <- compiler::cmpfun(p)

# prepare data and function
fun = function(data, job, N, ...){
   list(target=data$target, parms=data$parms, N=N)
}


# statistics to compute
stats = function(X){
   # expects X a coda::mcmc object
   Pr5 = mean(X[,1]>5)
   stat = aMTM::stats.aMTM(X)
   mean5 = apply(X[X[,1]>5,,drop=FALSE], 2, mean)
   quantiles = quantile(X[X[,1]<=5, 2], c(0.01, 0.99))
   names(quantiles) = c(1, 99)
   accrate = mean(apply(abs(diff(X)), 1, sum) > 0)
   results=c(
      Pr5=mean(X[,1]>5),
      mean_g5.=mean5,
      quantile_l5=quantiles,
      stat,
      accrate=accrate
   )
}