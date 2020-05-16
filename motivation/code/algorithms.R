IID.wrapper = function(data, job, instance,  ...){
   N = instance$N
   components = sample.int(instance$parms$M, N, T, instance$parms$w)
   X = matrix(rnorm(N*instance$parms$M), nrow=N, ncol=instance$parms$M)
   for(i in seq(instance$parms$M)){
      X[components == i, ] = X[components == i, ] %*% chol(solve(instance$parms$Sinv[, , i]))
   }
   X = X + instance$parms$mu[components, ]
   X = coda::mcmc(X)
   list(
      results=stats(X)
   )
}

aMTM.wrapper = function(data, job, instance, scales, adapt, ...){
   N = instance$N
   d = nrow(instance$parms$mu)
   K = length(scales)
   x0 = rnorm(ncol(instance$parms$mu), 5, 5)
   sig0 = array(0, dim=c(d,d,K))
   for(k in seq(K)) sig0 [, , k] = diag(d) * scales[[k]]
   mcmc = aMTM::aMTM(
      target=instance$target, 
      x0=x0,
      N=N,
      K=K,
      sig0=sig0,
      adapt=adapt,
      parms=instance$parms,
      burnin=0.5,
      accrate=0.5, global=F, scale=F, local=F, gamma=0.7, weight=-1
   )
   list(
      results=stats(mcmc$X)
   )
}

