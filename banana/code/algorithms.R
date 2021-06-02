aMTM.wrapper = function(data, job, instance, K, ...){
   N = instance$N
   d = 5
   x0 = rnorm(d, 0, 1)
   sig0 = array(0, dim=c(d,d,K))
   scales = 10^seq(-1, 1, length.out=K)
   if (K==1) scales = c(1.)
   for(k in seq(K)) sig0 [, , k] = oracle_cov * scales[k]
   t0 = proc.time()
   mcmc = aMTM::aMTM(
      target=instance$target, 
      x0=x0,
      N=N,
      K=K,
      sig0=sig0,
      parms=instance$parms,
      ...
   )
   time = (proc.time() - t0)[[3]]
   out = stats(mcmc$X)
   out$time = time
   out$esscpu = out$ess / out$time
   out$essnbeval = out$ess / (2*K-1)
   return(unlist(out))
}

