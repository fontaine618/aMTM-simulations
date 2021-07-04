aMTM.wrapper = function(data, job, instance, K, gamma, adapt, ...){
   N = instance$N
   d = 5
   x0 = rnorm(d, 0, 1)
   sig0 = array(0, dim=c(d,d,K))
   scales = 10^seq(-2, 2, length.out=K)
   if (K==1) scales = c(1.)
   for(k in seq(K)) sig0 [, , k] = diag(c(1, 3, 1, 1, 3)) * scales[k]
   gamma = ifelse(adapt==3, 0.5, 0.7)
   mcmc = aMTM::aMTM(
      target=instance$target, 
      x0=x0,
      N=N,
      K=K,
      sig0=sig0,
      parms=instance$parms,
      gamma=gamma,
      adapt=adapt,
      ...
   )
   time = mcmc$time
   out = stats(mcmc$X)
   out$time = time
   out$esscpu = out$ess / out$time
   out$essnbeval = out$ess / (2*K-1)
   out$msjdcpu = out$msjd / out$time
   out$msjdnbeval = out$msjd / (2*K-1)
   return(unlist(out))
}

