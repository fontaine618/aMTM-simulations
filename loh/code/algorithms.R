aMTM.wrapper = function(data, job, instance, K, init, ...){
   N = instance$N
   d = 4
   x0 = c(runif(3), rnorm(1,0,5))
   
   if(init=="oracle") sig0 = Sig_oracle
   if(init=="constant") sig0 = Sig_constant
   if(init=="scale") sig0 = Sig_scale
   
   if(K<3) sig0 = sig0[,,seq(K), drop=FALSE] # take K first if K>3
   if(K>3) K=3
      
   mcmc = aMTM::aMTM(
      target=instance$target, 
      x0=x0,
      N=N,
      K=K,
      sig0=sig0,
      parms=instance$parms,
      burnin=1/11,
      ...
   )
   
   stats(mcmc$X)
}


