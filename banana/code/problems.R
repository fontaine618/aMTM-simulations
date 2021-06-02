#############################################
# Definition of the target density

# dimension
d = 5
# lo-gdensity
logpi = function(x, p) apply(x,1,function(x) {
   -0.5 * (
      x[1]^2/(2*p$a1^2) +
         (p$B1*x[1]^2+x[2])^2 +
         x[3]^2 +
         x[4]^2/(2*p$a2^2) +
         (p$B2*x[4]^2+x[5])^2
   )
})
parms = list(a1=1, a2=1, B1=3, B2=1)

# precompile for faster evaluation
logp <- compiler::cmpfun(logpi)

# iid sample for TV
iid = function(n, p){
   x = matrix(rnorm(n*5), n, 5)
   x[, 1] = x[, 1] * p$a1
   x[, 2] = x[, 2] + p$B1 * (x[, 1]^2)
   x[, 4] = x[, 4] * p$a2
   x[, 5] = x[, 5] + p$B2 * (x[, 4]^2)
   return(x)
}
set.seed(1)
sample = iid(1000000, p)

# statistics to compute
regions = function(x){
   reg1 = (x[, 4]>=0) * (x[, 5]>=5)
   reg2 = (x[, 4]<0) * (x[, 5]>=5)
   reg3 = (x[, 1]>=0) * (x[, 2]>=10)
   reg4 = (x[, 1]<0) * (x[, 2]>=10)
   reg5 = (1-reg1)*(1-reg2)*(1-reg3)*(1-reg4)
   regions = rep(1, nrow(x))
   regions[reg1==1] = 2
   regions[reg2==1] = 3
   regions[reg3==1] = 4
   regions[reg4==1] = 5
   return(list(regions=regions, weights=c(mean(reg1), mean(reg2), mean(reg3), mean(reg4), mean(reg5))))
}
out = regions(sample)
iid_weights = out$weights
oracle_cov = cov(sample)

stats = function(X){
   # expects X a coda::mcmc object
   stat = aMTM::stats.aMTM(X)
   ess = stat[["ess"]]
   out = regions(X)
   weights = out$weights
   tv = sum(abs(weights-iid_weights))
   results=list(
      ess=ess,
      tv=tv
   )
}

# prepare data and function
data = list(
   parms = parms,
   target = logp
)
fun = function(data, job, N, ...){
   list(target=data$target, parms=data$parms, N=N)
}