#############################################
# Definition of the target density

# dimension
d = 5
# log-density
logpi = function(x, p) apply(x,1,function(x) {
   for(i in seq(20)){ # artificially increase computing time
      eval = -0.5 * (
         x[1]^2 / p$a1^2 +
         ( x[2] - p$B1 * x[1]^2 )^2 +
         x[3]^2 +
         x[4]^2 / p$a2^2 +
         ( x[5] - p$B2 * x[4]^2 )^2
      )
   }
   return(eval)
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
sample = iid(1000000, parms)

# statistics to compute
regions = function(x){
   reg1 = (x[, 4]>=0) * (x[, 5]>=3)
   reg2 = (x[, 4]<0) * (x[, 5]>=3)
   reg3 = (x[, 1]>=0) * (x[, 2]>=3)
   reg4 = (x[, 1]<0) * (x[, 2]>=3)
   reg5 = (1-reg1)*(1-reg2)*(1-reg3)*(1-reg4)
   regions = rep(5, nrow(x))
   regions[reg1==1] = 1
   regions[reg2==1] = 2
   regions[reg3==1] = 3
   regions[reg4==1] = 4
   return(list(regions=regions, weights=table(regions)/nrow(x)))
}
out = regions(sample)
iid_weights = out$weights
oracle_cov = cov(sample)

stats = function(X){
   # expects X a coda::mcmc object
   stat = aMTM::stats.aMTM(X, cov=oracle_cov)
   ess = stat[["ess"]]
   msjd = stat[["msjd"]]
   weights = regions(as.matrix(X))$weights
   tv = sum(abs(weights-iid_weights))
   results=list(
      msjd=msjd,
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