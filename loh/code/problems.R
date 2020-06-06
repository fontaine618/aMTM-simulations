#############################################
# Generating the problem 
setwd("~/Documents/aMTM-simulations/loh")
df = read.table("data/BarrettsLOH.dat", header=F)
colnames(df) = c("x", "n")


parms = list(
   x=df$x, 
   n=df$n,
   gam_bound=30,
   dbetabinom=rmutil::dbetabinom
)



logp = function(theta, parms){
   K = nrow(theta)
   eta = pmin(pmax(theta[,1], 1.0e-12),1-1.0e-12)
   pi1 = pmin(pmax(theta[,2], 1.0e-12),1-1.0e-12)
   pi2 = pmin(pmax(theta[,3], 1.0e-12),1-1.0e-12)
   gam = theta[,4]
   w = exp(gam) /(2*(1+exp(gam))) 
   cmp1 = sapply(pi1, function(pi) 
      dbinom(parms$x, parms$n, pi)
      )
   cmp2 = mapply(function(pi, s) 
      parms$dbetabinom(parms$x, parms$n, pi, 1/s),
      pi2, w
      )
   eta_mat = matrix(eta, length(parms$x), K, T)
   f = eta_mat * cmp1 + (1-eta_mat) * cmp2
   llk = apply(log(f), 2, sum)
   prior = rep(0, K)
   prior = prior - (eta < 1e-10)
   prior = prior - (eta > 1 - 1e-10)
   prior = prior - (pi1 < 1e-10)
   prior = prior - (pi1 > 1 - 1e-10)
   prior = prior - (pi2 < 1e-10)
   prior = prior - (pi2 > 1 - 1e-10)
   prior = prior - (abs(gam) > parms$gam_bound)
   return(llk + prior*1e5 + 100)
}





# precompile for faster evaluation
target <- compiler::cmpfun(logp)

# prepare data and function
data = list(
   parms = parms,
   target = logp
)

fun = function(data, job, N, ...){
   list(target=data$target, parms=data$parms, N=N)
}


# get pre-sample

grid = expand.grid(
   eta=seq(0.01, 0.99, length.out = 20),
   pi1=seq(0.01, 0.99, length.out = 20),
   pi2=seq(0.01, 0.99, length.out = 20),
   gam=seq(-29.99, 29.99, length.out = 20)
)
evals = logp(as.matrix(grid), parms)

ps = exp(evals)
ids1 = grid$pi1 > 0.4
sum(ps[ids1]) / sum(ps)


# global covariance
ids = evals > -10
grid_good = as.matrix(grid)[ids, ]
sig = round(cov(grid_good), 3)

# covariance mode 1
ids1 = grid$pi1 > 0.4
grid_good = as.matrix(grid)[ids * ids1 == 1, ]
sig1 = round(cov(grid_good), 3)

# covariance mode 2
ids2 = grid$pi1 < 0.4
grid_good = as.matrix(grid)[ids * ids2 == 1, ]
sig2 = round(cov(grid_good), 3)

Sig_oracle = array(c(sig, sig1, sig2), dim=c(4, 4, 3))

Sig_constant = array(c(sig, sig, sig), dim=c(4, 4, 3))

Sig_scale = array(c(sig, sig*0.1, sig*0.01), dim=c(4, 4, 3))

# statistics to compute
stats = function(X){
   N = nrow(X)
   # expects X a coda::mcmc object
   # get two modes ids
   ids1 = X[, 2] > 0.4
   ids2 = !ids1
   # proportion of smaller mode
   prop1 = mean(ids1)
   # parameter estimates
   ests = sapply(seq(4), function(i){
      m = mcmcse::mcse(X[, i])
      q = quantile(X[, i], c(0.025, 0.975), names=FALSE)
      sd = sd(X[, i])
      if(sum(ids1) > 0){
         m1 = mcmcse::mcse(X[ids1, i])
         q1 = quantile(X[ids1, i], c(0.025, 0.975), names=FALSE)
         sd1 = sd(X[ids1, i])
      }else{m1=list(est=NA, se=NA); q1=c(NA, NA); sd1=NA}
      if(sum(ids2) > 0){
         q2 = quantile(X[ids2, i], c(0.025, 0.975), names=FALSE)
         sd2 = sd(X[ids2, i])
         m2 = mcmcse::mcse(X[ids2, i])
      }else{m2=list(est=NA, se=NA); q2=c(NA, NA); sd2=NA}
      c(
         m=m$est, mse=m$se/sqrt(N), sd=sd, qL=q[1], qU=q[2], 
         m1=m1$est, m1se=m1$se/sqrt(N), sd1=sd1, q1L=q1[1], q1U=q1[2], 
         m2=m2$est, m2se=m2$se/sqrt(N), sd2=sd2, q2L=q2[1], q2U=q2[2],
         acf=mean(abs(stats::acf(X[, i], lag.max=40*40, plot=F)$acf))
      )
   })
   # chain statistics
   dif <- diff(X)
   stats_mcmc = c(
      msjd=mean(sqrt(apply(dif^2,1,sum))),
      accrate=mean(apply(abs(dif), 1, sum) > 0),
      prop1=mean(ids1)
   )
   colnames(ests) = c("eta", "pi1", "pi2", "gamma")
   list(stats=stats_mcmc, estimates=ests)
}





