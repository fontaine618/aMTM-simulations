# Setup
setwd("~/Documents/aMTM-simulations/motivation")
source('code/problems.R')
N = 1e4
d = 2
res = list()
target = logp

# IID Sample
set.seed(1)
components = sample.int(parms$M, N, T, parms$w)
X = matrix(rnorm(N*parms$M), nrow=N, ncol=parms$M)
for(i in seq(parms$M)){
   X[components == i, ] = X[components == i, ] %*% chol(solve(parms$Sinv[, , i]))
}
X = X + parms$mu[components, ]
X = coda::mcmc(X)
results = stats(X)
plot(as.matrix(X))
res[["IID"]] = list(X=X, results=results, sel=1)

# Metropolis Sample
set.seed(1)
scales = list(100)
K = length(scales)
adapt = 0
x0 = rnorm(d, 5, 10)
sig0 = array(0, dim=c(d,d,K))
for(k in seq(K)) sig0 [, , k] = diag(d) * scales[[k]]
mcmc = aMTM::aMTM(
   target=target, 
   x0=x0,
   N=N,
   K=K,
   sig0=sig0,
   adapt=adapt,
   parms=parms,
   burnin=0.2,
   accrate=0.5, global=F, scale=F, local=F, gamma=1.
)
X = mcmc$X
results = stats(X)
plot(as.matrix(X))
res[["Metropolis"]] = list(X=X, results=results, sig=mcmc$Sig, sel=mcmc$sel, lam=mcmc$lam)

# AM one modes
set.seed(1)
scales = list(100)
K = length(scales)
adapt = 1
x0 = rnorm(d, 5, 10)
sig0 = array(0, dim=c(d,d,K))
for(k in seq(K)) sig0 [, , k] = diag(d) * scales[[k]]
mcmc = aMTM::aMTM(
   target=target, 
   x0=x0,
   N=N,
   K=K,
   sig0=sig0,
   adapt=adapt,
   parms=parms,
   burnin=0.2,
   accrate=0.5, global=F, scale=F, local=F, gamma=1.
)
X = mcmc$X
results = stats(X)
plot(as.matrix(X))
res[["AM (one mode)"]] = list(X=X, results=results, sig=mcmc$Sig, sel=mcmc$sel, lam=mcmc$lam)

# AM two modes
set.seed(4)
scales = list(100)
K = length(scales)
adapt = 1
x0 = rnorm(d, 5, 10)
sig0 = array(0, dim=c(d,d,K))
for(k in seq(K)) sig0 [, , k] = diag(d) * scales[[k]]
mcmc = aMTM::aMTM(
   target=target, 
   x0=x0,
   N=N,
   K=K,
   sig0=sig0,
   adapt=adapt,
   parms=parms,
   burnin=0.2,
   accrate=0.5, global=F, scale=F, local=F, gamma=1.
)
X = mcmc$X
results = stats(X)
plot(as.matrix(X))
res[["AM (two modes)"]] = list(X=X, results=results, sig=mcmc$Sig, sel=mcmc$sel, lam=mcmc$lam)

# MTM
set.seed(1)
scales = list(400, 100, 9)
K = length(scales)
adapt = 0
x0 = rnorm(d, 5, 10)
sig0 = array(0, dim=c(d,d,K))
for(k in seq(K)) sig0 [, , k] = diag(d) * scales[[k]]
mcmc = aMTM::aMTM(
   target=target, 
   x0=x0,
   N=N,
   K=K,
   sig0=sig0,
   adapt=adapt,
   parms=parms,
   burnin=0.2,
   accrate=0.5, global=F, scale=F, local=F, gamma=1.
)
X = mcmc$X
results = stats(X)
plot(as.matrix(X), col=mcmc$sel)
res[["MTM"]] = list(X=X, results=results, sig=mcmc$Sig, sel=mcmc$sel, lam=mcmc$lam)


# aMTM
set.seed(1)
scales = list(400, 50, 50)
K = length(scales)
adapt = 3
x0 = rnorm(d, 5, 10)
sig0 = array(0, dim=c(d,d,K))
for(k in seq(K)) sig0 [, , k] = diag(d) * scales[[k]]
mcmc = aMTM::aMTM(
   target=target, 
   x0=x0,
   N=N,
   K=K,
   sig0=sig0,
   adapt=adapt,
   parms=parms,
   burnin=0.5,
   accrate=0.5, global=F, scale=F, local=T, gamma=0.7, weight=-1
)
X = mcmc$X
results = stats(X)
plot(as.matrix(X), col=mcmc$sel)
res[["aMTM"]] = list(X=X, results=results, sig=mcmc$Sig, sel=mcmc$sel, lam=mcmc$lam)


# save results
save(res, file="~/Documents/aMTM-simulations/motivation/results/samples.Rdata")
