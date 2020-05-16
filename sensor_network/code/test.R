setwd("~/git/aMTM/simulations/sensor_network")
source('code/problems.R')


# test evaluation
K = 3
x = Xu + rnorm(2*n_unknown*K, 0, 0.5)
x = matrix(x, K, 2*n_unknown, byrow=T)
data$target(x,data$parms)

# plot the sensors and detections
plot_target = function(){   
   par(mfrow=c(1,1))
   plot(X, xlim=c(-1,1.5), ylim=c(-0.5,1.5), 
        col=c(1, 2, 3, 4, 5, 5), pch=19)
   text(x=X[,1]+0.05, y=X[,2]+0.05, labels=1:n_sensors)
   segments(segs$x0, segs$y0, segs$x1, segs$y1)
   abline(h=0, lty=3)
   abline(v=0, lty=3)
}
plot_target()

n = 1e5
XX = matrix(rnorm(2*n_unknown*n, 0.5, 0.4), n, 2*n_unknown, byrow=T)
evals = data$target(XX,data$parms)
hist(evals)

ids = evals > -100
plot_target()
points(XX[ids,1:2], col=1)
points(XX[ids,3:4], col=2)
points(XX[ids,5:6], col=3)
points(XX[ids,7:8], col=4)

sig1 = round(cov(XX[ids, ]), 2)


# aMTM tests
library(aMTM)

K = 3

log_seq = function(start, end, n){
   l = seq(log(start), log(end), length.out = n)
   exp(l)
}

sig0 = sapply(log_seq(0.1, 0.0001, K), function(sig2){
   diag(rep(1,n_unknown*2)^2) * sig2
}, simplify="array")

sig0[, , 1] = sig1


N = 1e5
x0 = rnorm(2*n_unknown, 0.5, 0.2)

mcmc = aMTM::aMTM(
   target=data$target,
   x0=x0,
   parms=data$parms,
   sig0=sig0,
   K=K,
   N=N,
   global=T,
   local=F,
   scale=F,
   accrate=0.3,
   proposal=3,
   gamma=0.7,
   adapt=3,
   burnin=0.1,
   weight=-1
)

ids = seq(1,N,by=N/1000)
plot_target()
points(as.matrix(mcmc$X[ids,1:2]), col=1)
points(as.matrix(mcmc$X[ids,3:4]), col=2)
points(as.matrix(mcmc$X[ids,5:6]), col=3)
points(as.matrix(mcmc$X[ids,7:8]), col=4)


mcmc$sel.prop
mcmc$acc.rate

round(mcmc$Sig, 5)


data$target(mcmc$X[1:10, ],data$parms)


plot.aMTM(mcmc, vars=c(1,2), type="b")
plot.aMTM(mcmc, vars=c(3,4), type="b")
plot.aMTM(mcmc, vars=c(5,6), type="b")
plot.aMTM(mcmc, vars=c(7,8), type="b")
