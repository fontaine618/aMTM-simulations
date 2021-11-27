library("ggplot2")
library("GGally")

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

# iid = function(n, p){
#     x = matrix(rnorm(n*5), n, 5)
#     x[, 1] = x[, 1] * p$a1
#     x[, 2] = x[, 2] + p$B1 * (x[, 1]^2)
#     x[, 4] = x[, 4] * p$a2
#     x[, 5] = x[, 5] + p$B2 * (x[, 4]^2)
#     return(x)
# }
# 
# sample = iid(1e5, parms)


#
set.seed(1)
N = 1e5
K = 5
x0 = rnorm(d, 0, 1)
sig0 = array(0, dim=c(d,d,K))
scales = 10^seq(-2, 2, length.out=K)
if (K==1) scales = c(1.)
for(k in seq(K)) sig0 [, , k] = diag(c(1, 3, 1, 1, 3)) * scales[k]

mcmc = aMTM::aMTM(
    target=logpi, 
    x0=x0,
    N=N,
    K=K,
    sig0=sig0,
    parms=parms,
    proposal=2,
    adapt=2,
    gamma=0.5,
    weight=0,
    accrate=0.2,
    burnin=1/2,
    local=T
)

# final covariances
round(mcmc$Sig,1)

# plot
sample = mcmc$X
colors = regions(sample)$regions
data = data.frame(sample, Region=colors)[seq(1, N, 100), ]
data$Region = as.factor(data$Region)
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#000000", "#0072B2", "#D55E00", "#CC79A7")
lowerFn = function(data, mapping){
    ggplot(data=data, mapping=mapping) + geom_point(alpha=0.3) +
        theme(panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank())
}
g = ggpairs(
    data, 
    columns=1:5, 
    upper=list(continuous=lowerFn, mapping=aes(colour=Region)),
    diag="blank",
    lower=list(continuous=lowerFn, mapping=aes(colour=Region)),
    legend=6
)
g = g + scale_colour_manual(values=cbPalette)
g


regions(sample)$weights

