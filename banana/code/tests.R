library("ggplot2")
library("GGally")

logpi = function(x, p) apply(x,1,function(x) {
    -0.5 * (
        x[1]^2/(p$a1^2) +
        (x[2]-p$B1*x[1]^2)^2 +
        x[3]^2 +
        x[4]^2/(p$a2^2) +
        (x[5]-p$B2*x[4]^2)^2
    )
})

iid = function(n, p){
    x = matrix(rnorm(n*5), n, 5)
    x[, 1] = x[, 1] * p$a1
    x[, 2] = x[, 2] + p$B1 * (x[, 1]^2)
    x[, 4] = x[, 4] * p$a2
    x[, 5] = x[, 5] + p$B2 * (x[, 4]^2)
    return(x)
}


parms = list(a1=1, a2=1, B1=5, B2=1)


sample = iid(1e5, parms)


#
N = 1e5
d = 5
K = 3
x0 = rnorm(d, 0, 1)
sig0 = array(0, dim=c(d,d,K))
scales = 10^seq(2, 0, length.out=K)
if (K==1) scales = c(1.)
for(k in seq(K)) sig0 [, , k] = diag(c(1, 4, 1, 1, 4)) * scales[k]

mcmc = aMTM::aMTM(
    target=logpi, 
    x0=x0,
    N=N,
    K=K,
    sig0=sig0,
    parms=parms,
    adapt=3,
    gamma=0.7,
    weight=0,
    local=F,
    global=F, 
    accrate=0.2
)
sample = mcmc$X

# plot
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

