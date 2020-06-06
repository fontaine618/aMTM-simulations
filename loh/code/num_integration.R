library(batchtools)
library(snow)
setwd("~/Documents/aMTM-simulations/loh")
# add problem
source('code/problems.R')

res = 100

grid_pi = expand.grid(
    pi1=seq(0.001, 0.999, length.out = res),
    pi2=seq(0.001, 0.999, length.out = res)
)

grid_eta_gam = expand.grid(
    eta=seq(0.001, 0.999, length.out = res),
    gam=seq(-29.999, 29.999, length.out = res)
)

out = apply(grid_eta_gam, 1, function(x){
    eta = x[1]
    gam = x[2]
    grid = cbind(eta, grid_pi, gam)
    mode2 = grid[, 2] > 0.4
    evals = logp(as.matrix(grid), parms)
    ps = exp(evals)
    weight = sum(ps)
    weight2 = sum(ps[mode2])
    means = apply(grid * matrix(ps, ncol=1), 2, sum)
    means1 = apply(grid[!mode2, ] * matrix(ps[!mode2], ncol=1), 2, sum)
    means2 = apply(grid[mode2, ] * matrix(ps[mode2], ncol=1), 2, sum)
    c(weight=weight, weight2=weight2, means, means1, means2)
})
summed = apply(out, 1, sum)
results = t(t(summed / summed[1]))
write.csv(results, "./results/ses.csv")