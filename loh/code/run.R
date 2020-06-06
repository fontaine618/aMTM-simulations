library(batchtools)
library(snow)
setwd("~/Documents/aMTM-simulations/loh")

reg = makeExperimentRegistry(file.dir = NA, seed = 1)
reg$cluster.functions = makeClusterFunctionsSocket(ncpus = 15)
# add problem
source('code/problems.R')
addProblem(name="loh", data=list(target=logp, parms=parms), fun=fun)
batchExport(reg=reg, export=list(
    stats=stats, Sig_oracle=Sig_oracle,
    Sig_constant=Sig_constant, Sig_scale=Sig_scale
))
# add algorithms
source("code/algorithms.R")
addAlgorithm(name="aMTM", fun=aMTM.wrapper)
# add experiments
source("code/experiments.R")
addExperiments(pdes, ades, repls = 100)

summarizeExperiments(by = c("problem", "algorithm"))

testJob(1, external = FALSE)
testJob(1, external = TRUE)

submitJobs()
getStatus()

save(reg, file="./results/batchtools/reg.Rdata")
