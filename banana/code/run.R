library(batchtools)
library(snow)
setwd("~/Documents/aMTM-simulations/banana")

reg = makeExperimentRegistry(
   file.dir = NA, 
   seed = 1
)
reg$cluster.functions = makeClusterFunctionsSocket(ncpus = 15)
# add problem
source('code/problems.R')
addProblem(name="banana", data=data, fun=fun)
batchExport(reg=reg, export=list(
    stats=stats, regions=regions, 
    oracle_cov=oracle_cov, iid_weights=iid_weights
))
# add algorithms
source("code/algorithms.R")
addAlgorithm(name = "aMTM", fun = aMTM.wrapper)
# add experiments
source("code/experiments.R")
addExperiments(pdes, ades, repls = 100)

summarizeExperiments(by = c("problem", "algorithm"))

# testJob(1, external = FALSE)
# testJob(1, external = TRUE)
# testJob(10, external = FALSE)
# testJob(10, external = TRUE)
# testJob(18, external = FALSE)
# testJob(18, external = TRUE)

submitJobs()
getStatus()


save(reg, file="~/Documents/aMTM-simulations/banana/results/batchtools/reg.Rdata")

