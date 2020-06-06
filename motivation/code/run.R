library(batchtools)
library(snow)
setwd("~/Documents/aMTM-simulations/motivation")

reg = makeExperimentRegistry(
   file.dir = NA, 
   seed = 1
)
reg$cluster.functions = makeClusterFunctionsSocket(ncpus = 15)
# add problem
source('code/problems.R')
addProblem(name = "motivation", data = list(target=logp, parms=parms), fun = fun)
batchExport(reg=reg, export=list(stats=stats))
# add algorithms
source("code/algorithms.R")
addAlgorithm(name = "IID", fun = IID.wrapper)
addAlgorithm(name = "Metropolis", fun = aMTM.wrapper)
addAlgorithm(name = "AM", fun = aMTM.wrapper)
addAlgorithm(name = "MTM", fun = aMTM.wrapper)
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


save(reg, file="~/Documents/aMTM-simulations/motivation/results/batchtools/reg.Rdata")

