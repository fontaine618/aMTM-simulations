library(batchtools)
library(xtable)
# load reg
load("~/Documents/aMTM-simulations/banana/results/batchtools/reg.Rdata")
# check errors
findErrors(reg=reg)
getErrorMessages(reg=reg)
# load into data table
id = function(res) res
results = unwrap(reduceResultsDataTable(reg=reg, fun=id))
jobs = getJobPars(ids=results$job.id, reg=reg)
results$adapt = sapply(jobs$algo.pars, function(x) x$adapt)
results$K = sapply(jobs$algo.pars, function(x) x$K)
# show results
results
# table
means = results[
   ,
   lapply(.SD, mean, na.rm=TRUE),
   by=.(adapt, K)
]
means[order(adapt, K)]

ses = results[
   ,
   lapply(.SD, function(x) sd(x, na.rm=TRUE) / sqrt(length(x))),
   by=.(adapt, K)
]


write.csv(means, "./results/means.csv")
write.csv(ses, "./results/ses.csv")
