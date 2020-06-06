library(batchtools)
library(data.table)
# load reg
setwd("~/Documents/aMTM-simulations/loh")
load("./results/batchtools/reg.Rdata")
# check errors
findErrors(reg=reg)
getErrorMessages(reg=reg)

# summary statistics
id_stats = function(res) res$stats
results = unwrap(reduceResultsDataTable(reg=reg, fun=id_stats))
jobs = getJobPars(ids=results$job.id, reg=reg)

results$N = sapply(jobs$prob.pars, function(x) x$N)
results$init = sapply(jobs$algo.pars, function(x) x$init)
results$adapt = sapply(jobs$algo.pars, function(x) x$adapt)
results$K = sapply(jobs$algo.pars, function(x) x$K)
results$gamma = sapply(jobs$algo.pars, function(x) x$gamma)
results$target_accrate = sapply(jobs$algo.pars, function(x) x$accrate)
results$weight = sapply(jobs$algo.pars, function(x) x$weight)


# estimates
id_estimates = function(res) res$estimates
estimates = unwrap(reduceResultsDataTable(reg=reg, fun=id_estimates))
est = loadResult(1, reg=reg)$estimates
rows = rownames(est)
cols = colnames(est)
name = paste(
    cols[floor((seq(length(cols)*length(rows))-1)/length(rows))+1],
    rep(rows, length(cols)),
    sep="_"
)
colnames(estimates) = c("job.id", name)

results = cbind(results, estimates[, names(estimates) %in% name, with=FALSE])


# average over replications

means = results[
    ,
    lapply(.SD, mean, na.rm=TRUE),
    by=list(init, adapt, K, gamma, weight, target_accrate)
]

ses = results[
    ,
    lapply(.SD, function(x) sd(x, na.rm=TRUE) / sqrt(length(x))),
    by=list(init, adapt, K, gamma, weight, target_accrate)
]


write.csv(means, "./results/means.csv")
write.csv(ses, "./results/ses.csv")



