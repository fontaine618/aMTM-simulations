library(batchtools)
library(xtable)
# load reg
load("~/Documents/aMTM-simulations/motivation/results/batchtools/reg.Rdata")
# check errors
findErrors(reg=reg)
getErrorMessages(reg=reg)
# load into data table
id = function(res) res$results
results = unwrap(reduceResultsDataTable(reg=reg, fun = id))
results$algorithm = getJobPars(ids=results$job.id, reg=reg)$algorithm
# show results
results
# table
cols = c("accrate", "msjd", "Pr5", "quantile_l5.99", "mean_g5.2")
out = results[, 
   sapply(.SD, function(x) list(
      mean=mean(x, na.rm=T), 
      sd=sd(x, na.rm=T)/sqrt(sum(1-is.na(x)))
      )),
   .SDcols = cols, by = list(algorithm)
   ][, 
     setnames(.SD, paste("V", seq(length(cols)*2), sep=""), 
              paste(rep(cols, each=2), c("mean", "se"), sep = "_"))
     ]
# to latex
out = data.frame(out)
tab = data.frame(row.names = out$algorithm)
for(i in seq(length(cols))){
   mean = format(round(out[, paste(cols[i], "mean", sep="_")], 3), digits=3)
   sd = format(round(out[, paste(cols[i], "se", sep="_")], 3), digits=3)
   tab[, 2*i-1] = mean
   tab[, 2*i] = paste(" (", sd, ")", sep="")
}

xtab = xtable(
   tab,
   caption="",
   label="tab:motiv.stats",
   align="lcccccccccc"
)

print(xtab)
