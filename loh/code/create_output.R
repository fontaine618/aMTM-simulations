library(data.table)
library(xtable)
setwd("~/Documents/aMTM-simulations/loh")

# import
means = data.table(read.csv("./results/means.csv"), stringsAsFactors=TRUE)
means[, c("X", "job.id"):=NULL]
ses = data.table(read.csv("./results/ses.csv"), stringsAsFactors=TRUE)
ses[, c("X", "job.id"):=NULL]

# 
to_parms_matrix = function(row){
    mat = matrix(row[,11:74], nrow=16, ncol=4, byrow=FALSE)
    rownames(mat) = c(
        "m", "mse", "sd", "qL", "qU",
        "m1", "m1se", "sd1", "q1L", "q1U",
        "m2", "m2se", "sd2", "q2L", "q2U",
        "acf"
    )
    colnames(mat) = c("eta", "pi1", "pi2", "gamma")
    storage.mode(mat) = "numeric"
    data.frame(mat)
}
get_row = function(res, i, a, k, g, w, ta){
    res[init==i & adapt==a & K==k & 
        gamma==g & weight==w & target_accrate==ta
    ]
}

ram3m = to_parms_matrix(get_row(means, "scale", 3, 3, 1.0, -1, 0.2))
ram3s = to_parms_matrix(get_row(ses, "scale", 3, 3, 1.0, -1, 0.2))

# K=1, adapt
am = means[, adapt==1 & K==1 & init=="constant"]
aswam = means[, adapt==2 & K==1 & init=="constant"]
ram = means[, adapt==3 & K==1 & init=="constant"]

# best AM
cols = colnames(means)[seq(9)]
am_best = means[, am & prop1 < 0.3]
# best ASWAM
aswam_best = means[, aswam & prop1 < 0.28]
# best RAM
ram_best = means[, ram & prop1 < 0.1]

# K=3, no adapt (MTM)
mtm = means[, adapt==0 & K==3 & target_accrate==0.2 & gamma==1.0]

# selected aMTM
aMTM_selected = means[, init=="scale" & K==3 & adapt==3 & gamma==1.0 & weight==-1& target_accrate==0.2]

means[ram_best | am_best | aswam_best, ..cols]
ses[ram_best | am_best | aswam_best, ..cols]

means[mtm, ..cols]
ses[mtm, ..cols]

means[aMTM_selected, ..cols]
ses[aMTM_selected, ..cols]


