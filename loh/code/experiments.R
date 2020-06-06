library(data.table)
pdes = list("loh" = CJ(N=c(1e4)))
ades = list(
   aMTM = CJ(
       init=c("oracle", "scale", "constant"),
       adapt=0:3,
       K=1:3,
       gamma=c(0.5, 0.7, 1.0),
       accrate=c(0.2, 0.3, 0.4, 0.5),
       weight=c(-1, 0)
   )
)
